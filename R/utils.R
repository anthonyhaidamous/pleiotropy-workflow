# =============================================================================
# R/utils.R
# Shared helpers used across the pipeline. Source this at the top of every
# stage script:
#     source("R/utils.R")
#
# Conventions
# -----------
#   - All functions are pure where possible (no implicit side effects).
#   - Functions that need config receive it as an argument, not via globals.
#   - I/O helpers (load_intermediate, save_intermediate) standardise the
#     RDS naming so stages don't have to remember exact filenames.
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(gtools)
})


# -----------------------------------------------------------------------------
# Config loader
# -----------------------------------------------------------------------------
load_config <- function(path) {
  cfg <- yaml::read_yaml(path)

  # Derive a filesystem-safe slug from the phenotype name, used to
  # namespace output directories so multiple phenotypes can coexist
  # without stomping on each other's files.
  if (is.null(cfg$phenotype$name) || !nzchar(cfg$phenotype$name)) {
    stop("config must define phenotype.name (used to namespace outputs)")
  }
  slug <- tolower(cfg$phenotype$name)
  slug <- gsub("[^a-z0-9]+", "_", slug)   # spaces, punctuation -> underscore
  slug <- gsub("^_+|_+$",  "", slug)      # trim leading/trailing underscores
  slug <- gsub("_{2,}",    "_", slug)     # collapse repeated underscores
  if (!nzchar(slug)) {
    stop("phenotype.name '", cfg$phenotype$name,
         "' produces an empty slug after normalisation")
  }
  cfg$phenotype$slug <- slug

  # For each output directory, strip any trailing slashes from the YAML
  # value, then append the phenotype slug. So the YAML can specify a
  # base path like "results/" and the actual paths become e.g.
  # "results/leptin/intermediate", "results/leptin/review", etc.
  for (k in c("intermediate", "review", "figures")) {
    if (!is.null(cfg$paths[[k]])) {
      base <- sub("/+$", "", cfg$paths[[k]])
      cfg$paths[[k]] <- file.path(base, slug)
    }
  }

  # Make sure working directories exist
  for (d in c(cfg$paths$intermediate, cfg$paths$review, cfg$paths$figures)) {
    if (!is.null(d) && !dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  # Carried into the run manifest so a result can be tied to the exact config
  # bytes that produced it, not merely to a filename.
  attr(cfg, "config_path") <- path
  attr(cfg, "config_md5")  <- unname(tools::md5sum(path))
  cfg
}


# -----------------------------------------------------------------------------
# Standardised intermediate I/O
# Stage outputs live in cfg$paths$intermediate as <stage_label>.rds
# -----------------------------------------------------------------------------
save_intermediate <- function(obj, label, cfg) {
  path <- file.path(cfg$paths$intermediate, paste0(label, ".rds"))
  saveRDS(obj, path)
  message("  saved: ", path)
  invisible(path)
}

load_intermediate <- function(label, cfg) {
  path <- file.path(cfg$paths$intermediate, paste0(label, ".rds"))
  if (!file.exists(path)) stop("Missing intermediate: ", path)
  readRDS(path)
}


# -----------------------------------------------------------------------------
# Reviewed-file convention
# Stages with manual checkpoints write <label>_review.tsv. Analyst saves
# their edits as <label>_reviewed.tsv. This helper reads the reviewed file
# if present, else the auto file, else errors.
# -----------------------------------------------------------------------------
load_review <- function(label, cfg, required = TRUE) {
  reviewed <- file.path(cfg$paths$review, paste0(label, "_reviewed.tsv"))
  auto     <- file.path(cfg$paths$review, paste0(label, "_review.tsv"))
  path <- if (file.exists(reviewed)) reviewed else auto
  if (!file.exists(path)) {
    if (required) stop("No review file found at: ", reviewed, " or ", auto)
    return(NULL)
  }
  message("  reading: ", path)
  read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE,
             check.names = FALSE)
}

save_review <- function(df, label, cfg) {
  path <- file.path(cfg$paths$review, paste0(label, "_review.tsv"))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message("  review file: ", path)
  invisible(path)
}


# -----------------------------------------------------------------------------
# chr:pos -> rsID via Ensembl REST
# -----------------------------------------------------------------------------
rsid_from_chrpos <- function(x, build = c("GRCh38", "GRCh37")) {
  build <- match.arg(build)
  base  <- if (build == "GRCh37") "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"
  chr   <- sub("^chr", "", sub(":.*$", "", x), ignore.case = TRUE)
  pos   <- sub("^.*:", "", x)
  url <- sprintf("%s/overlap/region/human/%s:%s-%s?feature=variation;content-type=application/json",
                 base, chr, pos, pos)
  out <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  if (is.null(out) || length(out) == 0) return(x)
  rs <- unique(out$id[grepl("^rs\\d+$", out$id)])
  if (length(rs) == 0) x else rs[1]
}

resolve_chrpos_to_rsid <- function(snp_vec, build = "GRCh38") {
  is_chrpos <- grepl("^chr[0-9XYM]+:[0-9]+$", snp_vec, ignore.case = TRUE)
  if (!any(is_chrpos)) return(snp_vec)
  u <- unique(snp_vec[is_chrpos])
  message("  resolving ", length(u), " chr:pos coordinates via Ensembl ...")
  map <- setNames(vapply(u, rsid_from_chrpos, character(1), build = build), u)
  snp_vec[is_chrpos] <- unname(map[snp_vec[is_chrpos]])
  snp_vec
}


# -----------------------------------------------------------------------------
# Ancestry encoding via cohort-proportion rule
# -----------------------------------------------------------------------------
# Methodology (a convention adopted for this analysis at the analyst's
# discretion, not a GWAS Catalog standard — the 90% threshold in step 4 and the
# descriptor dictionary in step 2 are both adjustable below):
#   For each sample-size string (INITIAL or REPLICATION), separately:
#     1. Parse comma-separated cohorts into (count, descriptor) pairs.
#     2. Map each descriptor to a 1000G superpopulation via the hardcoded
#        dictionary below. Unmappable descriptors become "NR".
#     3. Sum counts per superpop within the field.
#     4. If the largest superpop is >=90% of the field total, emit ONLY that
#        superpop. Otherwise, emit all superpops present in the field.
#   Then union the INITIAL and REPLICATION sets to get the study's final
#   ancestry list.
#
# encode_ancestries_for_rows(df, cfg) returns a character vector of comma-
# joined superpops parallel to df, one entry per row. Works on any data frame
# with `INITIAL SAMPLE SIZE` and `REPLICATION SAMPLE SIZE` columns — no need
# to touch the GWAS Catalog ancestry file.
# -----------------------------------------------------------------------------

# Descriptor -> superpop mapping (verbatim port of the OG hardcoded vectors;
# adjust here if a phenotype-specific descriptor needs reassigning)
.ancestry_descriptor_map <- function() {
  list(
    EUR = c("European", "European -", "British", "European I", "European II",
            "European III", "Finnish", "Old Order Amish", "Erasmus Rucphen",
            "Icelandic", "European Crohn'", "European Erasmus Rucphen",
            "Orcadian", "Northern Finnish", "Caucasian", "Swedish",
            "Turkish", "Greek", "European Crohn' Mexican"),
    AMR = c("Native Hawaiian", "Native American", "Hispanic/Latino",
            "Hispanic -", "Hispanic", "Hispanic Latin American",
            "Hispanic or Latin American", "Mexican", 
            "Admixed African African", "Samoan", "Latino Amercian",
            "Mixed", "Afro-Caribbean", "Mexican-American",
            "Greenlandic Inuit", "Inuit", "Greenlandic"),
    AFR = c("African American", "African American/Afro-Caribbean",
            "African American or Afro-Caribbean",
            "African American -", "African", "Black", "Moroccan",
            "African American Afro-Caribbean", "Jamaican", "Nigerian",
            "Sub-Saharan African", "Zulu", "European African American",
            "Surinamese", "Black British"),
    EAS = c("East Asian", "Japanese", "Asian", "Chinese", "Taiwanese",
            "Asian Native Hawaiian", "Han Chinese", "East Asian American",
            "Korean", "Thai", "Filipino", "Asian Pacific Islander"),
    SAS = c("Indian Asian", "South Asian", "Indian",
            "Central Asian", "South Asian American", "Pakistani",
            "Saudi Arabian", "Bangladeshi", "Asian Oceanian",
            "Central Asian South Asian", "Qatari"),
    NR  = c("NR")
  )
}

# Parse one sample-size string into a tibble of (count, superpop).
# Each comma-separated cohort phrase is classified as either:
#   (a) a NEW cohort: starts with a leading number (optionally prefixed by
#       "up to"). Pull the count, find the descriptor, attach to that cohort.
#   (b) a CONTINUATION cohort: no leading number. Treated as an additional
#       ancestry label for the most recent NEW cohort, with that cohort's
#       count distributed equally across all listed ancestries.
# Example:
#   "56,708 European ancestry, African American or Afro-Caribbean,
#    East Asian ancestry, Hispanic or Latin American individuals"
#   -> 1 cohort of 56,708 listing 4 ancestries -> 14,177 each into
#      EUR, AFR, EAS, AMR.
.parse_cohorts <- function(str, descriptor_map) {
  if (is.na(str) || str == "" || str == "NA") {
    return(tibble::tibble(count = integer(0), superpop = character(0)))
  }
  # Split on commas, but NOT commas that are thousands separators (i.e. comma
  # followed by exactly 3 digits then either end-of-string or non-digit).
  cohorts <- unlist(strsplit(str, ",(?!\\d{3}(?:\\D|$))", perl = TRUE))
  cohorts <- trimws(cohorts)
  cohorts <- cohorts[nzchar(cohorts)]
  if (length(cohorts) == 0) {
    return(tibble::tibble(count = integer(0), superpop = character(0)))
  }

  # Build descriptor -> superpop reverse lookup (longest descriptor first
  # so e.g. "African American" beats "African" when both could match).
  rev_lookup <- unlist(lapply(names(descriptor_map), function(sp) {
    setNames(rep(sp, length(descriptor_map[[sp]])), descriptor_map[[sp]])
  }))
  rev_lookup <- rev_lookup[order(-nchar(names(rev_lookup)))]

  match_descriptor <- function(text) {
    for (descriptor in names(rev_lookup)) {
      pattern <- paste0("\\b", descriptor, "\\b")
      if (grepl(pattern, text, perl = TRUE, ignore.case = FALSE)) {
        return(rev_lookup[[descriptor]])
      }
    }
    NA_character_
  }

  # First pass: classify each chunk as NEW (has leading count) or CONTINUATION.
  # Build a list of cohort groups, each with a count and a vector of superpops.
  groups <- list()
  current <- NULL  # list(count = integer, superpops = character)

  for (chunk in cohorts) {
    # Leading-count regex: tolerates a "softener" prefix (approximately,
    # about, around, up to, upto, at least, at most, ~, >=, etc.) before
    # the number. Captures the number itself.
    softener <- "(?:approximately|about|around|up\\s*to|upto|at\\s+least|at\\s+most|nearly|over|more\\s+than|~|>=|>|<)"
    m <- regmatches(chunk,
                    regexpr(paste0("^(?:", softener, "\\s+)?",
                                   "\\d{1,3}(?:,\\d{3})*|^\\d+"),
                            chunk, perl = TRUE, ignore.case = TRUE))
    has_count <- length(m) > 0 && nchar(m) > 0

    if (has_count) {
      # Flush previous group if any
      if (!is.null(current)) groups[[length(groups) + 1]] <- current
      n <- as.numeric(gsub("[^0-9]", "", m))
      sp <- match_descriptor(chunk)
      current <- list(count = n,
                      superpops = if (is.na(sp)) "NR" else sp)
    } else {
      # Continuation cohort: just add its descriptor to the running group.
      if (is.null(current)) {
        # Continuation with no preceding NEW cohort — anomalous, skip.
        next
      }
      sp <- match_descriptor(chunk)
      if (!is.na(sp)) current$superpops <- unique(c(current$superpops, sp))
    }
  }
  if (!is.null(current)) groups[[length(groups) + 1]] <- current
  if (length(groups) == 0) {
    return(tibble::tibble(count = integer(0), superpop = character(0)))
  }

  # Distribute each group's count equally across its ancestry list.
  out <- lapply(groups, function(g) {
    sps  <- g$superpops
    # If multi-ancestry but every label is NR, keep as one NR row, not divided.
    sps_real <- sps[sps != "NR"]
    if (length(sps_real) == 0) {
      return(tibble::tibble(count = g$count, superpop = "NR"))
    }
    per <- g$count / length(sps_real)
    tibble::tibble(count = rep(per, length(sps_real)), superpop = sps_real)
  })
  bind_rows(out)
}

# Apply the >=90% rule to a single field's parsed cohorts.
# Returns a character vector of superpops (the field's contribution).
.dominant_or_all <- function(parsed, threshold = 0.90) {
  if (nrow(parsed) == 0) return(character(0))
  agg <- parsed %>%
    group_by(superpop) %>%
    summarise(n = sum(count), .groups = "drop")
  total <- sum(agg$n)
  if (total <= 0) return(character(0))
  if (max(agg$n) / total >= threshold) {
    return(agg$superpop[which.max(agg$n)])
  }
  agg$superpop
}

# Vectorised public entry point: takes a data frame with INITIAL/REPLICATION
# columns and returns one comma-joined superpop string per row.
encode_ancestries_for_rows <- function(df, threshold = 0.90) {
  descriptor_map <- .ancestry_descriptor_map()
  init <- df$`INITIAL SAMPLE SIZE`
  repl <- df$`REPLICATION SAMPLE SIZE`
  if (is.null(init)) stop("encode_ancestries_for_rows: missing 'INITIAL SAMPLE SIZE' column")
  if (is.null(repl)) repl <- rep(NA_character_, length(init))

  vapply(seq_along(init), function(i) {
    init_sp <- .dominant_or_all(.parse_cohorts(init[i], descriptor_map), threshold)
    repl_sp <- .dominant_or_all(.parse_cohorts(repl[i], descriptor_map), threshold)
    union_sp <- unique(c(init_sp, repl_sp))
    # NR is "Not Reported" — it can't be predominant. Drop it if anything
    # concrete is also present. (Keeps NR only when it's the sole signal.)
    if (length(union_sp) > 1) union_sp <- union_sp[union_sp != "NR"]
    if (length(union_sp) == 0) return("Check study")
    paste(sort(union_sp), collapse = ",")
  }, character(1))
}


# -----------------------------------------------------------------------------
# Block aggregation helpers
# -----------------------------------------------------------------------------
aggregate_blocks <- function(leadSNPs_LDs) {
  # Aggregate SNPs per lead SNP into block descriptors.
  # Input: data frame with columns query_snp, RS_Number, Coord (chr:pos).
  blocks <- aggregate(leadSNPs_LDs[, c("RS_Number", "Coord")],
                      list(leadSNPs_LDs$query_snp), paste)
  colnames(blocks) <- c("lSNP", "LDblock_SNPs", "allPos")
  blocks <- cbind(blocks[, 1:2], t(sapply(blocks$allPos, function(x) {
    y <- matrix(unlist(strsplit(x, ":")), ncol = length(x), nrow = 2)
    m <- min(as.numeric(y[2, ]))
    n <- max(as.numeric(y[2, ]))
    c(chr = unique(y[1, ]),
      LDblock_pos_min = m,
      LDblock_pos_max = n,
      LDblock_span    = n - m)
  })))
  blocks$chr             <- as.character(blocks$chr)
  blocks$LDblock_pos_min <- as.numeric(as.character(blocks$LDblock_pos_min))
  blocks$LDblock_pos_max <- as.numeric(as.character(blocks$LDblock_pos_max))
  blocks$LDblock_span    <- as.numeric(as.character(blocks$LDblock_span))
  blocks$LDblock_nbSNPs  <- sapply(blocks$LDblock_SNPs, length)
  blocks
}

find_block_overlaps <- function(blocks_df) {
  # Returns a list of overlapping block-pairs by name.
  block_min <- as.numeric(sapply(blocks_df$LDblock_pos_min, function(x) min(as.numeric(x))))
  block_max <- as.numeric(sapply(blocks_df$LDblock_pos_max, function(x) max(as.numeric(x))))
  block_chr <- as.character(sapply(blocks_df$chr, function(x) x[1]))
  block_names <- blocks_df$blockname
  overlaps <- list()
  if (nrow(blocks_df) < 2) return(overlaps)
  for (i in 1:(nrow(blocks_df) - 1)) {
    for (j in (i + 1):nrow(blocks_df)) {
      if (block_chr[i] != block_chr[j]) next
      if (block_min[i] <= block_max[j] && block_min[j] <= block_max[i]) {
        pair <- paste(block_names[i], block_names[j], sep = "_")
        overlaps[[pair]] <- c(block_names[i], block_names[j])
      }
    }
  }
  overlaps
}


# -----------------------------------------------------------------------------
# GxE / multi-trait filtering — default regex patterns
# Returns a list of named patterns. Stage 6 builds the auto-exclude regex
# from these. Override via cfg$gxe_exclusions$extra_patterns if needed.
# -----------------------------------------------------------------------------
default_gxe_patterns <- function() {
  list(
    pleiotropy   = "\\bpleiotrop(?:y|ic)?\\b",
    # The suffix list was exhaustive-by-enumeration and missed the PLURAL:
    # after "interaction" the trailing "s" is a word character, so \b fails and
    # "interactions" did not match. Study titles overwhelmingly use the plural
    # ("Genome-wide physical activity interactions in adiposity",
    # "Genotype-by-environment interactions", "Gene-Alcohol Interactions"), so
    # the pattern missed 19 leptin and 1 CLL rows -- 20 GxE associations
    # surviving the filter the stage is named after. Word-stem instead of a
    # suffix list; verified zero false positives on both retained sets.
    interaction  = "\\binteract\\w*\\b",
    pairwise     = "\\bpairw(?:ize|ise)?\\b",
    variate      = "\\b(?:bi|tri|multi)[ -]?variate\\b",
    multitrait   = "\\bmulti[ -]?trait\\b",
    txt          = "\\b[^,;:]+\\s+[x\u00d7]\\s+[^,;:]+\\b",
    # A ratio trait is a composite phenotype: an association with "X:Y ratio"
    # cannot be attributed to either component, and its direction is ambiguous
    # (a higher ratio can mean numerator up OR denominator down), so it is
    # actively misleading in the effect-direction and concordance analyses.
    #
    # The previous pattern was "\\b\\w+\\s*(?:to|:|/|per)\\s*\\w+\\s*ratio\\b",
    # which required BOTH operands to be a single word sitting immediately
    # either side of "ratio". Real catalog traits are almost never written that
    # way, so it fired on essentially nothing: "phospholipids:total lipids
    # ratio", "triglyceride:HDL cholesterol ratio", "polyunsaturated fatty
    # acids to monounsaturated fatty acids ratio", "waist-hip ratio" and
    # "CD58/ERBB3 protein level ratio in blood" all missed. 66 ratio traits
    # carrying 192 associations survived into the leptin output as a result.
    #
    # Matching the standalone word is sufficient and auditable: every trait
    # naming a ratio IS a ratio of two quantities. Verified against all 66
    # surviving leptin traits with no false positives. Note the word boundary
    # is load-bearing -- a bare "ratio" substring also matches "concentration",
    # "filtration" and "unsaturation".
    ratio        = "\\bratios?\\b",
    # The same composite shape as `ratio`, written the way the Nightingale NMR
    # panel and the blood-count panels write it. "cholesterol to total lipids
    # in large LDL percentage" IS cholesterol / total lipids x 100; the
    # direction ambiguity is identical and only the catalog's wording differs.
    #
    # Three distinct surface forms, all composite:
    #   pct_of_total : "X to Y [in Z] percentage"    (58 leptin, 24 CLL)
    #   pct_of       : "X percentage of Y"           ( 3 leptin,  7 CLL)
    # A BARE "<X> percentage" with no named denominator (e.g. "body fat
    # percentage") is NOT matched here: it is a standalone phenotype with a
    # conventional direction, like waist-hip ratio, and is ruled on via
    # ratio_keep_list rather than by regex.
    #
    # `proportion` and `fraction` carry word boundaries deliberately. A bare
    # "fraction" substring also matches "refraction" ("Abnormality of
    # refraction" is a real retained trait) and "ejection fraction" is a
    # standalone cardiac phenotype, not a composite — same class of failure
    # the `ratio` comment above warns about. Neither word occurs in either
    # retained set today; these are guards for future phenotypes.
    pct_of_total = "\\bto\\b[^,;]*\\b(?:percentage|proportions?|fractions?)\\b",
    pct_of       = "\\b(?:percentage|proportions?|fractions?)\\s+of\\b",
    # Covers three forms the original missed:
    #   "Y-adjusted X"            hyphenated (BMI-adjusted hip circumference)
    #   "X (adjusted)"            parenthetical
    #   the original "adjusted for/by/on"
    # Named composite indices whose adjustment is constitutive of the phenotype
    # definition (WHRadjBMI, ABSI, Hip Index) are exempted by name via
    # ratio_keep_list -- see .context/decisions.md C5.
    conditional  = paste0("\\b(?:adjusted|conditioned|stratified)\\s+(?:for|by|on)\\b",
                          "|\\b\\w+[- ]adjusted\\b",
                          "|\\(\\s*adjusted\\s*\\)"),
    # --- Stratification family (.context/decisions.md, S-series) -------------
    # Stated rule: analyses stratified by a MODIFIABLE EXPOSURE are excluded --
    # a stratified estimate is a gene-environment interaction presented in
    # another form. Analyses stratified by an INVARIANT demographic
    # characteristic (sex, age) are RETAINED: the within-stratum estimate is
    # still a marginal effect, so there is deliberately no pattern for those.
    # Conditioning on genotype is excluded as collider-conditioned.
    strat_genotype = paste0("\\bin\\s+rs[0-9]+\\s+[ACGT]{1,2}\\s+genotype\\b",
                            "|\\b(?:non[- ])?carriers?\\b"),
    strat_smoking  = "\\bin\\s+(?:non[- ]|ever[- ]|never[- ])?smokers?\\b",
    strat_activity = "\\bin\\s+(?:in)?active\\s+individuals\\b",
    strat_adiposity= "\\bin\\s+(?:obese|non[- ]obese|lean|overweight)\\s+individuals\\b",
    # --- Estimand guards -----------------------------------------------------
    # Only MEAN-EFFECT estimates enter the direction-concordance analysis. A
    # variance (vQTL) effect and a principal-component loading have no
    # interpretable direction against a mean effect.
    estimand_variance  = "\\b(?:variance|vQTL|dispersion)\\b",
    estimand_component = paste0("\\bprincipal\\s+component\\s*[0-9]*\\b",
                                "|\\bMOSTest\\b|\\bmultivariate\\b"),
    crosstrait   = "\\bcross[ -]?trait\\b",
    mr           = paste0("\\b(?:mendelian\\s+random(?:iz|is)ation|",
                          "genetically\\s+determined|",
                          "genetic(?:ally)?\\s+(?:instrument|proxy|proxies)|",
                          "drug\\s+target\\s+mr|two[ -]sample\\s+mr|",
                          "instrumental\\s+variable|causal\\s+effect\\s+of)\\b"),
    drugresponse = paste0("\\b(?:in response to|drug[ -]induced|drug response|",
                          "pharmacogen|in treated|following treatment)\\b")
  )
}


# -----------------------------------------------------------------------------
# Patterns that route a row to MANUAL REVIEW rather than auto-exclusion.
# These take precedence over default_gxe_patterns() in Stage 6.
#
# Anthropometric ratios (waist-hip, visceral:abdominal adipose tissue, ...) are
# structurally composite and so match the `ratio` rule, but unlike a lipid or
# pQTL ratio they are long-established standalone anthropometric phenotypes
# with a conventional direction. Whether to keep them is a scientific judgement
# rather than a mechanical one, so they are surfaced in 06_gxe_review.tsv for
# the analyst to decide instead of being dropped silently.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Multi-trait (composite) detection
#
# The methods already state that multi-trait associations are excluded: an
# association reported against "trait A, trait B" cannot be attributed to
# either, exactly like a ratio. The previous implementation treated ANY comma
# in MAPPED_TRAIT as multi-trait and sent it to manual review, which both
# over-reached and under-delivered — a comma is often part of ONE trait's name.
#
# Two kinds of comma:
#   (a) separating genuinely distinct phenotypes
#       "total cholesterol measurement, blood VLDL cholesterol amount"   -> exclude
#   (b) internal to a single trait's name
#       "diabetes mellitus, family history"          (qualifier)
#       "monounsaturated fatty acids; 16:1, 18:1 measurement" (nomenclature)
#                                                                        -> keep
#
# The discriminator is how many PHENOTYPE CORES the string contains, not how
# many commas. A core is a fragment ending in a head-noun that names a measured
# thing. A fragment that is only a qualifier ("family history", "hip") is not a
# core, so "osteoarthritis, hip" is one phenotype while
# "osteoarthritis, hip, osteoarthritis, knee" is two.
#
# Verified against all 38 comma-bearing traits Sébastien flagged `composite`
# in the leptin reference: 35 (a), 3 (b). See
# results/review/leptin/06_composite_split.tsv.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Trait-name normalisation — THE single point of normalisation
#
# Every join, group-by, uniqueness count and keep-list lookup on a trait name
# must go through this function, on BOTH sides of the comparison. The previous
# state of play was that the reference table and the pipeline output matched
# only because NEITHER was normalised: 117 leptin trait names and 37 CLL trait
# names carry leading/trailing whitespace, and the two sides happened to carry
# the same whitespace. That is a latent break, not a resolution — normalising
# either side alone silently loses 117 joins.
#
# Normalisation is deliberately conservative:
#   - trim leading/trailing whitespace
#   - collapse internal whitespace runs to a single space (7 leptin traits
#     carry double spaces, e.g. "fibroblast growth factor 21  level")
#   - fold case
# It does NOT strip punctuation or lowercase for DISPLAY — the returned value
# is a join key, not a label. Keep the original string for output.
#
# Documented in .context/conventions.md §2.
# -----------------------------------------------------------------------------
normalize_trait_name <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  tolower(trimws(gsub("[[:space:]]+", " ", x)))
}


# -----------------------------------------------------------------------------
# Trait-reference reader — THE single point of reference parsing
#
# The reference is a TSV whose writer QUOTES any field containing a comma.
# Reading it with quote = "" (the previous behaviour) treats those quote marks
# as literal data, so a field parses as
#     "diabetes mellitus, family history"      <- quote characters included
# and never matches the catalog's unquoted string. Five leptin traits fell
# through to `new` that way.
#
# The lookup failures were the smaller half. 327 `Hierarchized Trait` values
# are quoted in the leptin reference, and those are not join keys -- they are
# sub-category LABELS, carried untouched into 07_trait_table and everything
# built from it. They would have reached the supplement with literal quote
# marks around them:
#     T0001  Lipids  "Triglycerides and other glycerides (TAG, DAG, MAG)"  ...
#
# The file is also CRLF, so the final column of every row carries a trailing
# \r -- the same class of parse hazard in the same file.
#
# The assertion is the point: a re-export from Excel can silently reintroduce
# either, and nothing downstream would notice. Fail loudly at read time.
# -----------------------------------------------------------------------------
read_trait_reference <- function(path) {
  ref <- read.delim(path, sep = "\t", quote = "\"",
                    stringsAsFactors = FALSE, check.names = FALSE)

  # Strip CR and other control characters left by a CRLF / Excel round-trip.
  chr_cols <- names(ref)[vapply(ref, is.character, logical(1))]
  for (cl in chr_cols) {
    ref[[cl]] <- gsub("[[:cntrl:]]", "", ref[[cl]])
  }

  # Assert: nothing may carry a leading/trailing quote or a control character.
  # The CSV artefact is a MATCHED leading AND trailing quote wrapping the whole
  # value -- that is what an unparsed quoted field looks like. Testing for
  # either alone is too loose: it fires on legitimate prose that happens to end
  # in a quoted phrase (a Flag note ending `...renamed it 'zinc ribbon domain
  # containing 2'` tripped it). Require the pair.
  offenders <- list()
  for (cl in chr_cols) {
    v <- trimws(ref[[cl]])
    wrapped <- !is.na(v) & nchar(v) >= 2 &
      ((startsWith(v, '"')  & endsWith(v, '"')) |
       (startsWith(v, "'")  & endsWith(v, "'")))
    bad <- wrapped | (!is.na(v) & grepl("[[:cntrl:]]", v))
    if (any(bad)) {
      offenders[[cl]] <- unique(v[bad])
    }
  }
  if (length(offenders) > 0) {
    msg <- paste(vapply(names(offenders), function(cl)
      sprintf("  column '%s': %d value(s), e.g. %s", cl, length(offenders[[cl]]),
              encodeString(offenders[[cl]][1], quote = "`")),
      character(1)), collapse = "\n")
    stop("Trait reference at ", path, " has values carrying a leading/trailing ",
         "quote or a control character after parsing:\n", msg,
         "\nThis means the file was written with a quoting convention this ",
         "reader does not expect. Re-export it as plain TSV, or fix ",
         "read_trait_reference() -- do not work around it downstream.")
  }
  ref
}


# -----------------------------------------------------------------------------
# Variant identity — THE single point of variant normalisation
#
# rsIDs are NOT stable across the sources this pipeline joins. LDlink returns
# variants under current/merged dbSNP accessions while the GWAS Catalog holds
# their associations under the accession current at curation time. Joining on
# the rsID string therefore returns nothing, silently, and the association
# never enters the table -- a FALSE NEGATIVE, invisible by construction.
#
# Measured at the FTO locus, leptin:
#   catalog rs1421085 (T/C), 212 associations  <-> block SNP rs17817021 (T/C)
#   catalog rs8050136 (C/A),  23 associations  <-> block SNP rs17817527 (C/A)
#   catalog rs9939609 (T/A),   9 associations  <-> block SNP rs58374961 (T/A)
# Same position, same alleles, different accession, zero retrieved.
# Across the whole set: +941 associations in leptin, +156 in CLL.
#
# Position is the identity; the rsID is a display label. See
# .context/decisions.md B4 and .context/conventions.md.
# -----------------------------------------------------------------------------

# Build "chr:pos" from separate chromosome and position vectors. Fails LOUDLY
# on a malformed or missing coordinate rather than emitting a key that silently
# matches nothing -- the failure mode that let a literal "." into a block SNP
# list and contaminated a diagnostic.
make_position_key <- function(chr, pos, context = "position key") {
  chr <- toupper(trimws(sub("^chr", "", as.character(chr), ignore.case = TRUE)))
  pos <- suppressWarnings(as.numeric(as.character(pos)))
  bad <- is.na(chr) | !nzchar(chr) | chr %in% c(".", "NA") |
         is.na(pos) | pos <= 0
  key <- ifelse(bad, NA_character_,
                paste0(chr, ":", format(pos, scientific = FALSE, trim = TRUE)))
  attr(key, "n_bad") <- sum(bad)
  if (any(bad)) {
    message("  ", context, ": ", sum(bad), " / ", length(key),
            " row(s) have no usable coordinate and get NA (they can still ",
            "match on the rsID arm).")
  }
  key
}

# Sentinels that must NEVER be used as a join key on the rsID arm. A literal
# "." is a live value in LDlink output (a variant with no rsID); if it reached
# the join, every no-rsID block variant would match every no-rsID catalog row.
# This is an assertion, not a filter: a valid rsID is ^rs[0-9]+$ and anything
# else is rejected, so a new sentinel form cannot slip through.
is_valid_rsid <- function(x) {
  x <- as.character(x)
  !is.na(x) & grepl("^rs[0-9]+$", trimws(x))
}

# Allele-membership constraint for a positional match.
#
# A full ref/alt key is NOT constructible: the catalog carries only
# `STRONGEST SNP-RISK ALLELE` (one allele, e.g. "rs4929923-T"), and 21% of rows
# report it as "?". So this is a one-sided CONSTRAINT rather than part of the
# key: where the risk allele is known it must be one of the two alleles LDlink
# reports at that coordinate, on either strand.
#
# Measured on both retained sets: of the positional matches with a known risk
# allele, 761/762 (leptin) and 136/136 (CLL) matched directly, 1 matched only
# after a strand flip, and ZERO matched neither. The constraint therefore
# rejects nothing today; it exists to stop a SNP/indel or multi-allelic
# collision at a shared coordinate in future data.
alleles_compatible <- function(risk_allele, allele_pair) {
  ra <- toupper(sub("^.*-", "", as.character(risk_allele)))
  ap <- toupper(gsub("[() ]", "", as.character(allele_pair)))
  a1 <- sub("/.*$", "", ap)
  a2 <- sub("^.*/", "", ap)
  comp <- function(a) chartr("ACGT", "TGCA", a)
  known <- ra %in% c("A", "C", "G", "T")
  pair_known <- !is.na(a1) & !is.na(a2) & nzchar(a1) & nzchar(a2)
  # Unknown on either side -> pass. Never lose a row for missing metadata.
  !known | !pair_known |
    ra == a1 | ra == a2 | comp(ra) == a1 | comp(ra) == a2
}


trait_qualifier_terms <- function() {
  c("family history", "hip", "knee", "female", "male", "females", "males",
    "total", "chronic", "acute", "left", "right", "adjusted for bmi",
    "bmi-adjusted", "unspecified", "other", "nos")
}

trait_core_pattern <- function() {
  paste0("(measurement|amount|level|levels|count|counts|ratio|rate|trait|",
         "disease|disorder|mellitus|behaviour|behavior|dependence|quality|",
         "arthroplasty|osteoarthritis|hypertension|cancer|carcinoma|anaemia|",
         "anemia|abscess|cellulitis|gout|inhibitor|syndrome|infection|",
         "use measurement|status)\\b")
}

# The head-noun approach FAILED OPEN, and that is why it kept missing traits.
# It required >= 2 fragments to match a hardcoded vocabulary of head-nouns
# ("measurement|amount|level|...|osteoarthritis|cancer|..."). Any phenotype
# whose name is a bare disease term not in that list scored zero cores, so
# "diffuse large B-cell lymphoma, multiple sclerosis" and "Hodgkins lymphoma,
# multiple myeloma, non-Hodgkins lymphoma" were classified single-trait. Over
# the two retained sets the vocabulary caught 0 of 15 comma-bearing traits.
# Extending the word list would only move the boundary, not remove it: every
# future phenotype brings disease names nobody enumerated.
#
# This version FAILS CLOSED. A comma makes a trait composite UNLESS the comma
# is structurally explained as internal to one name. That is the safe
# direction for an exclusion rule: an unrecognised construction is excluded
# and visible in the log, rather than silently retained.
#
# The four structural explanations, all mechanical and closed-form:
#   locant       "5,6-dihydrothymine measurement"        (numeric locant)
#   nomenclature "monounsaturated fatty acids; 16:1, 18:1 measurement"
#   uniprot      "...serine hydroxymethyltransferase, cytosolic in blood"
#   qualifier    "diabetes mellitus, family history"     (qualifier-only tail)
#
# Anything else with a comma is composite by rule. Genuine single-trait names
# that are NOT mechanically explainable — appositions such as "brain volume,
# neuroimaging measurement" or "signal peptide, CUB and EGF-like
# domain-containing protein 1 measurement", where the second fragment renames
# the first rather than adding a phenotype — cannot be told from an
# enumeration by regex, because the distinction is semantic. Those are named
# individually in cfg$gxe_exclusions$multitrait_keep_list, exactly as
# composite ratios are named in ratio_keep_list. The rule is then complete by
# inspection: every comma-bearing trait is either excluded by rule or listed
# in config as a reviewed exception.
is_composite_multitrait <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  qual  <- trait_qualifier_terms()
  parts <- strsplit(x, ",\\s*")

  has_comma <- grepl(",", x)

  # (1) Numeric locant: "5,6-dihydrothymine", "1,5-anhydroglucitol".
  locant <- grepl("(^|[^0-9])[0-9]+,[0-9]+[-[:alpha:]]", x, perl = TRUE)
  # (2) Chemical / lipid nomenclature list: "16:1, 18:1".
  nomenclature <- grepl("[0-9]+:[0-9]+\\s*,\\s*[0-9]+:[0-9]+", x)
  # (3) UniProt subcellular / isoform suffix: ", mitochondrial", ", cytosolic".
  uniprot <- grepl(paste0(",\\s*(mitochondrial|cytosolic|secreted|membrane|",
                          "nuclear|precursor|isoform\\b[^,]*)$|",
                          ",\\s*(mitochondrial|cytosolic|secreted|membrane|",
                          "nuclear|precursor)\\s"),
                   x, ignore.case = TRUE, perl = TRUE)
  # (4) Every post-comma fragment is a bare qualifier.
  qualifier_only <- vapply(parts, function(p) {
    tail <- tolower(trimws(p[-1]))
    length(tail) > 0 && all(tail %in% qual)
  }, logical(1))

  has_comma & !locant & !nomenclature & !uniprot & !qualifier_only
}


default_manual_review_patterns <- function() {
  list(
    anthropometric_ratio = paste0(
      "\\b(?:waist|hip|abdominal|gluteofemoral|visceral|subcutaneous|adipose)\\b",
      "[^,;]*\\bratios?\\b")
  )
}


# -----------------------------------------------------------------------------
# Sample size: parse "X,XXX, Y,YYY" style strings to a numeric total
# -----------------------------------------------------------------------------
sum_sample_sizes <- function(column) {
  # GWAS Catalog INITIAL/REPLICATION SAMPLE SIZE strings list cohorts
  # separated by commas, each cohort starting with a count that may itself
  # contain comma thousands-separators, e.g.:
  #   "86,142 European ancestry female individuals, 557 Old Order Amish ...,
  #    1,103 Erasmus Rucphen ..."
  # Old approach — gsub("[^0-9,]","") then strsplit(",") — broke "86,142" into
  # 86 + 142, severely under-counting samples and dropping ~50% of cross-traits
  # at the >=1000 filter.
  #
  # New regex: match either (a) thousands-grouped numbers (\d{1,3}(?:,\d{3})+)
  # or (b) plain runs of digits, anywhere in the string. Sum all matches.
  # Safe because the SAMPLE SIZE columns contain only cohort counts —
  # no rsIDs, platform IDs, years, or other incidental numbers.
  sapply(column, function(x) {
    if (is.na(x) || x == "") return(0)
    matches <- regmatches(
      x,
      gregexpr("\\d{1,3}(?:,\\d{3})+|\\d+", x, perl = TRUE)
    )[[1]]
    if (length(matches) == 0) return(0)
    sum(as.numeric(gsub(",", "", matches)), na.rm = TRUE)
  }, USE.NAMES = FALSE)
}


# -----------------------------------------------------------------------------
# OR -> beta conversion (Stage 9)
# -----------------------------------------------------------------------------
or_to_beta <- function(or_value) log(as.numeric(or_value))

ci_or_to_beta <- function(ci_text) {
  if (is.na(ci_text)) return(ci_text)
  if (grepl("increase|decrease", ci_text, ignore.case = TRUE)) return(ci_text)
  bounds <- as.numeric(unlist(regmatches(ci_text, gregexpr("[0-9\\.]+", ci_text))))
  if (length(bounds) == 2) {
    sprintf("[%f-%f]", or_to_beta(bounds[1]), or_to_beta(bounds[2]))
  } else {
    ci_text
  }
}


# -----------------------------------------------------------------------------
# p-value pretty-printer for messages
# -----------------------------------------------------------------------------
format_pval <- function(p) {
  if (is.na(p) || p >= 0.001) return(formatC(p, digits = 2, format = "f"))
  parts <- strsplit(formatC(p, digits = 2, format = "e"), "e")[[1]]
  paste0(parts[1], "\u00d710^", as.integer(parts[2]))
}


# -----------------------------------------------------------------------------
# Run manifest
#
# Appends one JSON record per stage execution to
#   <intermediate>/run_manifest.jsonl
# This is the provenance record behind every number the analysis reports: given
# a manifest you can say exactly which commit, config, package set and input
# bytes produced a given row count.
#
# JSON Lines (one object per line, appended) rather than one JSON document, so
# concurrent stages cannot corrupt each other and the file is greppable.
#
# Never fails a stage: manifest writing is wrapped so that a provenance problem
# can never destroy an analysis run.
# -----------------------------------------------------------------------------
.git_sha <- function() {
  sha <- suppressWarnings(tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_))
  dirty <- suppressWarnings(tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0,
    error = function(e) NA))
  if (length(sha) == 0 || is.na(sha[1])) return(list(sha = NA_character_, dirty = NA))
  list(sha = sha[1], dirty = isTRUE(dirty))
}

.file_record <- function(paths) {
  paths <- paths[!is.na(paths)]
  lapply(paths, function(p) {
    if (!file.exists(p)) return(list(path = p, exists = FALSE))
    list(path = p, exists = TRUE,
         bytes = as.numeric(file.info(p)$size),
         md5   = unname(tools::md5sum(p)))
  })
}

write_stage_manifest <- function(stage, cfg, inputs = character(0),
                                 outputs = character(0), counts = list(),
                                 notes = NULL) {
  tryCatch({
    g <- .git_sha()
    pkgs <- c("dplyr", "tidyr", "stringr", "yaml", "jsonlite", "gtools",
              "tibble", "readr", "biomaRt", "LDlinkR", "igraph")
    pkg_versions <- vapply(pkgs, function(p) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
      if (is.na(v)) "" else v
    }, character(1))
    pkg_versions <- pkg_versions[nzchar(pkg_versions)]

    rec <- list(
      stage          = stage,
      phenotype      = cfg$phenotype$name,
      phenotype_slug = cfg$phenotype$slug,
      timestamp      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      git_sha        = g$sha,
      git_dirty      = g$dirty,
      config_file    = attr(cfg, "config_path"),
      config_md5     = attr(cfg, "config_md5"),
      r_version      = R.version.string,
      platform       = R.version$platform,
      packages       = as.list(pkg_versions),
      counts         = counts,
      inputs         = .file_record(inputs),
      outputs        = .file_record(outputs)
    )
    if (!is.null(notes)) rec$notes <- notes

    path <- file.path(cfg$paths$intermediate, "run_manifest.jsonl")
    cat(jsonlite::toJSON(rec, auto_unbox = TRUE, null = "null", digits = NA),
        "\n", file = path, sep = "", append = TRUE)
    message("  manifest: ", path, " (stage ", stage, ")")
    invisible(path)
  }, error = function(e) {
    warning("run manifest not written for stage ", stage, ": ", conditionMessage(e))
    invisible(NULL)
  })
}


# -----------------------------------------------------------------------------
# Convenience: stage banner for log readability
# -----------------------------------------------------------------------------
stage_banner <- function(stage_num, title) {
  bar <- paste(rep("=", 70), collapse = "")
  message("\n", bar, "\nStage ", stage_num, " — ", title, "\n", bar)
}


# -----------------------------------------------------------------------------
# LDlinkR helpers (Stage 2)
# -----------------------------------------------------------------------------

# Map cfg$ensembl$build to LDproxy's genome_build argument.
ldlink_build <- function(cfg) {
  switch(toupper(cfg$ensembl$build),
    "GRCH38" = "grch38",
    "GRCH37" = "grch37",
    stop("Unknown genome build: ", cfg$ensembl$build))
}

# Build the (SNP, population) pair table from the lSNPs data frame.
# Splits "Encoded Populations" on commas and drops NR (no LDlink data) and
# "Check study" placeholders.
build_snp_pop_pairs <- function(lSNPs) {
  ldlink_pops <- c("EUR", "AFR", "EAS", "SAS", "AMR")
  lSNPs %>%
    select(snp = SNPS, pops = `Encoded Populations`) %>%
    mutate(pop = strsplit(pops, ",")) %>%
    tidyr::unnest(pop) %>%
    mutate(pop = trimws(pop)) %>%
    filter(pop %in% ldlink_pops) %>%
    select(snp, pop) %>%
    distinct() %>%
    arrange(snp, pop)
}

# Post-call check: did LDproxy return a usable data frame?
# Monoallelic variants and some API errors come back as malformed frames
# rather than as thrown errors, so we inspect the structure directly.
ldproxy_result_is_valid <- function(res) {
  is.data.frame(res) &&
    nrow(res) > 0 &&
    "RS_Number" %in% names(res) &&
    "R2"        %in% names(res)
}