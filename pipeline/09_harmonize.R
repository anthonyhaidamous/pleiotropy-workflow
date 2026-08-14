# =============================================================================
# Stage 9 — Effect direction and allele harmonization
#
# For each filtered cross-trait, ensure the reported effect (beta) is on a
# common allele scale and a common direction, so that betas across SNPs and
# studies are comparable.
#
# Steps:
#   (a) OR -> beta : if the 95% CI text doesn't mention "increase"/"decrease",
#                    the OR or BETA value is an OR -> log-transform.
#   (b) Sign       : if the CI text mentions "decrease", apply a negative sign
#                    to beta and flip CI bounds.
#   (c) Allele     : match each SNP against the LDproxy table to retrieve
#                    its allele pair (CA1, CA2). If risk allele == CA1, flip
#                    beta + CI to reference CA2 instead. (Convention: emit
#                    effect relative to CA2.)
#   (d) Block dir  : per block, the lead SNP's beta sign defines the
#                    phenotype direction in that block (Positive / Negative).
#
# Input  : results/intermediate/06_crosstraits_filtered.rds
#          results/intermediate/04_blocks.rds
#          results/intermediate/03_leadSNPs_LDs.rds
#          results/intermediate/07_trait_table.rds
#          results/intermediate/08_trait_enrichment.rds   (for fold/p_adj/Stars)
#
# Output : results/intermediate/09_harmonized_trait_table.rds
#          results/intermediate/09_harmonized_trait_table.tsv
#          results/intermediate/09_block_directions.rds
#          results/intermediate/09_block_directions.tsv
#
# Run from the project root:
#     Rscript pipeline/09_harmonize.R config/<trait>.yml
# =============================================================================

source("R/utils.R")
source("R/figure_utils.R")   # extract_ci_bounds, derive_signed_beta

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(9, "Effect direction and allele harmonization")

# -----------------------------------------------------------------------------
# 9.1 — Load
# -----------------------------------------------------------------------------
ct           <- load_intermediate("06_crosstraits_filtered", cfg)
blocks       <- load_intermediate("04_blocks",               cfg)
leadSNPs_LDs <- load_intermediate("03_leadSNPs_LDs",         cfg)
trait_table  <- load_intermediate("07_trait_table",          cfg)
enrichment   <- load_intermediate("08_trait_enrichment",     cfg)

message("Loaded ", nrow(ct), " cross-trait associations to harmonize.")

# Same guard as stage 8: a placeholder category must never reach the master
# table, because every figure downstream groups by Category and would silently
# render "new" as though it were a finding.
.placeholders <- c("new", "", "NA", "na", "TBD", "tbd", "?", "-",
                   "unknown", "Unknown", "TODO", "todo")
.bad <- is.na(trait_table$Category) |
        trimws(trait_table$Category) %in% .placeholders
if (any(.bad)) {
  stop("Stage 9 refuses to run: ", sum(.bad), " of ", nrow(trait_table),
       " traits in 07_trait_table carry a placeholder category. ",
       "Complete the stage-7 review first — see stage 8's message for the ",
       "exact file.")
}

# -----------------------------------------------------------------------------
# 9.2 — Drop rows with unknown risk allele ("?" / "-")
# -----------------------------------------------------------------------------
ct <- ct %>%
  tidyr::separate(`STRONGEST SNP-RISK ALLELE`,
                  into   = c("snp_for_RA", "RiskAllele"),
                  sep    = "-",
                  remove = FALSE,
                  fill   = "right",
                  extra  = "merge") %>%
  filter(!RiskAllele %in% c("?", "-", "")) %>%
  filter(!is.na(RiskAllele))

message("After dropping unknown-risk-allele rows: ", nrow(ct), " remaining.")

# -----------------------------------------------------------------------------
# 9.3 — OR -> beta conversion (driven by 95% CI text)
# -----------------------------------------------------------------------------
# Helper: a CI text mentioning "increase" or "decrease" indicates beta scale;
# anything else (or empty) is OR scale.
ci_is_beta <- function(ci) grepl("increase|decrease", ci, ignore.case = TRUE)
ci_is_decrease <- function(ci) grepl("decrease", ci, ignore.case = TRUE)

ct$beta <- as.numeric(ct$`OR or BETA`)

# (a) OR -> log(OR) where CI text doesn't say increase/decrease
or_rows <- !ci_is_beta(ct$`95% CI (TEXT)`) & !is.na(ct$beta)
ct$beta[or_rows] <- log(ct$beta[or_rows])

# (b) Apply negative sign where CI text says "decrease"
dec_rows <- ci_is_decrease(ct$`95% CI (TEXT)`) & !is.na(ct$beta)
ct$beta[dec_rows] <- -abs(ct$beta[dec_rows])

# Extract numeric CI bounds from the CI text (regardless of OR/beta scale).
# extract_ci_bounds() lives in R/figure_utils.R so that this stage and the
# forest-plot figures cannot drift apart again.
ci_bounds <- t(vapply(ct$`95% CI (TEXT)`, extract_ci_bounds, numeric(2)))
ct$CI_lower <- ci_bounds[, 1]
ct$CI_upper <- ci_bounds[, 2]

# Order matters: log-transform OR-scale CIs FIRST (must be positive),
# then apply the decrease-sign for beta-scale CIs.
# Guard against malformed CI text producing 0/negative/NaN bounds.
or_lower_ok <- or_rows & !is.na(ct$CI_lower) & ct$CI_lower > 0
or_upper_ok <- or_rows & !is.na(ct$CI_upper) & ct$CI_upper > 0
ct$CI_lower[or_lower_ok] <- log(ct$CI_lower[or_lower_ok])
ct$CI_upper[or_upper_ok] <- log(ct$CI_upper[or_upper_ok])
# Anything left as a non-log-able OR-scale CI gets NA so it doesn't
# contaminate downstream comparisons.
ct$CI_lower[or_rows & !or_lower_ok] <- NA_real_
ct$CI_upper[or_rows & !or_upper_ok] <- NA_real_

# Apply negative sign to CIs where CI text mentioned "decrease"
# (these are by definition beta-scale, so or_rows == FALSE for them)
ct$CI_lower[dec_rows] <- -abs(ct$CI_lower[dec_rows])
ct$CI_upper[dec_rows] <- -abs(ct$CI_upper[dec_rows])

# -----------------------------------------------------------------------------
# 9.4 — Allele harmonization: pull (CA1, CA2) from LDproxy, flip if needed
#
# Each cross-trait row has a SNP and a risk allele. We match the SNP against
# its row in leadSNPs_LDs (which carries the alleles in "Alleles" as A/B).
# Then if the risk allele matches CA1, we flip the beta and CI to reference
# CA2; if it matches CA2, we leave it. Convention: report effect relative
# to CA2 (the second allele in the LDlink Alleles column).
# -----------------------------------------------------------------------------
allele_lookup <- leadSNPs_LDs %>%
  distinct(RS_Number, Alleles) %>%
  filter(!is.na(Alleles), Alleles != "") %>%
  tidyr::separate(Alleles, into = c("CA1_raw", "CA2_raw"),
                  sep = "/", remove = FALSE, extra = "merge") %>%
  mutate(CA1 = gsub("[^A-Za-z]", "", CA1_raw),
         CA2 = gsub("[^A-Za-z]", "", CA2_raw)) %>%
  select(RS_Number, CA1, CA2)

ct <- ct %>%
  left_join(allele_lookup, by = c("SNP" = "RS_Number"))

# Diagnostic: how many rows lack CA1/CA2?
n_no_alleles <- sum(is.na(ct$CA1) | is.na(ct$CA2))
if (n_no_alleles > 0) {
  message("Note: ", n_no_alleles,
          " cross-trait rows have no allele info from LDproxy ",
          "(harmonization will pass beta through unchanged for these).")
}

# Flip beta + CI when risk allele == CA1; emit effect relative to CA2
flip_idx <- !is.na(ct$CA1) & !is.na(ct$RiskAllele) & ct$RiskAllele == ct$CA1
message("Flipping beta/CI for ", sum(flip_idx),
        " rows where risk allele == CA1.")

ct$harmonized_beta     <- ct$beta
ct$harmonized_CI_lower <- ct$CI_lower
ct$harmonized_CI_upper <- ct$CI_upper

ct$harmonized_beta[flip_idx]     <- -ct$beta[flip_idx]
# When flipping the beta sign, CI bounds also flip sign AND swap order
new_lower <- -ct$CI_upper[flip_idx]
new_upper <- -ct$CI_lower[flip_idx]
ct$harmonized_CI_lower[flip_idx] <- new_lower
ct$harmonized_CI_upper[flip_idx] <- new_upper

# Reference allele (the allele the harmonized beta is now relative to)
ct$ref_allele <- ifelse(!is.na(ct$CA2), ct$CA2, NA_character_)

# Drop rows without allele info — their betas can't be harmonized to a
# common reference, so leaving them in would mix harmonized and unharmonized
# effects in downstream figures.
n_before_drop <- nrow(ct)
ct <- ct[!is.na(ct$CA1) & !is.na(ct$CA2), ]
n_dropped <- n_before_drop - nrow(ct)
if (n_dropped > 0) {
  message("Dropped ", n_dropped,
          " cross-trait rows lacking allele info from LDproxy.")
}

# Safety: ensure CI lower <= upper
swap_idx <- !is.na(ct$harmonized_CI_lower) & !is.na(ct$harmonized_CI_upper) &
            ct$harmonized_CI_lower > ct$harmonized_CI_upper
if (any(swap_idx)) {
  tmp <- ct$harmonized_CI_lower[swap_idx]
  ct$harmonized_CI_lower[swap_idx] <- ct$harmonized_CI_upper[swap_idx]
  ct$harmonized_CI_upper[swap_idx] <- tmp
}

# -----------------------------------------------------------------------------
# 9.5 — Per-block phenotype direction (from top lead SNP's effect)
#
# blocks$Top_lead_SNP encodes "rs..., pos, pval, OR-or-beta". The 4th token
# is the raw catalog "OR or BETA" value — but ORs are always positive, so we
# can't infer direction from sign alone. We look the SNP up in the original
# lSNPs table to read its CI text, then apply the same OR -> beta logic
# used for cross-traits.
# -----------------------------------------------------------------------------
lSNPs <- load_intermediate("01_lSNPs", cfg)

# derive_signed_beta() lives in R/figure_utils.R, shared with
# analysis/diagnose_allele_alignment.R.

parse_top_snp <- function(s) {
  if (is.na(s) || s == "") return(NA_character_)
  trimws(unlist(strsplit(s, ","))[1])
}

block_dir <- blocks %>%
  mutate(top_lSNP_rsid = vapply(Top_lead_SNP, parse_top_snp, character(1)))

# Vectorized lookups (avoid rowwise() + $ interaction issues)
match_idx <- match(block_dir$top_lSNP_rsid, lSNPs$SNPS)
block_dir$top_lSNP_OR_beta <- lSNPs$`OR or BETA`[match_idx]
block_dir$top_lSNP_CI_text <- lSNPs$`95% CI (TEXT)`[match_idx]

# Per-row signed-beta derivation
block_dir$top_lSNP_beta <- mapply(
  derive_signed_beta,
  block_dir$top_lSNP_OR_beta,
  block_dir$top_lSNP_CI_text,
  USE.NAMES = FALSE
)

block_dir$phenotype_direction <- with(block_dir, case_when(
  is.na(top_lSNP_beta) ~ NA_character_,
  top_lSNP_beta < 0    ~ "Negative",
  TRUE                 ~ "Positive"
))

block_dir <- block_dir %>%
  select(Block, top_lSNP_rsid, top_lSNP_OR_beta, top_lSNP_CI_text,
         top_lSNP_beta, phenotype_direction)

n_pos <- sum(block_dir$phenotype_direction == "Positive", na.rm = TRUE)
n_neg <- sum(block_dir$phenotype_direction == "Negative", na.rm = TRUE)
n_unk <- sum(is.na(block_dir$phenotype_direction))
message("Block phenotype direction (from top lead SNP): ",
        n_pos, " positive, ", n_neg, " negative, ", n_unk, " indeterminate.")

# -----------------------------------------------------------------------------
# 9.6 — Build the harmonized_trait_table (the master table for figures)
# -----------------------------------------------------------------------------
# The join key has to be normalised first, or the join quietly loses rows.
#
# Stage 5 keeps MAPPED_TRAIT exactly as the GWAS Catalog writes it, and the
# Catalog's export carries a trailing space on some traits: 13 of the 707 CLL
# cross-trait rows, 323 of leptin's 3290. Stage 7 trims when it aggregates, so
# the trait table holds the clean form — and an exact-string join therefore
# matches neither, handing those rows a blank Category instead of failing. In
# CLL that silently dropped 13 lipid associations out of the category figures
# (block B23 showed 27 lipid traits where it has 38); in leptin it stranded 60
# traits.
#
# Trimming is safe as well as necessary: no trait occurs in both forms in
# either phenotype, so this merges no rows and changes no count. It only lets
# keys that were always meant to be equal meet.
#
# The trimmed value is what goes on into the output, so 09 and 07 agree on how
# a trait is spelled. The only genuinely uncategorised traits left after this
# are reported by the check below.
ct$MAPPED_TRAIT           <- trimws(ct$MAPPED_TRAIT)
trait_table$MAPPED_TRAIT  <- trimws(trait_table$MAPPED_TRAIT)
enrichment$MAPPED_TRAIT   <- trimws(enrichment$MAPPED_TRAIT)

# Join: cross-traits + categories + block direction + enrichment.
ht <- ct %>%
  # Trait categories
  left_join(trait_table %>%
              select(MAPPED_TRAIT, Category, `Hierarchized Trait`),
            by = "MAPPED_TRAIT") %>%
  # Block direction
  left_join(block_dir %>% select(Block, phenotype_direction),
            by = "Block") %>%
  # Trait-level enrichment
  left_join(enrichment %>%
              select(MAPPED_TRAIT, `Fold Enrichment`,
                     `Normalized Enrichment`,
                     P_Value, P_Adjusted, Stars),
            by = "MAPPED_TRAIT") %>%
  # Tag SNPs as lead vs proxy (already on ct via SNP_type from Stage 5)
  rename(snp_type = SNP_type)

# Final column order (drop intermediate cols that are no longer needed)
keep_cols <- c(
  "Block", "phenotype_direction",
  "SNP", "snp_type", "RiskAllele", "ref_allele",
  "CA1", "CA2",
  "Category", "Hierarchized Trait", "MAPPED_TRAIT",
  "DISEASE/TRAIT", "STUDY",
  "P-VALUE",
  "OR or BETA", "beta",
  "harmonized_beta", "harmonized_CI_lower", "harmonized_CI_upper",
  "95% CI (TEXT)",
  "Encoded Populations",
  "INITIAL SAMPLE SIZE", "REPLICATION SAMPLE SIZE", "total_sample_size",
  "MAPPED_GENE", "PUBMEDID", "FIRST AUTHOR", "DATE", "STUDY ACCESSION",
  "Fold Enrichment", "Normalized Enrichment",
  "P_Value", "P_Adjusted", "Stars"
)
keep_cols <- intersect(keep_cols, names(ht))
ht <- ht[, keep_cols]

message("Built harmonized trait table: ", nrow(ht), " rows, ",
        ncol(ht), " columns.")

# A row with no Category is a cross-trait that never reached the trait table,
# and every category figure downstream drops it. That is a reasonable outcome
# for a trait the analyst chose not to categorise, and an invisible data loss
# otherwise, so it is named here rather than discovered in a figure.
.uncat <- is.na(ht$Category) | !nzchar(trimws(ht$Category))
if (any(.uncat)) {
  .traits <- sort(unique(ht$MAPPED_TRAIT[.uncat]))
  message("NOTE: ", sum(.uncat), " row(s) have no Category — ", length(.traits),
          " trait(s) absent from 07_trait_table. Category figures will drop them:")
  message("      ", paste(.traits, collapse = " | "))
}

# -----------------------------------------------------------------------------
# 9.7 — Save
# -----------------------------------------------------------------------------
save_intermediate(ht,        "09_harmonized_trait_table", cfg)
save_intermediate(block_dir, "09_block_directions",       cfg)

ht_path  <- file.path(cfg$paths$intermediate, "09_harmonized_trait_table.tsv")
bd_path  <- file.path(cfg$paths$intermediate, "09_block_directions.tsv")
write.table(ht,        ht_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
write.table(block_dir, bd_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
message("  saved: ", ht_path)
message("  saved: ", bd_path)

# -----------------------------------------------------------------------------
# 9.8 — Summary
# -----------------------------------------------------------------------------
bar <- paste(rep("\u2500", 60), collapse = "")
message("\n", bar)
message("Harmonization summary")
message(bar)
message("Cross-traits in (post Stage 6) : ", nrow(ct))
message("Cross-traits with valid beta   : ",
        sum(!is.na(ht$harmonized_beta)))
message("Cross-traits with allele info  : ",
        sum(!is.na(ht$ref_allele)))
message("Blocks with direction assigned : ",
        sum(!is.na(block_dir$phenotype_direction)),
        " of ", nrow(block_dir))
message(bar)

write_stage_manifest(
  stage   = 9,
  cfg     = cfg,
  inputs  = file.path(cfg$paths$intermediate,
                      c("06_crosstraits_filtered.rds", "04_blocks.rds",
                        "03_leadSNPs_LDs.rds", "07_trait_table.rds",
                        "08_trait_enrichment.rds", "01_lSNPs.rds")),
  outputs = file.path(cfg$paths$intermediate,
                      c("09_harmonized_trait_table.rds", "09_harmonized_trait_table.tsv",
                        "09_block_directions.rds", "09_block_directions.tsv")),
  counts  = list(
    n_rows            = nrow(ht),
    n_cols            = ncol(ht),
    n_valid_beta      = sum(!is.na(ht$harmonized_beta)),
    n_with_alleles    = sum(!is.na(ht$ref_allele)),
    n_blocks          = nrow(block_dir),
    n_blocks_directed = sum(!is.na(block_dir$phenotype_direction)),
    n_flipped         = sum(flip_idx)
  ))

message("\nStage 9 complete.")