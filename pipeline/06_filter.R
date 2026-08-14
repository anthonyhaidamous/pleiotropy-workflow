# =============================================================================
# Stage 6 — Cross-trait filtering
#
# Four sequential filters, in order:
#   (1) Ancestry concordance — cross-trait population must overlap block pop
#   (2) Total sample size    — initial + replication >= cfg$filters$min_sample_size
#   (3) CNV exclusion        — keep only CNV == "N"
#   (4) GxE / multi-trait    — two-tier regex: auto-exclude + manual review
#
# Manual checkpoint
# -----------------
# Filter (4) flags rows for either auto_exclude or manual_check. The manual
# review file lists the manual_check rows with disease/trait, mapped trait,
# and study fields. After review, save the file as 06_gxe_reviewed.tsv with
# the gxe_decision column set to either "keep" or "auto_exclude".
#
# Re-running this stage:
#   - First run produces 06_gxe_review.tsv (auto-flagged manual_check rows).
#   - Save edits as 06_gxe_reviewed.tsv. Re-run picks the reviewed file up
#     automatically and applies your decisions.
#
# Input  : results/intermediate/05_crosstraits_raw.rds
#          results/intermediate/04_blocks.rds        (for block ancestry)
#
# Output : results/intermediate/06_crosstraits_filtered.rds
#          results/intermediate/06_crosstraits_filtered.tsv
#          results/review/06_gxe_review.tsv          (manual_check rows)
#          results/intermediate/06_filtering_log.tsv (per-filter row counts)
#
# Run from the project root:
#     Rscript pipeline/06_filter.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(6, "Cross-trait filtering")

# -----------------------------------------------------------------------------
# 6.1 — Load
# -----------------------------------------------------------------------------
ct     <- load_intermediate("05_crosstraits_raw", cfg)
blocks <- load_intermediate("04_blocks",          cfg)

# Make sure we have block ancestry. Stage 5 attaches it as `Block Population`,
# but if that's missing we recover from blocks$Ethnicities.
if (!"Block Population" %in% names(ct)) {
  block_anc_lookup <- setNames(blocks$Ethnicities, blocks$Block)
  ct$`Block Population` <- block_anc_lookup[ct$Block]
}

initial_n <- nrow(ct)
log_rows <- list(
  list(filter = "0_input", n_before = NA_integer_, n_after = initial_n,
       n_dropped = 0L)
)
message("Starting from ", initial_n, " cross-trait associations.")

# -----------------------------------------------------------------------------
# 6.2 — Filter 1: ancestry concordance
# A row is kept if any superpop in the cross-trait's `Encoded Populations`
# is ALSO in the block's `Block Population`.
# -----------------------------------------------------------------------------
split_pops <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  parts <- trimws(unlist(strsplit(x, ",")))
  parts <- toupper(parts[parts != ""])
  unique(parts)
}

ancestry_overlap <- mapply(
  function(ct_pops, blk_pops) length(intersect(split_pops(ct_pops),
                                               split_pops(blk_pops))) > 0,
  ct$`Encoded Populations`, ct$`Block Population`
)

n_before <- nrow(ct)
ct <- ct[ancestry_overlap, ]
n_after <- nrow(ct)

message("(1) Ancestry concordance: kept ", n_after, " / ", n_before,
        " (dropped ", n_before - n_after, ").")
log_rows[[length(log_rows) + 1]] <- list(
  filter = "1_ancestry", n_before = n_before, n_after = n_after,
  n_dropped = n_before - n_after
)

# -----------------------------------------------------------------------------
# 6.3 — Filter 2: total sample size >= cfg$filters$min_sample_size
# -----------------------------------------------------------------------------
ct$total_sample_size <- sum_sample_sizes(ct$`INITIAL SAMPLE SIZE`) +
                        sum_sample_sizes(ct$`REPLICATION SAMPLE SIZE`)

n_before <- nrow(ct)
ct <- ct[ct$total_sample_size >= cfg$filters$min_sample_size, ]
n_after <- nrow(ct)

message("(2) Sample size \u2265 ", cfg$filters$min_sample_size, ": kept ",
        n_after, " / ", n_before, " (dropped ", n_before - n_after, ").")
log_rows[[length(log_rows) + 1]] <- list(
  filter = "2_sample_size", n_before = n_before, n_after = n_after,
  n_dropped = n_before - n_after
)

# -----------------------------------------------------------------------------
# 6.4 — Filter 3: CNV exclusion (keep CNV == "N")
# -----------------------------------------------------------------------------
n_before <- nrow(ct)
ct <- ct[!is.na(ct$CNV) & ct$CNV == "N", ]
n_after <- nrow(ct)

message("(3) CNV == \"N\": kept ", n_after, " / ", n_before,
        " (dropped ", n_before - n_after, ").")
log_rows[[length(log_rows) + 1]] <- list(
  filter = "3_cnv", n_before = n_before, n_after = n_after,
  n_dropped = n_before - n_after
)

# -----------------------------------------------------------------------------
# 6.5 — Filter 4: GxE / multi-trait exclusion
# Two-tier system:
#   Tier 1 (auto_exclude): pleiotropy / interaction / variate / multitrait /
#                          ratio / conditional / cross_trait / MR /
#                          drug_response patterns
#   Tier 2 (manual_check): comma in MAPPED_TRAIT but no Tier 1 hit
#                          (i.e. study lists multiple mapped traits)
# -----------------------------------------------------------------------------
patterns <- default_gxe_patterns()
extra_patterns <- cfg$gxe_exclusions$extra_patterns
if (!is.null(extra_patterns)) {
  for (i in seq_along(extra_patterns)) {
    patterns[[paste0("custom_", i)]] <- extra_patterns[[i]]
  }
}

na0 <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
disease_trait <- na0(ct$`DISEASE/TRAIT`)
mapped_trait  <- na0(ct$MAPPED_TRAIT)
study         <- na0(ct$STUDY)

# Build a matrix of which patterns triggered, on any of the three text fields
flag_matrix <- vapply(patterns, function(p) {
  grepl(p, disease_trait, ignore.case = TRUE, perl = TRUE) |
  grepl(p, mapped_trait,  ignore.case = TRUE, perl = TRUE) |
  grepl(p, study,         ignore.case = TRUE, perl = TRUE)
}, logical(nrow(ct)))

flag_auto <- if (is.matrix(flag_matrix)) {
  apply(flag_matrix, 1, any)
} else {
  any(flag_matrix)   # single-row guard
}

# Multi-trait: a comma separating genuinely distinct phenotypes. A comma that
# is internal to one trait's name (a qualifier, or a nomenclature list) is NOT
# multi-trait. See is_composite_multitrait() in R/utils.R.
flag_multi_mapped <- is_composite_multitrait(mapped_trait)

# Anthropometric ratios are tracked separately for REPORTING, not as a tier.
# The keep-list below is the decision: an anthropometric ratio on it is kept,
# one not on it is auto-excluded like any other composite ratio. Routing them
# to manual_check instead would mean an unreviewed row gets dropped by the
# "no reviewed decision" path — an exclusion by omission rather than by rule,
# which is the failure mode this whole pass exists to remove. They are named
# individually in the log so the exclusion is never invisible.
manual_patterns <- default_manual_review_patterns()
manual_matrix <- vapply(manual_patterns, function(p) {
  grepl(p, disease_trait, ignore.case = TRUE, perl = TRUE) |
  grepl(p, mapped_trait,  ignore.case = TRUE, perl = TRUE)
}, logical(nrow(ct)))
flag_anthro_ratio <- if (is.matrix(manual_matrix)) {
  apply(manual_matrix, 1, any)
} else {
  any(manual_matrix)
}

# Multi-trait associations are excluded by policy -- the same stated rule the
# ratio filter implements, applied to the other composite shape. They are no
# longer routed to manual review: an exclusion the methods already declare
# should not depend on a review file being present.
flag_manual <- rep(FALSE, nrow(ct))

# Config-declared exceptions to the composite-ratio rule. These take
# precedence over every exclusion flag: the analyst has already ruled on them,
# so they are neither dropped nor re-surfaced for review. Exact match
# (case-insensitive, whitespace-trimmed) on MAPPED_TRAIT or DISEASE/TRAIT.
# Both keep-lists use normalize_trait_name() on BOTH sides — see
# .context/conventions.md §2. Matching on the raw string would miss any keep-list
# entry whose whitespace or capitalisation differs from the catalog's.
# Matching is EXACT on the normalised whole string, never substring. A
# substring keep-list let "body fat percentage" sweep in "Body fat percentage
# variance" -- a vQTL rather than a mean effect -- which is the same failure
# shape as everything else this pass removed: a mechanism doing more than its
# name says. Regional adiposity phenotypes are wanted, but they are LISTED
# individually rather than swept in. See .context/decisions.md D6.
apply_keep_list <- function(entries, label) {
  flag <- rep(FALSE, nrow(ct))
  if (is.null(entries) || length(entries) == 0) return(flag)
  kl  <- normalize_trait_name(unlist(entries))
  nmt <- normalize_trait_name(mapped_trait)
  ndt <- normalize_trait_name(disease_trait)
  flag <- nmt %in% kl | ndt %in% kl          # exact, whole-string
  matched <- unique(c(mapped_trait[nmt %in% kl], disease_trait[ndt %in% kl]))
  message("    ", label, ": ", length(kl), " declared exception(s), ",
          sum(flag), " row(s) retained",
          if (length(matched)) paste0(" (", paste(matched, collapse = "; "), ")") else "")
  unmatched <- setdiff(kl, c(nmt, ndt))
  if (length(unmatched) > 0) {
    message("    ", label, ": ", length(unmatched),
            " entry matched nothing in this phenotype: ",
            paste(unmatched, collapse = "; "))
  }
  flag
}

flag_ratio_keep <- apply_keep_list(cfg$gxe_exclusions$ratio_keep_list,
                                   "ratio_keep_list")
# Composite-multitrait exceptions. is_composite_multitrait() now fails CLOSED:
# any comma that is not mechanically explained makes a trait composite. That
# is correct for enumerations but over-reaches on APPOSITIONS, where the second
# fragment renames the first rather than adding a phenotype ("brain volume,
# neuroimaging measurement"). Apposition cannot be told from enumeration by
# regex — the distinction is semantic — so those traits are named here.
# Between the rule and this list, every comma-bearing trait has an explicit
# disposition and the exclusion rule is complete by inspection.
flag_multi_keep <- apply_keep_list(cfg$gxe_exclusions$multitrait_keep_list,
                                   "multitrait_keep_list")

flag_keep_listed <- flag_ratio_keep | flag_multi_keep

# -----------------------------------------------------------------------------
# Keep-lists are SCOPED TO THE RULE THAT CREATED THEM.
#
# Previously a keep-list took precedence over EVERY exclusion flag, so
# `ratio_keep_list` -- whose only job is to exempt an established composite
# index from the ratio/percentage rule -- was also rescuing genuine GxE and
# stratified analyses: 5 rows of "Waist-to-hip ratio adjusted for BMI in active
# individuals" from a physical-activity interaction meta-analysis, and 4
# age-stratified WHRadjBMI sub-analyses. Same class of failure as the
# reviewed-file fall-through: a mechanism doing more than its name says.
#
#   ratio_keep_list      exempts ONLY: ratio, pct_of_total, pct_of, conditional
#                        (conditional because a named index like WHRadjBMI or
#                         ABSI carries its adjustment in the phenotype
#                         definition -- see .context/decisions.md C5)
#   multitrait_keep_list exempts ONLY: the composite-multitrait rule
#
# Everything else -- pleiotropy, interaction, MR, drug response, the
# stratification family, the estimand guards -- is unexemptable.
# -----------------------------------------------------------------------------
RATIO_SCOPE <- c("ratio", "pct_of_total", "pct_of", "conditional")

in_scope <- function(nms) intersect(nms, names(patterns))
ratio_cols <- in_scope(RATIO_SCOPE)

flag_ratio_rules <- if (is.matrix(flag_matrix) && length(ratio_cols)) {
  apply(flag_matrix[, ratio_cols, drop = FALSE], 1, any)
} else rep(FALSE, nrow(ct))

# Exclusion flags that NO keep-list may override.
other_cols <- setdiff(colnames(flag_matrix), ratio_cols)
flag_unexemptable <- if (is.matrix(flag_matrix) && length(other_cols)) {
  apply(flag_matrix[, other_cols, drop = FALSE], 1, any)
} else rep(FALSE, nrow(ct))

# A row survives only if every rule that fired on it is one its keep-list
# actually covers.
rescued_ratio <- flag_ratio_keep & flag_ratio_rules & !flag_unexemptable
rescued_multi <- flag_multi_keep & flag_multi_mapped & !flag_unexemptable

overreach <- (flag_ratio_keep | flag_multi_keep) & flag_unexemptable
if (any(overreach)) {
  message("    keep-list scoping: ", sum(overreach),
          " keep-listed row(s) ALSO match an unexemptable rule and are ",
          "excluded (previously they were rescued):")
  for (s in sort(unique(disease_trait[overreach]))) {
    hits <- names(patterns)[flag_matrix[which(disease_trait == s & overreach)[1], ]]
    message("      - ", s, "  [", paste(hits, collapse = ","), "]")
  }
}

flag_keep_listed <- rescued_ratio | rescued_multi

ct$gxe_decision <- ifelse(flag_keep_listed, "keep",
                   ifelse(flag_auto | flag_multi_mapped, "auto_exclude",
                   ifelse(flag_manual, "manual_check", "keep")))

# Name every anthropometric ratio that is being excluded, so the decision is
# auditable from the run log alone rather than only from the config.
anthro_dropped <- flag_anthro_ratio & !flag_keep_listed
if (any(anthro_dropped)) {
  message("    anthropometric ratios NOT on the keep-list, auto-excluded (",
          sum(anthro_dropped), " rows):")
  for (s in sort(unique(mapped_trait[anthro_dropped]))) message("      - ", s)
  for (s in setdiff(sort(unique(disease_trait[anthro_dropped])), "")) {
    message("      - (DISEASE/TRAIT) ", s)
  }
}

# Reasons (for the review file)
ct$gxe_reason <- vapply(seq_len(nrow(ct)), function(i) {
  if (flag_keep_listed[i]) return("ratio_keep_list")
  hits <- if (is.matrix(flag_matrix)) names(patterns)[flag_matrix[i, ]] else character(0)
  if (is.matrix(manual_matrix)) {
    hits <- c(hits, names(manual_patterns)[manual_matrix[i, ]])
  }
  if (flag_multi_mapped[i]) hits <- c(hits, "multiple_MAPPED_TRAIT")
  if (length(hits) == 0) "" else paste(hits, collapse = ";")
}, character(1))

message("(4) GxE / multi-trait classification:")
message("       auto_exclude : ", sum(ct$gxe_decision == "auto_exclude"))
message("       manual_check : ", sum(ct$gxe_decision == "manual_check"))
message("         of auto_exclude, multi-trait: ", sum(flag_multi_mapped & !flag_keep_listed))
message("       keep         : ", sum(ct$gxe_decision == "keep"))

# -----------------------------------------------------------------------------
# 6.6 — Manual review handling
# -----------------------------------------------------------------------------
review_label <- "06_gxe"
reviewed_path <- file.path(cfg$paths$review, paste0(review_label, "_reviewed.tsv"))

if (file.exists(reviewed_path)) {
  # Load reviewed decisions and apply them, joining on a stable row identity
  message("Reviewed GxE file found; applying analyst decisions.")
  reviewed <- load_review(review_label, cfg, required = TRUE)

  # Use SNP + STUDY ACCESSION + Block as join key (unique enough in practice)
  key_cols <- c("SNP", "STUDY ACCESSION", "Block")
  missing_keys <- setdiff(key_cols, names(reviewed))
  if (length(missing_keys) > 0) {
    stop("Reviewed file missing key columns: ",
         paste(missing_keys, collapse = ", "))
  }
  if (!"gxe_decision" %in% names(reviewed)) {
    stop("Reviewed file missing 'gxe_decision' column.")
  }

  # Build lookup
  rev_key <- do.call(paste, c(reviewed[, key_cols], sep = "__"))
  ct_key  <- do.call(paste, c(ct[,       key_cols], sep = "__"))
  hits <- match(ct_key, rev_key)

  # Apply: only update rows that were originally manual_check AND appear in reviewed file
  apply_idx <- !is.na(hits) & ct$gxe_decision == "manual_check"
  ct$gxe_decision[apply_idx] <- reviewed$gxe_decision[hits[apply_idx]]
  message("Updated ", sum(apply_idx), " manual_check rows from reviewed file.")

  # Anything still labelled manual_check at this point gets dropped with a warning
  remaining_manual <- sum(ct$gxe_decision == "manual_check")
  if (remaining_manual > 0) {
    warning(remaining_manual, " manual_check rows had no reviewed decision; ",
            "they will be excluded from downstream analysis.")
  }

} else {
  # First run: write a review file and treat manual_check as keep-pending
  manual_rows <- ct[ct$gxe_decision == "manual_check",
                    c("Block", "SNP", "SNP_type", "STUDY ACCESSION",
                      "MAPPED_TRAIT", "DISEASE/TRAIT", "STUDY",
                      "PUBMEDID", "FIRST AUTHOR", "gxe_reason",
                      "gxe_decision")]
  if (nrow(manual_rows) > 0) {
    save_review(manual_rows, review_label, cfg)
    message("First-run review file written. After editing gxe_decision values ",
            "(keep / auto_exclude), save as ", reviewed_path, " and re-run.")
    message("For now, manual_check rows are treated as 'keep' (provisional).")
    ct$gxe_decision[ct$gxe_decision == "manual_check"] <- "keep"
  } else {
    message("No manual_check rows; skipping review file.")
  }
}

# -----------------------------------------------------------------------------
# 6.7 — Apply final GxE filter
# -----------------------------------------------------------------------------
n_before <- nrow(ct)
ct <- ct[ct$gxe_decision == "keep", ]
n_after <- nrow(ct)

message("(4) GxE / multi-trait: kept ", n_after, " / ", n_before,
        " (dropped ", n_before - n_after, ").")
log_rows[[length(log_rows) + 1]] <- list(
  filter = "4_gxe", n_before = n_before, n_after = n_after,
  n_dropped = n_before - n_after
)

# -----------------------------------------------------------------------------
# 6.8 — Save outputs
# -----------------------------------------------------------------------------
save_intermediate(ct, "06_crosstraits_filtered", cfg)

tsv_path <- file.path(cfg$paths$intermediate, "06_crosstraits_filtered.tsv")
write.table(ct, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("  saved: ", tsv_path)

log_df <- bind_rows(lapply(log_rows, as_tibble))
log_path <- file.path(cfg$paths$intermediate, "06_filtering_log.tsv")
write.table(log_df, log_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("  saved: ", log_path)

# -----------------------------------------------------------------------------
# 6.9 — Summary
# -----------------------------------------------------------------------------
bar <- paste(rep("\u2500", 60), collapse = "")
message("\n", bar)
message("Cross-trait filtering summary")
message(bar)
message(sprintf("%-25s  %10s  %10s  %10s",
                "filter", "in", "out", "dropped"))
for (r in log_rows[-1]) {
  message(sprintf("%-25s  %10d  %10d  %10d",
                  r$filter, r$n_before, r$n_after, r$n_dropped))
}
message(bar)
message(sprintf("Initial -> final : %d -> %d (%d dropped, %.1f%%)",
                initial_n, n_after, initial_n - n_after,
                100 * (initial_n - n_after) / initial_n))
message(bar)

write_stage_manifest(
  stage   = 6,
  cfg     = cfg,
  inputs  = file.path(cfg$paths$intermediate,
                      c("05_crosstraits_raw.rds", "04_blocks.rds")),
  outputs = c(file.path(cfg$paths$intermediate,
                        c("06_crosstraits_filtered.rds", "06_crosstraits_filtered.tsv",
                          "06_filtering_log.tsv"))),
  counts  = list(
    n_input          = initial_n,
    n_after_ancestry = log_rows[[2]]$n_after,
    n_after_samplesize = log_rows[[3]]$n_after,
    n_after_cnv      = log_rows[[4]]$n_after,
    n_after_gxe      = n_after,
    n_output_rows    = nrow(ct),
    n_output_cols    = ncol(ct),
    n_keep_listed    = sum(flag_keep_listed),
    n_anthro_ratio_excluded = sum(flag_anthro_ratio & !flag_keep_listed),
    n_multitrait_excluded   = sum(flag_multi_mapped & !flag_keep_listed)
  ))

message("\nStage 6 complete.")