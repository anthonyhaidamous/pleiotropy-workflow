# =============================================================================
# Stage 7 — Trait categorization
#
# Build a unique-trait frequency table from the filtered cross-traits, match
# against an optional reference table of categories and hierarchized labels,
# and flag unmatched traits as "new" for analyst review.
#
# Manual checkpoint
# -----------------
# First run: writes 07_trait_table_review.tsv with Category and
# `Hierarchized Trait` set to "new" for traits not found in the reference.
# Analyst fills these in, saves as 07_trait_table_reviewed.tsv, re-runs.
#
# If cfg$paths$trait_ref is not set or doesn't exist, every trait is "new"
# on first run — fine for analyzing a new phenotype from scratch.
#
# Self-trait categories (cfg$self_trait_categories) are removed AFTER the
# review step, so they don't pollute the enrichment computation.
#
# Input  : results/intermediate/06_crosstraits_filtered.rds
#          [optional] cfg$paths$trait_ref  — prior phenotype reference table
#
# Output : results/intermediate/07_trait_table.rds
#          results/intermediate/07_trait_table.tsv
#          results/review/07_trait_table_review.tsv  (first run only)
#
# Run from the project root:
#     Rscript pipeline/07_categorize.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(7, "Trait categorization")

# -----------------------------------------------------------------------------
# 7.1 — Load cross-traits and build the unique-trait frequency table
# -----------------------------------------------------------------------------
ct <- load_intermediate("06_crosstraits_filtered", cfg)

# Trait names are grouped and joined on a NORMALISED key, never on the raw
# string. 117 leptin trait names and 37 CLL trait names carry leading or
# trailing whitespace, 7 carry double internal spaces, and one pair differs
# only by case ("Hyperlipidemia" / "hyperlipidemia"). Grouping on the raw
# string counted that pair as two distinct traits and inflated the unique-trait
# total by one.
#
# The previous arrangement — reference and pipeline agreeing because NEITHER
# side was normalised — held only while both carried identical whitespace.
# Normalising either side alone would have silently lost 117 joins. Both sides
# now go through normalize_trait_name() (R/utils.R) at this single point.
# See .context/conventions.md §2.
ct$trait_key <- normalize_trait_name(ct$MAPPED_TRAIT)

n_raw_forms <- length(unique(ct$MAPPED_TRAIT))

trait_table <- ct %>%
  group_by(trait_key) %>%
  summarise(
    Frequency = n(),
    # NOTE: summarise() evaluates these expressions in order, and later ones
    # see the columns created by earlier ones. n_raw_variants MUST be computed
    # before MAPPED_TRAIT is overwritten, or it reads the new scalar and is
    # always 1 — which silently disables the collision report below.
    n_raw_variants = length(unique(MAPPED_TRAIT)),
    # Display label: the most frequent raw spelling, ties broken alphabetically
    # so the choice is deterministic across runs.
    MAPPED_TRAIT = names(sort(table(trimws(MAPPED_TRAIT)),
                              decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  arrange(desc(Frequency))

n_collisions <- n_raw_forms - nrow(trait_table)
if (n_collisions > 0) {
  message("Trait-name normalisation merged ", n_collisions,
          " raw spelling(s) into an existing trait:")
  for (k in trait_table$trait_key[trait_table$n_raw_variants > 1]) {
    variants <- unique(ct$MAPPED_TRAIT[ct$trait_key == k])
    message("    ", k, "  <-  ",
            paste(sprintf("%s", encodeString(variants, quote = "\"")),
                  collapse = " | "))
  }
} else {
  message("Trait-name normalisation: no collisions (",
          n_raw_forms, " raw spellings, ", nrow(trait_table), " traits).")
}

# Add a stable trait_id so the review-file join can't be broken by Excel
# mangling long MAPPED_TRAIT strings on save (truncation, quote handling,
# whitespace normalisation, etc.). The id is derived deterministically
# from the trait string, so re-running on the same data yields the same id.
trait_table$trait_id <- sprintf("T%04d", seq_len(nrow(trait_table)))

message("Built frequency table: ", nrow(trait_table), " unique traits, ",
        sum(trait_table$Frequency), " total associations.")

# -----------------------------------------------------------------------------
# 7.2 — Reference table join (optional)
# -----------------------------------------------------------------------------
ref_path <- cfg$paths$trait_ref
ref_available <- !is.null(ref_path) && nzchar(ref_path) && file.exists(ref_path)

if (ref_available) {
  message("Joining against reference table: ", ref_path)
  # Parses with proper quoting, strips CRLF/control characters, and ASSERTS
  # that nothing survives with a leading/trailing quote. See
  # read_trait_reference() in R/utils.R.
  ref <- read_trait_reference(ref_path)

  # Be flexible about reference column naming
  trait_col_candidates <- c("Trait Name in GWAS Catalog", "Trait", "MAPPED_TRAIT")
  ref_trait_col <- intersect(trait_col_candidates, names(ref))[1]
  if (is.na(ref_trait_col)) {
    stop("Reference trait file at ", ref_path,
         " has none of the expected columns: ",
         paste(trait_col_candidates, collapse = ", "))
  }

  # Both sides normalised, at the same point, with the same function.
  ref_key <- normalize_trait_name(ref[[ref_trait_col]])
  m <- match(trait_table$trait_key, ref_key)

  if (anyDuplicated(ref_key)) {
    dup <- unique(ref_key[duplicated(ref_key)])
    warning("Reference table has ", length(dup),
            " trait name(s) that are duplicates after normalisation; ",
            "the first occurrence wins: ", paste(head(dup, 10), collapse = "; "))
  }

  trait_table$Category <- ref$Category[m]
  trait_table$`Hierarchized Trait` <- ref$`Hierarchized Trait`[m]
  # Carry explicit exclusions through from the reference, if it declares them.
  trait_table$exclude <- if ("exclude" %in% names(ref)) ref$exclude[m] else ""
  trait_table$exclude_reason <-
    if ("exclude_reason" %in% names(ref)) ref$exclude_reason[m] else ""
  trait_table$exclude[is.na(trait_table$exclude)] <- ""
  trait_table$exclude_reason[is.na(trait_table$exclude_reason)] <- ""
} else {
  message("No reference trait file configured (cfg$paths$trait_ref). ",
          "All traits will be flagged 'new'.")
  trait_table$Category <- NA_character_
  trait_table$`Hierarchized Trait` <- NA_character_
}

trait_table$Category[is.na(trait_table$Category)]                       <- "new"
trait_table$`Hierarchized Trait`[is.na(trait_table$`Hierarchized Trait`)] <- "new"

n_new <- sum(trait_table$Category == "new")
message("Traits flagged 'new' (no reference match): ", n_new,
        " / ", nrow(trait_table))

# Reorder columns: trait_id first so it sits in column A in Excel
if (!"exclude" %in% names(trait_table)) trait_table$exclude <- ""
if (!"exclude_reason" %in% names(trait_table)) trait_table$exclude_reason <- ""
trait_table <- trait_table[, c("trait_id", "Category", "Hierarchized Trait",
                               "MAPPED_TRAIT", "Frequency",
                               "exclude", "exclude_reason")]

# -----------------------------------------------------------------------------
# 7.3 — Manual checkpoint: handle reviewed file or write review file
# -----------------------------------------------------------------------------
review_label <- "07_trait_table"
reviewed_path <- file.path(cfg$paths$review,
                           paste0(review_label, "_reviewed.tsv"))

if (file.exists(reviewed_path)) {
  message("Reviewed trait table found; applying analyst decisions.")
  reviewed <- load_review(review_label, cfg, required = TRUE)

  required_cols <- c("trait_id", "Category", "Hierarchized Trait")
  missing_cols <- setdiff(required_cols, names(reviewed))
  if (length(missing_cols) > 0) {
    stop("Reviewed file missing columns: ",
         paste(missing_cols, collapse = ", "),
         ".\nIf you reviewed an older review file (without trait_id), ",
         "re-run Stage 7 once to regenerate the review file with trait_id ",
         "and re-do the review against the new file.")
  }

  hits <- match(trait_table$trait_id, reviewed$trait_id)
  apply_idx <- !is.na(hits) & trait_table$Category == "new"
  trait_table$Category[apply_idx] <- reviewed$Category[hits[apply_idx]]
  trait_table$`Hierarchized Trait`[apply_idx] <-
    reviewed$`Hierarchized Trait`[hits[apply_idx]]
  message("Updated ", sum(apply_idx), " 'new' traits from reviewed file.")

  # Exclusions are read for EVERY row, not just the 'new' ones: marking an
  # already-categorised trait for exclusion is a legitimate later decision.
  if ("exclude" %in% names(reviewed)) {
    got <- !is.na(hits)
    trait_table$exclude[got] <- reviewed$exclude[hits[got]]
    if ("exclude_reason" %in% names(reviewed)) {
      trait_table$exclude_reason[got] <- reviewed$exclude_reason[hits[got]]
    }
    trait_table$exclude[is.na(trait_table$exclude)] <- ""
    trait_table$exclude_reason[is.na(trait_table$exclude_reason)] <- ""
  }

  # NOTE: absence from the reviewed file is NO LONGER an exclusion signal.
  # Deleting a row leaves no record of the decision, so it evaporates the
  # moment the file is superseded. Exclude by setting exclude = TRUE instead;
  # see section 7.3b. A 'new' trait missing from the reviewed file is now
  # reported as unfinished work rather than silently dropped.
  missing_idx <- trait_table$Category == "new" &
                 !trait_table$trait_id %in% reviewed$trait_id
  if (any(missing_idx)) {
    warning(sum(missing_idx), " 'new' trait(s) are absent from the reviewed ",
            "file. They are RETAINED as 'new' (stage 8 will refuse to run). ",
            "Add them with a Category, or with exclude=TRUE.")
    for (s in trait_table$MAPPED_TRAIT[missing_idx]) message("  ! unresolved: ", s)
  }

  remaining_new <- sum(trait_table$Category == "new")
  if (remaining_new > 0) {
    warning(remaining_new, " traits in the reviewed file are still ",
            "labelled 'new'. Fix their Category/Hierarchized Trait values.")
  }

} else if (n_new > 0) {
  save_review(trait_table, review_label, cfg)
  message("First-run review file written. After editing Category and ",
          "`Hierarchized Trait` for the 'new' rows, save as ",
          basename(reviewed_path), " and re-run.")
  message("For now, 'new' traits are passed through provisionally.")
} else {
  message("All traits matched the reference; no manual review needed.")
}

# -----------------------------------------------------------------------------
# 7.3b — Explicit exclusions
#
# CONVENTION: no analytical decision is ever encoded as the absence of a row.
#
# Traits are excluded by marking `exclude = TRUE`, never by deleting them from
# the reviewed file. Deleting a row leaves no record, so the decision evaporates
# the moment the file is superseded — which is exactly how four analyst
# exclusions from April 2026 (X-19438, X-23639, X-24295, YKL40) silently came
# back when a stale reviewed file was retired.
#
# The excluded set is written out as a first-class artifact so it can be
# tabulated in the supplement rather than inferred from what is missing.
# -----------------------------------------------------------------------------
truthy <- function(x) {
  x <- trimws(tolower(as.character(x)))
  !is.na(x) & x %in% c("true", "t", "yes", "y", "1", "x", "exclude")
}

if (!"exclude" %in% names(trait_table)) trait_table$exclude <- ""
if (!"exclude_reason" %in% names(trait_table)) trait_table$exclude_reason <- ""

is_excluded <- truthy(trait_table$exclude)
if (any(is_excluded)) {
  excluded <- trait_table[is_excluded,
                          c("trait_id", "Category", "Hierarchized Trait",
                            "MAPPED_TRAIT", "Frequency", "exclude_reason")]
  message("Explicit exclusions (exclude = TRUE): ", nrow(excluded),
          " trait(s), ", sum(excluded$Frequency), " association(s):")
  for (i in seq_len(nrow(excluded))) {
    message("  - ", excluded$MAPPED_TRAIT[i],
            if (nzchar(excluded$exclude_reason[i]))
              paste0("   [", excluded$exclude_reason[i], "]") else "")
  }
  excl_path <- file.path(cfg$paths$intermediate, "07_excluded_traits.tsv")
  write.table(excluded, excl_path, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
  message("  saved: ", excl_path)
  trait_table <- trait_table[!is_excluded, ]
} else {
  message("No traits carry exclude = TRUE.")
}

# -----------------------------------------------------------------------------
# 7.4 — Self-trait removal (Category in cfg$self_trait_categories)
# -----------------------------------------------------------------------------
self_cats <- cfg$self_trait_categories
if (!is.null(self_cats) && length(self_cats) > 0) {
  is_self <- trait_table$Category %in% self_cats
  n_self <- sum(is_self)
  if (n_self > 0) {
    message("Removing ", n_self,
            " self-phenotype traits (categories: ",
            paste(self_cats, collapse = ", "), ").")
    trait_table <- trait_table[!is_self, ]
  } else {
    message("No traits matched self-phenotype categories.")
  }
}

# -----------------------------------------------------------------------------
# 7.5 — Save
# -----------------------------------------------------------------------------
save_intermediate(trait_table, "07_trait_table", cfg)

tsv_path <- file.path(cfg$paths$intermediate, "07_trait_table.tsv")
write.table(trait_table, tsv_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
message("  saved: ", tsv_path)

# -----------------------------------------------------------------------------
# 7.6 — Summary
# -----------------------------------------------------------------------------
bar <- paste(rep("\u2500", 60), collapse = "")
message("\n", bar)
message("Trait categorization summary")
message(bar)
message("Unique traits                   : ", nrow(trait_table))
message("Total associations covered      : ", sum(trait_table$Frequency))
message("Distinct categories             : ",
        length(unique(trait_table$Category)))
message("Distinct hierarchized labels    : ",
        length(unique(trait_table$`Hierarchized Trait`)))
message("Traits still flagged 'new'      : ",
        sum(trait_table$Category == "new"))
message(bar)

cat_summary <- trait_table %>%
  group_by(Category) %>%
  summarise(traits = n(), associations = sum(Frequency),
            .groups = "drop") %>%
  arrange(desc(associations))

message("\nCategory breakdown:")
for (i in seq_len(nrow(cat_summary))) {
  message(sprintf("  %-30s  %4d traits  %5d assocs",
                  cat_summary$Category[i],
                  cat_summary$traits[i],
                  cat_summary$associations[i]))
}

write_stage_manifest(
  stage   = 7,
  cfg     = cfg,
  inputs  = c(file.path(cfg$paths$intermediate, "06_crosstraits_filtered.rds"),
              if (ref_available) ref_path else NULL),
  outputs = file.path(cfg$paths$intermediate,
                      c("07_trait_table.rds", "07_trait_table.tsv")),
  counts  = list(
    n_traits          = nrow(trait_table),
    n_cols            = ncol(trait_table),
    n_associations    = sum(trait_table$Frequency),
    n_categories      = length(unique(trait_table$Category)),
    n_hierarchized    = length(unique(trait_table$`Hierarchized Trait`)),
    n_still_new       = sum(trait_table$Category == "new"),
    reference_used    = if (ref_available) basename(ref_path) else "(none)"
  ))

message("\nStage 7 complete.")