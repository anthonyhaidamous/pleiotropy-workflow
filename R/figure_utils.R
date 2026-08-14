# =============================================================================
# R/figure_utils.R
# Helpers shared by the CLI-driven scripts (figure_scripts/*, analysis/*,
# pipeline/10_extract_genes.R) and by the effect-harmonization code in
# pipeline/09_harmonize.R.
#
#     source("R/figure_utils.R")     # from the project root
#
# Scope and rules
# ---------------
#   - Base R only. This file loads NO packages and has NO side effects, so it
#     is safe to source from any script regardless of what that script has
#     already attached. Do not add library() calls here.
#   - These functions were previously copy-pasted across 9 files (read_any x8,
#     safe_dir_create x3, extract_ci_bounds x3, derive_signed_beta x2). The
#     copies had drifted apart; see the notes on each function below for which
#     variant is canonical and why.
# =============================================================================


# -----------------------------------------------------------------------------
# Read a stage artifact, dispatching on file extension
#
# Stage handoffs exist as both .rds (authoritative, preserves list-columns and
# types) and .tsv (human-readable mirror). Scripts accept either.
#
# `as_data_frame` reproduces a genuine difference between the previous copies:
# five call sites coerced the RDS to a base data.frame, four used it as-is.
# That is not cosmetic — the stage RDS files are tibbles, and coercion changes
# `[` semantics (tibble always returns a tibble; data.frame drops to a vector
# for a single column). Each call site therefore passes the value that matches
# the copy it used to carry.
# -----------------------------------------------------------------------------
read_any <- function(path, as_data_frame = FALSE) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (as_data_frame) {
      as.data.frame(obj, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      obj
    }
  } else {
    read.delim(path, sep = "\t", quote = "", header = TRUE,
               check.names = FALSE, stringsAsFactors = FALSE)
  }
}


# -----------------------------------------------------------------------------
# Create a directory if it isn't already there
# All three previous copies were byte-identical.
# -----------------------------------------------------------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


# -----------------------------------------------------------------------------
# Parse the two numeric bounds out of a GWAS Catalog "95% CI (TEXT)" string
#
# The catalog writes intervals as "[lower-upper]", optionally followed by a
# unit phrase: "[0.0039-0.0051] unit increase", "[1.45-2.12]".
#
# The separator is a HYPHEN, which is also the minus sign. A global numeric
# match such as gregexpr("-?\\d+\\.?\\d*") therefore reads the separator as
# part of the second number and returns (0.0039, -0.0051) — silently negating
# the upper bound. Downstream, the allele flip and the lower<=upper "safety"
# swap in stage 9 turn that into an interval STRADDLING ZERO: 1,429 of 3,040
# leptin rows and 176 of 445 CLL rows, i.e. intervals that read as
# non-significant in a forest plot when the true interval excludes zero.
# See AUDIT.md 3.1(A).
#
# Effect estimates were never affected: beta comes from the catalog's
# "OR or BETA" column, and the CI text is used only via
# grepl("increase|decrease") for sign detection, which the defect cannot
# reach. Only the interval bounds were wrong.
#
# The fix matches the interval STRUCTURALLY — two numbers on either side of an
# explicit separator token, captured as regex groups — instead of scavenging
# numbers out of a global match. A minus sign is consumed only where it can
# belong to a number rather than to the separator.
#
# Returns c(lower, upper), or c(NA, NA) when no interval is present.
# -----------------------------------------------------------------------------
extract_ci_bounds <- function(ci) {
  if (is.na(ci) || ci == "") return(c(NA_real_, NA_real_))

  # A number: optional sign, optional integer part, optional decimal point,
  # required digits, optional scientific notation ("1.2e-4").
  num_pat <- "(-?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?)"

  # Two numbers separated by hyphen, en-dash, "to", comma or semicolon.
  pattern <- paste0(num_pat, "\\s*(?:-|to|,|;|\u2013)\\s*", num_pat)

  m <- regmatches(ci, regexec(pattern, ci, ignore.case = TRUE))[[1]]

  # m[1] is the whole match; m[2] and m[3] are the two captured numbers.
  if (length(m) >= 3) {
    c(as.numeric(m[2]), as.numeric(m[3]))
  } else {
    c(NA_real_, NA_real_)
  }
}


# -----------------------------------------------------------------------------
# Derive a signed beta from the catalog's "OR or BETA" value plus its CI text
#
# The catalog stores odds ratios and betas in one column and does not say which
# is which. The 95% CI text is the discriminator: a phrase containing
# "increase" or "decrease" means the value is already on the beta scale;
# anything else means it is an odds ratio and needs log-transforming. ORs are
# always positive, so "decrease" is the only signal that the effect is negative.
#
# suppressWarnings() on the coercion is deliberate: a non-numeric value must
# return NA quietly rather than emitting a warning per row.
# -----------------------------------------------------------------------------
derive_signed_beta <- function(or_or_beta, ci_text) {
  if (is.na(or_or_beta)) return(NA_real_)
  if (is.na(ci_text)) ci_text <- ""
  is_beta_scale <- grepl("increase|decrease", ci_text, ignore.case = TRUE)
  val <- suppressWarnings(as.numeric(or_or_beta))
  if (is.na(val)) return(NA_real_)
  beta <- if (is_beta_scale) val else log(val)
  if (grepl("decrease", ci_text, ignore.case = TRUE)) beta <- -abs(beta)
  beta
}
