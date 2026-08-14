# =============================================================================
# Stage 1 — Load and prepare lead SNPs (lSNPs)
#
# Input  : GWAS catalog studies, associations, and ancestry files
# Output : results/intermediate/01_lSNPs.rds         — lead SNPs with ancestry
#          results/intermediate/01_gwassos_sigP.rds  — full filtered associations
#          results/review/01_lSNPs_review.tsv        — ancestry encoding failures
#          results/review/01_pvalue_review.tsv       — lSNPs at the p-value
#                                                       threshold (likely
#                                                       rounded); analyst
#                                                       can drop or correct
#
# Run from the project root:
#     Rscript pipeline/01_load_snps.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(1, paste0("Lead SNPs for phenotype: ", cfg$phenotype$name))

# -----------------------------------------------------------------------------
# 1.1 — Studies: identify the studies that match our trait terms
# -----------------------------------------------------------------------------
studies <- read.delim(cfg$paths$studies, sep = "\t", quote = "",
                      check.names = FALSE, stringsAsFactors = FALSE)
studies <- subset(studies, MAPPED_TRAIT %in% cfg$phenotype$trait_terms)

message("Found ", nrow(studies), " studies covering ",
        sum(studies$`ASSOCIATION COUNT`, na.rm = TRUE),
        " associations for phenotype '", cfg$phenotype$name, "'.")

# -----------------------------------------------------------------------------
# 1.2 — Associations: filter by p-value, filter by trait, dedupe by lowest p
# -----------------------------------------------------------------------------
gwassos_sigP <- read.delim(cfg$paths$associations,
                           sep = "\t", quote = "",
                           check.names = FALSE, stringsAsFactors = FALSE)
gwassos_sigP$`P-VALUE` <- as.numeric(gwassos_sigP$`P-VALUE`)
gwassos_sigP <- gwassos_sigP[!is.na(gwassos_sigP$`P-VALUE`) &
                               gwassos_sigP$`P-VALUE` <= cfg$filters$pvalue, ]

lSNPs <- gwassos_sigP[gwassos_sigP$MAPPED_TRAIT %in% cfg$phenotype$trait_terms, ]
lSNPs <- lSNPs[!is.na(lSNPs$SNPS) & lSNPs$SNPS != "", ]

# Resolve any chr:pos coordinates to rsIDs
lSNPs$SNPS <- resolve_chrpos_to_rsid(lSNPs$SNPS, build = cfg$ensembl$build)

# Sort by SNP, then p-value (smallest first), and keep one row per SNP
lSNPs <- lSNPs[order(lSNPs$SNPS, lSNPs$`P-VALUE`, na.last = TRUE), ]
lSNPs <- lSNPs[!duplicated(lSNPs$SNPS), ]

message("Retained ", nrow(lSNPs),
        " unique lead SNPs after p-value filter (P \u2264 ",
        format(cfg$filters$pvalue, scientific = TRUE), ").")

# -----------------------------------------------------------------------------
# 1.3 — Ancestry: cohort-proportion rule on INITIAL + REPLICATION sample sizes
#
# In each field, parse comma-separated cohorts -> superpop counts. If one
# superpop >= 90% of the field total, emit just that one; otherwise emit all.
# Union the discovery and replication results.
# -----------------------------------------------------------------------------
lSNPs$`Encoded Populations` <- encode_ancestries_for_rows(lSNPs)

# -----------------------------------------------------------------------------
# 1.4 — Manual checkpoint: rows whose ancestry could not be determined
# -----------------------------------------------------------------------------
needs_review <- lSNPs[lSNPs$`Encoded Populations` == "Check study",
                      c("SNPS", "STUDY ACCESSION", "INITIAL SAMPLE SIZE",
                        "REPLICATION SAMPLE SIZE",
                        "PUBMEDID", "FIRST AUTHOR", "STUDY")]

if (nrow(needs_review) > 0) {
  warning(nrow(needs_review),
          " lead SNPs could not be ancestry-encoded automatically. ",
          "See review file.")
  save_review(needs_review, "01_lSNPs", cfg)
} else {
  message("All lead SNPs successfully ancestry-encoded.")
}

# -----------------------------------------------------------------------------
# 1.5 — Manual checkpoint: p-values exactly at the significance threshold
#
# The GWAS Catalog stores p-values at varying precision; some studies report
# a value like p = 5e-8 that has actually been rounded down from a sub-
# threshold value (e.g. true p = 5.4e-8 reported as 5e-8). Such rows pass
# our P <= cfg$filters$pvalue filter when they shouldn't. This checkpoint
# flags every lSNP whose p-value equals the filter threshold exactly so the
# analyst can verify against the original publication.
#
# Review file columns the analyst fills in:
#   pvalue_decision  : 'keep' or 'exclude'
#   pvalue_corrected : numeric, optional. If provided AND decision is 'keep',
#                      the corrected value is written back into both lSNPs
#                      and gwassos_sigP (matched on SNPS + STUDY ACCESSION).
#
# First run  : writes 01_pvalue_review.tsv (kept rows pass through provisionally).
# Re-run     : reads 01_pvalue_reviewed.tsv, applies decisions and corrections.
# -----------------------------------------------------------------------------
pval_review_label  <- "01_pvalue"
pval_reviewed_path <- file.path(cfg$paths$review,
                                paste0(pval_review_label, "_reviewed.tsv"))

threshold    <- cfg$filters$pvalue
pval_tol     <- threshold * 1e-6
at_threshold <- !is.na(lSNPs$`P-VALUE`) &
                abs(lSNPs$`P-VALUE` - threshold) < pval_tol

if (file.exists(pval_reviewed_path)) {
  message("Reviewed p-value file found; applying analyst decisions.")
  reviewed <- load_review(pval_review_label, cfg, required = TRUE)

  required_cols <- c("SNPS", "STUDY ACCESSION",
                     "pvalue_decision", "pvalue_corrected")
  missing_cols <- setdiff(required_cols, names(reviewed))
  if (length(missing_cols) > 0) {
    stop("Reviewed p-value file missing column(s): ",
         paste(missing_cols, collapse = ", "),
         ". Expected columns: SNPS, STUDY ACCESSION, ",
         "pvalue_decision (keep|exclude), pvalue_corrected (numeric or blank).")
  }

  # Normalise: trim whitespace, drop empty SNPs rows
  reviewed$pvalue_decision  <- trimws(reviewed$pvalue_decision)
  reviewed$pvalue_corrected <- trimws(as.character(reviewed$pvalue_corrected))
  reviewed <- reviewed[reviewed$SNPS != "" & !is.na(reviewed$SNPS), ]

  # Validate decisions
  empty_decisions <- reviewed$pvalue_decision == "" | is.na(reviewed$pvalue_decision)
  if (any(empty_decisions)) {
    stop(sum(empty_decisions),
         " row(s) in p-value review have an empty 'pvalue_decision'. ",
         "Set each to 'keep' or 'exclude' before re-running. ",
         "Offending SNPs: ",
         paste(reviewed$SNPS[empty_decisions], collapse = ", "))
  }
  bad_decisions <- !reviewed$pvalue_decision %in% c("keep", "exclude")
  if (any(bad_decisions)) {
    stop("Unrecognised pvalue_decision value(s): ",
         paste(unique(reviewed$pvalue_decision[bad_decisions]),
               collapse = ", "),
         ". Expected 'keep' or 'exclude'.")
  }

  # Parse corrections (blank/NA stays NA; otherwise must be numeric)
  reviewed$pvalue_corrected_num <- suppressWarnings(
    as.numeric(reviewed$pvalue_corrected)
  )
  bad_numeric <- reviewed$pvalue_corrected != "" &
                 !is.na(reviewed$pvalue_corrected) &
                 is.na(reviewed$pvalue_corrected_num)
  if (any(bad_numeric)) {
    stop("Non-numeric pvalue_corrected value(s) in review: ",
         paste(reviewed$pvalue_corrected[bad_numeric], collapse = ", "))
  }

  # Catch contradictions: keep + corrected > threshold
  bad_contradictions <- reviewed$pvalue_decision == "keep" &
                        !is.na(reviewed$pvalue_corrected_num) &
                        reviewed$pvalue_corrected_num > threshold
  if (any(bad_contradictions)) {
    stop("Contradictory review row(s): pvalue_decision='keep' but ",
         "pvalue_corrected exceeds the significance threshold (",
         format(threshold, scientific = TRUE), "). ",
         "Either change decision to 'exclude' or revise the corrected value. ",
         "Offending SNPs: ",
         paste(reviewed$SNPS[bad_contradictions], collapse = ", "))
  }

  # ---- Apply 'exclude' decisions -------------------------------------------
  excludes <- reviewed[reviewed$pvalue_decision == "exclude",
                       c("SNPS", "STUDY ACCESSION")]
  if (nrow(excludes) > 0) {
    drop_lSNP <- paste(lSNPs$SNPS, lSNPs$`STUDY ACCESSION`) %in%
                 paste(excludes$SNPS, excludes$`STUDY ACCESSION`)
    drop_gwas <- paste(gwassos_sigP$SNPS, gwassos_sigP$`STUDY ACCESSION`) %in%
                 paste(excludes$SNPS, excludes$`STUDY ACCESSION`)

    message("Dropping ", sum(drop_lSNP),
            " lead SNP(s) marked 'exclude' (and ",
            sum(drop_gwas),
            " matching row(s) in gwassos_sigP):")
    for (s in lSNPs$SNPS[drop_lSNP]) message("  - ", s)

    lSNPs        <- lSNPs[!drop_lSNP, ]
    gwassos_sigP <- gwassos_sigP[!drop_gwas, ]
  }

  # ---- Apply 'keep' + correction overrides ---------------------------------
  corrections <- reviewed[reviewed$pvalue_decision == "keep" &
                          !is.na(reviewed$pvalue_corrected_num),
                          c("SNPS", "STUDY ACCESSION", "pvalue_corrected_num")]
  if (nrow(corrections) > 0) {
    n_lSNP <- 0L; n_gwas <- 0L
    for (i in seq_len(nrow(corrections))) {
      key_snp   <- corrections$SNPS[i]
      key_study <- corrections$`STUDY ACCESSION`[i]
      new_p     <- corrections$pvalue_corrected_num[i]

      hit_lSNP <- lSNPs$SNPS == key_snp &
                  lSNPs$`STUDY ACCESSION` == key_study
      lSNPs$`P-VALUE`[hit_lSNP] <- new_p
      n_lSNP <- n_lSNP + sum(hit_lSNP)

      hit_gwas <- gwassos_sigP$SNPS == key_snp &
                  gwassos_sigP$`STUDY ACCESSION` == key_study
      gwassos_sigP$`P-VALUE`[hit_gwas] <- new_p
      n_gwas <- n_gwas + sum(hit_gwas)
    }
    message("Applied ", nrow(corrections),
            " p-value correction(s): updated ", n_lSNP,
            " lSNP row(s) and ", n_gwas,
            " gwassos_sigP row(s).")
  }

  if (nrow(excludes) == 0 && nrow(corrections) == 0) {
    message("Review file present but no exclusions or corrections; ",
            "all flagged SNPs retained as-is.")
  }

} else if (any(at_threshold)) {
  pval_review <- lSNPs[at_threshold,
                       c("SNPS", "P-VALUE", "STUDY ACCESSION",
                         "PUBMEDID", "FIRST AUTHOR", "STUDY",
                         "MAPPED_GENE", "STRONGEST SNP-RISK ALLELE")]
  pval_review$pvalue_decision  <- ""   # 'keep' or 'exclude'
  pval_review$pvalue_corrected <- ""   # numeric, optional

  message(sum(at_threshold),
          " lead SNPs have a p-value exactly at the significance threshold (",
          format(threshold, scientific = TRUE),
          "). These may be rounded values from the catalog.")
  message("Edit the review file:")
  message("  - set 'pvalue_decision' to 'keep' or 'exclude' (required)")
  message("  - optionally fill in 'pvalue_corrected' with the true value ",
          "from the publication")
  message("Save as ", pval_review_label, "_reviewed.tsv and re-run.")
  message("For now, all flagged SNPs are kept (provisional).")

  save_review(pval_review, pval_review_label, cfg)

} else {
  message("No lead SNPs at the p-value threshold; ",
          "no p-value review needed.")
}

# -----------------------------------------------------------------------------
# 1.6 — Save
# -----------------------------------------------------------------------------
save_intermediate(lSNPs,         "01_lSNPs",         cfg)
save_intermediate(gwassos_sigP,  "01_gwassos_sigP",  cfg)  # reused downstream

message("\nStage 1 complete.")