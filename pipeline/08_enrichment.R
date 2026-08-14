# =============================================================================
# Stage 8 — Enrichment analysis
#
# For each trait and each category, compare the frequency observed within the
# phenotype's pleiotropy set against the frequency observed across the entire
# GWAS Catalog. Compute fold enrichment, log2-normalized enrichment, Fisher's
# exact p-value, and BH-FDR adjusted p-value.
#
# What this stage does NOT do:
#   - Filter or reshape the data — it just adds enrichment columns
#   - Plot — Stage 10 handles figures
#
# Catalog-frequency counting note
# -------------------------------
# The catalog frequency for a trait is the number of GWAS Catalog rows whose
# MAPPED_TRAIT field contains the trait name as a substring. This matches the
# original pipeline's behavior. It will count compound rows like
#   "X, Y" -> contributes 1 to both X and Y's catalog frequency
# which is consistent with how the cross-trait counts were built upstream.
#
# Input  : results/intermediate/07_trait_table.rds         — categorized traits
#          results/intermediate/01_gwassos_sigP.rds        — full GWAS Catalog
#                                                            (post p-value filter)
#
# Output : results/intermediate/08_trait_enrichment.rds    — per-trait stats
#          results/intermediate/08_trait_enrichment.tsv
#          results/intermediate/08_category_enrichment.rds — per-category stats
#          results/intermediate/08_category_enrichment.tsv
#
# Run from the project root:
#     Rscript pipeline/08_enrichment.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(8, "Enrichment analysis")

# -----------------------------------------------------------------------------
# 8.1 — Load
# -----------------------------------------------------------------------------
trait_table  <- load_intermediate("07_trait_table",  cfg)
gwassos_sigP <- load_intermediate("01_gwassos_sigP", cfg)

message("Trait table : ", nrow(trait_table), " rows.")
message("GWAS Catalog: ", nrow(gwassos_sigP), " p-value-filtered associations.")

# -----------------------------------------------------------------------------
# 8.1b — HARD GUARD: refuse to enrich uncategorised traits
#
# A placeholder is not a category. Left unchecked, stage 8 will happily compute
# a fold enrichment and a Fisher p-value for a category literally named "new"
# and return a plausible-looking, entirely meaningless result — which is more
# dangerous than a crash, because it survives review. Aborting is the point.
# -----------------------------------------------------------------------------
placeholders <- c("new", "", "NA", "na", "TBD", "tbd", "?", "-", "unknown",
                  "Unknown", "TODO", "todo", "check", "Check")
bad <- is.na(trait_table$Category) |
       trimws(trait_table$Category) %in% placeholders
if (any(bad)) {
  offenders <- sort(unique(trait_table$Category[bad]))
  stop("Stage 8 refuses to run: ", sum(bad), " of ", nrow(trait_table),
       " traits carry a placeholder category (",
       paste(sprintf("'%s'", offenders), collapse = ", "), ").\n",
       "A placeholder is not a category — enriching one produces a ",
       "plausible-looking but meaningless result.\n",
       "Complete ", file.path(cfg$paths$review, "07_trait_table_review.tsv"),
       " (set Category, or mark exclude=TRUE), save it as ",
       "07_trait_table_reviewed.tsv, and re-run stage 7.\n",
       "First few offending traits:\n  - ",
       paste(utils::head(trait_table$MAPPED_TRAIT[bad], 5), collapse = "\n  - "))
}
message("Category guard: all ", nrow(trait_table),
        " traits carry a real category.")

# -----------------------------------------------------------------------------
# 8.2 — Catalog frequency per trait (substring match on MAPPED_TRAIT)
# -----------------------------------------------------------------------------
clean_text <- function(x) {
  x <- gsub("[^[:alnum:],: ]+", "", x)
  x <- trimws(x)
  tolower(x)
}

mapped_clean <- clean_text(gwassos_sigP$MAPPED_TRAIT)

message("Computing catalog frequency for ", nrow(trait_table),
        " traits (substring match against ", length(mapped_clean),
        " catalog rows) ...")

trait_table$`Frequency in Catalog` <- vapply(
  seq_len(nrow(trait_table)),
  function(i) {
    target <- clean_text(trait_table$MAPPED_TRAIT[i])
    sum(grepl(target, mapped_clean, fixed = TRUE))
  },
  integer(1)
)

# -----------------------------------------------------------------------------
# 8.3 — Per-trait enrichment
# -----------------------------------------------------------------------------
total_pheno   <- sum(trait_table$Frequency)
total_catalog <- sum(trait_table$`Frequency in Catalog`)

# Drop traits that don't appear in the catalog at all (would force /0 enrichment)
zero_in_cat <- sum(trait_table$`Frequency in Catalog` == 0)
if (zero_in_cat > 0) {
  message("Dropping ", zero_in_cat,
          " traits with 0 catalog hits (cannot enrich).")
  trait_table <- trait_table[trait_table$`Frequency in Catalog` > 0, ]
  total_pheno   <- sum(trait_table$Frequency)
  total_catalog <- sum(trait_table$`Frequency in Catalog`)
}

trait_table$`Fold Enrichment` <-
  (trait_table$Frequency / trait_table$`Frequency in Catalog`) /
  (total_pheno / total_catalog)

trait_table$`Normalized Enrichment` <- log2(trait_table$`Fold Enrichment`)

# -----------------------------------------------------------------------------
# 8.4 — Per-category aggregation + enrichment
# -----------------------------------------------------------------------------
catagg <- trait_table %>%
  group_by(Category) %>%
  summarise(
    n_traits             = n(),
    Frequency            = sum(Frequency),
    `Frequency in Catalog` = sum(`Frequency in Catalog`),
    .groups = "drop"
  )

cat_total_pheno   <- sum(catagg$Frequency)
cat_total_catalog <- sum(catagg$`Frequency in Catalog`)

catagg$`Fold Enrichment` <-
  (catagg$Frequency / catagg$`Frequency in Catalog`) /
  (cat_total_pheno / cat_total_catalog)

catagg$`Normalized Enrichment` <- log2(catagg$`Fold Enrichment`)

# -----------------------------------------------------------------------------
# 8.5 — Per-trait Fisher's exact test
# -----------------------------------------------------------------------------
fisher_p <- function(in_pheno, in_catalog, total_pheno, total_catalog) {
  contingency <- matrix(c(
    in_pheno,                                          # this trait, in pheno
    in_catalog - in_pheno,                             # this trait, in catalog \ pheno
    total_pheno  - in_pheno,                           # other traits, in pheno
    total_catalog - in_catalog + in_pheno              # other traits, in catalog \ pheno
  ), nrow = 2)
  tryCatch(fisher.test(contingency)$p.value, error = function(e) NA_real_)
}

trait_table$P_Value <- vapply(
  seq_len(nrow(trait_table)),
  function(i) fisher_p(trait_table$Frequency[i],
                       trait_table$`Frequency in Catalog`[i],
                       total_pheno, total_catalog),
  numeric(1)
)

trait_table$P_Adjusted <- p.adjust(trait_table$P_Value, method = "BH")

# -----------------------------------------------------------------------------
# 8.6 — Per-category Fisher's exact test
# -----------------------------------------------------------------------------
catagg$P_Value <- vapply(
  seq_len(nrow(catagg)),
  function(i) fisher_p(catagg$Frequency[i],
                       catagg$`Frequency in Catalog`[i],
                       cat_total_pheno, cat_total_catalog),
  numeric(1)
)

catagg$P_Adjusted <- p.adjust(catagg$P_Value, method = "BH")

# Star annotation (used by Stage 10's lollipop plot)
star <- function(p) {
  ifelse(is.na(p), "-",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",  "-"))))
}
trait_table$Stars <- star(trait_table$P_Adjusted)
catagg$Stars      <- star(catagg$P_Adjusted)

# -----------------------------------------------------------------------------
# 8.7 — Order results
# -----------------------------------------------------------------------------
trait_table <- trait_table %>%
  arrange(desc(`Fold Enrichment`), Category)

catagg <- catagg %>%
  arrange(desc(`Fold Enrichment`))

# -----------------------------------------------------------------------------
# 8.8 — Save
# -----------------------------------------------------------------------------
save_intermediate(trait_table, "08_trait_enrichment",    cfg)
save_intermediate(catagg,      "08_category_enrichment", cfg)

trait_path <- file.path(cfg$paths$intermediate, "08_trait_enrichment.tsv")
cat_path   <- file.path(cfg$paths$intermediate, "08_category_enrichment.tsv")
write.table(trait_table, trait_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(catagg, cat_path,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("  saved: ", trait_path)
message("  saved: ", cat_path)

# -----------------------------------------------------------------------------
# 8.9 — Summary
# -----------------------------------------------------------------------------
bar <- paste(rep("\u2500", 70), collapse = "")
message("\n", bar)
message("Category enrichment summary (sorted by fold enrichment)")
message(bar)
message(sprintf("%-30s %6s %8s %10s %5s",
                "Category", "n_tr", "fold", "p_adj", ""))
for (i in seq_len(nrow(catagg))) {
  message(sprintf("  %-30s %6d %8.2f %10s %5s",
                  catagg$Category[i],
                  catagg$n_traits[i],
                  catagg$`Fold Enrichment`[i],
                  format_pval(catagg$P_Adjusted[i]),
                  catagg$Stars[i]))
}
message(bar)
message("Significance: *** padj<0.001  ** padj<0.01  * padj<0.05  - n.s.")
message("(BH-FDR-corrected Fisher's exact test)")

write_stage_manifest(
  stage   = 8,
  cfg     = cfg,
  inputs  = file.path(cfg$paths$intermediate,
                      c("07_trait_table.rds", "01_gwassos_sigP.rds")),
  outputs = file.path(cfg$paths$intermediate,
                      c("08_trait_enrichment.rds", "08_trait_enrichment.tsv",
                        "08_category_enrichment.rds", "08_category_enrichment.tsv")),
  counts  = list(
    n_traits            = nrow(trait_table),
    n_categories        = nrow(catagg),
    n_traits_dropped_zero_catalog = zero_in_cat,
    total_pheno         = total_pheno,
    total_catalog       = total_catalog,
    n_sig_traits_padj05 = sum(trait_table$P_Adjusted < 0.05, na.rm = TRUE),
    n_sig_cats_padj05   = sum(catagg$P_Adjusted < 0.05, na.rm = TRUE)
  ))

message("\nStage 8 complete.")