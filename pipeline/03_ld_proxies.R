# =============================================================================
# Stage 3 — LD proxy retrieval (filtering)
#
# Splits the LDproxy results into two views needed downstream:
#   - leadSNPs_LDs : full LD set at R² > ld_block_r2, used by Stage 4 to
#                    define block boundaries (lead + proxies + correlated SNPs)
#   - pSNPs        : high-LD proxies at R² ≥ ld_proxy_r2, self-rows removed,
#                    aggregated to one row per lead SNP with a list-column
#                    of unique pSNPs (union across queried populations)
#
# Input  : results/intermediate/02_proxy_table.rds
#          results/intermediate/02_lSNPs.rds
#
# Output : results/intermediate/03_leadSNPs_LDs.rds
#          results/intermediate/03_pSNPs_per_lSNP.rds
#
# Run from the project root:
#     Rscript pipeline/03_ld_proxies.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(3, "LD proxy retrieval")

# -----------------------------------------------------------------------------
# 3.1 — Load
# -----------------------------------------------------------------------------
proxy_table <- load_intermediate("02_proxy_table", cfg)
lSNPs       <- load_intermediate("02_lSNPs",       cfg)

message("Loaded ", nrow(proxy_table), " LDproxy rows across ",
        n_distinct(proxy_table$query_snp), " lead SNPs.")

# -----------------------------------------------------------------------------
# 3.2 — leadSNPs_LDs : full LD set for block construction (R² > ld_block_r2)
# -----------------------------------------------------------------------------
leadSNPs_LDs <- proxy_table %>%
  filter(R2 > cfg$filters$ld_block_r2)

message("Retained ", nrow(leadSNPs_LDs),
        " rows at R\u00b2 > ", cfg$filters$ld_block_r2,
        " (block construction set).")

# -----------------------------------------------------------------------------
# 3.3 — pSNPs : high-LD proxies (R² ≥ ld_proxy_r2), self-rows excluded
# A self-row is one where the proxy IS the queried SNP (R²=1, Distance=0).
# Filter explicitly on RS_Number != query_snp rather than positional drops.
# -----------------------------------------------------------------------------
pSNPs_rows <- proxy_table %>%
  filter(R2 >= cfg$filters$ld_proxy_r2,
         RS_Number != query_snp)

# Aggregate to one row per lead SNP, dedup across populations
pSNPs_per_lSNP <- pSNPs_rows %>%
  group_by(query_snp) %>%
  summarise(pSNPs = list(unique(RS_Number)), .groups = "drop") %>%
  mutate(n_pSNPs = lengths(pSNPs))

# Make sure every surviving lead SNP appears, even if it has zero pSNPs
missing_lSNPs <- setdiff(lSNPs$SNPS, pSNPs_per_lSNP$query_snp)
if (length(missing_lSNPs) > 0) {
  pSNPs_per_lSNP <- bind_rows(
    pSNPs_per_lSNP,
    tibble::tibble(query_snp = missing_lSNPs,
                   pSNPs     = vector("list", length(missing_lSNPs)),
                   n_pSNPs   = 0L)
  )
}

n_with_proxies <- sum(pSNPs_per_lSNP$n_pSNPs > 0)
n_total_pSNPs  <- length(unique(unlist(pSNPs_per_lSNP$pSNPs)))

message("pSNPs at R\u00b2 \u2265 ", cfg$filters$ld_proxy_r2,
        " (self-rows dropped):")
message("  lead SNPs with \u22651 pSNP : ", n_with_proxies,
        " of ", nrow(pSNPs_per_lSNP))
message("  unique pSNPs total    : ", n_total_pSNPs)

# -----------------------------------------------------------------------------
# 3.4 — Save
# -----------------------------------------------------------------------------
save_intermediate(leadSNPs_LDs,   "03_leadSNPs_LDs",   cfg)
save_intermediate(pSNPs_per_lSNP, "03_pSNPs_per_lSNP", cfg)

message("\nStage 3 complete.")
