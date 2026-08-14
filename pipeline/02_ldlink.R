# =============================================================================
# Stage 2 — LDproxy retrieval via LDlinkR
#
# Input  : results/intermediate/01_lSNPs.rds
#          LDLINK_TOKEN environment variable
#          (request one at https://ldlink.nih.gov/?tab=apiaccess)
#
# Output : results/intermediate/02_proxy_table.rds   — combined LDproxy results
#          results/intermediate/02_proxy_table.tsv   — TSV mirror for archival
#          results/intermediate/02_lSNPs.rds         — lSNPs after dropping any
#                                                      SNP that failed in EVERY
#                                                      queried population
#          results/review/02_failures_review.tsv     — failed SNP-population pairs
#
# Run from the project root:
#     export LDLINK_TOKEN=<your_token>
#     Rscript pipeline/02_ldlink.R config/<trait>.yml
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages(library(LDlinkR))

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(2, "LDproxy retrieval via LDlinkR")

# -----------------------------------------------------------------------------
# 2.1 — Token
# -----------------------------------------------------------------------------
token <- Sys.getenv("LDLINK_TOKEN")
if (token == "") {
  stop("LDLINK_TOKEN environment variable is not set.\n",
       "Get one at https://ldlink.nih.gov/?tab=apiaccess and run:\n",
       "    export LDLINK_TOKEN=<your_token>\n",
       "or add it to ~/.Renviron")
}

genome_build <- tolower(cfg$ensembl$build)   # "grch38" / "grch37"

# -----------------------------------------------------------------------------
# 2.2 — Load Stage 1 lead SNPs
# -----------------------------------------------------------------------------
lSNPs <- load_intermediate("01_lSNPs", cfg)

# -----------------------------------------------------------------------------
# 2.3 — Build SNP–population map (one row per SNP-pop pair)
# -----------------------------------------------------------------------------
valid_pops <- c("EUR", "AFR", "EAS", "SAS", "AMR")

snp_pop_map <- lSNPs %>%
  select(SNPS, `Encoded Populations`) %>%
  distinct() %>%
  mutate(`Encoded Populations` = strsplit(`Encoded Populations`, ",\\s*")) %>%
  tidyr::unnest(`Encoded Populations`) %>%
  mutate(`Encoded Populations` = trimws(`Encoded Populations`)) %>%
  filter(`Encoded Populations` %in% valid_pops)

message("Querying ", nrow(snp_pop_map), " SNP-population pairs (",
        n_distinct(snp_pop_map$SNPS), " unique SNPs).")

# -----------------------------------------------------------------------------
# 2.4 — LDproxy loop
# Each iteration returns a populated data frame on success, or NULL on failure
# (with a row appended to failed_log via <<- in the parent scope).
# -----------------------------------------------------------------------------
failed_log <- list()

results_list <- lapply(seq_len(nrow(snp_pop_map)), function(i) {
  snp <- snp_pop_map$SNPS[i]
  pop <- snp_pop_map$`Encoded Populations`[i]

  message("  [", i, "/", nrow(snp_pop_map), "] ", snp, " — ", pop)

  res <- tryCatch(
    LDproxy(snp = snp, pop = pop, r2d = "r2",
            token = token, genome_build = genome_build),
    error = function(e) {
      failed_log[[paste(snp, pop, sep = "_")]] <<- data.frame(
        query_snp = snp, population = pop, reason = e$message
      )
      NULL
    }
  )
  if (is.null(res)) return(NULL)

  # LDproxy returns a 1-row df with an error string instead of throwing
  if (ncol(res) == 1 ||
      (nrow(res) == 1 && grepl("error", res[1, 1], ignore.case = TRUE))) {
    failed_log[[paste(snp, pop, sep = "_")]] <<- data.frame(
      query_snp = snp, population = pop, reason = as.character(res[1, 1])
    )
    return(NULL)
  }

  res$query_snp  <- snp
  res$population <- pop
  res
})

# -----------------------------------------------------------------------------
# 2.5 — Combine results
# -----------------------------------------------------------------------------
proxy_table <- bind_rows(results_list) %>%
  filter(!is.na(query_snp)) %>%
  select(query_snp, RS_Number, Coord, Alleles, MAF,
         Distance, Dprime, R2, Correlated_Alleles,
         FORGEdb, RegulomeDB, Function)

failed_table <- bind_rows(failed_log)

# -----------------------------------------------------------------------------
# 2.6 — Drop SNPs that failed in EVERY queried population
# -----------------------------------------------------------------------------
attempted_snps    <- unique(snp_pop_map$SNPS)
successful_snps   <- unique(proxy_table$query_snp)
fully_failed_snps <- setdiff(attempted_snps, successful_snps)

if (length(fully_failed_snps) > 0) {
  message("Dropping ", length(fully_failed_snps),
          " SNPs absent from LDproxy in every queried population:")
  for (s in fully_failed_snps) message("  - ", s)
  lSNPs <- lSNPs[!lSNPs$SNPS %in% fully_failed_snps, ]
}

# -----------------------------------------------------------------------------
# 2.7 — Save
# -----------------------------------------------------------------------------
save_intermediate(proxy_table, "02_proxy_table", cfg)
save_intermediate(lSNPs,       "02_lSNPs",       cfg)

# TSV mirror of the proxy table for archival / sharing
proxy_tsv <- file.path(cfg$paths$intermediate, "02_proxy_table.tsv")
write.table(proxy_table, proxy_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("  saved: ", proxy_tsv)

if (nrow(failed_table) > 0) save_review(failed_table, "02_failures", cfg)

# -----------------------------------------------------------------------------
# 2.8 — Summary
# -----------------------------------------------------------------------------
total_attempted  <- nrow(snp_pop_map)
failed_pairs     <- nrow(failed_table)
successful_pairs <- total_attempted - failed_pairs

bar <- paste(rep("\u2500", 41), collapse = "")
message("\n", bar)
message("LDproxy processing summary")
message(bar)
message("Total SNP-population pairs attempted : ", total_attempted)
message("Successfully processed pairs         : ", successful_pairs)
message("Failed / excluded pairs              : ", failed_pairs)
if (nrow(failed_table) > 0) {
  message("  Failed pairs:")
  for (i in seq_len(nrow(failed_table))) {
    message("  \u2192 ", failed_table$query_snp[i], " / ",
            failed_table$population[i], " : ", failed_table$reason[i])
  }
}
message("")
message("Unique SNPs attempted                : ", length(attempted_snps))
message("Unique SNPs with \u22651 result          : ", length(successful_snps))
message("SNPs fully excluded (all pops failed): ", length(fully_failed_snps))
if (length(fully_failed_snps) > 0) {
  message("  \u2192 ", paste(fully_failed_snps, collapse = ", "))
}
message(bar)

message("\nStage 2 complete.")
