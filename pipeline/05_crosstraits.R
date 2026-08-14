# =============================================================================
# Stage 5 — Cross-trait retrieval
#
# For every SNP in every block (both lSNPs and pSNPs), pull every other
# significant association from the GWAS Catalog. This is the raw set of
# pleiotropic relationships, before any quality filtering.
#
# What this stage does NOT do:
#   - Quality filtering (sample size, CNV, GxE) — that's Stage 6
#   - Ancestry concordance against the block's own ancestry — also Stage 6
#
# Manual checkpoint
# -----------------
# If any cross-trait studies can't be ancestry-mapped, they're flagged in
# 05_unmapped_ancestry_review.tsv. Update cfg$ancestry_mapping in the YAML
# if you see legitimate categories you want to add, then re-run.
#
# Input  : results/intermediate/04_blocks.rds
#          results/intermediate/01_gwassos_sigP.rds
#          GWAS Catalog ancestry file (path from cfg)
#
# Output : results/intermediate/05_crosstraits_raw.rds
#          results/intermediate/05_crosstraits_raw.tsv
#          results/review/05_unmapped_ancestry_review.tsv  (if any)
#
# Run from the project root:
#     Rscript pipeline/05_crosstraits.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(5, "Cross-trait retrieval")

# -----------------------------------------------------------------------------
# 5.1 — Load
# -----------------------------------------------------------------------------
blocks       <- load_intermediate("04_blocks", cfg)
gwassos_sigP <- load_intermediate("01_gwassos_sigP", cfg)
# Needed for the positional arm of the union join (5.3): supplies each proxy's
# coordinate and allele pair.
proxy_table  <- load_intermediate("02_proxy_table", cfg)

# Defensive checks: stage will fail clearly if the columns it expects are missing
required_cols <- c("Block", "lSNPs", "pSNPs")
missing <- setdiff(required_cols, names(blocks))
if (length(missing) > 0) {
  stop("Stage 5 expects columns ", paste(required_cols, collapse = ", "),
       " in 04_blocks.rds. Missing: ", paste(missing, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 5.2 — Explode block -> SNP, tagging lSNP vs pSNP
# -----------------------------------------------------------------------------
explode_to_snps <- function(blocks_df, snp_col, snp_type) {
  rows <- lapply(seq_len(nrow(blocks_df)), function(i) {
    snps <- unlist(blocks_df[[snp_col]][[i]])
    if (length(snps) == 0) return(NULL)
    tibble::tibble(
      Block    = blocks_df$Block[i],
      SNP      = snps,
      SNP_type = snp_type
    )
  })
  bind_rows(rows)
}

block_snp_long <- bind_rows(
  explode_to_snps(blocks, "lSNPs", "lSNP"),
  explode_to_snps(blocks, "pSNPs", "pSNP")
)

message("Block-SNP table: ",
        nrow(block_snp_long), " rows (",
        sum(block_snp_long$SNP_type == "lSNP"), " lSNP, ",
        sum(block_snp_long$SNP_type == "pSNP"), " pSNP).")

# -----------------------------------------------------------------------------
# 5.3 — Join against the GWAS Catalog (significant associations)
# -----------------------------------------------------------------------------
# Columns we want from the catalog. Anything missing is silently skipped.
catalog_cols <- c(
  "SNPS", "P-VALUE", "OR or BETA", "MAPPED_TRAIT", "DISEASE/TRAIT",
  "STUDY", "95% CI (TEXT)", "STRONGEST SNP-RISK ALLELE",
  "RISK ALLELE FREQUENCY", "INITIAL SAMPLE SIZE",
  "REPLICATION SAMPLE SIZE", "PUBMEDID", "MAPPED_GENE",
  "FIRST AUTHOR", "DATE", "JOURNAL", "GENOTYPING TECHNOLOGY",
  "CNV", "STUDY ACCESSION",
  # Coordinates drive the positional arm of the union join below.
  "CHR_ID", "CHR_POS"
)
catalog_cols <- intersect(catalog_cols, names(gwassos_sigP))
for (rc in c("CHR_ID", "CHR_POS")) {
  if (!rc %in% catalog_cols) {
    stop("Stage 5 needs '", rc, "' in 01_gwassos_sigP.rds for the positional ",
         "join arm, but it is absent. Re-run stage 1.")
  }
}

# UNION JOIN: position OR rsID.
#
# rsIDs are not stable across sources (see make_position_key() in R/utils.R and
# .context/decisions.md B4). Joining on the rsID string alone silently missed
# 941 leptin and 156 CLL associations, including the entire FTO locus.
#
# The union is deliberate rather than a straight switch to positional matching:
#   - the positional arm ADDS recall (merged / re-accessioned variants);
#   - the rsID arm is retained as a fallback, so nothing previously retrieved
#     can be lost -- the join is MONOTONIC and has no loss to account for.
# Positional matching alone would have dropped 2 CLL rows: one with a +/-1
# indel-alignment offset, one lead SNP with no proxy-table coordinate at all.

# --- coordinates for the block SNPs -------------------------------------------
# Lead SNPs are sourced from the CATALOG's own coordinates, not the proxy
# table: 8 of 29 leptin leads never appear in their own LDproxy result, which
# is precisely where the missing-coordinate gap came from.
proxy_coords <- proxy_table %>%
  dplyr::filter(is_valid_rsid(RS_Number)) %>%
  dplyr::transmute(
    SNP     = RS_Number,
    chr     = sub(":.*$", "", sub("^chr", "", Coord)),
    pos     = sub("^.*:", "", Coord),
    Alleles = Alleles) %>%
  dplyr::distinct()

catalog_coords <- gwassos_sigP %>%
  dplyr::filter(is_valid_rsid(SNPS)) %>%
  dplyr::transmute(SNP = SNPS, chr = CHR_ID, pos = CHR_POS,
                   Alleles = NA_character_) %>%
  dplyr::distinct(SNP, .keep_all = TRUE)

snp_coords <- dplyr::bind_rows(
  proxy_coords,
  dplyr::anti_join(catalog_coords, proxy_coords, by = "SNP")
) %>% dplyr::distinct(SNP, .keep_all = TRUE)

block_keyed <- block_snp_long %>%
  dplyr::left_join(snp_coords, by = "SNP") %>%
  dplyr::mutate(pos_key = make_position_key(chr, pos, "block SNP coordinates"))

catalog_keyed <- gwassos_sigP[, catalog_cols] %>%
  dplyr::mutate(cat_row = dplyr::row_number(),
                pos_key = make_position_key(CHR_ID, CHR_POS, "catalog coordinates"))

# --- arm 1: rsID. Sentinels are REJECTED, not filtered ------------------------
stopifnot(!any(block_keyed$SNP[is_valid_rsid(block_keyed$SNP)] == "."))
n_sentinel <- sum(!is_valid_rsid(block_keyed$SNP))
if (n_sentinel > 0) {
  message("  rsID arm: rejecting ", n_sentinel,
          " block entry/entries that are not a valid rsID (e.g. \".\"); ",
          "they remain reachable on the positional arm.")
}
hit_rsid <- block_keyed %>%
  dplyr::filter(is_valid_rsid(SNP)) %>%
  dplyr::inner_join(catalog_keyed %>% dplyr::filter(is_valid_rsid(SNPS)),
                    by = c("SNP" = "SNPS"), suffix = c("", ".cat"))

# --- arm 2: position, with the allele-membership constraint -------------------
hit_pos <- block_keyed %>%
  dplyr::filter(!is.na(pos_key)) %>%
  dplyr::inner_join(catalog_keyed %>% dplyr::filter(!is.na(pos_key)),
                    by = "pos_key", suffix = c("", ".cat")) %>%
  dplyr::filter(alleles_compatible(`STRONGEST SNP-RISK ALLELE`, Alleles))

n_allele_rejected <- nrow(
  block_keyed %>% dplyr::filter(!is.na(pos_key)) %>%
    dplyr::inner_join(catalog_keyed %>% dplyr::filter(!is.na(pos_key)),
                      by = "pos_key", suffix = c("", ".cat"))) - nrow(hit_pos)

# --- union, then DEDUPLICATE at the association level -------------------------
# A catalog row can reach the same block through BOTH arms, or through two
# different block SNPs that are the same physical variant under two accessions.
# Counting it twice would inflate association counts, trait frequencies, and
# therefore the enrichment and concordance -- silently, and in the direction
# that flatters the new numbers. Identity is (Block, catalog row).
combined <- dplyr::bind_rows(
  hit_rsid %>% dplyr::mutate(match_route = "rsid"),
  hit_pos  %>% dplyr::mutate(match_route = "position"))

n_before_dedup <- nrow(combined)
crosstraits <- combined %>%
  dplyr::arrange(Block, cat_row, match_route) %>%
  dplyr::distinct(Block, cat_row, .keep_all = TRUE)
n_dupes <- n_before_dedup - nrow(crosstraits)

message("Joining the GWAS Catalog (union: position OR rsID)")
message("  rsID arm      : ", nrow(hit_rsid), " (block, association) pairs")
message("  position arm  : ", nrow(hit_pos),  " (block, association) pairs")
message("  allele constraint rejected: ", n_allele_rejected)
message("  union before dedup        : ", n_before_dedup)
message("  duplicates removed        : ", n_dupes,
        sprintf(" (%.1f%%)", 100 * n_dupes / max(n_before_dedup, 1)))
message("After joining the GWAS Catalog: ", nrow(crosstraits), " associations.")

route <- table(crosstraits$match_route)
message("  retained by route: ",
        paste(sprintf("%s=%d", names(route), as.integer(route)), collapse = ", "))

crosstraits <- crosstraits %>%
  dplyr::select(-dplyr::any_of(c("chr", "pos", "Alleles", "pos_key",
                                 "CHR_ID", "CHR_POS", "cat_row")))

# -----------------------------------------------------------------------------
# 5.4 — Drop self-phenotype rows
# -----------------------------------------------------------------------------
self_pattern <- paste0("\\b(", paste(cfg$phenotype$trait_terms, collapse = "|"), ")\\b")
self_idx <- grepl(self_pattern, crosstraits$MAPPED_TRAIT, ignore.case = TRUE)

n_self <- sum(self_idx)
crosstraits <- crosstraits[!self_idx, ]

message("Dropped ", n_self, " self-phenotype rows for '",
        cfg$phenotype$name, "'.")
message("Remaining cross-trait associations: ", nrow(crosstraits), ".")

# -----------------------------------------------------------------------------
# 5.5 — Ancestry encoding (same cohort-proportion rule as Stage 1)
# -----------------------------------------------------------------------------
crosstraits$`Encoded Populations` <- encode_ancestries_for_rows(crosstraits)

# -----------------------------------------------------------------------------
# 5.6 — Apply reviewed ancestry decisions, if a reviewed file exists
# -----------------------------------------------------------------------------
review_label <- "05_unmapped_ancestry"
reviewed_path <- file.path(cfg$paths$review,
                           paste0(review_label, "_reviewed.tsv"))

if (file.exists(reviewed_path)) {
  reviewed <- load_review(review_label, cfg, required = TRUE)
  if (!"STUDY ACCESSION" %in% names(reviewed)) {
    stop("Reviewed file missing 'STUDY ACCESSION' column.")
  }
  if (!"Encoded Populations" %in% names(reviewed)) {
    stop("Reviewed file must include an 'Encoded Populations' column ",
         "filled in for the studies you want to keep.")
  }

  # Build accession -> encoding lookup, keep only rows the analyst filled in
  fixed <- reviewed[!is.na(reviewed$`Encoded Populations`) &
                    reviewed$`Encoded Populations` != "" &
                    reviewed$`Encoded Populations` != "Check study", ]
  fix_lookup <- setNames(fixed$`Encoded Populations`,
                         fixed$`STUDY ACCESSION`)

  apply_idx <- crosstraits$`Encoded Populations` == "Check study" &
               crosstraits$`STUDY ACCESSION` %in% names(fix_lookup)
  crosstraits$`Encoded Populations`[apply_idx] <-
    fix_lookup[crosstraits$`STUDY ACCESSION`[apply_idx]]
  message("Applied reviewed ancestry decisions to ",
          sum(apply_idx), " cross-trait rows from ",
          length(unique(crosstraits$`STUDY ACCESSION`[apply_idx])),
          " studies.")
}

# -----------------------------------------------------------------------------
# 5.7 — Manual checkpoint: rows whose ancestry could not be determined
# -----------------------------------------------------------------------------
needs_review <- crosstraits %>%
  filter(`Encoded Populations` == "Check study") %>%
  distinct(`STUDY ACCESSION`, PUBMEDID, `FIRST AUTHOR`,
           `INITIAL SAMPLE SIZE`, `REPLICATION SAMPLE SIZE`, STUDY) %>%
  mutate(`Encoded Populations` = "")   # blank column for analyst to fill in

if (nrow(needs_review) > 0) {
  message("\u26a0  ", nrow(needs_review),
          " cross-trait studies still couldn't be ancestry-encoded.")
  message("   Edit results/review/", review_label, "_review.tsv: fill in the ",
          "'Encoded Populations' column with comma-joined superpops ",
          "(e.g. EUR or EUR,EAS).")
  message("   Save as ", review_label, "_reviewed.tsv and re-run.")
  save_review(needs_review, review_label, cfg)
} else {
  message("\u2713  All cross-trait studies are ancestry-encoded.")
}

# -----------------------------------------------------------------------------
# 5.7 — Add block-level ancestry for downstream concordance check (Stage 6)
# -----------------------------------------------------------------------------
if ("Ethnicities" %in% names(blocks)) {
  block_anc_lookup <- setNames(blocks$Ethnicities, blocks$Block)
  crosstraits$`Block Population` <- block_anc_lookup[crosstraits$Block]
} else {
  warning("'Ethnicities' column missing from blocks table; ",
          "Stage 6 will need to derive block ancestry on the fly.")
}

# -----------------------------------------------------------------------------
# 5.8 — Sort by block (natural order) and save
# -----------------------------------------------------------------------------
crosstraits <- crosstraits[mixedorder(crosstraits$Block), ]
rownames(crosstraits) <- NULL

save_intermediate(crosstraits, "05_crosstraits_raw", cfg)

tsv_path <- file.path(cfg$paths$intermediate, "05_crosstraits_raw.tsv")
write.table(crosstraits, tsv_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
message("  saved: ", tsv_path)

# -----------------------------------------------------------------------------
# 5.9 — Summary
# -----------------------------------------------------------------------------
bar <- paste(rep("\u2500", 41), collapse = "")
message("\n", bar)
message("Cross-trait retrieval summary")
message(bar)
message("Cross-trait associations          : ", nrow(crosstraits))
message("  from lSNPs                      : ", sum(crosstraits$SNP_type == "lSNP"))
message("  from pSNPs                      : ", sum(crosstraits$SNP_type == "pSNP"))
message("Unique cross-trait SNPs           : ", n_distinct(crosstraits$SNP))
message("Unique mapped traits              : ", n_distinct(crosstraits$MAPPED_TRAIT))
message("Unique blocks with cross-traits   : ", n_distinct(crosstraits$Block),
        " of ", nrow(blocks))
message(bar)

message("\nStage 5 complete.")