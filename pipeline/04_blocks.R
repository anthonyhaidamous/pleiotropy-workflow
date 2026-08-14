# =============================================================================
# Stage 4 — Block construction and merging
#
# Aggregates lead SNPs into LD-defined blocks (LD-fusion via connected
# components), then optionally merges blocks specified in cfg$block_merges
# (filled in by the analyst after reviewing the overlap report on first run).
#
# Input  : results/intermediate/02_lSNPs.rds
#          results/intermediate/03_leadSNPs_LDs.rds
#          results/intermediate/03_pSNPs_per_lSNP.rds
#
# Output : results/intermediate/04_blocks.rds                — final blocks (rich)
#          results/intermediate/04_blocks.tsv                — flattened TSV mirror
#          results/review/04_block_overlaps_review.tsv       — overlap pairs to inspect
#
# Run from the project root:
#     Rscript pipeline/04_blocks.R config/<trait>.yml
# =============================================================================

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else "config/leptin.yml"
cfg <- load_config(cfg_path)

stage_banner(4, "Block construction and merging")

# -----------------------------------------------------------------------------
# 4.1 — Load
# -----------------------------------------------------------------------------
lSNPs           <- load_intermediate("02_lSNPs",          cfg)
leadSNPs_LDs    <- load_intermediate("03_leadSNPs_LDs",   cfg)
pSNPs_per_lSNP  <- load_intermediate("03_pSNPs_per_lSNP", cfg)

# -----------------------------------------------------------------------------
# 4.2 — Per-lSNP descriptors: chr, position bounds, SNP set, pSNP set
# -----------------------------------------------------------------------------
parse_coord <- function(coords) {
  parts <- strsplit(coords, ":", fixed = TRUE)
  list(chr = vapply(parts, `[`, character(1), 1),
       pos = as.numeric(vapply(parts, `[`, character(1), 2)))
}

per_lSNP <- leadSNPs_LDs %>%
  group_by(query_snp) %>%
  summarise(
    LDblock_SNPs   = list(unique(RS_Number)),
    LDblock_coords = list(unique(Coord)),
    .groups = "drop"
  )

parsed <- lapply(per_lSNP$LDblock_coords, parse_coord)
per_lSNP$chr             <- vapply(parsed, function(p) unique(p$chr)[1], character(1))
per_lSNP$LDblock_pos_min <- vapply(parsed, function(p) min(p$pos), numeric(1))
per_lSNP$LDblock_pos_max <- vapply(parsed, function(p) max(p$pos), numeric(1))
per_lSNP$LDblock_span    <- per_lSNP$LDblock_pos_max - per_lSNP$LDblock_pos_min

# Lead SNP's own genomic position (Distance == 0 row)
self_rows <- leadSNPs_LDs[leadSNPs_LDs$Distance == 0,
                          c("query_snp", "Coord")]
self_rows <- self_rows[!duplicated(self_rows$query_snp), ]
self_pos <- as.numeric(sub("^.*:", "", self_rows$Coord))
names(self_pos) <- self_rows$query_snp
per_lSNP$lSNP_pos <- self_pos[per_lSNP$query_snp]

# Attach pSNPs from Stage 3
per_lSNP <- per_lSNP %>%
  left_join(
    pSNPs_per_lSNP %>%
      select(query_snp,
             LDblock_pSNPs    = pSNPs,
             LDblock_nbpSNPs  = n_pSNPs),
    by = "query_snp"
  ) %>%
  rename(lSNP = query_snp) %>%
  select(-LDblock_coords)

# Sort by chromosome (natural order) then position so block names follow
# genomic order
chr_order <- gtools::mixedsort(unique(per_lSNP$chr))
per_lSNP <- per_lSNP %>%
  arrange(match(chr, chr_order), lSNP_pos)

# -----------------------------------------------------------------------------
# 4.3 — LD-fusion via connected components
# Two lSNPs are in the same component if either appears in the other's
# LDblock_SNPs list. Compute components by iteratively merging any sets
# that share an element.
# -----------------------------------------------------------------------------
n <- nrow(per_lSNP)

fusion_sets <- lapply(seq_len(n), function(i) {
  hits <- which(vapply(per_lSNP$LDblock_SNPs,
                       function(snps) per_lSNP$lSNP[i] %in% snps,
                       logical(1)))
  unique(c(i, hits))
})

components <- list()
for (s in fusion_sets) {
  hit <- which(vapply(components, function(g) any(g %in% s), logical(1)))
  if (length(hit) == 0) {
    components <- c(components, list(s))
  } else {
    components <- c(components[-hit],
                    list(unique(c(unlist(components[hit]), s))))
  }
}

# Order components so block names go in genomic order (rows are pre-sorted)
components <- components[order(vapply(components, min, integer(1)))]

per_lSNP$blockname <- character(n)
for (k in seq_along(components)) {
  per_lSNP$blockname[components[[k]]] <- paste0("B", k)
}

n_components <- length(components)
n_fused      <- sum(lengths(components) > 1)
message("Auto-fused ", n, " lead SNPs into ", n_components,
        " LD-defined blocks (", n_fused, " block",
        ifelse(n_fused == 1, "", "s"), " contain multiple lead SNPs).")

# -----------------------------------------------------------------------------
# 4.4 — Aggregate to one row per block
# -----------------------------------------------------------------------------
aggregate_to_blocks <- function(rows_df) {
  agg <- rows_df %>%
    group_by(blockname) %>%
    summarise(
      lSNP            = list(lSNP),
      lSNP_pos        = list(lSNP_pos),
      chr             = unique(chr)[1],
      LDblock_pos_min = min(LDblock_pos_min),
      LDblock_pos_max = max(LDblock_pos_max),
      LDblock_SNPs    = list(unique(unlist(LDblock_SNPs))),
      LDblock_pSNPs   = list(unique(unlist(LDblock_pSNPs))),
      .groups = "drop"
    )
  agg %>%
    mutate(
      LDblock_pSNPs   = lapply(seq_along(lSNP),
                               function(i) setdiff(LDblock_pSNPs[[i]],
                                                   lSNP[[i]])),
      LDblock_span    = LDblock_pos_max - LDblock_pos_min,
      LDblock_nblSNPs = lengths(lSNP),
      LDblock_nbSNPs  = lengths(LDblock_SNPs),
      LDblock_nbpSNPs = lengths(LDblock_pSNPs)
    )
}

blocks <- aggregate_to_blocks(per_lSNP) %>%
  arrange(match(chr, chr_order), LDblock_pos_min)

# -----------------------------------------------------------------------------
# 4.5 — Find genomically overlapping blocks → review file
# -----------------------------------------------------------------------------
overlap_pairs <- tibble::tibble(blockA = character(), blockB = character(),
                                chr = character(),
                                A_range = character(), B_range = character())
if (nrow(blocks) > 1) {
  for (i in 1:(nrow(blocks) - 1)) {
    for (j in (i + 1):nrow(blocks)) {
      if (blocks$chr[i] != blocks$chr[j]) next
      if (blocks$LDblock_pos_min[i] <= blocks$LDblock_pos_max[j] &&
          blocks$LDblock_pos_min[j] <= blocks$LDblock_pos_max[i]) {
        overlap_pairs <- bind_rows(overlap_pairs, tibble::tibble(
          blockA  = blocks$blockname[i],
          blockB  = blocks$blockname[j],
          chr     = blocks$chr[i],
          A_range = paste0(blocks$LDblock_pos_min[i], "-",
                           blocks$LDblock_pos_max[i]),
          B_range = paste0(blocks$LDblock_pos_min[j], "-",
                           blocks$LDblock_pos_max[j])
        ))
      }
    }
  }
}

if (nrow(overlap_pairs) > 0) {
  message("Found ", nrow(overlap_pairs),
          " pair(s) of genomically overlapping blocks. ",
          "Review and add to cfg$block_merges if you want them merged.")
  save_review(overlap_pairs, "04_block_overlaps", cfg)
} else {
  message("No genomically overlapping blocks detected.")
}

# -----------------------------------------------------------------------------
# 4.6 — Apply YAML-driven merges
# -----------------------------------------------------------------------------
merge_groups <- cfg$block_merges
if (!is.null(merge_groups) && length(merge_groups) > 0) {
  referenced <- unique(unlist(merge_groups))
  missing    <- setdiff(referenced, blocks$blockname)
  if (length(missing) > 0) {
    warning("cfg$block_merges references unknown blocks (skipped): ",
            paste(missing, collapse = ", "))
  }

  applied <- 0
  for (group in merge_groups) {
    group <- intersect(group, blocks$blockname)
    if (length(group) < 2) next
    new_name <- group[1]
    per_lSNP$blockname[per_lSNP$blockname %in% group] <- new_name
    applied <- applied + 1
  }

  if (applied > 0) {
    blocks <- aggregate_to_blocks(per_lSNP) %>%
      arrange(match(chr, chr_order), LDblock_pos_min)
    message("Applied ", applied, " merge group(s) from cfg$block_merges.")
  }
}

# Reindex blocknames sequentially in chr+position order
old_names <- blocks$blockname
blocks$blockname <- paste0("B", seq_len(nrow(blocks)))
rename_map <- setNames(blocks$blockname, old_names)
per_lSNP$blockname <- unname(rename_map[per_lSNP$blockname])

message("Final block count: ", nrow(blocks))

# -----------------------------------------------------------------------------
# 4.7 — Block-level metadata + top lead SNP per block (min p-value)
# -----------------------------------------------------------------------------
lsnp_meta_cols <- c("STRONGEST SNP-RISK ALLELE", "P-VALUE", "OR or BETA",
                    "95% CI (TEXT)", "Encoded Populations",
                    "INITIAL SAMPLE SIZE", "MAPPED_GENE",
                    "GENOTYPING TECHNOLOGY")
lSNP_meta <- lSNPs[, c("SNPS", lsnp_meta_cols)]
names(lSNP_meta)[1] <- "lSNP"

lsnp_lookup <- function(snps, col) {
  unname(lSNP_meta[[col]][match(snps, lSNP_meta$lSNP)])
}

blocks$lSNPs_alleles  <- lapply(blocks$lSNP, lsnp_lookup, col = "STRONGEST SNP-RISK ALLELE")
blocks$lSNPs_pval     <- lapply(blocks$lSNP, lsnp_lookup, col = "P-VALUE")
blocks$lSNPs_ORbeta   <- lapply(blocks$lSNP, lsnp_lookup, col = "OR or BETA")
blocks$lSNPs_effectCI <- lapply(blocks$lSNP, lsnp_lookup, col = "95% CI (TEXT)")
blocks$Genes          <- lapply(blocks$lSNP, lsnp_lookup, col = "MAPPED_GENE")
blocks$Genotyping     <- lapply(blocks$lSNP, lsnp_lookup, col = "GENOTYPING TECHNOLOGY")
blocks$Cohorts        <- lapply(blocks$lSNP, lsnp_lookup, col = "INITIAL SAMPLE SIZE")

blocks$Ethnicities <- vapply(seq_len(nrow(blocks)), function(i) {
  pops <- lsnp_lookup(blocks$lSNP[[i]], "Encoded Populations")
  pops <- unlist(strsplit(pops, ",\\s*"))
  paste(sort(unique(trimws(pops[pops != ""]))), collapse = ",")
}, character(1))

blocks$Top_lead_SNP <- vapply(seq_len(nrow(blocks)), function(i) {
  pvals <- as.numeric(blocks$lSNPs_pval[[i]])
  idx   <- which.min(pvals)
  paste(blocks$lSNP[[i]][idx],
        blocks$lSNP_pos[[i]][idx],
        format(pvals[idx], scientific = TRUE),
        blocks$lSNPs_ORbeta[[i]][idx],
        sep = ", ")
}, character(1))

blocks$Genomic_position <- paste0(blocks$chr, ":",
                                   blocks$LDblock_pos_min, "-",
                                   blocks$LDblock_pos_max)

# -----------------------------------------------------------------------------
# 4.8 — Final shape (matches existing leptin15_blocks columns)
# -----------------------------------------------------------------------------
blocks_table <- blocks %>%
  transmute(
    Block                   = blockname,
    Genomic_position        = Genomic_position,
    Genomic_span            = LDblock_span,
    Top_lead_SNP            = Top_lead_SNP,
    Genes                   = Genes,
    Ethnicities             = Ethnicities,
    nb_lSNPs                = LDblock_nblSNPs,
    nb_pSNPs                = LDblock_nbpSNPs,
    nb_SNPs                 = LDblock_nbSNPs,
    lSNPs                   = lSNP,
    lSNPs_positions         = lSNP_pos,
    lSNPs_alleles           = lSNPs_alleles,
    lSNPs_pval              = lSNPs_pval,
    lSNPs_ORbeta            = lSNPs_ORbeta,
    lSNPs_effectCIdirection = lSNPs_effectCI,
    pSNPs                   = LDblock_pSNPs,
    SNPs                    = LDblock_SNPs,
    Genotyping              = Genotyping,
    Cohorts                 = Cohorts
  )

# -----------------------------------------------------------------------------
# 4.9 — Save
# -----------------------------------------------------------------------------
save_intermediate(blocks_table, "04_blocks", cfg)

flat <- as.data.frame(lapply(blocks_table, function(col) {
  if (is.list(col)) sapply(col, paste, collapse = ", ") else col
}), check.names = FALSE)
flat_path <- file.path(cfg$paths$intermediate, "04_blocks.tsv")
write.table(flat, flat_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
message("  saved: ", flat_path)

# -----------------------------------------------------------------------------
# 4.10 — Summary
# -----------------------------------------------------------------------------
bar <- paste(rep("\u2500", 41), collapse = "")
message("\n", bar)
message("Block construction summary")
message(bar)
message("Lead SNPs in:                 ", n)
message("LD-fused blocks (auto):       ", n_components)
message("Final block count:            ", nrow(blocks_table))
message("Blocks with >1 lead SNP:      ", sum(blocks_table$nb_lSNPs > 1))
message("Total LD-set SNPs in blocks:  ", sum(blocks_table$nb_SNPs))
message("Total pSNPs in blocks:        ", sum(blocks_table$nb_pSNPs))
message(bar)

# Re-check overlaps on the FINAL (post-merge) blocks
final_overlaps <- tibble::tibble(blockA = character(), blockB = character(),
                                 chr = character(),
                                 A_range = character(), B_range = character())
if (nrow(blocks_table) > 1) {
  for (i in 1:(nrow(blocks_table) - 1)) {
    for (j in (i + 1):nrow(blocks_table)) {
      bi <- blocks_table[i, ]; bj <- blocks_table[j, ]
      chr_i <- sub(":.*", "", bi$Genomic_position)
      chr_j <- sub(":.*", "", bj$Genomic_position)
      min_i <- as.numeric(sub(".*:(\\d+)-.*", "\\1", bi$Genomic_position))
      max_i <- as.numeric(sub(".*-",            "",      bi$Genomic_position))
      min_j <- as.numeric(sub(".*:(\\d+)-.*", "\\1", bj$Genomic_position))
      max_j <- as.numeric(sub(".*-",            "",      bj$Genomic_position))
      if (chr_i != chr_j) next
      if (min_i <= max_j && min_j <= max_i) {
        final_overlaps <- bind_rows(final_overlaps, tibble::tibble(
          blockA = bi$Block, blockB = bj$Block, chr = chr_i,
          A_range = paste0(min_i, "-", max_i),
          B_range = paste0(min_j, "-", max_j)
        ))
      }
    }
  }
}
if (nrow(final_overlaps) > 0) {
  message("\n\u26a0  ", nrow(final_overlaps),
          " overlap(s) remain after merging. ",
          "Add them to cfg$block_merges and re-run.")
  save_review(final_overlaps, "04_block_overlaps", cfg)
} else {
  message("\n\u2713  No overlapping blocks in final output.")
  review_path <- file.path(cfg$paths$review, "04_block_overlaps_review.tsv")
  if (file.exists(review_path)) file.remove(review_path)
}

message("\nStage 4 complete.")