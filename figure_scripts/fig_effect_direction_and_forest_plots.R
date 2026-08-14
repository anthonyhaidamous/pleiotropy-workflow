#!/usr/bin/env Rscript
# =============================================================================
# Effect direction harmonization + concordance analysis + forest plots
#
# Self-contained: re-derives the effect-direction harmonization from the
# pipeline's stage outputs (it does NOT consume 09_harmonized_trait_table).
#
# Inputs (each accepts .tsv or .rds):
#   --trait-table   results/intermediate/<trait>/07_trait_table.{tsv,rds}
#   --blocks        results/intermediate/<trait>/04_blocks.{tsv,rds}
#   --ldlink        results/intermediate/<trait>/02_proxy_table.{tsv,rds}
#   --crosstraits   results/intermediate/<trait>/06_crosstraits_filtered.{tsv,rds}
#
# Outputs (all under --outdir):
#   direction_summary_by_trait.txt
#   direction_summary_by_category.txt
#   direction_summary_by_block.txt
#   direction_summary_by_hierarchized_trait.txt
#   direction_concordance_by_category.{png,svg}
#   direction_concordance_by_block.png
#   forest_plots/<trait>_forest_plot.png        (unless --no-forest-plots)
#   concordance_per_cat/direction_concordance_<cat>.{png,svg}
#
# Usage:
#   Rscript figure_scripts/fig_effect_direction_and_forest_plots.R \
#       --trait-table results/intermediate/leptin/07_trait_table.rds \
#       --blocks      results/intermediate/leptin/04_blocks.rds \
#       --ldlink      results/intermediate/leptin/02_proxy_table.rds \
#       --crosstraits results/intermediate/leptin/06_crosstraits_filtered.rds \
#       --outdir      figures/leptin
#
#   # skip the per-trait forest plots (slow):
#   Rscript figure_scripts/fig_effect_direction_and_forest_plots.R ... --no-forest-plots
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tools)
  library(stringr)
})

rename <- dplyr::rename

# SVG device: prefer svglite (cleaner output) but fall back to the built-in
# cairo svg() device so the script runs without svglite installed.
# fix_text_size = FALSE keeps text editable (no textLength/lengthAdjust) in
# downstream vector editors; the cairo fallback has no such attributes anyway.
svg_device <- if (requireNamespace("svglite", quietly = TRUE)) {
  function(...) svglite::svglite(..., fix_text_size = FALSE)
} else grDevices::svg

option_list <- list(
  make_option(c("--trait-table"), type = "character",
              help = "07_trait_table file (.tsv or .rds) [required]"),
  make_option(c("--blocks"), type = "character",
              help = "04_blocks file (.tsv or .rds) [required]"),
  make_option(c("--ldlink"), type = "character",
              help = "02_proxy_table file (.tsv or .rds) [required]"),
  make_option(c("--crosstraits"), type = "character",
              help = "06_crosstraits_filtered file (.tsv or .rds) [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = "figures",
              help = "Output directory [default: %default]"),
  make_option(c("--no-forest-plots"), action = "store_true", default = FALSE,
              help = "Skip per-trait forest plot generation"),
  make_option(c("--drop-others"), action = "store_true", default = FALSE,
              help = "Exclude the 'Others' catch-all category from the concordance/discordance summaries (and its multiple-testing correction)"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Effect direction analysis and forest plots")
opt <- parse_args(parser)

required_args <- c("trait-table", "blocks", "ldlink", "crosstraits")
for (arg in required_args) {
  val <- opt[[arg]]
  if (is.null(val)) {
    cat("ERROR: --", arg, " is required\n\n", sep = "", file = stderr())
    print_help(parser)
    quit(status = 2)
  }
  if (!file.exists(val)) {
    cat("ERROR: File not found:", val, "\n", file = stderr())
    quit(status = 2)
  }
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================

# Read a pipeline output that may be supplied as .rds or .tsv
source("R/figure_utils.R")
source("R/palettes.R")

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged. The existing
# convention is cool = concordant, warm = discordant; `neg` and `pos` are every
# palette's cool and warm poles, so that reading survives a palette change.
pal      <- get_palette(opt$palette)
col_ink  <- pal_col(pal, "ink",     "black")
col_grid <- pal_col(pal, "grid",    "gray90")
col_conc <- pal_col(pal, "neg",     "steelblue")
col_disc <- pal_col(pal, "pos",     "firebrick")
col_ns   <- pal_col(pal, "neutral", "grey70")

# Normalise a single 04_blocks list-field cell to a character vector.
# Handles both .rds (list column -> vector element) and .tsv (comma string).
as_items <- function(cell) {
  if (is.list(cell)) cell <- cell[[1]]
  cell <- as.character(cell)
  if (length(cell) == 1) cell <- strsplit(cell, ",\\s*")[[1]]
  trimws(cell)
}

# Extract the two numeric bounds from a 95% CI text string
ci_is_beta     <- function(ci) grepl("increase|decrease", ci, ignore.case = TRUE)
ci_is_decrease <- function(ci) grepl("decrease", ci, ignore.case = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================

trait_table <- read_any(opt[["trait-table"]])
blocks      <- read_any(opt$blocks)
proxy       <- read_any(opt$ldlink)
ct          <- read_any(opt$crosstraits)

message("Loaded: ", nrow(ct), " cross-traits, ", nrow(blocks), " blocks, ",
        nrow(trait_table), " trait-table rows.")

# =============================================================================
# 1. PER-BLOCK PHENOTYPE DIRECTION  (from 04_blocks lead-SNP lists)
#
# For each block, the representative lead SNP is the one with the lowest
# p-value; its signed effect (negative if the CI text says "decrease")
# defines the block's phenotype direction (Positive / Negative).
# =============================================================================

block_dir <- do.call(rbind, lapply(seq_len(nrow(blocks)), function(i) {
  rs  <- as_items(blocks$lSNPs[[i]])
  pv  <- suppressWarnings(as.numeric(as_items(blocks$lSNPs_pval[[i]])))
  orb <- suppressWarnings(as.numeric(as_items(blocks$lSNPs_ORbeta[[i]])))
  dir <- as_items(blocks$lSNPs_effectCIdirection[[i]])

  k <- if (all(is.na(pv))) 1L else which.min(pv)
  d <- if (k <= length(dir)) dir[k] else ""
  b <- if (k <= length(orb)) orb[k] else NA_real_
  signed <- if (is.na(b)) NA_real_ else if (grepl("decrease", d, ignore.case = TRUE)) -abs(b) else abs(b)

  data.frame(
    Block = blocks$Block[i],
    phenotype_direction = if (is.na(signed)) NA_character_ else if (signed < 0) "Negative" else "Positive",
    stringsAsFactors = FALSE
  )
}))

n_pos <- sum(block_dir$phenotype_direction == "Positive", na.rm = TRUE)
n_neg <- sum(block_dir$phenotype_direction == "Negative", na.rm = TRUE)
message("Block phenotype direction: ", n_pos, " positive, ", n_neg, " negative, ",
        sum(is.na(block_dir$phenotype_direction)), " indeterminate.")

# =============================================================================
# 2. ALLELE LOOKUP  (from 02_proxy_table)
#
# The LDproxy "Alleles" pair (A/B) is a property of the variant itself, so we
# key on RS_Number alone. CA1 = first allele, CA2 = second allele.
# =============================================================================

allele_lookup <- proxy %>%
  distinct(RS_Number, Alleles) %>%
  filter(!is.na(Alleles), Alleles != "") %>%
  separate(Alleles, into = c("CA1_raw", "CA2_raw"), sep = "/",
           remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(CA1 = gsub("[^A-Za-z]", "", CA1_raw),
         CA2 = gsub("[^A-Za-z]", "", CA2_raw)) %>%
  select(RS_Number, CA1, CA2) %>%
  filter(!is.na(CA1), !is.na(CA2), CA1 != "", CA2 != "") %>%
  distinct(RS_Number, .keep_all = TRUE)

# =============================================================================
# 3. HARMONIZE CROSS-TRAIT EFFECTS  (06_crosstraits_filtered)
# =============================================================================

# Parse SNP + risk allele; drop unknown risk alleles
ct <- ct %>%
  separate(`STRONGEST SNP-RISK ALLELE`,
           into = c("snp_for_RA", "RiskAllele"),
           sep = "-", remove = FALSE, fill = "right", extra = "merge") %>%
  filter(!RiskAllele %in% c("?", "-", ""), !is.na(RiskAllele))

# (a) OR -> beta where CI text is not on a beta scale
ct$beta <- suppressWarnings(as.numeric(ct$`OR or BETA`))
or_rows  <- !ci_is_beta(ct$`95% CI (TEXT)`) & !is.na(ct$beta)
ct$beta[or_rows] <- log(ct$beta[or_rows])
# (b) negative sign where CI text says "decrease"
dec_rows <- ci_is_decrease(ct$`95% CI (TEXT)`) & !is.na(ct$beta)
ct$beta[dec_rows] <- -abs(ct$beta[dec_rows])

# CI bounds
ci_bounds   <- t(vapply(ct$`95% CI (TEXT)`, extract_ci_bounds, numeric(2)))
ct$CI_lower <- ci_bounds[, 1]
ct$CI_upper <- ci_bounds[, 2]
# log OR-scale bounds (positive only); NA out non-log-able ones
or_lo_ok <- or_rows & !is.na(ct$CI_lower) & ct$CI_lower > 0
or_hi_ok <- or_rows & !is.na(ct$CI_upper) & ct$CI_upper > 0
ct$CI_lower[or_lo_ok] <- log(ct$CI_lower[or_lo_ok])
ct$CI_upper[or_hi_ok] <- log(ct$CI_upper[or_hi_ok])
ct$CI_lower[or_rows & !or_lo_ok] <- NA_real_
ct$CI_upper[or_rows & !or_hi_ok] <- NA_real_
# sign for decrease (beta-scale rows)
ct$CI_lower[dec_rows] <- -abs(ct$CI_lower[dec_rows])
ct$CI_upper[dec_rows] <- -abs(ct$CI_upper[dec_rows])

# Allele harmonization: flip beta/CI when risk allele == CA1, so the effect
# is always reported relative to CA2 (= reference / effect allele).
ct <- ct %>% left_join(allele_lookup, by = c("SNP" = "RS_Number"))

flip_idx <- !is.na(ct$CA1) & !is.na(ct$RiskAllele) & ct$RiskAllele == ct$CA1
ct$new_beta     <- ct$beta
ct$new_CI_lower <- ct$CI_lower
ct$new_CI_upper <- ct$CI_upper
ct$new_beta[flip_idx]     <- -ct$beta[flip_idx]
ct$new_CI_lower[flip_idx] <- -ct$CI_upper[flip_idx]
ct$new_CI_upper[flip_idx] <- -ct$CI_lower[flip_idx]
ct$new_RA <- ct$CA2   # effect allele convention: relative to CA2

# Drop rows without allele info (betas can't be harmonized to a common ref)
n_before <- nrow(ct)
ct <- ct[!is.na(ct$CA1) & !is.na(ct$CA2), ]
if (n_before - nrow(ct) > 0)
  message("Dropped ", n_before - nrow(ct), " rows lacking allele info from LDproxy.")

# Ensure CI lower <= upper
swap_idx <- !is.na(ct$new_CI_lower) & !is.na(ct$new_CI_upper) &
            ct$new_CI_lower > ct$new_CI_upper
if (any(swap_idx)) {
  tmp <- ct$new_CI_lower[swap_idx]
  ct$new_CI_lower[swap_idx] <- ct$new_CI_upper[swap_idx]
  ct$new_CI_upper[swap_idx] <- tmp
}

# =============================================================================
# 4. ATTACH CATEGORIES + BLOCK DIRECTION, COLLAPSE TO ONE ROW PER SNP/TRAIT/BLOCK
# =============================================================================

ct$`P-VALUE` <- suppressWarnings(as.numeric(ct$`P-VALUE`))

updated_trait_table <- ct %>%
  left_join(trait_table %>% select(MAPPED_TRAIT, Category, `Hierarchized Trait`),
            by = "MAPPED_TRAIT") %>%
  left_join(block_dir, by = "Block") %>%
  filter(!is.na(Category)) %>%
  group_by(SNP, MAPPED_TRAIT, Block) %>%
  summarise(
    Trait               = first(MAPPED_TRAIT),
    Category            = first(Category),
    `Hierarchized Trait`= first(`Hierarchized Trait`),
    MAPPED_GENE         = paste(unique(na.omit(MAPPED_GENE)[nzchar(MAPPED_GENE)]), collapse = ", "),
    new_RA              = first(new_RA),
    new_beta            = mean(new_beta),
    new_CI_lower        = mean(new_CI_lower),
    new_CI_upper        = mean(new_CI_upper),
    `P-VALUE`           = mean(`P-VALUE`, na.rm = TRUE),
    phenotype_direction = first(phenotype_direction),
    SNP_type            = first(SNP_type),
    .groups             = "drop"
  )

cat("Final table:", nrow(updated_trait_table), "rows,",
    length(unique(updated_trait_table$Trait)), "unique traits\n")

# =============================================================================
# 5. FOREST PLOTS
# =============================================================================

create_forest_plot <- function(data, trait_name) {
  plot_data <- data %>%
    dplyr::filter(tolower(Trait) == tolower(trait_name)) %>%
    dplyr::ungroup() %>%
    filter(!is.na(new_beta)) %>%
    arrange(desc(new_beta)) %>%
    mutate(
      Beta_CI = sprintf("%.3f (%.3f to %.3f)", new_beta, new_CI_lower, new_CI_upper),
      SNP_RA  = paste(SNP, new_RA, sep = "-"),
      pheno_arrow = case_when(
        phenotype_direction == "Negative" ~ "◄",
        phenotype_direction == "Positive" ~ "►",
        TRUE ~ as.character(phenotype_direction)
      ),
      SNP_RA = factor(SNP_RA, levels = unique(SNP_RA))
    )
  if (nrow(plot_data) == 0) stop("No valid entries for: ", trait_name)

  p_info <- plot_data %>%
    ggplot(aes(y = SNP_RA)) +
    geom_text(aes(x = 0, label = SNP_RA), hjust = 0,
              fontface = if_else(plot_data$SNP_type == "lSNP", "bold", "plain")) +
    geom_text(aes(x = 0.15, label = Block),       hjust = 0) +
    geom_text(aes(x = 0.25, label = Beta_CI),     hjust = 0) +
    geom_text(aes(x = 0.5,  label = MAPPED_GENE), hjust = 0) +
    geom_text(aes(x = 0.7,  label = pheno_arrow), hjust = -3.5) +
    annotate("text", x = 0,    y = Inf, label = "SNP",            vjust = 1, hjust = 0, fontface = "bold", size = 4.5) +
    annotate("text", x = 0.15, y = Inf, label = "Block",          vjust = 1, hjust = 0, fontface = "bold", size = 4.5) +
    annotate("text", x = 0.25, y = Inf, label = "Beta (95 % CI)", vjust = 1, hjust = 0, fontface = "bold", size = 4.5) +
    annotate("text", x = 0.5,  y = Inf, label = "Mapped Gene",    vjust = 1, hjust = 0, fontface = "bold", size = 4.5) +
    annotate("text", x = 0.7,  y = Inf, label = "Block direction",vjust = 1, hjust = 0, fontface = "bold", size = 4.5) +
    theme_void() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, nrow(plot_data) + 1))

  p_plot <- ggplot(plot_data, aes(y = SNP_RA)) +
    theme_classic() +
    geom_point(aes(x = new_beta), shape = 15, size = 3) +
    geom_linerange(aes(xmin = new_CI_lower, xmax = new_CI_upper)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = seq(nrow(plot_data)), color = col_ink, linetype = "dotted") +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), axis.title.y = element_blank()) +
    labs(x = "Beta Coefficient") +
    coord_cartesian(ylim = c(0, nrow(plot_data) + 1))

  layout <- c(area(t = 0, l = 0, b = 30, r = 6), area(t = 0, l = 6, b = 30, r = 9))
  wrap_plots(p_info, p_plot, design = layout) +
    plot_annotation(title = toTitleCase(str_to_sentence(trait_name))) &
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
}

if (!opt[["no-forest-plots"]]) {
  forest_dir <- file.path(opt$outdir, "forest_plots")
  dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)
  unique_traits <- unique(updated_trait_table$Trait)
  cat("Generating forest plots for", length(unique_traits), "traits\n")

  for (trait_name in unique_traits) {
    trait_data <- updated_trait_table %>%
      filter(tolower(Trait) == tolower(trait_name), !is.na(new_beta))
    if (nrow(trait_data) == 0) next

    plot_height <- max(2.5, nrow(trait_data) * 0.4)
    file_path   <- file.path(forest_dir, paste0(gsub("[^A-Za-z0-9]", "_", trait_name), "_forest_plot.png"))

    forest_plot <- create_forest_plot(updated_trait_table, trait_name) +
      plot_annotation(title = str_wrap(str_to_sentence(trait_name), width = 50)) &
      theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            text       = element_text(family = "Arial Unicode MS"))

    ggsave(file_path, plot = forest_plot, device = png,
           width = 16, height = plot_height, dpi = 300)
    message("Saved: ", file_path)
  }
} else {
  message("--no-forest-plots set; skipping forest plot generation.")
}

# =============================================================================
# 6. DIRECTION CONCORDANCE  (concordant if trait beta sign matches block dir)
# =============================================================================

updated_trait_table <- updated_trait_table %>%
  mutate(
    pheno_sign      = ifelse(phenotype_direction == "Positive", 1, -1),
    trait_sign      = sign(new_beta),
    direction_match = case_when(
      is.na(new_beta) | is.na(phenotype_direction) ~ NA_character_,
      pheno_sign == trait_sign ~ "Concordant",
      pheno_sign != trait_sign ~ "Discordant",
      TRUE                     ~ NA_character_
    )
  )

# Optionally drop the 'Others' catch-all so it is neither tested nor counted
# in the BH correction. Scoped to the concordance summaries below; the forest
# plots above are unaffected.
if (opt[["drop-others"]]) {
  n_before <- nrow(updated_trait_table)
  updated_trait_table <- updated_trait_table %>% filter(tolower(Category) != "others")
  message("--drop-others: removed ", n_before - nrow(updated_trait_table),
          " 'Others' row(s) from the concordance/discordance analysis.")
}

cat("Direction match distribution:\n")
print(table(updated_trait_table$direction_match, useNA = "always"))

binom_p_vec <- function(k, n) mapply(function(kk, nn) {
  if (nn < 3) return(NA_real_)
  binom.test(kk, nn, p = 0.5, alternative = "two.sided")$p.value
}, k, n)

# --- Trait-level ---
direction_trait <- updated_trait_table %>%
  filter(!is.na(direction_match)) %>%
  group_by(Trait, Category) %>%
  summarise(
    n_snps = n(),
    n_concordant = sum(direction_match == "Concordant"),
    n_discordant = sum(direction_match == "Discordant"),
    pct_concordant = round(100 * n_concordant / n_snps, 1),
    .groups = "drop"
  ) %>%
  mutate(
    binom_p = binom_p_vec(n_concordant, n_snps),
    binom_p_adj = p.adjust(binom_p, method = "BH"),
    predominant = case_when(
      is.na(binom_p) ~ "Too few SNPs",
      binom_p_adj < 0.05 & pct_concordant > 50 ~ "Significantly concordant",
      binom_p_adj < 0.05 & pct_concordant < 50 ~ "Significantly discordant",
      pct_concordant > 50 ~ "Nominally concordant",
      pct_concordant < 50 ~ "Nominally discordant",
      TRUE ~ "Mixed"
    )
  ) %>%
  arrange(desc(pct_concordant))

# --- Category-level ---
direction_category <- updated_trait_table %>%
  filter(!is.na(direction_match)) %>%
  group_by(Category) %>%
  summarise(
    n_traits = n_distinct(Trait),
    n_snps = n(),
    n_concordant = sum(direction_match == "Concordant"),
    n_discordant = sum(direction_match == "Discordant"),
    pct_concordant = round(100 * n_concordant / n_snps, 1),
    .groups = "drop"
  ) %>%
  mutate(
    binom_p = binom_p_vec(n_concordant, n_snps),
    binom_p_adj = p.adjust(binom_p, method = "BH"),
    predominant = case_when(
      binom_p_adj < 0.05 & pct_concordant > 50 ~ "Significantly concordant",
      binom_p_adj < 0.05 & pct_concordant < 50 ~ "Significantly discordant",
      pct_concordant > 50 ~ "Nominally concordant",
      pct_concordant < 50 ~ "Nominally discordant",
      TRUE ~ "Mixed"
    )
  ) %>%
  arrange(desc(pct_concordant))

# --- Block-level ---
direction_block <- updated_trait_table %>%
  filter(!is.na(direction_match)) %>%
  group_by(Block, phenotype_direction) %>%
  summarise(
    n_traits = n_distinct(Trait),
    n_snps = n(),
    n_concordant = sum(direction_match == "Concordant"),
    n_discordant = sum(direction_match == "Discordant"),
    pct_concordant = round(100 * n_concordant / n_snps, 1),
    .groups = "drop"
  ) %>%
  mutate(
    binom_p = binom_p_vec(n_concordant, n_snps),
    binom_p_adj = p.adjust(binom_p, method = "BH")
  ) %>%
  arrange(Block)

write.table(direction_trait,    file.path(opt$outdir, "direction_summary_by_trait.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(direction_category, file.path(opt$outdir, "direction_summary_by_category.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(direction_block,    file.path(opt$outdir, "direction_summary_by_block.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Direction summaries saved to ", opt$outdir)

# =============================================================================
# 7. CONCORDANCE BAR CHART — by category
# =============================================================================

cat_concordance <- direction_category %>%
  filter(n_snps >= 3) %>%
  mutate(
    Category   = stringr::str_to_sentence(Category),
    Category   = factor(Category, levels = Category[order(pct_concordant)]),
    fill_color = case_when(
      binom_p_adj < 0.05 & pct_concordant > 50 ~ "Sig. concordant",
      binom_p_adj < 0.05 & pct_concordant < 50 ~ "Sig. discordant",
      TRUE ~ "Non-significant"
    ),
    pval_label = case_when(
      is.na(binom_p_adj) ~ "NA",
      binom_p_adj == 1 ~ "1",
      TRUE ~ formatC(binom_p_adj, format = "e", digits = 1)
    )
  )

p_pval <- ggplot(cat_concordance, aes(x = Category, y = 0, label = pval_label)) +
  geom_text(size = 3.5, color = col_ink, hjust = 0.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  theme_void() +
  theme(axis.text.y = element_text(size = 11, face = "bold", hjust = 1),
        plot.margin = unit(c(1, 0, 1, 1), "cm"))

p_bars <- ggplot(cat_concordance, aes(x = Category, y = pct_concordant, fill = fill_color)) +
  geom_col(width = 0.7, color = col_ink, linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dashed", color = col_ink, linewidth = 0.25) +
  scale_fill_manual(values = c("Sig. concordant" = col_conc,
                               "Sig. discordant" = col_disc,
                               "Non-significant" = col_ns)) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "% Concordant with phenotype direction", fill = NULL) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 12),
        axis.ticks.x = element_line(linewidth = 0.6),
        panel.grid.major = element_line(color = col_grid, linewidth = 0.5),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        legend.position = "top", legend.text = element_text(size = 11),
        plot.margin = unit(c(1, 3, 1, 0), "cm"))

concordance_plot <- p_pval + p_bars + plot_layout(widths = c(0.25, 3))

ggsave(file.path(opt$outdir, "direction_concordance_by_category.png"),
       plot = concordance_plot, width = 12, height = 6, dpi = 300)
ggsave(file.path(opt$outdir, "direction_concordance_by_category.svg"),
       plot = concordance_plot, device = svg_device, width = 12, height = 6)
message("Saved: direction_concordance_by_category.{png,svg}")

# =============================================================================
# 8. CONCORDANCE BAR CHART — by block
# =============================================================================

block_plot <- direction_block %>%
  mutate(
    fill_color = case_when(
      binom_p_adj < 0.05 & pct_concordant > 50 ~ "Sig. concordant",
      binom_p_adj < 0.05 & pct_concordant < 50 ~ "Sig. discordant",
      TRUE ~ "Non-significant"
    ),
    bar_label = paste0(n_concordant, "/", n_snps,
                       " (p=", ifelse(is.na(binom_p_adj), "NA",
                                      formatC(binom_p_adj, format = "e", digits = 1)), ")"),
    block_label = paste0(Block, " (", ifelse(phenotype_direction == "Positive", "►", "◄"), ")"),
    block_label = factor(block_label, levels = block_label[order(pct_concordant)])
  ) %>%
  ggplot(aes(x = block_label, y = pct_concordant, fill = fill_color)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = col_ink) +
  geom_text(aes(label = bar_label, y = pct_concordant + 1), hjust = 0, size = 3) +
  scale_fill_manual(values = c("Sig. concordant" = col_conc,
                               "Sig. discordant" = col_disc,
                               "Non-significant" = col_ns)) +
  scale_y_continuous(limits = c(0, 115)) +
  coord_flip() +
  labs(x = NULL, y = "% Concordant with phenotype direction", fill = NULL,
       title = "Effect direction concordance by block") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "top")

ggsave(file.path(opt$outdir, "direction_concordance_by_block.png"),
       plot = block_plot, width = 10, height = 6, dpi = 300)
message("Saved: direction_concordance_by_block.png")

# =============================================================================
# 9. CONCORDANCE BAR CHARTS — per category (hierarchized trait level)
# =============================================================================

direction_hierarchized <- updated_trait_table %>%
  filter(!is.na(direction_match)) %>%
  group_by(Category, `Hierarchized Trait`) %>%
  summarise(
    n_traits = n_distinct(Trait),
    n_snps = n(),
    n_concordant = sum(direction_match == "Concordant"),
    n_discordant = sum(direction_match == "Discordant"),
    pct_concordant = round(100 * n_concordant / n_snps, 1),
    .groups = "drop"
  ) %>%
  mutate(
    binom_p = binom_p_vec(n_concordant, n_snps),
    binom_p_adj = p.adjust(binom_p, method = "BH"),
    fill_color = case_when(
      binom_p_adj < 0.05 & pct_concordant > 50 ~ "Sig. concordant",
      binom_p_adj < 0.05 & pct_concordant < 50 ~ "Sig. discordant",
      TRUE ~ "Non-significant"
    ),
    bar_label = paste0(n_concordant, "/", n_snps,
                       " (p=", ifelse(is.na(binom_p_adj), "NA",
                                      formatC(binom_p_adj, format = "e", digits = 1)), ")")
  )

write.table(direction_hierarchized,
            file.path(opt$outdir, "direction_summary_by_hierarchized_trait.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

per_cat_dir <- file.path(opt$outdir, "concordance_per_cat")
dir.create(per_cat_dir, recursive = TRUE, showWarnings = FALSE)

for (cat in unique(direction_hierarchized$Category)) {
  plot_data <- direction_hierarchized %>%
    filter(Category == cat) %>%
    mutate(`Hierarchized Trait` = factor(`Hierarchized Trait`,
                                         levels = `Hierarchized Trait`[order(pct_concordant)]))
  if (nrow(plot_data) == 0) next

  p <- ggplot(plot_data, aes(x = `Hierarchized Trait`, y = pct_concordant, fill = fill_color)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 50, linetype = "dashed", color = col_ink, linewidth = 0.25) +
    geom_text(aes(label = bar_label, y = pct_concordant + 1), hjust = 0, size = 3) +
    scale_fill_manual(values = c("Sig. concordant" = col_conc,
                                 "Sig. discordant" = col_disc,
                                 "Non-significant" = col_ns), drop = FALSE) +
    scale_y_continuous(limits = c(0, 115)) +
    coord_flip() +
    labs(x = NULL, y = "% Concordant with phenotype direction", fill = NULL,
         title = paste0("Effect direction concordance — ", stringr::str_to_sentence(cat))) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "top")

  plot_height <- max(3, nrow(plot_data) * 0.4)
  safe_cat <- gsub("[^A-Za-z0-9]", "_", cat)
  png_path <- file.path(per_cat_dir, paste0("direction_concordance_", safe_cat, ".png"))
  svg_path <- file.path(per_cat_dir, paste0("direction_concordance_", safe_cat, ".svg"))

  ggsave(png_path, plot = p, width = 10, height = plot_height, dpi = 300)
  message("Saved: ", png_path)
  ggsave(svg_path, plot = p, device = svg_device, width = 10, height = plot_height)
  message("Saved: ", svg_path)
}

message("All done.")
