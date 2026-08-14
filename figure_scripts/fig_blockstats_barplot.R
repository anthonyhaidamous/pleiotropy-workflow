#!/usr/bin/env Rscript
# ============================================================
# fig_blockstats_barplot.R
# Purpose:
#   Figure showing the pleiotropy distribution among LD blocks.
#   Cross-traits per block are drawn as bars; block span, SNP count and
#   overlapping-gene count are overlaid as dashed lines. All metrics are
#   log2(count + 1) transformed so they share one axis.
#
# Derives the per-block statistics directly from the pipeline outputs
# (no separate "blockstats.txt" step needed):
#   --blocks       04_blocks.{tsv,rds}  -> Block, span (kb), nb lSNPs,
#                                          nb pSNPs, overlapping-gene count
#   --crosstraits  06_crosstraits_filtered.{tsv,rds} -> cross-traits per block
#
# Outputs:
#   <outdir>/blockstats_plot.png
#
# Usage:
#   Rscript figure_scripts/fig_blockstats_barplot.R \
#       --blocks      results/intermediate/leptin/04_blocks.rds \
#       --crosstraits results/intermediate/leptin/06_crosstraits_filtered.rds \
#       --outdir      figures/leptin
#
#   # count cross-trait associations instead of distinct traits:
#   Rscript figure_scripts/fig_blockstats_barplot.R ... --crosstrait-metric associations
# ============================================================

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))  # unit()

# -----------------------------
# Helpers
# -----------------------------
source("R/figure_utils.R")
source("R/palettes.R")

# Count distinct gene symbols in one 04_blocks "Genes" cell.
# Handles .rds (list column -> character vector) and .tsv (comma string).
# Splits on commas, slashes, and " - " flanking-gene separators.
count_genes <- function(cell) {
  if (is.list(cell)) cell <- cell[[1]]
  toks <- unlist(strsplit(as.character(cell), "\\s*[,/]\\s*|\\s+-\\s+"))
  toks <- trimws(toks)
  toks <- toks[toks != "" & toks != "." & !is.na(toks)]
  length(unique(toks))
}

# -----------------------------
# CLI
# -----------------------------
option_list <- list(
  make_option(c("-b", "--blocks"), type = "character",
              help = "04_blocks file (.tsv or .rds) [required]"),
  make_option(c("-c", "--crosstraits"), type = "character",
              help = "06_crosstraits_filtered file (.tsv or .rds) [required]"),
  make_option(c("--crosstrait-metric"), type = "character", default = "traits",
              help = "Count 'traits' (distinct MAPPED_TRAIT) or 'associations' (rows) per block [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option(c("--outfile"), type = "character", default = "blockstats_plot.png",
              help = "Output figure filename [default: %default]"),
  make_option(c("--width"), type = "double", default = 8,
              help = "Figure width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 4.5,
              help = "Figure height in inches [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "PNG DPI [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Per-block pleiotropy distribution figure")
opt <- parse_args(parser)

if (is.null(opt$blocks)) {
  cat("ERROR: --blocks is required\n\n", file = stderr())
  print_help(parser); quit(status = 2)
}
if (is.null(opt$crosstraits)) {
  cat("ERROR: --crosstraits is required (cross-trait counts are not in 04_blocks;\n",
      "       point this at 06_crosstraits_filtered)\n\n", sep = "", file = stderr())
  print_help(parser); quit(status = 2)
}
for (f in c(opt$blocks, opt$crosstraits)) {
  if (!file.exists(f)) { cat("ERROR: File not found:", f, "\n", file = stderr()); quit(status = 2) }
}
ct_metric <- opt[["crosstrait-metric"]]
if (!ct_metric %in% c("traits", "associations")) {
  cat("ERROR: --crosstrait-metric must be 'traits' or 'associations'\n", file = stderr())
  quit(status = 2)
}

safe_dir_create(opt$outdir)

message("blocks      = ", opt$blocks)
message("crosstraits = ", opt$crosstraits, " (metric: ", ct_metric, ")")
message("output      = ", file.path(opt$outdir, opt$outfile))

# -----------------------------
# Derive per-block statistics from 04_blocks + 06_crosstraits
# -----------------------------
blocks <- read_any(opt$blocks)
ct     <- read_any(opt$crosstraits)

req <- c("Block", "Genomic_span", "nb_lSNPs", "nb_pSNPs", "Genes")
miss <- setdiff(req, names(blocks))
if (length(miss) > 0) {
  cat("ERROR: 04_blocks missing column(s): ", paste(miss, collapse = ", "), "\n",
      sep = "", file = stderr())
  quit(status = 2)
}

# Cross-traits per block (distinct traits or associations)
ct_per_block <- ct %>%
  group_by(Block) %>%
  summarise(`Cross-traits` = if (ct_metric == "traits") n_distinct(MAPPED_TRAIT) else n(),
            .groups = "drop")

blockstats <- blocks %>%
  transmute(
    Block,
    Size                = as.numeric(Genomic_span) / 1000,          # bp -> kb
    lSNPs               = as.numeric(nb_lSNPs),
    pSNPs               = as.numeric(nb_pSNPs),
    `Overlapping genes` = vapply(Genes, count_genes, integer(1)),
    SNPs                = as.numeric(nb_lSNPs) + as.numeric(nb_pSNPs)
  ) %>%
  left_join(ct_per_block, by = "Block") %>%
  mutate(`Cross-traits` = tidyr::replace_na(`Cross-traits`, 0))

# Order blocks by descending cross-traits (NA last)
block_order <- blockstats %>%
  arrange(desc(replace_na(`Cross-traits`, -Inf))) %>%
  pull(Block)

# Reshape & log2-transform (log2(x + 1))
df_long <- blockstats %>%
  select(Block, Size, SNPs, `Cross-traits`, `Overlapping genes`) %>%
  mutate(across(-Block, ~ log2(.x + 1))) %>%
  pivot_longer(-Block, names_to = "Metric", values_to = "Log2Value") %>%
  mutate(
    Block  = factor(Block, levels = block_order),
    Metric = factor(Metric, levels = c("SNPs", "Cross-traits", "Overlapping genes", "Size"))
  )

# -----------------------------
# Plot
# -----------------------------
use_chicklet <- requireNamespace("ggchicklet", quietly = TRUE)

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal      <- get_palette(opt$palette)
col_bar  <- pal_col(pal, "accent", "steelblue")
# The three overlaid metric lines read as a categorical set.
line_cols <- pal_cat(pal, 3, c("green4", "red2", "black"))
names(line_cols) <- c("SNPs", "Overlapping genes", "Size")

p <- ggplot()

if (use_chicklet) {
  p <- p +
    ggchicklet::geom_chicklet(
      data   = df_long %>% filter(Metric == "Cross-traits"),
      aes(x = Block, y = Log2Value),
      fill   = col_bar, width = 0.2, radius = grid::unit(3, "pt")
    )
} else {
  message("NOTE: Package 'ggchicklet' not installed; using standard geom_col (no rounded corners).")
  p <- p +
    geom_col(
      data  = df_long %>% filter(Metric == "Cross-traits"),
      aes(x = Block, y = Log2Value),
      fill  = col_bar, width = 0.2
    )
}

p <- p +
  geom_line(
    data = df_long %>% filter(Metric != "Cross-traits"),
    aes(x = Block, y = Log2Value, color = Metric, linetype = Metric,
        linewidth = Metric, group = Metric)
  ) +
  scale_color_manual(values = line_cols) +
  scale_linetype_manual(values = c(
    SNPs = "twodash", `Overlapping genes` = "twodash", Size = "twodash"
  )) +
  scale_linewidth_manual(values = c(
    SNPs = 0.6, `Overlapping genes` = 0.6, Size = 0.6
  )) +
  guides(linewidth = "none") +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
    axis.text.y      = element_text(size = 15),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.key.size  = unit(1.2, "lines"),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 0),
    legend.spacing.x = unit(0.5, "cm")
  ) +
  labs(x = NULL, y = "log₂(count + 1)", color = NULL, linetype = NULL)

# -----------------------------
# Save
# -----------------------------
out_path <- file.path(opt$outdir, opt$outfile)
ggsave(filename = out_path, plot = p, device = "png",
       width = opt$width, height = opt$height, units = "in", dpi = opt$dpi)
message("Figure written to: ", out_path)
