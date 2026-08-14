#!/usr/bin/env Rscript
# =============================================================================
# Block x Category pleiotropy heatmap
#
# Tile heatmap of LD blocks (rows) against cross-trait categories (columns),
# each cell coloured by how much that block is shared with that category.
# Reveals pleiotropic hotspot loci and which categories co-localise.
#
# Input  : 09_harmonized_trait_table.tsv (or .rds) — needs Block + Category
#          (and a trait column for --metric traits).
# Output : <outdir>/<outfile>.png   (and .svg with --svg)
#
# Metric (--metric):
#   associations  cell = number of harmonized SNP-trait associations  [default]
#   traits        cell = number of distinct traits (--trait-col)
#
# Usage:
#   Rscript figure_scripts/fig_block_category_heatmap.R \
#       -i results/intermediate/leptin/09_harmonized_trait_table.tsv \
#       -o figures/leptin
#   # CLL has 49 blocks — give it more height:
#   Rscript figure_scripts/fig_block_category_heatmap.R \
#       -i .../chronic_lymphocytic_leukemia/09_harmonized_trait_table.tsv \
#       -o figures/chronic_lymphocytic_leukemia --height 14 --show-counts
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(gtools)   # mixedsort, for natural B1, B2, ... B10 ordering
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Harmonized trait table (09_harmonized_trait_table.tsv or .rds) [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--outfile"), type = "character", default = "block_category_heatmap",
              help = "Output base filename, no extension [default: %default]"),
  make_option(c("--metric"), type = "character", default = "associations",
              help = "Cell value: 'associations' or 'traits' [default: %default]"),
  make_option(c("--trait-col"), type = "character", default = "Hierarchized Trait",
              help = "Trait column counted when --metric traits [default: %default]"),
  make_option(c("--scale"), type = "character", default = "log",
              help = "Fill scale: 'log', 'linear', or 'binary' (presence/absence) [default: %default]"),
  make_option(c("--show-counts"), action = "store_true", default = FALSE,
              help = "Print the cell value inside each tile"),
  make_option(c("--title"), type = "character", default = NULL,
              help = "Optional plot title"),
  make_option(c("--width"), type = "double", default = 11.0,
              help = "Figure width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 9.0,
              help = "Figure height in inches [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300L,
              help = "PNG resolution [default: %default]"),
  make_option(c("--svg"), action = "store_true", default = FALSE,
              help = "Also write an .svg (requires the svglite package)"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Block x Category pleiotropy heatmap")
opt <- parse_args(parser)

if (is.null(opt$input)) {
  cat("ERROR: --input is required\n\n", file = stderr()); print_help(parser); quit(status = 2)
}
if (!file.exists(opt$input)) {
  cat("ERROR: File not found:", opt$input, "\n", file = stderr()); quit(status = 2)
}
if (is.null(opt$outdir)) {
  cat("ERROR: --outdir is required\n\n", file = stderr()); print_help(parser); quit(status = 2)
}
if (!opt$metric %in% c("associations", "traits")) {
  cat("ERROR: --metric must be 'associations' or 'traits'\n", file = stderr()); quit(status = 2)
}
if (!opt$scale %in% c("log", "linear", "binary")) {
  cat("ERROR: --scale must be 'log', 'linear', or 'binary'\n", file = stderr()); quit(status = 2)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

source("R/figure_utils.R")
source("R/palettes.R")

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal         <- get_palette(opt$palette)
col_binary  <- pal_col(pal, "hi",    "#2c3e50")   # filled cell, binary mode
col_cell    <- pal_col(pal, "grid",  "grey85")    # tile separators
col_cell2   <- pal_col(pal, "grid",  "grey88")
grad_lo     <- pal_col(pal, "lo",    "#deebf7")
grad_hi     <- pal_col(pal, "hi",    "#08306b")
txt_on_hi   <- pal_col(pal, "paper", "white")     # label on a dark cell
txt_on_lo   <- pal_col(pal, "ink",   "grey20")    # label on a light cell
col_paper   <- pal_col(pal, "paper", "white")

dat <- read_any(opt$input, as_data_frame = TRUE)
need <- c("Block", "Category")
if (opt$metric == "traits") need <- c(need, opt[["trait-col"]])
missing <- setdiff(need, names(dat))
if (length(missing) > 0) {
  cat("ERROR: input missing column(s): ", paste(missing, collapse = ", "), "\n",
      sep = "", file = stderr())
  quit(status = 2)
}
dat <- dat[!is.na(dat$Block) & nzchar(dat$Block) &
           !is.na(dat$Category) & nzchar(dat$Category), ]
# Sentence-case category names for display
dat$Category <- stringr::str_to_sentence(dat$Category)

# Aggregate to one value per (Block, Category)
if (opt$metric == "associations") {
  agg <- dat %>% count(Block, Category, name = "value")
  fill_lab <- "Associations"
} else {
  tcol <- opt[["trait-col"]]
  agg <- dat %>%
    group_by(Block, Category) %>%
    summarise(value = n_distinct(.data[[tcol]]), .groups = "drop")
  fill_lab <- "Distinct traits"
}

# Natural block order (B1..B10..) on y, busiest categories first on x
block_levels <- mixedsort(unique(agg$Block))
cat_order <- agg %>% group_by(Category) %>%
  summarise(tot = sum(value), .groups = "drop") %>%
  arrange(desc(tot)) %>% pull(Category)

agg$Block    <- factor(agg$Block, levels = rev(block_levels))  # B1 at top
agg$Category <- factor(agg$Category, levels = cat_order)

p <- ggplot(agg, aes(x = Category, y = Block)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 0, vjust = 0, size = 10),
    axis.text.y      = element_text(size = 10),
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.margin      = unit(c(0.5, 1.2, 0.5, 0.5), "cm")
  )

if (opt$scale == "binary") {
  # Presence/absence: every (block, category) with >=1 hit is a filled tile.
  p <- p + geom_tile(fill = col_binary, color = col_cell, linewidth = 0.3)
} else {
  p <- p + geom_tile(aes(fill = value), color = col_cell2, linewidth = 0.3)
  if (opt$scale == "log") {
    p <- p + scale_fill_gradient(low = grad_lo, high = grad_hi,
                                 trans = "log10", name = paste0(fill_lab, "\n(log10)"))
  } else {
    p <- p + scale_fill_gradient(low = grad_lo, high = grad_hi, name = fill_lab)
  }
}

if (opt[["show-counts"]]) {
  txt_col <- if (opt$scale == "binary") txt_on_hi
             else ifelse(agg$value > max(agg$value) * 0.6, txt_on_hi, txt_on_lo)
  p <- p + geom_text(aes(label = value), size = 3, color = txt_col)
}
if (!is.null(opt$title)) p <- p + ggtitle(opt$title)

png_path <- file.path(opt$outdir, paste0(opt$outfile, ".png"))
ggsave(png_path, plot = p, device = "png",
       width = opt$width, height = opt$height, units = "in", dpi = opt$dpi, bg = col_paper)
message("PNG saved -> ", png_path)

if (opt$svg) {
  svg_path <- file.path(opt$outdir, paste0(opt$outfile, ".svg"))
  ggsave(svg_path, plot = p, device = "svg", width = opt$width, height = opt$height,
         fix_text_size = FALSE)
  message("SVG saved -> ", svg_path)
}

message("Heatmap: ", length(block_levels), " blocks x ", length(cat_order),
        " categories (", opt$metric, ", ", opt$scale, " scale).")
