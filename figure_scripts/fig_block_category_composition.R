#!/usr/bin/env Rscript
# ============================================================
# fig_block_category_composition.R
# Purpose:
#   Companion panel to fig_blockstats_barplot.R. Where blockstats shows
#   the *amount* of pleiotropy per LD block (cross-traits, SNPs, genes,
#   size), this shows its *categorical make-up*: a stacked bar per block
#   coloured by cross-trait category.
#
#   Designed to sit side-by-side with the blockstats figure: blocks are
#   placed on x in the same order (descending cross-traits) so the two
#   panels' columns line up. Only blocks that actually carry cross-traits
#   get a bar, and every bar is full -- so unlike the Block x Category
#   heatmap there is no empty whitespace (the CLL "80% empty" problem).
#
#   To keep the legend readable, only the top-N categories are coloured;
#   the rest are lumped into "Other".
#
#   Counts the SAME cross-trait population blockstats does, so the stacked
#   bars *decompose* the blockstats cross-trait bars exactly:
#     - source is 06_crosstraits_filtered (the pre-harmonization cross-traits),
#       NOT 09 (09 drops rows whose alleles can't be harmonized, which is an
#       effect-direction concern unrelated to counting category membership);
#     - default metric is 'traits' (distinct --trait-col), matching the
#       blockstats default. Because each trait maps to exactly one category,
#       summing distinct traits over categories == distinct traits per block
#       == the blockstats cross-trait count for that block.
#   Category labels are attached from the stage-07 trait table (--categories),
#   the same map 09 itself uses.
#
# Input  : --input      06_crosstraits_filtered.{tsv,rds} -> Block, --trait-col
#          --categories  07_trait_table.{tsv,rds}         -> --trait-col, Category
# Output : <outdir>/<outfile>.png  (and .svg with --svg)
#
# Metric (--metric):
#   traits        bar segment = number of distinct traits (--trait-col)  [default]
#   associations  bar segment = number of SNP-trait association rows
#
# Usage:
#   Rscript figure_scripts/fig_block_category_composition.R \
#       -i results/intermediate/chronic_lymphocytic_leukemia/06_crosstraits_filtered.tsv \
#       --categories results/intermediate/chronic_lymphocytic_leukemia/07_trait_table.tsv \
#       -o figures/chronic_lymphocytic_leukemia
#
#   # 100%-stacked composition instead of raw counts, keep top 8 categories:
#   Rscript figure_scripts/fig_block_category_composition.R ... --proportional --top-n 8
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(grid)   # unit()
})

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------
source("R/figure_utils.R")
source("R/palettes.R")

# Distinguishable qualitative palette; "Other" is always grey.
# First 12 are the original Tableau-10/20 set (kept so existing top-10 figures
# are unchanged); the next 8 extend it to 20 distinct hues for --top-n 0.
cat_palette <- c(
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2",
  "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
  "#393B79", "#843C39", "#BCBD22", "#17BECF", "#E377C2", "#637939",
  "#8C6D31", "#7B4173"
)

# ----------------------------------------------------------------
# CLI
# ----------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Cross-traits table (06_crosstraits_filtered.tsv or .rds) -> Block, --trait-col [required]"),
  make_option(c("--categories"), type = "character", default = NULL,
              help = "Trait->category map (07_trait_table.tsv or .rds) -> --trait-col, --category-col [required]"),
  make_option(c("--blocks"), type = "character", default = NULL,
              help = "04_blocks file (.tsv/.rds); show ALL its blocks as columns, incl. zero-cross-trait ones as empty columns (aligns with the blockstats panel)"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option(c("--outfile"), type = "character", default = "block_category_composition.png",
              help = "Output figure filename [default: %default]"),
  make_option(c("--metric"), type = "character", default = "traits",
              help = "Bar segment value: 'traits' (distinct --trait-col) or 'associations' (rows) [default: %default]"),
  make_option(c("--trait-col"), type = "character", default = "MAPPED_TRAIT",
              help = "Trait column: join key to --categories, and unit counted when --metric traits [default: %default]"),
  make_option(c("--category-col"), type = "character", default = "Category",
              help = "Column in --categories holding the category names [default: %default]"),
  make_option(c("--exclude"), type = "character", default = NULL,
              help = "Comma-separated category names to exclude, applied before --top-n (case-insensitive)"),
  make_option(c("--top-n"), type = "integer", default = 10L,
              help = "Colour the N most frequent categories; lump the rest into 'Other'. Use 0 to keep ALL categories (no 'Other' bucket) [default: %default]"),
  make_option(c("--proportional"), action = "store_true", default = FALSE,
              help = "100%-stacked bars (per-block composition) instead of raw counts"),
  make_option(c("--title"), type = "character", default = NULL,
              help = "Optional plot title"),
  make_option(c("--width"), type = "double", default = 8.0,
              help = "Figure width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 4.5,
              help = "Figure height in inches [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300L,
              help = "PNG DPI [default: %default]"),
  make_option(c("--svg"), action = "store_true", default = FALSE,
              help = "Also write an .svg (requires the svglite package)"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Per-block cross-trait category composition (companion to blockstats)")
opt <- parse_args(parser)

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal         <- get_palette(opt$palette)
cat_palette <- pal_cat(pal, length(cat_palette), cat_palette)
col_other   <- pal_col(pal, "neutral", "grey75")   # the "Other" catch-all
col_paper   <- pal_col(pal, "paper",   "white")

if (is.null(opt$input)) {
  cat("ERROR: --input is required\n\n", file = stderr()); print_help(parser); quit(status = 2)
}
if (!file.exists(opt$input)) {
  cat("ERROR: File not found:", opt$input, "\n", file = stderr()); quit(status = 2)
}
if (is.null(opt$categories)) {
  cat("ERROR: --categories is required (the 07_trait_table that maps traits to categories)\n\n",
      file = stderr()); print_help(parser); quit(status = 2)
}
if (!file.exists(opt$categories)) {
  cat("ERROR: File not found:", opt$categories, "\n", file = stderr()); quit(status = 2)
}
if (!opt$metric %in% c("associations", "traits")) {
  cat("ERROR: --metric must be 'associations' or 'traits'\n", file = stderr()); quit(status = 2)
}
if (opt[["top-n"]] < 0) {
  cat("ERROR: --top-n must be >= 0 (0 = keep all categories)\n", file = stderr()); quit(status = 2)
}
safe_dir_create(opt$outdir)

cat_col   <- opt[["category-col"]]
trait_col <- opt[["trait-col"]]

# ----------------------------------------------------------------
# Read & validate: cross-traits (06) + trait->category map (07)
# ----------------------------------------------------------------
dat <- read_any(opt$input, as_data_frame = TRUE)
miss_in <- setdiff(c("Block", trait_col), names(dat))
if (length(miss_in) > 0) {
  cat("ERROR: --input missing column(s): ", paste(miss_in, collapse = ", "), "\n",
      sep = "", file = stderr())
  quit(status = 2)
}
dat <- dat[!is.na(dat$Block) & nzchar(dat$Block) &
           !is.na(dat[[trait_col]]) & nzchar(dat[[trait_col]]), ]
if (nrow(dat) == 0) {
  cat("ERROR: no --input rows with both a Block and a ", trait_col, "\n",
      sep = "", file = stderr()); quit(status = 2)
}

cats <- read_any(opt$categories, as_data_frame = TRUE)
miss_cat <- setdiff(c(trait_col, cat_col), names(cats))
if (length(miss_cat) > 0) {
  cat("ERROR: --categories missing column(s): ", paste(miss_cat, collapse = ", "), "\n",
      sep = "", file = stderr())
  quit(status = 2)
}
# One category per trait (verified for this pipeline). Guard against a trait that
# maps to >1 category, which would multiply association rows and break the
# decomposition; keep distinct (trait, category) pairs and warn if any conflict.
cmap <- cats[!is.na(cats[[trait_col]]) & nzchar(cats[[trait_col]]) &
             !is.na(cats[[cat_col]])   & nzchar(cats[[cat_col]]),
             c(trait_col, cat_col)]
cmap <- cmap[!duplicated(cmap), ]
multi <- names(which(table(cmap[[trait_col]]) > 1))
if (length(multi) > 0)
  message("WARNING: ", length(multi), " trait(s) map to >1 category; bars will no longer ",
          "sum exactly to the blockstats cross-trait count. e.g. ", multi[1])

# Attach categories to the cross-traits (same map 09 uses)
#
# Both keys are trimmed first, for the reason documented at the join in
# pipeline/09_harmonize.R: the Catalog's MAPPED_TRAIT carries a trailing space
# on some traits while 07_trait_table holds the trimmed form, so an exact-string
# join drops them — 13 lipid rows in CLL, 60 traits in leptin — and they arrive
# here as "no category" rather than as an error. No trait exists in both forms,
# so this merges nothing and changes no count.
dat[[trait_col]]  <- trimws(dat[[trait_col]])
cmap[[trait_col]] <- trimws(cmap[[trait_col]])
cmap <- cmap[!duplicated(cmap), ]

n_before <- nrow(dat)
dat <- merge(dat, cmap, by = trait_col, all.x = TRUE)
orphan <- is.na(dat[[cat_col]]) | !nzchar(dat[[cat_col]])
if (any(orphan)) {
  message("NOTE: ", sum(orphan), " of ", n_before, " cross-trait rows have no category in ",
          "--categories and are dropped (", length(unique(dat[[trait_col]][orphan])), " trait(s)).")
  dat <- dat[!orphan, ]
}
if (nrow(dat) == 0) {
  cat("ERROR: no cross-trait rows could be categorised\n", file = stderr()); quit(status = 2)
}

# Sentence-case category names for display
dat[[cat_col]] <- stringr::str_to_sentence(dat[[cat_col]])

# Block order reference: descending cross-traits (distinct trait_col) per block,
# captured BEFORE any --exclude so the x-axis still matches the blockstats panel.
block_ct <- dat %>% group_by(Block) %>%
  summarise(ct = n_distinct(.data[[trait_col]]), .groups = "drop") %>%
  arrange(desc(ct)) %>% pull(Block)

# Drop excluded categories (matched case-insensitively), before the top-N cut
if (!is.null(opt$exclude)) {
  drop <- trimws(strsplit(opt$exclude, ",")[[1]])
  drop <- drop[nzchar(drop)]
  unmatched <- drop[!tolower(drop) %in% tolower(dat[[cat_col]])]
  if (length(unmatched))
    message("WARNING: --exclude value(s) not found: ", paste(unmatched, collapse = ", "))
  dat <- dat[!tolower(dat[[cat_col]]) %in% tolower(drop), ]
  if (nrow(dat) == 0) {
    cat("ERROR: every category was excluded; nothing to plot\n", file = stderr()); quit(status = 2)
  }
}

# ----------------------------------------------------------------
# Aggregate to one value per (Block, Category)
# ----------------------------------------------------------------
if (opt$metric == "associations") {
  agg <- dat %>% count(Block, !!sym(cat_col), name = "value")
  y_lab <- "Associations"
} else {
  agg <- dat %>%
    group_by(Block, !!sym(cat_col)) %>%
    summarise(value = n_distinct(.data[[trait_col]]), .groups = "drop")
  y_lab <- "Cross-traits"
}
names(agg)[names(agg) == cat_col] <- "Category"

# ----------------------------------------------------------------
# Keep the top-N categories, lump the rest into "Other"
# ----------------------------------------------------------------
cat_totals <- agg %>% group_by(Category) %>%
  summarise(tot = sum(value), .groups = "drop") %>%
  arrange(desc(tot))
# --top-n 0 means "keep every category" (no 'Other' bucket)
n_keep    <- if (opt[["top-n"]] == 0) nrow(cat_totals) else min(opt[["top-n"]], nrow(cat_totals))
top_cats  <- cat_totals$Category[seq_len(n_keep)]
n_lumped  <- nrow(cat_totals) - n_keep
if (n_keep > length(cat_palette))
  message("NOTE: ", n_keep, " categories exceed the ", length(cat_palette),
          "-colour palette; some fills will repeat. Consider --top-n or --exclude.")

agg$Category <- ifelse(agg$Category %in% top_cats, agg$Category, "Other")
# Re-aggregate so each block has a single "Other" segment
agg <- agg %>% group_by(Block, Category) %>%
  summarise(value = sum(value), .groups = "drop")

# Stack order: biggest category at the bottom, "Other" always last (top, grey)
lvl <- c(top_cats, if (n_lumped > 0 || any(agg$Category == "Other")) "Other")
lvl <- lvl[lvl %in% unique(agg$Category)]
agg$Category <- factor(agg$Category, levels = lvl)

# ----------------------------------------------------------------
# Block order: descending cross-traits, matching the blockstats panel.
# block_ct was computed from the full (pre-exclude) cross-trait set above;
# keep only blocks that still have a bar, any stragglers last.
# ----------------------------------------------------------------
block_order <- block_ct[block_ct %in% agg$Block]
block_order <- c(block_order, setdiff(unique(agg$Block), block_order))

# With --blocks, append every block from 04_blocks that has no bar here
# (zero cross-traits, or all its categories excluded) as an empty column, so
# the x-axis covers the full block set like blockstats does.
if (!is.null(opt$blocks)) {
  if (!file.exists(opt$blocks)) {
    cat("ERROR: --blocks file not found:", opt$blocks, "\n", file = stderr()); quit(status = 2)
  }
  allb <- read_any(opt$blocks, as_data_frame = TRUE)
  if (!"Block" %in% names(allb)) {
    cat("ERROR: --blocks file has no 'Block' column\n", file = stderr()); quit(status = 2)
  }
  full <- unique(allb$Block[!is.na(allb$Block) & nzchar(allb$Block)])
  empties <- setdiff(full, block_order)
  block_order <- c(block_order, empties)
  message("Including ", length(empties), " block(s) with no cross-traits as empty columns.")
}
agg$Block <- factor(agg$Block, levels = block_order)

# ----------------------------------------------------------------
# Plot
# ----------------------------------------------------------------
fill_vals <- setNames(cat_palette[((seq_along(lvl) - 1) %% length(cat_palette)) + 1], lvl)
if ("Other" %in% lvl) fill_vals[["Other"]] <- col_other

pos <- if (opt$proportional) "fill" else "stack"

p <- ggplot(agg, aes(x = Block, y = value, fill = Category)) +
  geom_col(position = pos, width = 0.8, color = col_paper, linewidth = 0.15) +
  scale_x_discrete(drop = FALSE) +     # keep empty (zero-cross-trait) blocks as columns
  scale_fill_manual(values = fill_vals, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
    axis.text.y      = element_text(size = 15),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 12),
    legend.title     = element_blank(),
    legend.key.size  = unit(0.9, "lines")
  ) +
  labs(x = NULL,
       y = if (opt$proportional) paste0(y_lab, " (proportion)") else y_lab)

if (opt$proportional) p <- p + scale_y_continuous(labels = scales::percent_format())
if (!is.null(opt$title)) p <- p + ggtitle(opt$title) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

# ----------------------------------------------------------------
# Save
# ----------------------------------------------------------------
png_path <- file.path(opt$outdir, opt$outfile)
ggsave(png_path, plot = p, device = "png",
       width = opt$width, height = opt$height, units = "in", dpi = opt$dpi, bg = col_paper)
message("PNG saved -> ", png_path)

if (opt$svg) {
  svg_path <- sub("\\.png$", ".svg", png_path)
  if (svg_path == png_path) svg_path <- paste0(png_path, ".svg")
  ggsave(svg_path, plot = p, device = "svg", width = opt$width, height = opt$height,
         fix_text_size = FALSE)
  message("SVG saved -> ", svg_path)
}

message(sprintf("Composition: %d blocks x %d categories shown (%s, %s)%s.",
                length(block_order), length(lvl), opt$metric,
                if (opt$proportional) "proportional" else "counts",
                if (n_lumped > 0) sprintf("; %d rare categories lumped into 'Other'", n_lumped) else ""))
