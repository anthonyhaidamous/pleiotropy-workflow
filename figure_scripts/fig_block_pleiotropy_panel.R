#!/usr/bin/env Rscript
# ============================================================
# fig_block_pleiotropy_panel.R
# The two per-block figures as one panel, on one block axis.
#
#   above   the pleiotropy distribution: cross-traits per block as bars, with
#           block span, SNP count and overlapping-gene count as dashed lines on
#           the same log2 axis (fig_blockstats_barplot.R)
#   below   what those cross-traits ARE: the same bar decomposed by trait
#           category (fig_block_category_composition.R)
#   between a rule, and one set of block labels under the lower panel
#
# WHY THE TWO BELONG TOGETHER
# ---------------------------
# Separately they answer half a question each. The upper panel says the
# pleiotropy is concentrated — B17, B23 and B10 carry a third of it — and that
# the concentration is not explained by the obvious confounders: the SNP,
# gene-count and span lines stay flat while the bars decay. The lower panel says
# what the concentrated loci are concentrated ON. The claim only lands with
# both, and it lands hardest when the same column carries both, which is what
# this script is for.
#
# ONE ORDER, ENFORCED
# -------------------
# The two sibling scripts each sort blocks by their own cross-trait count, and
# they disagree: blockstats counts every cross-trait, the composition figure can
# only count the ones that have a category, so a block can sit in a different
# place in each. Stacked, that would put two different blocks in one column and
# invite exactly the wrong reading.
#
# So the order is computed ONCE here, from the full cross-trait count (the
# blockstats rule), and imposed on both panels with drop = FALSE. Ties are
# broken by natural block order rather than by row order, so the figure is the
# same on every run. Because the upper bar counts uncategorised traits and the
# lower stack cannot, the two totals can differ by a trait or two; the count is
# reported on the console rather than left to be noticed.
#
# Inputs:
#   --blocks       04_blocks.{tsv,rds}                 span, SNP and gene counts
#   --crosstraits  06_crosstraits_filtered.{tsv,rds}   cross-traits per block
#   --categories   07_trait_table.{tsv,rds}            trait -> category
#
# Outputs (in --outdir):
#   <outfile>.png, and .svg with --svg
#
# Usage:
#   Rscript figure_scripts/fig_block_pleiotropy_panel.R \
#       --blocks      results/intermediate/chronic_lymphocytic_leukemia/04_blocks.rds \
#       --crosstraits results/intermediate/chronic_lymphocytic_leukemia/06_crosstraits_filtered.tsv \
#       --categories  results/intermediate/chronic_lymphocytic_leukemia/07_trait_table.tsv \
#       --outdir      figures/chronic_lymphocytic_leukemia --palette poster
#
#   --top-n N        colour the N most frequent categories, lump the rest into
#                    "Other". DEFAULT 0 — every category keeps its own colour,
#                    which is the opposite of the sibling script's default of 10.
#                    A panel this wide has room for the full legend.
#   --all-blocks     keep blocks with no cross-trait at all as empty columns
#                    (default: only blocks with at least one)
#   --no-lines       drop the SNP / gene / span lines from the upper panel,
#                    leaving the bare distribution
#   --heights A,B    relative panel heights [default 1,1.25]
#   --crosstrait-metric traits|associations   what a bar counts [traits]
#   --palette NM     article | poster | ppt, or a path to a YAML file. Sets the
#                    paper, ink, grid, bar and line colours — but NOT the
#                    category fills, which have their own flag:
#   --cat-palette C  contrast (default) | tableau | palette, or a comma-separated
#                    list of hex colours. See the note at CAT_RAMPS below for why
#                    the category fills do not follow --palette.
#   --base-size N    base type size in POINTS at the requested --width; every
#                    other size is relative to it. This is the knob for making
#                    the figure legible at a given PHYSICAL size: set --width to
#                    what the poster will print and raise --base-size until the
#                    type holds. Rendering wide and letting the layout program
#                    scale the file down is what turns 11 pt into 6.5 pt.
#   --legend-cols N  columns in the category legend; 0 (default) fits them to
#                    --width and --base-size
#   --width, --height, --dpi, --svg, --outfile
#
# For a 9.51-inch poster column:
#   … --width 9.51 --height 7 --base-size 10 --dpi 600 --svg
# ============================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(grid); library(gtools)
})

if (!file.exists("R/figure_utils.R")) stop("run from the repository root")
source("R/figure_utils.R")
source("R/palettes.R")

# ------------------------------------------------------------------
# Categorical ramps for the lower panel (--cat-palette)
#
# The category fills are deliberately NOT taken from --palette by default, and
# the reason is structural rather than a matter of taste. The `poster` palette's
# `cat` is twelve colours walking a hue circle at constant lightness: elegant
# for six or eight categories, and unusable for seventeen, because asking for
# more than twelve interpolates BETWEEN neighbouring hues and the four largest
# categories — which are adjacent in every stack, since the stack is ordered by
# size — come out as four shades of the same dusty red.
#
# `contrast` is the default instead. It is built from Paul Tol's qualitative
# schemes (colour-blind safe by construction) and ordered so that CONSECUTIVE
# entries differ in hue and in lightness, which is the property that matters
# here: consecutive entries are exactly the segments that touch each other in a
# stacked bar. Hue does get reused past a dozen categories — seventeen genuinely
# distinct colours do not exist — but the reuse is always far apart in the
# order, never between neighbours.
#
# --cat-palette palette restores the old behaviour for anyone who needs the
# figure to match the poster's own colour identity exactly.
# ------------------------------------------------------------------
CAT_RAMPS <- list(
  contrast = c(
    "#4477AA", "#EE7733", "#117733", "#EE3377", "#88CCEE", "#882255",
    "#DDCC77", "#009988", "#AA4499", "#999933", "#EE6677", "#332288",
    "#994F00", "#CCBB44", "#66CCEE", "#44AA99", "#771155", "#BBBBBB",
    "#0077BB", "#225522"
  ),
  # The fallback ramp of fig_block_category_composition.R, so the lower panel
  # can be made to match that figure's no-palette rendering exactly.
  tableau = c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2",
    "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
    "#393B79", "#843C39", "#BCBD22", "#17BECF", "#E377C2", "#637939",
    "#8C6D31", "#7B4173"
  )
)

# Count distinct gene symbols in one 04_blocks "Genes" cell. Kept identical to
# fig_blockstats_barplot.R: same separators, same de-duplication, so the gene
# line here is the gene line there.
count_genes <- function(cell) {
  if (is.list(cell)) cell <- cell[[1]]
  toks <- unlist(strsplit(as.character(cell), "\\s*[,/]\\s*|\\s+-\\s+"))
  toks <- trimws(toks)
  toks <- toks[toks != "" & toks != "." & !is.na(toks)]
  length(unique(toks))
}

option_list <- list(
  make_option(c("-b", "--blocks"),      type = "character", help = "04_blocks file [required]"),
  make_option(c("-c", "--crosstraits"), type = "character", help = "06_crosstraits_filtered file [required]"),
  make_option(c("--categories"),        type = "character", help = "07_trait_table file [required]"),
  make_option(c("-o", "--outdir"),      type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option(c("--outfile"),           type = "character", default = "block_pleiotropy_panel",
              help = "Output basename, no extension [default: %default]"),
  make_option(c("--trait-col"),         type = "character", default = "MAPPED_TRAIT"),
  make_option(c("--category-col"),      type = "character", default = "Category"),
  make_option(c("--crosstrait-metric"), type = "character", default = "traits",
              help = "'traits' (distinct) or 'associations' (rows) [default: %default]"),
  make_option(c("--top-n"),             type = "integer",   default = 0L,
              help = "Colour the N most frequent categories, 0 = all [default: %default]"),
  make_option(c("--all-blocks"),        action = "store_true", default = FALSE,
              help = "Keep blocks with zero cross-traits as empty columns"),
  make_option(c("--no-lines"),          action = "store_true", default = FALSE,
              help = "Drop the SNP / gene / span lines from the upper panel"),
  make_option(c("--heights"),           type = "character", default = "1,1.25",
              help = "Relative heights, upper,lower [default: %default]"),
  make_option(c("--bar-width"),         type = "double",    default = 0.62,
              help = "Bar width, 0-1 [default: %default]"),
  make_option(c("--width"),             type = "double",    default = 16),
  make_option(c("--height"),            type = "double",    default = 9),
  make_option(c("--dpi"),               type = "integer",   default = 300),
  make_option(c("--base-size"),         type = "double",    default = 13,
              help = paste("Base type size in POINTS at the requested --width.",
                           "Everything else is relative to it, so this is the one",
                           "number to raise when a figure has to stay legible at a",
                           "smaller physical size [default: %default]")),
  make_option(c("--legend-cols"),       type = "integer",   default = 0L,
              help = "Columns in the category legend; 0 = fit to --width [default: %default]"),
  make_option(c("--svg"),               action = "store_true", default = FALSE),
  make_option(c("--palette"),           type = "character", default = NULL,
              help = "article | poster | ppt, or a path to a YAML file"),
  make_option(c("--cat-palette"),       type = "character", default = "contrast",
              help = paste("Fills for the lower panel: contrast | tableau | palette,",
                           "or a comma-separated list of hex colours [default: %default]"))
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Per-block pleiotropy distribution over its category composition"))

for (a in c("blocks", "crosstraits", "categories")) {
  if (is.null(opt[[a]])) { cat("ERROR: --", a, " is required\n", sep = "", file = stderr()); quit(status = 2) }
  if (!file.exists(opt[[a]])) { cat("ERROR: file not found: ", opt[[a]], "\n", sep = "", file = stderr()); quit(status = 2) }
}
trait_col <- opt[["trait-col"]]
cat_col   <- opt[["category-col"]]
ct_metric <- opt[["crosstrait-metric"]]
if (!ct_metric %in% c("traits", "associations")) {
  cat("ERROR: --crosstrait-metric must be 'traits' or 'associations'\n", file = stderr()); quit(status = 2)
}
if (opt[["top-n"]] < 0) {
  cat("ERROR: --top-n must be >= 0 (0 = every category)\n", file = stderr()); quit(status = 2)
}
heights <- suppressWarnings(as.numeric(trimws(strsplit(opt$heights, ",")[[1]])))
if (length(heights) != 2 || anyNA(heights) || any(heights <= 0)) {
  cat("ERROR: --heights takes two positive numbers, e.g. 1,1.25\n", file = stderr()); quit(status = 2)
}
safe_dir_create(opt$outdir)

blocks <- read_any(opt$blocks)
ct     <- read_any(opt$crosstraits, as_data_frame = TRUE)
cats   <- read_any(opt$categories,  as_data_frame = TRUE)

req <- c("Block", "Genomic_span", "nb_lSNPs", "nb_pSNPs", "Genes")
miss <- setdiff(req, names(blocks))
if (length(miss)) {
  cat("ERROR: --blocks missing column(s): ", paste(miss, collapse = ", "), "\n", sep = "", file = stderr())
  quit(status = 2)
}
for (nm in c(trait_col)) if (!nm %in% names(ct)) {
  cat("ERROR: --crosstraits missing column: ", nm, "\n", sep = "", file = stderr()); quit(status = 2)
}
for (nm in c(trait_col, cat_col)) if (!nm %in% names(cats)) {
  cat("ERROR: --categories missing column: ", nm, "\n", sep = "", file = stderr()); quit(status = 2)
}

ct <- ct[!is.na(ct$Block) & nzchar(ct$Block) &
         !is.na(ct[[trait_col]]) & nzchar(ct[[trait_col]]), ]

# ------------------------------------------------------------------
# The block axis, decided once
# ------------------------------------------------------------------
count_ct <- function(d) {
  d %>% group_by(Block) %>%
    summarise(n = if (ct_metric == "traits") n_distinct(.data[[trait_col]]) else n(),
              .groups = "drop")
}
per_block <- count_ct(ct)

all_blocks <- unique(blocks$Block[!is.na(blocks$Block) & nzchar(blocks$Block)])
if (opt[["all-blocks"]]) {
  per_block <- bind_rows(per_block,
                         tibble(Block = setdiff(all_blocks, per_block$Block), n = 0L))
}
# Natural tiebreak: without it two blocks with the same count swap places
# between runs depending on row order, and the panel stops being comparable
# with itself.
per_block <- per_block[order(-per_block$n, gtools::mixedorder(per_block$Block)), ]
block_order <- per_block$Block
omitted <- setdiff(all_blocks, block_order)
if (length(omitted)) {
  message("blocks with no cross-trait, not drawn (--all-blocks keeps them): ",
          length(omitted), " (", paste(gtools::mixedsort(omitted), collapse = " "), ")")
}

# ------------------------------------------------------------------
# Upper panel: the distribution, and what does not explain it
# ------------------------------------------------------------------
blockstats <- blocks %>%
  transmute(
    Block,
    Size                = as.numeric(Genomic_span) / 1000,
    `Overlapping genes` = vapply(Genes, count_genes, integer(1)),
    SNPs                = as.numeric(nb_lSNPs) + as.numeric(nb_pSNPs)
  ) %>%
  left_join(per_block %>% rename(`Cross-traits` = n), by = "Block") %>%
  mutate(`Cross-traits` = tidyr::replace_na(`Cross-traits`, 0)) %>%
  filter(Block %in% block_order)

df_long <- blockstats %>%
  select(Block, Size, SNPs, `Cross-traits`, `Overlapping genes`) %>%
  mutate(across(-Block, ~ log2(.x + 1))) %>%
  pivot_longer(-Block, names_to = "Metric", values_to = "Log2Value") %>%
  mutate(Block  = factor(Block, levels = block_order),
         Metric = factor(Metric, levels = c("SNPs", "Cross-traits", "Overlapping genes", "Size")))

# ------------------------------------------------------------------
# Lower panel: the same bar, decomposed
#
# Keys are trimmed on both sides for the reason documented at the join in
# pipeline/09_harmonize.R — the Catalog leaves a trailing space on some traits
# while the trait table holds the trimmed form.
# ------------------------------------------------------------------
cmap <- cats[!is.na(cats[[trait_col]]) & nzchar(cats[[trait_col]]) &
             !is.na(cats[[cat_col]])   & nzchar(cats[[cat_col]]), c(trait_col, cat_col)]
cmap[[trait_col]] <- trimws(cmap[[trait_col]])
cmap <- cmap[!duplicated(cmap), ]
multi <- names(which(table(cmap[[trait_col]]) > 1))
if (length(multi))
  message("WARNING: ", length(multi), " trait(s) map to >1 category; the lower stack ",
          "will exceed the upper bar for their blocks. e.g. ", multi[1])

dat <- ct
dat[[trait_col]] <- trimws(dat[[trait_col]])
n_before <- nrow(dat)
dat <- merge(dat, cmap, by = trait_col, all.x = TRUE)
orphan <- is.na(dat[[cat_col]]) | !nzchar(dat[[cat_col]])
if (any(orphan)) {
  message("NOTE: ", sum(orphan), " of ", n_before, " cross-trait rows have no category (",
          length(unique(dat[[trait_col]][orphan])), " trait(s)): ",
          paste(sort(unique(dat[[trait_col]][orphan])), collapse = " | "))
  message("      they are counted in the upper bar and absent from the lower stack.")
  dat <- dat[!orphan, ]
}
dat[[cat_col]] <- stringr::str_to_sentence(dat[[cat_col]])

agg <- if (ct_metric == "associations") {
  dat %>% count(Block, .data[[cat_col]], name = "value")
} else {
  dat %>% group_by(Block, .data[[cat_col]]) %>%
    summarise(value = n_distinct(.data[[trait_col]]), .groups = "drop")
}
names(agg)[names(agg) == cat_col] <- "Category"
agg <- agg[agg$Block %in% block_order, ]

cat_totals <- agg %>% group_by(Category) %>% summarise(tot = sum(value), .groups = "drop") %>%
  arrange(desc(tot))
n_keep   <- if (opt[["top-n"]] == 0) nrow(cat_totals) else min(opt[["top-n"]], nrow(cat_totals))
top_cats <- cat_totals$Category[seq_len(n_keep)]
n_lumped <- nrow(cat_totals) - n_keep
agg$Category <- ifelse(agg$Category %in% top_cats, agg$Category, "Other")
agg <- agg %>% group_by(Block, Category) %>% summarise(value = sum(value), .groups = "drop")
lvl <- c(top_cats, if (n_lumped > 0) "Other")
agg$Category <- factor(agg$Category, levels = lvl)
agg$Block    <- factor(agg$Block,    levels = block_order)

# How far the two panels disagree, stated rather than left to be spotted.
top_tot <- sum(blockstats$`Cross-traits`)
bot_tot <- sum(agg$value)
message(sprintf("upper panel total %d %s, lower panel total %d (%s)",
                top_tot, ct_metric, bot_tot,
                if (top_tot == bot_tot) "identical" else
                  sprintf("%d uncategorised, drawn above and not below", top_tot - bot_tot)))
message(sprintf("panel: %d blocks x %d categories%s", length(block_order), n_keep,
                if (n_lumped) sprintf(", %d lumped into 'Other'", n_lumped) else ""))

# ------------------------------------------------------------------
# Colours
# ------------------------------------------------------------------
pal        <- get_palette(opt$palette)

cat_choice <- opt[["cat-palette"]]
cat_ramp <- if (cat_choice %in% names(CAT_RAMPS)) {
  CAT_RAMPS[[cat_choice]]
} else if (identical(cat_choice, "palette")) {
  pal_cat(pal, length(CAT_RAMPS$contrast), CAT_RAMPS$contrast)
} else {
  hex <- trimws(strsplit(cat_choice, ",")[[1]])
  hex <- hex[nzchar(hex)]
  bad <- hex[!grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{3})$", hex)]
  if (!length(hex) || length(bad)) {
    cat("ERROR: --cat-palette must be contrast, tableau, palette, or a comma-separated",
        " list of hex colours", if (length(bad)) paste0(" (bad: ", paste(bad, collapse = " "), ")"),
        "\n", sep = "", file = stderr())
    quit(status = 2)
  }
  hex
}
if (n_keep > length(cat_ramp)) {
  message("NOTE: ", n_keep, " categories exceed the ", length(cat_ramp),
          "-colour '", cat_choice, "' ramp; the extra fills are interpolated and ",
          "neighbouring hues will be close. Consider --top-n ", length(cat_ramp), ".")
  cat_ramp <- grDevices::colorRampPalette(cat_ramp)(n_keep)
}
col_other  <- pal_col(pal, "neutral", "grey75")
col_bar    <- pal_col(pal, "accent",  "steelblue")
INK        <- pal_col(pal, "ink",     "#0b0b0b")
AXIS       <- pal_col(pal, "axis",    "#9aa0ab")
GRIDC      <- pal_col(pal, "grid",    "#e8e8e4")
PAPER      <- pal_col(pal, "paper",   "white")
line_cols  <- pal_cat(pal, 3, c("green4", "red2", "black"))
names(line_cols) <- c("SNPs", "Overlapping genes", "Size")

fill_vals <- setNames(cat_ramp[((seq_along(lvl) - 1) %% length(cat_ramp)) + 1], lvl)
if ("Other" %in% lvl) fill_vals[["Other"]] <- col_other

# ------------------------------------------------------------------
# Draw
#
# The rule between the panels is the upper panel's own x axis line: it sits
# exactly on the shared block axis, which is what the divider is meant to mark.
# Both panels use drop = FALSE on the same levels, so patchwork aligns the two
# plot regions column for column.
# ------------------------------------------------------------------
# Every size below is relative to --base-size, and --base-size is in points at
# the requested --width. That is the whole trick to a figure that has to be
# legible at a given physical size: fix the width to what the poster will
# actually print, then choose the point size for it. Rendering wide and letting
# the layout program shrink the file is what makes 11 pt arrive as 6.5.
base_size <- opt[["base-size"]]
base_theme <- theme_bw(base_size = base_size) +
  theme(panel.background = element_rect(fill = PAPER, colour = NA),
        plot.background  = element_rect(fill = PAPER, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRIDC, linewidth = 0.35),
        panel.border     = element_blank(),
        axis.text.y      = element_text(size = rel(0.85), colour = INK),
        axis.title.y     = element_text(size = rel(0.95), colour = INK),
        legend.background = element_rect(fill = PAPER, colour = NA),
        legend.key        = element_rect(fill = PAPER, colour = NA),
        legend.key.size  = unit(base_size * 0.9, "pt"),
        legend.text      = element_text(size = rel(0.82), colour = INK),
        legend.title     = element_blank(),
        legend.margin    = margin(2, 2, 2, 2),
        legend.box.spacing = unit(base_size * 0.35, "pt"))

# Legend columns: a category key needs roughly 3.1 inches at base size 13 and
# scales with the type, so the number that fits is a function of the physical
# width, not a constant. Long labels ("Age-related biomarkers and traits") set
# the requirement, and at 9.5 inches only three fit where five did at 16.
legend_cols <- opt[["legend-cols"]]
if (legend_cols <= 0L) {
  legend_cols <- max(2L, round(opt$width / (3.1 * base_size / 13)))
}
legend_cols <- min(legend_cols, length(lvl))

# The bar gets a legend key of its own. In the standalone blockstats figure the
# bars are the only unlabelled mark on the plot and the reader has to infer them
# from the title; here the lower panel would answer it eventually, but naming it
# on the upper panel costs one key.
p_top <- ggplot() +
  geom_col(data = df_long %>% filter(Metric == "Cross-traits"),
           aes(x = Block, y = Log2Value, fill = "Cross-traits"),
           width = opt[["bar-width"]]) +
  scale_fill_manual(name = NULL, values = setNames(col_bar, "Cross-traits")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = expression(log[2] * "(count + 1)")) +
  base_theme +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x  = element_line(colour = AXIS, linewidth = 0.7),
        legend.position = "top",
        plot.margin  = margin(4, 6, 0, 6))

if (!opt[["no-lines"]]) {
  # Each context line gets its own dash pattern as well as its own hue. The
  # standalone figure dashes all three identically, which works when the
  # palette is green/red/black; the poster palette's three categorical slots
  # are close in hue by design, and the pattern is what keeps them apart.
  p_top <- p_top +
    geom_line(data = df_long %>% filter(Metric != "Cross-traits"),
              aes(x = Block, y = Log2Value, colour = Metric,
                  linetype = Metric, group = Metric), linewidth = 0.6) +
    scale_colour_manual(values = line_cols) +
    scale_linetype_manual(values = c(SNPs = "twodash",
                                     `Overlapping genes` = "dotted",
                                     Size = "longdash"))
}

p_bot <- ggplot(agg, aes(x = Block, y = value, fill = Category)) +
  geom_col(width = opt[["bar-width"]]) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = fill_vals, drop = FALSE) +
  labs(y = if (ct_metric == "traits") "Cross-traits" else "Associations") +
  base_theme +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                    size = rel(0.78), colour = INK),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin  = margin(0, 6, 4, 6)) +
  guides(fill = guide_legend(ncol = legend_cols, byrow = TRUE))

panel <- p_top / p_bot + plot_layout(heights = heights)

png_path <- file.path(opt$outdir, paste0(opt$outfile, ".png"))
ggsave(png_path, panel, width = opt$width, height = opt$height, dpi = opt$dpi,
       bg = PAPER)
message("PNG saved -> ", png_path)

if (opt$svg) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("NOTE: svglite not installed; skipping the SVG.")
  } else {
    svg_path <- file.path(opt$outdir, paste0(opt$outfile, ".svg"))
    ggsave(svg_path, panel, width = opt$width, height = opt$height,
           device = svglite::svglite, bg = PAPER)
    message("SVG saved -> ", svg_path)
  }
}
