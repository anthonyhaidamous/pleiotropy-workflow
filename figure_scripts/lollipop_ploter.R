#!/usr/bin/env Rscript
# =============================================================================
# Category enrichment lollipop plot
#
# Diverging lollipop of per-category enrichment (default: Normalized Enrichment,
# i.e. log2 fold enrichment around 0), one row per trait category, annotated
# with each category's frequency. Built from the stage-8 category enrichment
# table.
#
# Input  : 08_category_enrichment.tsv (or .rds) — needs the value, label and
#          Category columns (default: "Normalized Enrichment", "Frequency").
#          If a significance column is present (default "Stars", the */**/***
#          from stage 8's Fisher test), it is appended to the category labels.
# Output : <outdir>/<outfile>.png   (and .svg with --svg)
#
# Usage:
#   Rscript figure_scripts/lollipop_ploter.R \
#       -i results/intermediate/leptin/08_category_enrichment.tsv \
#       -o figures/leptin
#   # custom value column / title / bold threshold:
#   Rscript figure_scripts/lollipop_ploter.R -i .../08_category_enrichment.tsv \
#       -o figures/leptin --value "Fold Enrichment" --title "Leptin" --bold-threshold 100
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Category enrichment table (08_category_enrichment.tsv or .rds) [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--outfile"), type = "character", default = "category_lollipop",
              help = "Output base filename, no extension [default: %default]"),
  make_option(c("--value"), type = "character", default = "Normalized Enrichment",
              help = "Column to plot on the lollipop axis [default: %default]"),
  make_option(c("--label"), type = "character", default = "Frequency",
              help = "Column shown as the point annotation [default: %default]"),
  make_option(c("--category-col"), type = "character", default = "Category",
              help = "Column holding the category names [default: %default]"),
  make_option(c("--exclude"), type = "character", default = NULL,
              help = "Comma-separated category names to exclude (case-insensitive)"),
  make_option(c("--drop-others"), action = "store_true", default = FALSE,
              help = "Drop the 'Others' catch-all category (shortcut for --exclude Others)"),
  make_option(c("--stars-col"), type = "character", default = "Stars",
              help = "Column holding the significance stars (***/**/*, - = n.s.) from stage 8 [default: %default]"),
  make_option(c("--sig-display"), type = "character", default = "stars",
              help = "Significance panel content: 'stars' or 'pvalue' (the numeric value) [default: %default]"),
  make_option(c("--pvalue-col"), type = "character", default = "P_Adjusted",
              help = "Numeric significance column shown when --sig-display pvalue [default: %default]"),
  make_option(c("--no-significance"), action = "store_true", default = FALSE,
              help = "Draw a single-panel lollipop (no separate names/stars panels)"),
  make_option(c("--panel-widths"), type = "character", default = "3.2,0.5,8",
              help = "Relative widths of the names,stars,plot panels [default: %default]"),
  make_option(c("--bold-threshold"), type = "double", default = 50,
              help = "Annotations >= this value are drawn bold [default: %default]"),
  make_option(c("--title"), type = "character", default = NULL,
              help = "Optional plot title"),
  make_option(c("--width"), type = "double", default = 16.0,
              help = "Figure width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 7.5,
              help = "Figure height in inches [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300L,
              help = "PNG resolution [default: %default]"),
  make_option(c("--svg"), action = "store_true", default = FALSE,
              help = "Also write an .svg (requires the svglite package)"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Category enrichment lollipop plot")
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
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Read either the RDS or the TSV mirror.
source("R/figure_utils.R")
source("R/palettes.R")

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal        <- get_palette(opt$palette)
col_ink    <- pal_col(pal, "ink",    "#2c3e50")   # stems, labels, caption
col_point  <- pal_col(pal, "accent", "steelblue") # the lollipop head
col_grid   <- pal_col(pal, "grid",   "gray90")
col_paper  <- pal_col(pal, "paper",  "white")

dat <- read_any(opt$input, as_data_frame = TRUE)

val_col <- opt$value
lab_col <- opt$label
cat_col <- opt[["category-col"]]
for (col in c(cat_col, val_col, lab_col)) {
  if (!col %in% names(dat)) {
    cat("ERROR: column '", col, "' not found in ", opt$input,
        " (have: ", paste(names(dat), collapse = ", "), ")\n",
        sep = "", file = stderr())
    quit(status = 2)
  }
}

dat[[val_col]] <- as.numeric(dat[[val_col]])
dat <- dat[!is.na(dat[[val_col]]), ]
# Sentence-case category names
dat[[cat_col]] <- stringr::str_to_sentence(dat[[cat_col]])
# Drop excluded categories (matched case-insensitively against the names above)
drop <- character(0)
if (!is.null(opt$exclude)) {
  drop <- trimws(strsplit(opt$exclude, ",")[[1]])
  drop <- drop[nzchar(drop)]
}
if (opt[["drop-others"]]) drop <- unique(c(drop, "Others"))
if (length(drop)) {
  unmatched <- drop[!tolower(drop) %in% tolower(dat[[cat_col]])]
  if (length(unmatched))
    message("WARNING: excluded category(ies) not found: ", paste(unmatched, collapse = ", "))
  dat <- dat[!tolower(dat[[cat_col]]) %in% tolower(drop), ]
  if (nrow(dat) == 0) {
    cat("ERROR: every category was excluded; nothing to plot\n", file = stderr()); quit(status = 2)
  }
}
# Order categories by the plotted value (reorder + coord_flip => smallest at bottom)
dat[[cat_col]] <- factor(dat[[cat_col]], levels = dat[[cat_col]][order(dat[[val_col]])])

# Significance panel (from stage-8 Fisher's exact test). It shows either the
# ***/**/* stars ("-" when n.s.) or the actual numeric value, per --sig-display.
sig_mode <- tolower(opt[["sig-display"]])
if (!sig_mode %in% c("stars", "pvalue")) {
  message("NOTE: --sig-display must be 'stars' or 'pvalue'; using 'stars'.")
  sig_mode <- "stars"
}
star_col    <- opt[["stars-col"]]
pval_col    <- opt[["pvalue-col"]]
sig_src_col <- if (sig_mode == "pvalue") pval_col else star_col
have_sig    <- sig_src_col %in% names(dat)
use_panels  <- !opt[["no-significance"]] && have_sig
if (use_panels) {
  if (sig_mode == "pvalue") {
    pv <- suppressWarnings(as.numeric(dat[[pval_col]]))
    # Compact: scientific below 0.001, three decimals otherwise.
    dat$.sig <- ifelse(is.na(pv), "NA",
                ifelse(pv < 1e-3, formatC(pv, format = "e", digits = 1),
                       formatC(pv, format = "f", digits = 3)))
  } else {
    dat$.sig <- as.character(dat[[star_col]])
    dat$.sig[is.na(dat$.sig) | !nzchar(dat$.sig)] <- "-"   # blank/NA -> n.s.
  }
} else if (!opt[["no-significance"]] && !have_sig) {
  message("NOTE: no '", sig_src_col, "' column found; drawing a single-panel lollipop ",
          "without significance.")
}

# Bold the annotations at/above the threshold (only meaningful for numeric labels)
lab_num <- suppressWarnings(as.numeric(dat[[lab_col]]))
dat$.face <- ifelse(!is.na(lab_num) & lab_num >= opt[["bold-threshold"]], "bold", "plain")
# Nudge the annotation outward from the baseline, in the direction of the value
dat$.lab_y <- dat[[val_col]] + ifelse(dat[[val_col]] > 0, 0.3, -0.3)

# geom_lollipop() (ggalt) reproduced with base geoms: stem + point.
# The category axis labels live in their own panel (p_names) when use_panels,
# so the main plot hides them.
p_main <- ggplot(dat, aes(x = .data[[cat_col]], y = .data[[val_col]])) +
  geom_segment(aes(xend = .data[[cat_col]], y = 0, yend = .data[[val_col]]),
               color = col_ink, linewidth = 0.8) +
  geom_point(color = col_point, size = 6) +
  geom_text(aes(y = .lab_y, label = .data[[lab_col]], fontface = .face),
            size = 4.5, color = col_ink) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title       = element_blank(),
    axis.text.y      = if (use_panels) element_blank() else element_text(size = 15, face = "bold"),
    axis.ticks.y     = if (use_panels) element_blank() else element_line(linewidth = 0.6),
    axis.text.x      = element_text(size = 12),
    axis.ticks.x     = element_line(linewidth = 0.6),
    panel.grid.major = element_line(color = col_grid, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    # Follow the palette's paper, or theme_bw()'s white when none is active,
    # so a tinted palette does not leave a white panel on a tinted margin.
    panel.background = element_rect(fill = col_paper, color = NA),
    plot.background  = element_rect(fill = col_paper, color = NA),
    legend.position  = "none",
    plot.margin      = unit(c(1, 1, 1, 0.2), "cm")
  )

if (use_panels) {
  # Panel 1 (left): category names, right-aligned toward the stars.
  p_names <- ggplot(dat, aes(x = .data[[cat_col]], y = 1)) +
    geom_text(aes(label = .data[[cat_col]]), hjust = 1, size = 5,
              fontface = "bold", color = col_ink) +
    coord_flip(clip = "off") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    theme_void() +
    theme(plot.margin = unit(c(1, 0.1, 1, 1.4), "cm"))   # left room for long names

  # Panel 2 (middle): the significance mark for every category (stars or value).
  p_stars <- ggplot(dat, aes(x = .data[[cat_col]], y = 0)) +
    geom_text(aes(label = .sig), hjust = 0.5, size = 5, color = col_ink,
              fontface = "bold") +
    coord_flip() +
    scale_y_continuous(limits = c(-0.5, 0.5), expand = expansion(mult = c(0, 0))) +
    theme_void() +
    theme(plot.margin = unit(c(1, 0.1, 1, 0.1), "cm"))

  pw <- suppressWarnings(as.numeric(trimws(strsplit(opt[["panel-widths"]], ",")[[1]])))
  if (length(pw) != 3 || any(is.na(pw))) {
    message("NOTE: --panel-widths must be 3 numbers (names,stars,plot); using 2.4,0.5,8.")
    pw <- c(2.4, 0.5, 8)
  }
  # Numeric values need more room than a 1-3 char star mark: widen the middle
  # panel when showing p-values and the user kept the default widths.
  if (sig_mode == "pvalue" && identical(opt[["panel-widths"]], "3.2,0.5,8"))
    pw <- c(3.2, 1.3, 8)

  p <- p_names + p_stars + p_main +
    plot_layout(widths = pw) +
    plot_annotation(
      title   = opt$title,
      caption = if (sig_mode == "pvalue")
                  "Middle column: BH-FDR adjusted p-value (Fisher's exact test)"
                else
                  "BH-FDR Fisher's exact test:  * p<0.05   ** p<0.01   *** p<0.001   - n.s.",
      theme   = theme(plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
                      plot.caption = element_text(hjust = 1, size = 11, color = col_ink))
    )
} else {
  p <- p_main
  if (!is.null(opt$title)) p <- p + ggtitle(opt$title)
}

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
