#!/usr/bin/env Rscript
# =============================================================================
# Effect size per category — summary forest + violin/box, one of each per
# trait category.
#
# Within each category, the hierarchized traits are the sub-units (typically
# a handful each). For every category this writes:
#   <outdir>/effect_size_per_cat/<Category>_forest.{png,svg}
#       summary forest: one row per hierarchized trait, mean harmonized beta
#       with 95% CI of the mean (mean +/- 1.96*SE across that trait's
#       associations; no whisker when n = 1).
#   <outdir>/effect_size_per_cat/<Category>_violin.{png,svg}
#       violin + box + points of the harmonized betas, one violin per
#       hierarchized trait.
#
# Per-association CIs in the harmonized table are mostly empty, so the forest
# CI is derived from the spread of harmonized_beta within each trait rather
# than pooled from per-association CIs.
#
# Input  : 09_harmonized_trait_table.tsv (or .rds) — needs Category,
#          'Hierarchized Trait', harmonized_beta.
#
# Usage:
#   Rscript figure_scripts/fig_effect_size_per_category.R \
#       -i results/intermediate/leptin/09_harmonized_trait_table.tsv \
#       -o figures/leptin --label "Leptin"
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Harmonized trait table (09_harmonized_trait_table.tsv or .rds) [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory; figures go in <outdir>/effect_size_per_cat/ [required]"),
  make_option(c("--label"), type = "character", default = NULL,
              help = "Optional phenotype label prefixed to each title (e.g. 'Leptin')"),
  make_option(c("--svg"), action = "store_true", default = FALSE,
              help = "Also write .svg alongside each .png (requires svglite)"),
  make_option(c("--dpi"), type = "integer", default = 300L,
              help = "PNG resolution [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Per-category effect-size forest and violin plots")
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

source("R/figure_utils.R")
source("R/palettes.R")

dat <- read_any(opt$input, as_data_frame = TRUE)
need <- c("Category", "Hierarchized Trait", "harmonized_beta")
missing <- setdiff(need, names(dat))
if (length(missing) > 0) {
  cat("ERROR: input missing column(s): ", paste(missing, collapse = ", "), "\n",
      sep = "", file = stderr()); quit(status = 2)
}

# Hierarchized-trait (and category) labels can arrive wrapped in literal double
# quotes — embedded upstream when the name contains a comma. Strip a matched
# surrounding pair so axis labels read cleanly.
strip_quotes <- function(x) {
  x <- trimws(as.character(x))
  x <- sub('^"(.*)"$', "\\1", x)
  gsub('""', '"', x, fixed = TRUE)
}
dat$Category <- strip_quotes(dat$Category)
dat[["Hierarchized Trait"]] <- strip_quotes(dat[["Hierarchized Trait"]])

dat$harmonized_beta <- suppressWarnings(as.numeric(dat$harmonized_beta))
dat <- dat[!is.na(dat$harmonized_beta) &
           !is.na(dat$Category) & nzchar(dat$Category) &
           !is.na(dat$`Hierarchized Trait`) & nzchar(dat$`Hierarchized Trait`), ]
if (nrow(dat) == 0) {
  cat("ERROR: no rows with a numeric harmonized_beta and a Category/Trait.\n",
      file = stderr()); quit(status = 2)
}

out_dir <- file.path(opt$outdir, "effect_size_per_cat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal        <- get_palette(opt$palette)
DIR_COLORS <- c(Positive = pal_col(pal, "pos", "#d73027"),
                Negative = pal_col(pal, "neg", "#4575b4"))
col_paper  <- pal_col(pal, "paper",   "white")
col_ref    <- pal_col(pal, "axis",    "grey50")   # zero-reference lines
col_jitter <- pal_col(pal, "muted",   "grey25")

sign_lab <- function(x) ifelse(x >= 0, "Positive", "Negative")
title_for <- function(cat) {
  cat <- stringr::str_to_sentence(cat)
  if (is.null(opt$label)) cat else paste0(opt$label, " — ", cat)
}
safe <- function(x) gsub("[^A-Za-z0-9]", "_", x)

save_both <- function(plot, base, width, height) {
  ggsave(paste0(base, ".png"), plot = plot, device = "png",
         width = width, height = height, units = "in", dpi = opt$dpi, bg = col_paper)
  if (opt$svg) ggsave(paste0(base, ".svg"), plot = plot, device = "svg",
                      width = width, height = height, fix_text_size = FALSE)
}

categories <- sort(unique(dat$Category))
message("Generating effect-size plots for ", length(categories), " categories -> ", out_dir)

for (cat in categories) {
  sub <- dat[dat$Category == cat, ]

  # ---- summary forest: per hierarchized trait, mean beta +/- 95% CI of mean
  fdat <- sub %>%
    group_by(Trait = `Hierarchized Trait`) %>%
    summarise(
      n    = dplyr::n(),
      mean = mean(harmonized_beta),
      se   = if (dplyr::n() > 1) sd(harmonized_beta) / sqrt(dplyr::n()) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      lo  = mean - 1.96 * se,
      hi  = mean + 1.96 * se,
      dir = sign_lab(mean),
      Trait = factor(Trait, levels = Trait[order(mean)])
    )

  n_tr <- nrow(fdat)
  forest <- ggplot(fdat, aes(x = mean, y = Trait, color = dir)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = col_ref) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, na.rm = TRUE) +
    geom_point(aes(size = n)) +
    scale_color_manual(values = DIR_COLORS, name = "Direction") +
    scale_size_continuous(name = "N assoc.", range = c(2, 6)) +
    labs(x = "Harmonized effect (beta)", y = NULL, title = title_for(cat)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 10))
  save_both(forest, file.path(out_dir, paste0(safe(cat), "_forest")),
            width = 8, height = max(2.5, 0.55 * n_tr + 1.8))

  # ---- violin/box: harmonized beta distribution per hierarchized trait
  vdat <- sub
  trait_levels <- fdat$Trait                       # same order as forest
  vdat$Trait <- factor(vdat$`Hierarchized Trait`, levels = levels(trait_levels))
  med <- vdat %>% group_by(Trait) %>%
    summarise(dir = sign_lab(median(harmonized_beta)), .groups = "drop")
  vdat <- left_join(vdat, med, by = "Trait")

  violin <- ggplot(vdat, aes(x = Trait, y = harmonized_beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = col_ref) +
    geom_violin(aes(fill = dir), color = NA, alpha = 0.45, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = col_paper) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.35, color = col_jitter) +
    scale_fill_manual(values = DIR_COLORS, name = "Direction") +
    coord_flip() +
    labs(x = NULL, y = "Harmonized effect (beta)", title = title_for(cat)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 10))
  save_both(violin, file.path(out_dir, paste0(safe(cat), "_violin")),
            width = 8, height = max(2.5, 0.6 * n_tr + 1.8))
}

message("Done. ", length(categories), " categories x 2 plots written to ", out_dir)
