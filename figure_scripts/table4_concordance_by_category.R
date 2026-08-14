#!/usr/bin/env Rscript
# =============================================================================
# Table 4 — Effect direction concordance by category
# Shows the numeric detail behind the category-level concordance bar chart:
# N traits, N SNP-trait pairs, concordant/discordant counts, %, BH adj. p
#
# Output : <outdir>/table4_concordance_category.tsv
#          <outdir>/table4_concordance_category.svg
#
# Prerequisite: run fig_effect_direction_and_forest_plots.R first to
#               produce direction_summary_by_category.txt
#
# Usage:
#   Rscript figure_scripts/table4_concordance_by_category.R -i direction_summary_by_category.txt
#   Rscript figure_scripts/table4_concordance_by_category.R -i dir.txt -o Figures --svg-width 24
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(svglite)
  library(grid)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to direction_summary_by_category.txt [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--svg-width"), type = "double", default = 22.0,
              help = "SVG canvas width in cm [default: %default]"),
  make_option(c("--row-height"), type = "double", default = 0.52,
              help = "Row height in cm; raise to 0.65 if text clips [default: %default]"),
  make_option(c("--min-pairs"), type = "integer", default = 3L,
              help = "Min SNP-trait pairs to include a category [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Table 4: effect direction concordance by category")
opt <- parse_args(parser)

if (is.null(opt$input)) {
  cat("ERROR: --input is required\n\n", file = stderr())
  print_help(parser)
  quit(status = 2)
}
if (is.null(opt$outdir)) {
  cat("ERROR: --outdir is required\n\n", file = stderr())
  print_help(parser)
  quit(status = 2)
}
if (!file.exists(opt$input)) {
  cat("ERROR: File not found:", opt$input, "\n", file = stderr())
  quit(status = 2)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. BUILD DISPLAY TABLE
# =============================================================================

direction_category <- read.delim(opt$input, sep = "\t", check.names = FALSE)

concordance_display <- direction_category %>%
  filter(n_snps >= opt[["min-pairs"]]) %>%
  arrange(desc(pct_concordant)) %>%
  mutate(
    Category = stringr::str_to_sentence(Category),
    sig = case_when(
      is.na(binom_p_adj)  ~ "",
      binom_p_adj < 0.001 ~ "***",
      binom_p_adj < 0.01  ~ "**",
      binom_p_adj < 0.05  ~ "*",
      TRUE                ~ ""
    ),
    `BH adj. p` = case_when(
      is.na(binom_p_adj)  ~ "NA",
      binom_p_adj >= 0.01 ~ formatC(binom_p_adj, format = "f", digits = 3),
      TRUE                ~ formatC(binom_p_adj, format = "e", digits = 1)
    ),
    `% Concordant` = paste0(formatC(pct_concordant, format = "f", digits = 1), sig)
  ) %>%
  select(
    Category,
    `N traits`    = n_traits,
    `N pairs`     = n_snps,
    `Concordant`  = n_concordant,
    `Discordant`  = n_discordant,
    `% Concordant`,
    `BH adj. p`
  )

# =============================================================================
# 2. TSV EXPORT
# =============================================================================

tsv_path <- file.path(opt$outdir, "table4_concordance_category.tsv")
write.table(concordance_display, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("TSV saved -> ", tsv_path)

# =============================================================================
# 3. SVG FIGURE
# =============================================================================

FONT      <- "Arial"
FS_HEADER <- 8
FS_BODY   <- 8
LWD_THICK <- 1.5
LWD_THIN  <- 0.75

# Colours, resolved once. Each fallback is the literal this script used
# before it took a --palette, so the default rendering is unchanged.
source("R/palettes.R")
pal       <- get_palette(opt$palette)
COL_LINE  <- pal_col(pal, "axis", "black")
COL_TEXT  <- pal_col(pal, "ink",  "black")

svg_width     <- opt[["svg-width"]]
row_h         <- opt[["row-height"]]
header_h      <- 0.62
margin_top    <- 0.12
margin_bottom <- 0.20
pad_l         <- 0.25
pad_r         <- 0.25

n_rows <- nrow(concordance_display)
n_cols <- ncol(concordance_display)
W      <- svg_width
H      <- margin_top + header_h + n_rows * row_h + margin_bottom

# Column proportions (must sum to 1):
#   Category       ~25 chars
#   N traits       8 chars  [right]
#   N pairs        7 chars  [right]
#   Concordant    10 chars  [right]
#   Discordant    10 chars  [right]
#   % Concordant   8 chars  [right]
#   BH adj. p      8 chars  [right]
col_w <- c(0.34, 0.09, 0.09, 0.11, 0.11, 0.13, 0.13)
stopifnot("col_w must sum to 1" = abs(sum(col_w) - 1) < 1e-9)

right_cols <- c(2L, 3L, 4L, 5L, 6L, 7L)

x_left       <- cumsum(c(0, col_w[-n_cols]))
x_right      <- x_left + col_w
y_top_rule   <- 1 - margin_top / H
y_midrule    <- y_top_rule - header_h / H
y_header_mid <- (y_top_rule + y_midrule) / 2
y_bot_rule   <- margin_bottom / H

svg_path <- file.path(opt$outdir, "table4_concordance_category.svg")
svglite(svg_path, width = W / 2.54, height = H / 2.54,
        system_fonts = list(sans = FONT), fix_text_size = FALSE)
grid.newpage()

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_top_rule,y_top_rule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

for (j in seq_len(n_cols)) {
  is_r  <- j %in% right_cols
  x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
  grid.text(label = names(concordance_display)[j],
            x = unit(x_txt,"npc"), y = unit(y_header_mid,"npc"),
            hjust = if (is_r) 1 else 0, vjust = 0.5,
            gp = gpar(fontfamily = FONT, fontface = "bold",
                      fontsize = FS_HEADER, col = COL_TEXT))
}

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_midrule,y_midrule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THIN))

for (i in seq_len(n_rows)) {
  y_mid <- y_midrule - ((i - 0.5) * row_h) / H
  for (j in seq_len(n_cols)) {
    is_r  <- j %in% right_cols
    x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
    grid.text(label = as.character(concordance_display[[j]][i]),
              x = unit(x_txt,"npc"), y = unit(y_mid,"npc"),
              hjust = if (is_r) 1 else 0, vjust = 0.5,
              gp = gpar(fontfamily = FONT, fontface = "plain",
                        fontsize = FS_BODY, col = COL_TEXT))
  }
}

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_bot_rule,y_bot_rule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

dev.off()
message("SVG saved -> ", svg_path)
