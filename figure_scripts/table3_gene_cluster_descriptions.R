#!/usr/bin/env Rscript
# =============================================================================
# Table 3 — Gene cluster descriptions (panel C)
# Output : <outdir>/table3_gene_clusters.tsv   (data)
#          <outdir>/table3_gene_clusters.svg   (publication figure)
#
# Input  : a STRING "cluster descriptions" export (the k-means cluster file,
#          e.g. string_kmeans_cluster_descriptions_<trait>.tsv). Tab-separated,
#          header line begins with '#'. Expected columns (by name):
#            cluster number, hex color, gene count,
#            primary description, secondary description, tertiary description
#
# Column mapping:
#   Cluster            <- "Cluster <cluster number>"
#   Color (swatch)     <- hex color
#   Main pathway       <- primary description           ('-' shown blank)
#   Brief description  <- secondary (+ tertiary) joined  ('-' shown blank)
#   Genes              <- gene count
#
# By default, clusters with no annotation at all (primary/secondary/tertiary
# all '-') are dropped; pass --include-unannotated to keep them.
#
# Usage:
#   Rscript figure_scripts/table3_gene_cluster_descriptions.R \
#       -i results/intermediate/chronic_lymphocytic_leukemia/string/string_kmeans_cluster_descriptions_cll.tsv \
#       -o figures/chronic_lymphocytic_leukemia
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(svglite)
  library(grid)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "STRING cluster-descriptions TSV [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--outfile"), type = "character", default = "table3_gene_clusters.svg",
              help = "Output SVG filename [default: %default]"),
  make_option(c("--include-unannotated"), action = "store_true", default = FALSE,
              help = "Keep clusters whose descriptions are all '-' [default: drop them]"),
  make_option(c("--svg-width"), type = "double", default = 22.0,
              help = "SVG canvas width in cm [default: %default]"),
  make_option(c("--row-height"), type = "double", default = 0.60,
              help = "Row height in cm [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Table 3: gene cluster descriptions (panel C)")
opt <- parse_args(parser)

if (is.null(opt$input)) {
  cat("ERROR: --input is required\n\n", file = stderr())
  print_help(parser)
  quit(status = 2)
}
if (!file.exists(opt$input)) {
  cat("ERROR: File not found:", opt$input, "\n", file = stderr())
  quit(status = 2)
}
if (is.null(opt$outdir)) {
  cat("ERROR: --outdir is required\n\n", file = stderr())
  print_help(parser)
  quit(status = 2)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DATA
# =============================================================================

# read.delim with comment.char="" so the '#'-prefixed STRING header is kept.
clusters <- read.delim(opt$input, sep = "\t", quote = "", comment.char = "",
                       check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("cluster number", "hex color", "gene count",
                   "primary description", "secondary description",
                   "tertiary description")
missing <- setdiff(required_cols, names(clusters))
if (length(missing) > 0) {
  cat("ERROR: STRING input missing expected column(s): ",
      paste(missing, collapse = ", "), "\n", sep = "", file = stderr())
  quit(status = 2)
}

# STRING uses "-" for "no annotation". Blank those out.
dash_blank <- function(x) {
  x <- trimws(as.character(x))
  x[x == "-" | is.na(x)] <- ""
  x
}

main_desc <- dash_blank(clusters[["primary description"]])
sec_desc  <- dash_blank(clusters[["secondary description"]])
ter_desc  <- dash_blank(clusters[["tertiary description"]])

# Brief = secondary (+ tertiary), dropping empties, joined with "; ".
brief_desc <- mapply(function(a, b) {
  parts <- c(a, b)
  paste(parts[nzchar(parts)], collapse = "; ")
}, sec_desc, ter_desc, USE.NAMES = FALSE)

panel_c_data <- data.frame(
  Cluster = paste0("Cluster ", clusters[["cluster number"]]),
  color   = clusters[["hex color"]],
  Main    = main_desc,
  Brief   = brief_desc,
  Genes   = as.integer(clusters[["gene count"]]),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Drop fully-unannotated clusters unless asked to keep them.
if (!opt[["include-unannotated"]]) {
  keep <- nzchar(panel_c_data$Main) | nzchar(panel_c_data$Brief)
  message("Loaded ", nrow(panel_c_data), " clusters; keeping ", sum(keep),
          " annotated (dropped ", sum(!keep), " with no description).")
  panel_c_data <- panel_c_data[keep, ]
} else {
  message("Loaded ", nrow(panel_c_data), " clusters (keeping unannotated).")
}

if (nrow(panel_c_data) == 0) {
  cat("ERROR: no clusters left to tabulate.\n", file = stderr())
  quit(status = 2)
}

# TSV export (mirrors the figure data)
tsv_path <- file.path(opt$outdir, sub("\\.svg$", ".tsv", opt$outfile))
write.table(
  setNames(panel_c_data, c("Cluster", "Color", "Main pathway",
                           "Brief description", "Genes")),
  tsv_path, sep = "\t", quote = FALSE, row.names = FALSE
)
message("TSV saved -> ", tsv_path)

# =============================================================================
# STYLE & LAYOUT
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

n_rows  <- nrow(panel_c_data)
n_cols  <- 5L
W       <- svg_width
H       <- margin_top + header_h + n_rows * row_h + margin_bottom

# Cluster | Color | Main pathway | Brief description | Genes
col_w <- c(0.11, 0.07, 0.34, 0.40, 0.08)
stopifnot("col_w must sum to 1" = abs(sum(col_w) - 1) < 1e-9)

right_cols  <- c(5L)
center_cols <- c(2L)
col_labels  <- c("Cluster", "Color", "Main pathway", "Brief description", "Genes")

x_left       <- cumsum(c(0, col_w[-n_cols]))
x_right      <- x_left + col_w
x_center     <- (x_left + x_right) / 2
y_top_rule   <- 1 - margin_top / H
y_midrule    <- y_top_rule - header_h / H
y_header_mid <- (y_top_rule + y_midrule) / 2
y_bot_rule   <- margin_bottom / H

# =============================================================================
# SVG
# =============================================================================

svg_path <- file.path(opt$outdir, opt$outfile)
svglite(svg_path, width = W / 2.54, height = H / 2.54,
        system_fonts = list(sans = FONT), fix_text_size = FALSE)
grid.newpage()

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_top_rule,y_top_rule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

for (j in seq_len(n_cols)) {
  is_r  <- j %in% right_cols
  is_c  <- j %in% center_cols
  x_txt <- if (is_r) (x_right[j] - pad_r / W) else if (is_c) x_center[j] else (x_left[j] + pad_l / W)
  grid.text(label = col_labels[j],
            x = unit(x_txt,"npc"), y = unit(y_header_mid,"npc"),
            hjust = if (is_r) 1 else if (is_c) 0.5 else 0, vjust = 0.5,
            gp = gpar(fontfamily = FONT, fontface = "bold",
                      fontsize = FS_HEADER, col = COL_TEXT))
}

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_midrule,y_midrule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THIN))

for (i in seq_len(n_rows)) {
  y_mid     <- y_midrule - ((i - 0.5) * row_h) / H
  row_h_npc <- row_h / H

  for (j in seq_len(n_cols)) {
    is_r <- j %in% right_cols
    is_c <- j %in% center_cols

    if (j == 2L) {
      grid.rect(x = unit(x_center[j], "npc"), y = unit(y_mid, "npc"),
                width  = unit(col_w[j] * 0.55, "npc"),
                height = unit(row_h_npc * 0.65, "npc"),
                gp = gpar(fill = panel_c_data$color[i], col = COL_LINE, lwd = 0.4))
    } else {
      val <- switch(as.character(j),
        "1" = panel_c_data$Cluster[i],
        "3" = panel_c_data$Main[i],
        "4" = panel_c_data$Brief[i],
        "5" = as.character(panel_c_data$Genes[i])
      )
      x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
      grid.text(label = val,
                x = unit(x_txt,"npc"), y = unit(y_mid,"npc"),
                hjust = if (is_r) 1 else 0, vjust = 0.5,
                gp = gpar(fontfamily = FONT, fontface = "plain",
                          fontsize = FS_BODY, col = COL_TEXT))
    }
  }
}

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_bot_rule,y_bot_rule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

dev.off()
message("SVG saved -> ", svg_path)
