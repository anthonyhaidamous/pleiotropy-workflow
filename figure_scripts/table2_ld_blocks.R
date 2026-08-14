#!/usr/bin/env Rscript
# =============================================================================
# Table 2 — LD blocks identified in the pleiotropy analysis
# Style  : Nature/Cell 3-line booktabs, black & white, no row shading
# Output : <outdir>/table2_blocks.tsv   (data)
#          <outdir>/table2_blocks.svg   (publication figure)
#
# Input  : the stage-4 blocks table, either 04_blocks.rds (rich tibble with
#          list-columns) or 04_blocks.tsv (flattened mirror). The nine
#          display columns are taken by name from:
#            Block, Genomic_position, Genomic_span, Top_lead_SNP, Genes,
#            Ethnicities, nb_lSNPs, nb_pSNPs, nb_SNPs
#          Top_lead_SNP is "rsID, pos, pval, ORbeta"; only the rsID is shown.
#
# Usage:
#   Rscript figure_scripts/table2_ld_blocks.R \
#       -i results/intermediate/leptin/04_blocks.rds -o figures/leptin
#   Rscript figure_scripts/table2_ld_blocks.R \
#       -i results/intermediate/leptin/04_blocks.tsv -o figures/leptin --svg-width 28
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(svglite)
  library(grid)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to the stage-4 blocks table (04_blocks.rds or 04_blocks.tsv) [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--svg-width"), type = "double", default = 26.0,
              help = "SVG canvas width in cm [default: %default]"),
  make_option(c("--row-height"), type = "double", default = 0.52,
              help = "Row height in cm; raise to 0.65 if text clips [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Table 2: LD blocks for the pleiotropy analysis")
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
# 1. DATA PIPELINE
# =============================================================================

# Read either the rich RDS (list-columns) or the flattened TSV mirror.
ext <- tolower(tools::file_ext(opt$input))
if (ext == "rds") {
  blocks <- as.data.frame(readRDS(opt$input), check.names = FALSE)
} else if (ext %in% c("tsv", "txt")) {
  blocks <- read.delim(opt$input, sep = "\t", quote = "",
                       stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
} else {
  cat("ERROR: unsupported input extension '", ext,
      "' (expected .rds, .tsv or .txt)\n", sep = "", file = stderr())
  quit(status = 2)
}

needed <- c("Block", "Genomic_position", "Genomic_span", "Top_lead_SNP",
            "Genes", "Ethnicities", "nb_lSNPs", "nb_pSNPs", "nb_SNPs")
missing_cols <- setdiff(needed, names(blocks))
if (length(missing_cols) > 0) {
  cat("ERROR: input is missing expected 04_blocks column(s): ",
      paste(missing_cols, collapse = ", "), "\n", sep = "", file = stderr())
  quit(status = 2)
}

# Flatten a column to one string per block, whether it arrives as an RDS
# list-column or an already-collapsed TSV character column.
flatten_field <- function(col) {
  if (is.list(col)) vapply(col, function(x) paste(x, collapse = ", "), character(1))
  else as.character(col)
}

# Deduplicate a comma-separated field. Split on commas ONLY — gene labels use
# " - " for intergenic pairs (e.g. "MIR129-1 - LEP"), which must stay intact.
dedup_csv <- function(x) {
  vapply(x, function(s) {
    toks <- trimws(unlist(strsplit(s, ",")))
    toks <- toks[nzchar(toks) & toks != "NA"]
    paste(unique(toks), collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

# Thousands separators for whole-number display.
fmt_int <- function(x) format(as.numeric(x), big.mark = ",", trim = TRUE)

# "chr7:127239952-128764383" -> "chr7:127,239,952-128,764,383"
fmt_coord <- function(s) {
  vapply(s, function(one) {
    m <- regmatches(one, regexec("^([^:]+):(\\d+)-(\\d+)$", one))[[1]]
    if (length(m) == 4) {
      sprintf("%s:%s-%s", m[2],
              format(as.numeric(m[3]), big.mark = ",", trim = TRUE),
              format(as.numeric(m[4]), big.mark = ",", trim = TRUE))
    } else one
  }, character(1), USE.NAMES = FALSE)
}

blocks_final <- data.frame(
  Block         = flatten_field(blocks$Block),
  Coordinates   = fmt_coord(flatten_field(blocks$Genomic_position)),
  `Size (bp)`   = fmt_int(flatten_field(blocks$Genomic_span)),
  `Main lSNP`   = trimws(sub(",.*$", "", flatten_field(blocks$Top_lead_SNP))),
  Genes         = dedup_csv(flatten_field(blocks$Genes)),
  Ancestry      = dedup_csv(flatten_field(blocks$Ethnicities)),
  `N lSNPs`     = fmt_int(flatten_field(blocks$nb_lSNPs)),
  `N pSNPs`     = fmt_int(flatten_field(blocks$nb_pSNPs)),
  `N SNPs`      = fmt_int(flatten_field(blocks$nb_SNPs)),
  check.names = FALSE, stringsAsFactors = FALSE
)

# =============================================================================
# 2. TSV EXPORT
# =============================================================================

tsv_path <- file.path(opt$outdir, "table2_blocks.tsv")
write.table(blocks_final, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("TSV saved -> ", tsv_path)

# =============================================================================
# 3. SVG FIGURE
# =============================================================================

# --- 3a. Style ---
FONT      <- "Arial"
FS_HEADER <- 8       # pt - header text
FS_BODY   <- 8       # pt - body text
LWD_THICK <- 1.5     # pt - top rule and bottom rule
LWD_THIN  <- 0.75    # pt - midrule (under header)

# Colours, resolved once. Each fallback is the literal this script used
# before it took a --palette, so the default rendering is unchanged.
source("R/palettes.R")
pal       <- get_palette(opt$palette)
COL_LINE  <- pal_col(pal, "axis", "black")
COL_TEXT  <- pal_col(pal, "ink",  "black")

# --- 3b. Layout (cm) ---
svg_width     <- opt[["svg-width"]]
row_h         <- opt[["row-height"]]
header_h      <- 0.62
margin_top    <- 0.12
margin_bottom <- 0.20
pad_l         <- 0.22   # left text padding inside cell
pad_r         <- 0.22   # right text padding inside cell

n_rows <- nrow(blocks_final)
n_cols <- ncol(blocks_final)
W      <- svg_width
H      <- margin_top + header_h + n_rows * row_h + margin_bottom

# --- 3c. Column proportions (must sum exactly to 1) ---
# Calibrated to longest cell content at 8 pt Arial:
#   Block          "B16"                              3 chars
#   Coordinates    "chr7:127,239,952-128,764,383"     28 chars
#   Size (bp)      "1,524,431"                        9 chars  [right]
#   Main lSNP      "rs10487505 / rs780093178"          24 chars
#   Genes          "LEP, SLC35D2, ..."                varies
#   Ancestry       "EUR, AFR, EAS"                    13 chars
#   N lSNPs        "N lSNPs" / "3"                    7 chars  [right]
#   N pSNPs        "N pSNPs" / "596"                  7 chars  [right]
#   N SNPs         "N SNPs"  / "6,688"                6 chars  [right]
col_w <- c(0.05, 0.15, 0.09, 0.25, 0.16, 0.12, 0.06, 0.06, 0.06)
stopifnot("col_w must sum to 1" = abs(sum(col_w) - 1) < 1e-9)

# --- 3d. Right-aligned columns (1-based index) ---
right_cols <- c(3L, 7L, 8L, 9L)   # Size (bp), N lSNPs, N pSNPs, N SNPs

# --- 3e. Derived coordinates ---
x_left       <- cumsum(c(0, col_w[-n_cols]))
x_right      <- x_left + col_w
y_top_rule   <- 1 - margin_top / H
y_midrule    <- y_top_rule - header_h / H
y_header_mid <- (y_top_rule + y_midrule) / 2
y_bot_rule   <- margin_bottom / H

# --- 3f. Open SVG device ---
svg_path <- file.path(opt$outdir, "table2_blocks.svg")
svglite(svg_path, width = W / 2.54, height = H / 2.54,
        system_fonts = list(sans = FONT), fix_text_size = FALSE)
grid.newpage()

# --- 3g. Top rule (thick) ---
grid.lines(x  = unit(c(0, 1), "npc"), y = unit(c(y_top_rule, y_top_rule), "npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

# --- 3h. Header text ---
for (j in seq_len(n_cols)) {
  is_r  <- j %in% right_cols
  x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
  grid.text(label = names(blocks_final)[j],
            x = unit(x_txt, "npc"), y = unit(y_header_mid, "npc"),
            hjust = if (is_r) 1 else 0, vjust = 0.5,
            gp = gpar(fontfamily = FONT, fontface = "bold",
                      fontsize = FS_HEADER, col = COL_TEXT))
}

# --- 3i. Midrule (thin) ---
grid.lines(x  = unit(c(0, 1), "npc"), y = unit(c(y_midrule, y_midrule), "npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THIN))

# --- 3j. Body rows ---
for (i in seq_len(n_rows)) {
  y_mid <- y_midrule - ((i - 0.5) * row_h) / H
  for (j in seq_len(n_cols)) {
    is_r  <- j %in% right_cols
    x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
    grid.text(label = as.character(blocks_final[[j]][i]),
              x = unit(x_txt, "npc"), y = unit(y_mid, "npc"),
              hjust = if (is_r) 1 else 0, vjust = 0.5,
              gp = gpar(fontfamily = FONT, fontface = "plain",
                        fontsize = FS_BODY, col = COL_TEXT))
  }
}

# --- 3k. Bottom rule (thick) ---
grid.lines(x  = unit(c(0, 1), "npc"), y = unit(c(y_bot_rule, y_bot_rule), "npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

dev.off()
message("SVG saved -> ", svg_path)

# =============================================================================
# TUNING REFERENCE
# -----------------------------------------------------------------------------
# --svg-width   raw canvas width in cm (26 = comfortable for 9 cols;
#               shrink to ~17 cm in Inkscape for a Nature single column)
# --row-height  raise to 0.65 if any body text is clipped vertically
# col_w         9 values summing to 1 — edit in-script to adjust column widths
# =============================================================================
