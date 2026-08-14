#!/usr/bin/env Rscript
# ============================================================
# fig_string_network_by_module.R
#
# The MCL modules drawn on the same frozen layout as the block figure, so the
# two pictures can be compared node for node.
#
#   string_network_modules.png         survivors vs positional artifacts
#   string_network_module_facets.png   one panel per module
#
# Inputs:
#   <stage_dir>/string_02_graph_downstream.rds
#   <stage_dir>/string_03_modules.tsv, string_03_module_summary.tsv
#   <stage_dir>/string_02_layout.tsv          (written by fig_string_network_by_block.R)
#
# Usage:
#   Rscript figure_scripts/fig_string_network_by_module.R config/leptin.yml
#   Rscript figure_scripts/fig_string_network_by_module.R config/cll.yml \
#       --size 4872 --res 660
#
#   --size N   edge of the square canvas, in px, for the modules figure
#              (default: the 2400 x 1900 canvas this script has always used)
#   --res  N   nominal dpi. Text and node sizes are set in points, so they
#              scale with res, not with size: to enlarge the canvas without
#              shrinking the labels, raise both by the same factor.
#   --no-labels  drop the per-gene labels from the modules figure. Below about
#              1200 px the 171 labels are smaller than the type can carry and
#              read as texture rather than text; the legend still names the
#              modules, and the facet figure still names every gene.
# ============================================================

if (!file.exists("R/utils.R")) stop("run from the repository root")
source("R/utils.R")
source("R/string_utils.R")
source("R/palettes.R")
suppressPackageStartupMessages(library(igraph))

args <- commandArgs(trailingOnly = TRUE)

pal_opt <- NULL
if ("--palette" %in% args) {
  i <- which(args == "--palette"); pal_opt <- args[i + 1]; args <- args[-c(i, i + 1)]
}

# Canvas for the main modules figure. The defaults are the literals this script
# has always used, so an unflagged run is unchanged byte for byte; --size makes
# it square, which is also the honest shape, since the frozen layout is
# normalised into a square and igraph draws it at asp = 1 either way. The
# 2400 x 1900 canvas just padded the sides with paper.
fig_w <- 2400L; fig_h <- 1900L; fig_res <- 220L
do_labels <- !("--no-labels" %in% args)
args <- args[args != "--no-labels"]
if ("--size" %in% args) {
  i <- which(args == "--size")
  fig_w <- fig_h <- as.integer(args[i + 1]); args <- args[-c(i, i + 1)]
}
if ("--res" %in% args) {
  i <- which(args == "--res"); fig_res <- as.integer(args[i + 1]); args <- args[-c(i, i + 1)]
}
if (anyNA(c(fig_w, fig_h, fig_res)) || any(c(fig_w, fig_h, fig_res) <= 0)) {
  stop("--size and --res take a positive integer")
}
if (length(args) != 1L) {
  stop("usage: Rscript figure_scripts/fig_string_network_by_module.R config/<pheno>.yml",
       " [--size N] [--res N] [--palette NAME]")
}
cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

g    <- readRDS(file.path(paths$stage_dir, "string_02_graph_downstream.rds"))
mod  <- read.delim(file.path(paths$stage_dir, "string_03_modules.tsv"))
msum <- read.delim(file.path(paths$stage_dir, "string_03_module_summary.tsv"))

lay_path <- file.path(paths$stage_dir, "string_02_layout.tsv")
if (!file.exists(lay_path)) {
  stop("frozen layout not found: ", lay_path,
       "\n  Run figure_scripts/fig_string_network_by_block.R first.")
}
lt  <- read.delim(lay_path)
lay <- as.matrix(lt[match(igraph::V(g)$name, lt$string_id), c("x", "y")])
lay <- igraph::norm_coords(lay, -1, 1, -1, 1)

outdir <- file.path(cfg$paths$figures, "string")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Colours. --palette wins; with no flag get_palette() falls back to the
# FIGURE_PALETTE environment variable, and with neither, to the literals below,
# which are the ones this script used before it knew about palettes. So the
# default rendering is unchanged.
pal      <- get_palette(pal_opt)
CAT_FALLBACK <- c("#2a78d6", "#eb6834", "#1baf7a")
NEUTRAL  <- pal_col(pal, "neutral", "#898781")
INK      <- pal_col(pal, "ink",     "#0b0b0b")
SUBTLE   <- pal_col(pal, "muted",   "#c9c8c3")
COL_PAPER<- pal_col(pal, "paper",   "white")
COL_OFF  <- pal_col(pal, "grid",    "#ebebe8")   # node outside the highlighted set
COL_EDGE <- pal_col(pal, "grid",    "#e8e8e4")

mid  <- mod$module_id[match(igraph::V(g)$gene, mod$gene)]
art  <- msum$positional_artifact[match(mid, msum$module_id)]
keep <- sort(msum$module_id[!msum$positional_artifact])

# One colour per surviving module. This was fixed at three, on the argument
# that three is all a network form can carry — true, but leptin happened to
# have two survivors so the cap never bound. CLL has six, and SLOT[4:6] is NA,
# which igraph draws as an unfilled node: eleven genes silently vanish from the
# figure and from the legend. A hue that is hard to tell from its neighbour is
# a worse figure; a gene that isn't drawn at all is a wrong one. Beyond three
# survivors, read the panels (string_network_module_facets.png), which name
# every module individually.
SLOT <- pal_cat(pal, max(3L, length(keep)), CAT_FALLBACK)

# ---- survivors vs artifacts ------------------------------------------------
vcol <- ifelse(art, NEUTRAL, SLOT[match(mid, keep)])

png(file.path(outdir, "string_network_modules.png"),
    width = fig_w, height = fig_h, res = fig_res)
op <- par(mar = c(0, 0.8, 2.6, 0.8), bg = COL_PAPER)
deg <- igraph::degree(g)
plot(g, layout = lay, rescale = FALSE, xlim = c(-1, 1), ylim = c(-1, 1),
     vertex.color = vcol, vertex.frame.color = COL_PAPER, vertex.frame.width = 1.4,
     vertex.size = 7,
     vertex.label = if (do_labels) igraph::V(g)$gene else NA,
     vertex.label.cex = 0.46,
     vertex.label.color = INK, vertex.label.dist = ifelse(deg >= 6, 1.15, 0.85),
     vertex.label.degree = -pi/2, vertex.label.family = "sans",
     edge.color = SUBTLE, edge.width = 1.1)
title(sprintf("%s — MCL modules (inflation %s), positional artifacts greyed",
              cfg$phenotype$name, s$mcl$inflation),
      cex.main = 0.95, col.main = INK, font.main = 1, adj = 0.5, line = 0.8)
legend("bottomleft", bty = "n", cex = 0.72, pt.cex = 1.5, pch = 21, col = COL_PAPER,
       text.col = INK,
       legend = c(sprintf("M%d  (%d genes, %d blocks)", keep,
                          msum$size[match(keep, msum$module_id)],
                          msum$n_blocks[match(keep, msum$module_id)]),
                  sprintf("positional artifact  (%d genes)", sum(art))),
       pt.bg = c(SLOT[seq_along(keep)], NEUTRAL))
par(op); dev.off()
message("wrote: ", file.path(outdir, "string_network_modules.png"))

# ---- one panel per module --------------------------------------------------
ids  <- msum$module_id[order(-msum$size)]
n    <- length(ids)
ncol <- ceiling(sqrt(n)); nrow <- ceiling(n / ncol)

png(file.path(outdir, "string_network_module_facets.png"),
    width = 380 * ncol, height = 380 * nrow + 60, res = 150)
op <- par(mfrow = c(nrow, ncol), mar = c(0.4, 0.4, 2.2, 0.4),
          oma = c(0, 0, 2.6, 0), bg = COL_PAPER)
for (m in ids) {
  hit <- !is.na(mid) & mid == m
  r   <- msum[msum$module_id == m, ]
  plot(g, layout = lay, rescale = FALSE, xlim = c(-1, 1), ylim = c(-1, 1),
       vertex.color = ifelse(hit, if (r$positional_artifact) NEUTRAL else SLOT[1], COL_OFF),
       vertex.frame.color = ifelse(hit, COL_PAPER, COL_OFF),
       vertex.size = ifelse(hit, 9, 4), vertex.label = NA,
       edge.color = COL_EDGE, edge.width = 0.8)
  title(sprintf("M%d — %d genes, %d block%s%s", m, r$size, r$n_blocks,
                if (r$n_blocks == 1) "" else "s",
                if (r$positional_artifact) "  (artifact)" else ""),
        cex.main = 0.92, col.main = INK, font.main = 1, line = 0.4)
}
if (n < nrow * ncol) for (i in seq_len(nrow * ncol - n)) plot.new()
mtext(sprintf("%s — MCL modules on the STRING network (same layout in every panel)",
              cfg$phenotype$name),
      outer = TRUE, cex = 0.85, col = INK, adj = 0.5, line = 0.6)
par(op); dev.off()
message("wrote: ", file.path(outdir, "string_network_module_facets.png"))
