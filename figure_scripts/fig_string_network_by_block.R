#!/usr/bin/env Rscript
# ============================================================
# fig_string_network_by_block.R
#
# Draws the frozen STRING network with each gene coloured by its LD block
# of origin, so you can see how the blocks sit on the network.
#
# Two figures:
#   string_network_by_block.png      the network, blocks coloured
#   string_network_block_facets.png  one small panel per block, same layout
#
# The layout is computed once and written to disk, then reused. Freezing the
# coordinates is what makes the picture reproducible — set.seed() is not
# enough, because igraph's RNG use has changed between versions.
#
# Inputs (from pipeline/string_pipeline/stage_string_02_components.R):
#   results/intermediate/<slug>/string_modules/string_02_graph_full.rds
#
# Usage:
#   Rscript figure_scripts/fig_string_network_by_block.R config/leptin.yml
#   Rscript figure_scripts/fig_string_network_by_block.R config/leptin.yml --relayout
#
#   --relayout  recompute and overwrite the frozen layout coordinates
# ============================================================

if (!file.exists("R/utils.R")) stop("run from the repository root")
source("R/utils.R")
source("R/string_utils.R")
source("R/palettes.R")
suppressPackageStartupMessages(library(igraph))

args     <- commandArgs(trailingOnly = TRUE)
relayout <- "--relayout" %in% args
args     <- args[args != "--relayout"]
pal_opt <- NULL
if ("--palette" %in% args) {
  i <- which(args == "--palette"); pal_opt <- args[i + 1]; args <- args[-c(i, i + 1)]
}
if (length(args) != 1L) {
  stop("usage: Rscript figure_scripts/fig_string_network_by_block.R config/<pheno>.yml [--relayout] [--palette NAME]")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

g_path <- file.path(paths$stage_dir, "string_02_graph_full.rds")
if (!file.exists(g_path)) {
  stop("graph not found: ", g_path,
       "\n  Run pipeline/string_pipeline/stage_string_02_components.R first.")
}
g_full <- readRDS(g_path)

outdir <- file.path(cfg$paths$figures, "string")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- the network is the genes that actually have an association ------------
g <- igraph::induced_subgraph(g_full, igraph::V(g_full)[igraph::degree(g_full) > 0])
message("network: ", igraph::vcount(g), " genes, ", igraph::ecount(g), " associations")

# ---- frozen layout ---------------------------------------------------------
lay_path <- file.path(paths$stage_dir, "string_02_layout.tsv")
if (file.exists(lay_path) && !relayout) {
  lt  <- readr::read_tsv(lay_path, progress = FALSE, show_col_types = FALSE)
  lay <- as.matrix(lt[match(igraph::V(g)$name, lt$string_id), c("x", "y")])
  message("layout: reusing frozen coordinates -> ", lay_path)
} else {
  set.seed(1)                                   # only seeds THIS computation;
  lay <- igraph::layout_with_fr(g)              # the coordinates below are the freeze
  readr::write_tsv(
    tibble::tibble(string_id = igraph::V(g)$name, gene = igraph::V(g)$gene,
                   x = lay[, 1], y = lay[, 2]),
    lay_path)
  message("layout: computed and frozen -> ", lay_path)
}
lay <- igraph::norm_coords(lay, -1, 1, -1, 1)

# ---- colour by block -------------------------------------------------------
# Categorical hues are capped at three for a form where every pair can end up
# adjacent on screen. The three largest blocks get hues; the remainder share a
# neutral, and the facet panels below name every block individually.
# Colours. --palette wins; with no flag get_palette() falls back to the
# FIGURE_PALETTE environment variable, and with neither, to the literals below,
# which are the ones this script used before it knew about palettes. So the
# default rendering is unchanged.
pal      <- get_palette(pal_opt)
SLOT     <- pal_cat(pal, 3, c("#2a78d6", "#eb6834", "#1baf7a"))  # blue, orange, aqua
NEUTRAL  <- pal_col(pal, "neutral", "#898781")
INK      <- pal_col(pal, "ink",     "#0b0b0b")
SUBTLE   <- pal_col(pal, "muted",   "#c9c8c3")
COL_PAPER<- pal_col(pal, "paper",   "white")
COL_OFF  <- pal_col(pal, "grid",    "#ebebe8")   # node outside the highlighted set
COL_EDGE <- pal_col(pal, "grid",    "#e8e8e4")

blocks   <- igraph::V(g)$block
blk_size <- sort(table(blocks), decreasing = TRUE)
top3     <- names(blk_size)[seq_len(min(3L, length(blk_size)))]

vcol <- ifelse(blocks %in% top3, SLOT[match(blocks, top3)], NEUTRAL)
vcol[is.na(blocks)] <- NEUTRAL

png(file.path(outdir, "string_network_by_block.png"),
    width = 2400, height = 1900, res = 220)
op <- par(mar = c(0, 0.8, 2.6, 0.8), bg = COL_PAPER)

# Push labels away from the crowded centre: in the dense core a label sitting
# on top of its own node is unreadable, so label offset grows with local degree.
deg  <- igraph::degree(g)
ldeg <- ifelse(deg >= 6, 1.15, 0.85)

plot(g, layout = lay, rescale = FALSE, xlim = c(-1, 1), ylim = c(-1, 1),
     vertex.color       = vcol,
     vertex.frame.color = COL_PAPER,
     vertex.frame.width = 1.4,
     vertex.size        = 7,
     vertex.label       = igraph::V(g)$gene,
     vertex.label.cex   = 0.46,
     vertex.label.color = INK,
     vertex.label.dist  = ldeg,
     vertex.label.degree= -pi/2,
     vertex.label.family= "sans",
     edge.color         = SUBTLE,
     edge.width         = 1.1)

title(sprintf("%s — STRING v%s functional associations (score >= %d), by LD block",
              cfg$phenotype$name, s$version, s$score_threshold),
      cex.main = 0.95, col.main = INK, font.main = 1, adj = 0.5, line = 0.8)

legend("bottomleft", bty = "n", cex = 0.72, pt.cex = 1.5, pch = 21,
       text.col = INK, col = COL_PAPER,
       legend = c(sprintf("%s  (%d genes)", top3, as.integer(blk_size[top3])),
                  sprintf("other blocks  (%d genes)", sum(!blocks %in% top3))),
       pt.bg  = c(SLOT[seq_along(top3)], NEUTRAL))
par(op); dev.off()
message("wrote: ", file.path(outdir, "string_network_by_block.png"))

# ---- one panel per block, same layout --------------------------------------
# Colour alone cannot carry 12 blocks. Small multiples can: identical geometry
# in every panel, so the eye compares position rather than hue.
ublk <- names(blk_size)
n    <- length(ublk)
ncol <- ceiling(sqrt(n)); nrow <- ceiling(n / ncol)

png(file.path(outdir, "string_network_block_facets.png"),
    width = 380 * ncol, height = 380 * nrow + 60, res = 150)
op <- par(mfrow = c(nrow, ncol), mar = c(0.4, 0.4, 2.0, 0.4),
          oma = c(0, 0, 2.6, 0), bg = COL_PAPER)

for (b in ublk) {
  hit <- !is.na(blocks) & blocks == b
  plot(g, layout = lay, rescale = FALSE, xlim = c(-1, 1), ylim = c(-1, 1),
       vertex.color       = ifelse(hit, SLOT[1], COL_OFF),
       vertex.frame.color = ifelse(hit, COL_PAPER, COL_OFF),
       vertex.size        = ifelse(hit, 9, 4),
       vertex.label       = NA,
       edge.color         = COL_EDGE,
       edge.width         = 0.8)
  title(sprintf("%s  —  %d gene%s", b, sum(hit), if (sum(hit) == 1) "" else "s"),
        cex.main = 0.95, col.main = INK, font.main = 1, line = 0.4)
}
if (n < nrow * ncol) for (i in seq_len(nrow * ncol - n)) plot.new()
mtext(sprintf("%s — LD block of origin on the STRING network (same layout in every panel)",
              cfg$phenotype$name),
      outer = TRUE, cex = 0.85, col = INK, adj = 0.5, line = 0.6)
par(op); dev.off()
message("wrote: ", file.path(outdir, "string_network_block_facets.png"))
