#!/usr/bin/env Rscript
# ============================================================
# fig_string_chord_by_block.R
#
# The STRING network as a circular chord diagram, arranged by genomic position.
#
#   sectors   LD blocks, ordered around the circle by GRCh38 position
#             (chromosome, then block region start), sector width proportional
#             to the number of connected genes in the block
#   ribbons   PPI associations between genes of block i and block j
#   tracks    outward from the ribbons: block band, MCL module of each gene,
#             gene symbol, block name, chromosome ideogram
#
# The question it answers: are the CLL loci 49 separate genomic regions, or do
# they converge on shared biology? A chord diagram is the right form for that
# because it makes the between-block mass visible as area rather than asking
# the reader to infer it from a force-directed layout.
#
# WHAT THE RIBBONS ENCODE, AND WHY
# --------------------------------
# Raw edge counts let a big block fake connectivity: B21 carries 39 of the 136
# connected genes, so it collects inter-block edges simply by being large.
# The fix is to compare each pair against a null.
#
# This script carries three weightings (--weight):
#
#   observed    ribbon width = number of PPI associations between the pair.
#   oe_degree   width = observed / expected under a degree-preserving
#               (configuration-model) null: E_ij = K_i * K_j / 2m, where K_i is
#               the summed degree of block i's genes. This is the standard
#               modularity null and it is analytic, so it needs no simulation
#               and no seed — the value is reproducible exactly.
#   oe_size     width = observed / expected under a block-size null,
#               E_ij = n_i * n_j * density. Controls for how many genes a block
#               has, but not for how connected those genes are.
#
# The default is `observed`, and it is deliberate. In the CLL network 106 of
# the 139 linked block pairs share exactly ONE association, and those
# single-edge pairs carry 88.8% of the total O/E mass: a pair of two-gene
# blocks with one edge scores O/E = 70, while B21-B22, which actually share 16
# associations, scores 0.51. Under --weight oe_degree the widest ribbon in the
# figure is therefore a single association, and the densest real link is the
# thinnest line. O/E is the right statistic and the wrong width.
#
# So the null does its work through COLOUR instead: every ribbon is coloured by
# log2 of its degree-preserving O/E, diverging about 1. A thick blue ribbon is
# a pair that has many associations but fewer than its genes' degrees predict —
# exactly the "big block faking it" case, and now visible as such rather than
# silently inflating a width. Use --weight oe_degree when the claim is about
# enrichment per pair rather than about the mass of cross-block wiring, and
# raise --min-edges with it.
#
# Inputs:
#   <stage_dir>/string_02_graph_downstream.rds     (stage 02)
#   <stage_dir>/string_03_modules.tsv              (stage 03)
#   <stage_dir>/string_03_module_summary.tsv       (stage 03)
#   <intermediate>/protein_coding_genes_by_block.tsv   GRCh38 gene coordinates
#
# Outputs (in <figures>/string/):
#   string_chord_by_block.png
#   string_chord_by_block.svg
#   string_chord_by_block_pairs.tsv  every pair drawn, with obs, both expecteds
#
# Usage:
#   Rscript figure_scripts/fig_string_chord_by_block.R config/cll.yml
#   Rscript figure_scripts/fig_string_chord_by_block.R config/cll.yml \
#       --weight oe_degree --min-edges 2
#
#   --size N      edge of the square canvas in px            (default 812)
#   --res  N      nominal dpi; scales all type and rules     (default 104)
#   --weight W    observed | oe_degree | oe_size             (default observed)
#   --min-edges N draw only pairs with at least N observed associations (1)
#   --outfile NM  basename for the two image files  (default string_chord_by_block)
#   --palette NM  article | poster | ppt, or a path to a YAML file of
#                 role: "#hex" pairs. Falls back to $FIGURE_PALETTE, then to
#                 the script's own colours. See R/palettes.R.
#
# Type sizes — all independent, for tuning against a printed proof:
#   --cex-block  N  block names, e.g. B21                    (default 0.78)
#   --cex-gene   N  gene symbols                             (default 0.30)
#   --cex-chr    N  chromosome numbers in the outer band     (default 0.52)
#   --cex-legend N  multiplier on every legend text size     (default 1.0)
#
# Ring depths. The block and gene rings default to AUTO: each is measured to
# the widest label it actually has to carry, at the cex in force, plus a hair
# of padding. That is what stops the labels drifting away from the circle when
# they are small — a fixed depth leaves the slack as a gap, so raising a cex
# moves its ring outward rather than overflowing it. Override only to force a
# specific radial budget; the resolved values print on every run.
#   --track-block N  radial depth of the block-name ring     (default auto)
#   --track-gene  N  radial depth of the gene-symbol ring    (default auto)
#   --track-chr   N  radial depth of the chromosome band     (default 0.032)
#   --no-chr         drop the chromosome band entirely; --track-chr 0 also works
#
# A chromosome is labelled only where its arc can hold the text; any that
# cannot are named on the console, and widening --track-chr or lowering
# --cex-chr will fit them.
# ============================================================

if (!file.exists("R/utils.R")) stop("run from the repository root")
source("R/utils.R")
source("R/string_utils.R")
source("R/palettes.R")
suppressPackageStartupMessages({
  library(igraph); library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)
no_chr <- "--no-chr" %in% args
args   <- args[args != "--no-chr"]
take <- function(flag, default, cast = as.character) {
  if (!(flag %in% args)) return(default)
  i <- which(args == flag)
  v <- cast(args[i + 1]); args <<- args[-c(i, i + 1)]
  v
}
fig_size  <- take("--size",      812L, as.integer)
fig_res   <- take("--res",       104L, as.integer)
weight_by <- take("--weight",    "observed")
min_edges <- take("--min-edges", 1L,   as.integer)
outfile   <- take("--outfile",   "string_chord_by_block")
palette   <- take("--palette",   NULL)
cex_block <- take("--cex-block",  0.78, as.numeric)
cex_gene  <- take("--cex-gene",   0.30, as.numeric)
cex_chr   <- take("--cex-chr",    0.52, as.numeric)
cex_leg   <- take("--cex-legend", 1.00, as.numeric)
trk_block <- take("--track-block", NA_real_, as.numeric)
trk_gene  <- take("--track-gene",  NA_real_, as.numeric)
trk_chr   <- take("--track-chr",   0.032,    as.numeric)

if (any(c(cex_block, cex_gene, cex_chr, cex_leg) <= 0)) {
  stop("--cex-block, --cex-gene, --cex-chr and --cex-legend take positive numbers")
}
# --track-chr 0 is the obvious way to ask for no band, so accept it as one
# rather than handing circlize a zero-depth track.
if (!is.na(trk_chr) && trk_chr <= 0) no_chr <- TRUE

if (!weight_by %in% c("observed", "oe_degree", "oe_size")) {
  stop("--weight must be one of: observed, oe_degree, oe_size")
}
if (anyNA(c(fig_size, fig_res, min_edges)) || fig_size <= 0 || fig_res <= 0) {
  stop("--size, --res and --min-edges take positive integers")
}
if (length(args) != 1L) {
  stop("usage: Rscript figure_scripts/fig_string_chord_by_block.R config/<pheno>.yml",
       " [--size N] [--res N] [--weight W] [--min-edges N] [--outfile NAME]",
       " [--palette NAME]",
       " [--cex-block N] [--cex-gene N] [--cex-chr N] [--cex-legend N]",
       " [--track-block N] [--track-gene N] [--track-chr N] [--no-chr]")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

g_path <- file.path(paths$stage_dir, "string_02_graph_downstream.rds")
if (!file.exists(g_path)) {
  stop("graph not found: ", g_path,
       "\n  Run pipeline/string_pipeline/stage_string_02_components.R first.")
}
g    <- readRDS(g_path)
mod  <- read.delim(file.path(paths$stage_dir, "string_03_modules.tsv"))
msum <- read.delim(file.path(paths$stage_dir, "string_03_module_summary.tsv"))

outdir <- file.path(cfg$paths$figures, "string")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- genes without a block cannot sit in a sector --------------------------
# Named rather than dropped silently, matching stage 02's rule that nothing
# leaves the analysis anonymously.
no_blk <- igraph::V(g)$gene[is.na(igraph::V(g)$block)]
if (length(no_blk)) {
  message("dropping ", length(no_blk), " gene(s) with no LD block: ",
          paste(sort(no_blk), collapse = ", "))
  g <- igraph::induced_subgraph(g, igraph::V(g)[!is.na(igraph::V(g)$block)])
}
blk  <- igraph::V(g)$block
gene <- igraph::V(g)$gene
m    <- igraph::ecount(g)
message("chord input: ", igraph::vcount(g), " genes, ", m, " associations, ",
        length(unique(blk)), " LD blocks")

# ---- GRCh38 coordinates, for the order of everything -----------------------
# The circle is a genome: sectors run chromosome 1 -> X clockwise from twelve
# o'clock, and inside a sector the genes run in start-position order. Reading
# the perimeter is then reading the genome, and a ribbon is visibly a contact
# between two places on it.
bt_path <- file.path(cfg$paths$intermediate, s$input$block_table)
if (!file.exists(bt_path)) stop("block table not found: ", bt_path)
bt <- read.delim(bt_path)
gcol <- s$input$block_gene_col
bcol <- s$input$block_block_col

key   <- paste(bt[[gcol]], bt[[bcol]])
hit   <- match(paste(gene, blk), key)
g_chr <- bt$chromosome_name[hit]
g_pos <- bt$start_position[hit]
if (anyNA(g_pos)) {
  stop("no GRCh38 coordinate for: ",
       paste(gene[is.na(g_pos)], collapse = ", "),
       "\n  every connected gene must appear in ", bt_path)
}

chr_rank <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  ifelse(is.na(n), 100 + match(toupper(x), c("X", "Y", "MT")), n)
}

# one row per block: where it sits in the genome, and how many genes it holds
bo <- data.frame(block = blk, chr = g_chr, pos = g_pos, stringsAsFactors = FALSE)
bs <- do.call(rbind, lapply(split(bo, bo$block), function(d) data.frame(
  block = d$block[1], chr = d$chr[1], pos = min(d$pos), n = nrow(d),
  stringsAsFactors = FALSE)))
bs <- bs[order(chr_rank(bs$chr), bs$pos), ]
sectors <- bs$block
n_i     <- setNames(bs$n, bs$block)

# Chromosomes, for the outer ideogram band. Because the sectors are already in
# genomic order, every chromosome's blocks are contiguous, so a chromosome is
# just a run of adjacent sectors and the band needs no extra ordering logic.
chr_run    <- rle(as.character(bs$chr))
chr_groups <- split(sectors, rep(seq_along(chr_run$lengths), chr_run$lengths))
chr_labels <- chr_run$values
message("chromosomes represented: ", length(chr_labels), " (",
        paste(chr_labels, collapse = " "), ")")

# genes inside each sector, in position order
gene_of <- lapply(sectors, function(b) {
  d <- data.frame(gene = gene[blk == b], pos = g_pos[blk == b],
                  stringsAsFactors = FALSE)
  d[order(d$pos), ]
})
names(gene_of) <- sectors

# ---- observed and expected between-block association counts ----------------
el  <- igraph::as_edgelist(g, names = FALSE)
bi  <- blk[el[, 1]]; bj <- blk[el[, 2]]
obs <- matrix(0, length(sectors), length(sectors),
              dimnames = list(sectors, sectors))
for (k in seq_len(nrow(el))) {
  a <- bi[k]; b <- bj[k]
  if (a != b) { obs[a, b] <- obs[a, b] + 1; obs[b, a] <- obs[b, a] + 1 }
}
within <- sum(bi == bj)
between <- sum(bi != bj)

K   <- vapply(sectors, function(b) sum(igraph::degree(g)[blk == b]), numeric(1))
e_deg <- outer(K, K) / (2 * m)                       # configuration model
dens  <- m / choose(igraph::vcount(g), 2)
e_siz <- outer(as.numeric(n_i), as.numeric(n_i)) * dens
dimnames(e_deg) <- dimnames(e_siz) <- dimnames(obs)

safe_ratio <- function(o, e) { r <- o / e; r[!is.finite(r)] <- 0; r }
oe_deg <- safe_ratio(obs, e_deg)
oe_siz <- safe_ratio(obs, e_siz)

W <- switch(weight_by,
            observed  = obs,
            oe_degree = oe_deg,
            oe_size   = oe_siz)
W[obs < min_edges] <- 0

ut    <- which(upper.tri(obs) & W > 0, arr.ind = TRUE)
pairs <- data.frame(
  block_i  = sectors[ut[, 1]],
  block_j  = sectors[ut[, 2]],
  observed = obs[ut],
  exp_degree = round(e_deg[ut], 4),
  oe_degree  = round(oe_deg[ut], 4),
  exp_size   = round(e_siz[ut], 4),
  oe_size    = round(oe_siz[ut], 4),
  width_weight = W[ut],
  stringsAsFactors = FALSE
)
pairs <- pairs[order(-pairs$observed, -pairs$oe_degree), ]
# Keyed to --outfile: the table lists only the pairs actually drawn, so a
# variant run must not overwrite the default figure's table with a filtered one.
pairs_path <- file.path(outdir, paste0(outfile, "_pairs.tsv"))
write.table(pairs, pairs_path, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("associations: %d within-block, %d between-block (%.1f%% cross)",
                within, between, 100 * between / m))
message(sprintf("drawing %d block pairs (of %d linked) at --weight %s --min-edges %d",
                nrow(pairs), sum(obs[upper.tri(obs)] > 0), weight_by, min_edges))
# The figure carries no title, so the headline it used to print lives here and
# belongs in the poster caption.
message(sprintf(
  "\ncaption line: %d of %d associations (%.0f%%) cross a block boundary — all %d connected blocks link to at least one other.\n",
  between, m, 100 * between / m, sum(rowSums(obs) > 0)))

# ---- colours ---------------------------------------------------------------
# --palette wins; with no flag get_palette() falls back to the FIGURE_PALETTE
# environment variable, and with neither, to the literals below.
pal      <- get_palette(palette)
NEUTRAL  <- pal_col(pal, "neutral", "#898781")
INK      <- pal_col(pal, "ink",     "#0b0b0b")
MUTED    <- pal_col(pal, "muted",   "#4b505a")
COL_PAPER<- pal_col(pal, "paper",   "white")
GRID     <- pal_col(pal, "grid",    "#e8e8e4")
AXIS     <- pal_col(pal, "axis",    "#9aa0ab")
POS      <- pal_col(pal, "pos",     "#7a2738")
NEG      <- pal_col(pal, "neg",     "#00516f")

# Module colours match fig_string_network_by_module.R exactly, so the two
# figures can be read against each other: same module, same hue.
keep  <- sort(msum$module_id[!msum$positional_artifact])
SLOT  <- pal_cat(pal, max(3L, length(keep)), c("#2a78d6", "#eb6834", "#1baf7a"))
mid   <- mod$module_id[match(gene, mod$gene)]
art   <- msum$positional_artifact[match(mid, msum$module_id)]
gene_mod_col <- ifelse(is.na(mid), GRID,
                       ifelse(art, NEUTRAL, SLOT[match(mid, keep)]))
names(gene_mod_col) <- gene

# Ribbon colour: log2 of the degree-preserving O/E, diverging about 1. Capped
# at +/-2 (a 4x departure) so the tail of single-edge pairs cannot own the
# scale; the cap is stated in the legend rather than left implicit.
OE_CAP  <- 2
oe_ramp <- circlize::colorRamp2(c(-OE_CAP, 0, OE_CAP), c(NEG, GRID, POS))
ribbon_col <- function(oe, alpha = 0.62) {
  v <- log2(pmax(oe, 1e-6))
  circlize::add_transparency(oe_ramp(pmax(-OE_CAP, pmin(OE_CAP, v))), 1 - alpha)
}

# alternating band so adjacent sectors stay separable without 38 hues
band_col <- rep(c(MUTED, AXIS), length.out = length(sectors))
names(band_col) <- sectors

# ---- sector width and ribbon geometry --------------------------------------
# Ribbon feet are laid on ONE global scale (kscale x-units per unit weight), so
# a width means the same thing anywhere on the circle. The question is then how
# wide to make a sector.
#
# Gene count alone does not work: B44 holds a single gene carrying 23
# between-block associations, so pinning every sector to its gene count and
# shrinking the scale until B44 fits collapses the whole figure to hairlines
# (kscale 0.037 — B21 would spend 4% of its arc on ribbons).
#
# Total weight alone does not work either: it is what chordDiagram() does, and
# it would give B21 about a tenth of the circle to hold 39 gene labels.
#
# So a sector is as wide as it needs to be for BOTH:
#     L_i = max(n_i, kscale * W_i / FILL)
# and genes are spread evenly across whatever width results. Sector width is
# therefore "gene count, or connectivity where that is the binding constraint",
# which is reported per sector in the console output. Feet are ordered by the
# partner's angular position, which keeps the ribbons from braiding.
FILL  <- 0.86
Wsum  <- rowSums(W)
kscale <- sum(as.numeric(n_i)) / sum(Wsum)   # ribbons ask for as much arc as genes
L <- pmax(as.numeric(n_i), kscale * Wsum / FILL)
names(L) <- sectors
message(sprintf("ribbon scale: %.3f sector-units per unit weight; %d of %d sectors widened past their gene count (%s)",
                kscale, sum(L > as.numeric(n_i) + 1e-9), length(sectors),
                paste(sectors[L > as.numeric(n_i) + 1e-9], collapse = " ")))

ord  <- setNames(seq_along(sectors), sectors)
feet <- list()
for (b in sectors) {
  ptn <- sectors[W[b, sectors] > 0]
  if (!length(ptn)) next
  ptn <- ptn[order(ord[ptn])]
  w   <- kscale * W[b, ptn]
  st  <- (L[[b]] - sum(w)) / 2 + c(0, cumsum(w)[-length(w)])
  feet[[b]] <- data.frame(partner = ptn, x1 = st, x2 = st + w,
                          stringsAsFactors = FALSE)
}

# gene slot i of n, mapped onto the sector's actual width
gx <- function(b, i) (i - 0.5) * L[[b]] / length(gene_of[[b]]$gene)
ge <- function(b, i) c((i - 1), i) * L[[b]] / length(gene_of[[b]]$gene)

# ---- the drawing ------------------------------------------------------------
REPORTED <- NULL; SKIP_REPORTED <- NULL   # console lines print once, not per device

draw_chord <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2), bg = COL_PAPER, xpd = NA)
  circos.clear()
  circos.par(start.degree = 90, clock.wise = TRUE,
             gap.degree = 1.1, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0, 0.004), points.overflow.warning = FALSE,
             canvas.xlim = c(-1.04, 1.04), canvas.ylim = c(-1.22, 1.03))
  circos.initialize(factors = factor(sectors, levels = sectors),
                    xlim = cbind(0, as.numeric(L[sectors])))

  # circos.initialize has now established the coordinate system, so labels can
  # be measured. Both rings are radial, so a ring's depth is the width of the
  # longest string it carries — measured, not guessed, so changing a cex moves
  # the ring with it instead of leaving a gap.
  PAD <- 0.014
  tb <- if (is.na(trk_block)) max(strwidth(sectors, cex = cex_block, units = "user")) + PAD else trk_block
  tg <- if (is.na(trk_gene))  max(strwidth(gene,    cex = cex_gene,  units = "user")) + PAD else trk_gene
  if (is.null(REPORTED)) {
    message(sprintf("type: block cex %.2f in ring %.3f | gene cex %.2f in ring %.3f | chr cex %.2f | legend x%.2f",
                    cex_block, tb, cex_gene, tg, cex_chr, cex_leg))
    REPORTED <<- TRUE
  }

  # ---- chromosome ideogram, outermost -------------------------------------
  # A landmark band so the circle can be read as a genome rather than as 38
  # anonymous sectors. Drawn with highlight.sector() rather than per-sector
  # rectangles because it has to run across the inter-block gaps to read as one
  # chromosome. Alternating ink/muted, which is the ideogram convention and
  # keeps neighbouring chromosomes apart without spending a hue.
  # Skipping it costs nothing structurally: tracks are allocated outside in, so
  # with no band the remaining rings simply start further out.
  if (!no_chr) {
    circos.track(ylim = c(0, 1), track.height = trk_chr,
                 bg.col = NA, bg.border = NA, panel.fun = function(x, y) {})
    chr_fill <- rep(c(INK, MUTED), length.out = length(chr_labels))
    free_deg <- 360 - circos.par("gap.degree") * length(sectors)
    r_mid    <- 1 - trk_chr / 2
    skipped  <- character(0)
    for (i in seq_along(chr_groups)) {
      grp <- chr_groups[[i]]
      highlight.sector(grp, track.index = 1, col = chr_fill[i], border = NA)
      # Label only where the arc can actually hold it. A chromosome carrying one
      # narrow block has no room for tangential text, and an overflowing label
      # would collide with its neighbour; the band tone still marks the grouping,
      # and any skipped chromosome is named on the console.
      ang <- sum(as.numeric(L[grp])) / sum(as.numeric(L)) * free_deg +
             (length(grp) - 1) * circos.par("gap.degree")
      if (r_mid * ang * pi / 180 > strwidth(chr_labels[i], cex = cex_chr, units = "user") * 1.25) {
        highlight.sector(grp, track.index = 1, col = NA, border = NA,
                         text = chr_labels[i], text.col = COL_PAPER,
                         facing = "bending.inside", niceFacing = TRUE, cex = cex_chr)
      } else {
        skipped <- c(skipped, chr_labels[i])
      }
    }
    if (length(skipped) && is.null(SKIP_REPORTED)) {
      message("chromosome band: no room to label chr ", paste(skipped, collapse = ", "),
              " — raise --track-chr or lower --cex-chr to fit them")
      SKIP_REPORTED <<- TRUE
    }
  }

  # block name. Set larger than the gene symbols: the block is the
  # unit the figure is about, and the gene ring is its detail.
  circos.track(ylim = c(0, 1), track.height = tb, bg.border = NA,
               panel.fun = function(x, y) {
                 b <- CELL_META$sector.index
                 circos.text(mean(CELL_META$xlim), 0, b,
                             facing = "clockwise", niceFacing = TRUE,
                             adj = c(0, 0.5), cex = cex_block, col = INK)
               })

  # gene symbols
  circos.track(ylim = c(0, 1), track.height = tg, bg.border = NA,
               panel.fun = function(x, y) {
                 b  <- CELL_META$sector.index
                 gs <- gene_of[[b]]$gene
                 for (i in seq_along(gs)) {
                   circos.text(gx(b, i), 0.04, gs[i],
                               facing = "clockwise", niceFacing = TRUE,
                               adj = c(0, 0.5), cex = cex_gene, col = INK)
                 }
               })

  # MCL module of each gene
  circos.track(ylim = c(0, 1), track.height = 0.035, bg.border = NA,
               panel.fun = function(x, y) {
                 b  <- CELL_META$sector.index
                 gs <- gene_of[[b]]$gene
                 for (i in seq_along(gs)) {
                   e <- ge(b, i)
                   circos.rect(e[1], 0, e[2], 1,
                               col = gene_mod_col[[gs[i]]], border = NA)
                 }
               })

  # the block band itself
  circos.track(ylim = c(0, 1), track.height = 0.022, bg.border = NA,
               panel.fun = function(x, y) {
                 circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                             col = band_col[[CELL_META$sector.index]], border = NA)
               })

  # ribbons, widest first so the thin ones stay visible on top
  o <- order(-pairs$width_weight)
  for (r in o) {
    a <- pairs$block_i[r]; b <- pairs$block_j[r]
    fa <- feet[[a]][feet[[a]]$partner == b, ]
    fb <- feet[[b]][feet[[b]]$partner == a, ]
    if (!nrow(fa) || !nrow(fb)) next
    circos.link(a, c(fa$x1, fa$x2), b, c(fb$x1, fb$x2),
                col = ribbon_col(pairs$oe_degree[r]), border = NA)
  }

  # ---- legend band ---------------------------------------------------------
  # No title: the circle takes the space instead, and the poster caption carries
  # the headline (still printed to the console on every run, see below).
  #
  # Two keys, each with its own heading and a plain-language subtitle. Ribbon
  # width is deliberately NOT keyed — it is the one encoding a reader can guess
  # (thicker = more) and explaining it cost a third of the band.
  #
  # Both sit outboard of x = |0.26|, where the circle has already curved away
  # from the bottom of the canvas, so the band can start at y = -1.005 without
  # colliding with the block labels at six o'clock.
  hd   <- function(x, s) text(x, -1.005, s, cex = 0.42 * cex_leg, col = INK, adj = c(0, 0.5), font = 2)
  sub  <- function(x, s) text(x, -1.064, s, cex = 0.35 * cex_leg, col = MUTED, adj = c(0, 0.5))
  rule <- function(x0, x1) segments(x0, -1.032, x1, -1.032, col = AXIS, lwd = 1.2)

  # left — the gene ring
  lab <- c(sprintf("M%d", keep), "artifact")
  col <- c(SLOT[seq_along(keep)], NEUTRAL)
  x0  <- -1.02
  hd(x0, "GENE RING");  rule(x0, x0 + 0.60)
  sub(x0, "MCL module of each gene")
  for (i in seq_along(lab)) {
    xi <- x0 + ((i - 1) %% 4) * 0.152
    yi <- if (i <= 4) -1.112 else -1.166
    rect(xi, yi - 0.017, xi + 0.032, yi + 0.017, col = col[i], border = NA)
    text(xi + 0.043, yi, lab[i], cex = 0.35 * cex_leg, col = INK, adj = c(0, 0.5))
  }

  # right — ribbon colour
  x2 <- 0.26
  hd(x2, "RIBBON COLOUR"); rule(x2, x2 + 0.76)
  sub(x2, "how that count compares with chance")
  nseg <- 60; bw <- 0.50 / nseg
  for (i in seq_len(nseg)) {
    v <- -OE_CAP + 2 * OE_CAP * (i - 0.5) / nseg
    rect(x2 + (i - 1) * bw, -1.130, x2 + i * bw, -1.098,
         col = circlize::add_transparency(oe_ramp(v), 0.38), border = NA)
  }
  text(x2 + 0.25, -1.114, "1x", cex = 0.32 * cex_leg, col = MUTED, adj = c(0.5, 0.5))
  text(x2,        -1.162, "0.25x  fewer than chance", cex = 0.32 * cex_leg, col = NEG, adj = c(0, 0.5))
  text(x2 + 0.50, -1.162, "more than chance  4x",     cex = 0.32 * cex_leg, col = POS, adj = c(1, 0.5))
  text(x2, -1.198,
       "chance = degree-preserving null; every gene keeps its own number of partners",
       cex = 0.30 * cex_leg, col = MUTED, adj = c(0, 0.5))

  circos.clear()
}

png(file.path(outdir, paste0(outfile, ".png")),
    width = fig_size, height = fig_size, res = fig_res)
draw_chord(); dev.off()
message("wrote: ", file.path(outdir, paste0(outfile, ".png")))

# svglite sizes in inches at a fixed 72 px/in, and has no `res`. To land the
# same pixel box AND the same proportions as the PNG, the canvas is size/72
# inches and the base pointsize is scaled by res/72 — that keeps type the same
# fraction of the canvas in both files.
svg_path <- file.path(outdir, paste0(outfile, ".svg"))
svglite::svglite(svg_path,
                 width = fig_size / 72, height = fig_size / 72,
                 pointsize = 12 * fig_res / 72, bg = COL_PAPER)
draw_chord(); dev.off()

# svglite declares the root width/height in points, so an 812-unit viewBox is
# advertised as 812pt = 1083 px and the file no longer answers to "812 square".
# The viewBox is already right; only the declared unit is restated here.
svg <- readLines(svg_path, warn = FALSE)
svg[1:3] <- sub(sprintf("width='%.2fpt' height='%.2fpt'", fig_size, fig_size),
                sprintf("width='%d' height='%d'", fig_size, fig_size),
                svg[1:3], fixed = TRUE)
writeLines(svg, svg_path)
message("wrote: ", svg_path)
message("wrote: ", pairs_path)
