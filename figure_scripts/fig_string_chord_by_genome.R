#!/usr/bin/env Rscript
# ============================================================
# fig_string_chord_by_genome.R
#
# The STRING network as a circular chord diagram drawn on the WHOLE GENOME.
#
#   sectors   the 22 chromosomes, at their true GRCh38 lengths, in order
#             clockwise from twelve o'clock. Every chromosome is present
#             whether or not it carries an LD block.
#   marks     each LD block at its own genomic position inside its chromosome
#   ribbons   PPI associations between genes of block i and block j
#
# Sibling script, and how to choose
# ---------------------------------
# fig_string_chord_by_block.R gives every block its own sector, so the circle is
# a list of blocks in genomic order and the chromosome band is a grouping drawn
# over them. Chromosomes carrying no block are simply not on that circle, and
# neither is the distance between two blocks that happen to be 100 Mb apart.
#
# This script inverts that: the circle is the genome itself, at scale. The
# question it answers is where in the genome the signal is and, just as much,
# where it is NOT — the empty arc of chromosome 14 is a result, and the block
# view cannot show it. What it costs is per-gene legibility: at 3.1 Gb around a
# circle a 200 kb block is 0.02 degrees, so blocks cannot be drawn to scale (see
# MINIMUM ARC below) and gene symbols only fit on demand (--genes).
#
# Neither is a replacement for the other. Use the block view for the anatomy of
# the network, this one for its genomic distribution.
#
# MINIMUM ARC, AND WHY THE TICKS EXIST
# ------------------------------------
# The blocks in this analysis span 0.04-1.2 Mb. Drawn to scale, the widest of
# them is 0.14 degrees of arc — a hairline that can carry neither a ribbon foot
# nor a label. So each block is given a band of
#
#     span_b = clamp(kscale * W_b, --min-arc, --max-arc)
#
# centred on its true position, where W_b is the block's total ribbon weight and
# kscale is one global constant (--fill sets it: the block bands together take
# that fraction of the genome). Between the two bounds, ribbon width per
# association is the same everywhere on the circle, exactly as in the sibling
# script.
#
# The upper bound is not cosmetic. B21 and B22 are 1 Mb apart on chromosome 6
# and carry 59 and 52 cross-block associations between them; unbounded, they ask
# for 111 Mb of chromosome 6 to hold two loci that are neighbours, and
# everything else on the chromosome gets shoved out of the way — 119 Mb out of
# the way, in the run that produced this default. --max-arc caps a band, and the
# ribbons inside a capped band are scaled down together, so they stay
# proportional to each other but no longer to the rest of the circle. Blocks
# affected are named on the console.
#
# Bands that would still collide are pushed apart along the chromosome, and if a
# chromosome cannot hold its blocks at all its bands are scaled down together
# (both are reported on the console). That displacement is real and the figure
# admits it: a tick on the ideogram marks each block's TRUE position, with a
# leader line to the band it was moved to. Where a block was not moved the tick
# sits under its band and the leader is a plumb line.
#
# The one thing not to read off this figure is a distance between two blocks on
# the same chromosome: their bands are positioned, not scaled. Everything else —
# which chromosome, where along it, how much of the genome is untouched — is to
# scale.
#
# WHAT THE RIBBONS ENCODE
# -----------------------
# Unchanged from fig_string_chord_by_block.R, deliberately, so the two figures
# can be read against each other: width by --weight (observed | oe_degree |
# oe_size), colour always by log2 of the degree-preserving O/E ratio, diverging
# about 1 and capped at 4x. See the sibling's header for why width is observed
# counts by default and the null does its work through colour.
#
# NOTE ON DUPLICATION
# -------------------
# The association counting, the two null models and the module colours are
# copied from fig_string_chord_by_block.R rather than shared through R/, because
# that script is the frozen one that produced the circulated figure and is not
# to be touched. The two must stay in step: a change to the O/E maths or to the
# module colour slots belongs in both files or in neither.
#
# Inputs:
#   <stage_dir>/string_02_graph_downstream.rds     (stage 02)
#   <stage_dir>/string_03_modules.tsv              (stage 03)
#   <stage_dir>/string_03_module_summary.tsv       (stage 03)
#   <intermediate>/protein_coding_genes_by_block.tsv   GRCh38 block + gene coords
#   data/chromosomes/chrom_sizes.txt               GRCh38 chromosome lengths
#
# Outputs (in <figures>/string/):
#   string_chord_by_genome.png
#   string_chord_by_genome.svg
#   string_chord_by_genome_pairs.tsv   every pair drawn, with obs and both expecteds
#   string_chord_by_genome_blocks.tsv  every block, its true position, and the
#                                      band it was drawn at
#
# Usage:
#   Rscript figure_scripts/fig_string_chord_by_genome.R config/cll.yml
#   Rscript figure_scripts/fig_string_chord_by_genome.R config/cll.yml \
#       --genes --min-arc 1.3 --fill 0.2
#
#   --chrsizes P  chromosome length table  (default data/chromosomes/chrom_sizes.txt)
#   --size N      edge of the square canvas in px            (default 900)
#   --res  N      nominal dpi; scales all type and rules     (default 104)
#   --weight W    observed | oe_degree | oe_size             (default observed)
#   --min-edges N draw only pairs with at least N observed associations (1)
#   --outfile NM  basename for the output files   (default string_chord_by_genome)
#   --palette NM  article | poster | ppt, or a path to a YAML file of
#                 role: "#hex" pairs. Falls back to $FIGURE_PALETTE, then to
#                 the script's own colours. See R/palettes.R.
#
# Genome layout:
#   --min-arc D   smallest arc a block band may occupy, in degrees (default 0.9)
#   --max-arc D   largest, before a hub block starts displacing its neighbours
#                 off their own coordinates                        (default 3.0)
#   --fill F      fraction of the genome the block bands together take, which
#                 is what sets ribbon width per association        (default 0.18)
#   --gap D       gap between chromosome sectors, in degrees      (default 0.8)
#   --linked-only draw only blocks that carry a STRING association; the default
#                 draws every block in the block table, so a locus with no
#                 network partner is still visibly a locus
#
# Type sizes:
#   --cex-chr    N  chromosome numbers in the ideogram        (default 0.62)
#   --cex-block  N  block names, e.g. B21                     (default 0.50)
#   --cex-gene   N  gene symbols, with --genes                (default 0.26)
#   --cex-legend N  multiplier on every legend text size      (default 1.0)
#
# Rings (radial depth). The block and gene rings default to AUTO: each is
# measured to the widest label it must carry at the cex in force.
#   --track-chr   N  depth of the chromosome ideogram         (default 0.055)
#   --track-block N  depth of the block-name ring             (default auto)
#   --track-gene  N  depth of the gene-symbol ring            (default auto)
#   --genes          draw the gene-symbol ring. Off by default, and it needs a
#                    bigger canvas to be worth turning on: at --size 900 the 39
#                    symbols of B21 are unreadable inside a 3-degree band. Try
#                    --genes --size 1600 --cex-gene 0.34, or read the gene ring
#                    off fig_string_chord_by_block.R, which has room for it.
# ============================================================

if (!file.exists("R/utils.R")) stop("run from the repository root")
source("R/utils.R")
source("R/string_utils.R")
source("R/palettes.R")
suppressPackageStartupMessages({
  library(igraph); library(circlize)
})

args        <- commandArgs(trailingOnly = TRUE)
show_genes  <- "--genes"       %in% args
linked_only <- "--linked-only" %in% args
args <- args[!args %in% c("--genes", "--linked-only")]
take <- function(flag, default, cast = as.character) {
  if (!(flag %in% args)) return(default)
  i <- which(args == flag)
  v <- cast(args[i + 1]); args <<- args[-c(i, i + 1)]
  v
}
chrsizes  <- take("--chrsizes",  "data/chromosomes/chrom_sizes.txt")
fig_size  <- take("--size",      900L, as.integer)
fig_res   <- take("--res",       104L, as.integer)
weight_by <- take("--weight",    "observed")
min_edges <- take("--min-edges", 1L,   as.integer)
outfile   <- take("--outfile",   "string_chord_by_genome")
palette   <- take("--palette",   NULL)
min_arc   <- take("--min-arc",   0.90, as.numeric)
max_arc   <- take("--max-arc",   3.00, as.numeric)
fill_frac <- take("--fill",      0.18, as.numeric)
gap_deg   <- take("--gap",       0.80, as.numeric)
cex_chr   <- take("--cex-chr",    0.62, as.numeric)
cex_block <- take("--cex-block",  0.50, as.numeric)
cex_gene  <- take("--cex-gene",   0.26, as.numeric)
cex_leg   <- take("--cex-legend", 1.00, as.numeric)
trk_chr   <- take("--track-chr",   0.055,    as.numeric)
trk_block <- take("--track-block", NA_real_, as.numeric)
trk_gene  <- take("--track-gene",  NA_real_, as.numeric)

if (any(c(cex_chr, cex_block, cex_gene, cex_leg) <= 0)) {
  stop("--cex-chr, --cex-block, --cex-gene and --cex-legend take positive numbers")
}
if (!weight_by %in% c("observed", "oe_degree", "oe_size")) {
  stop("--weight must be one of: observed, oe_degree, oe_size")
}
if (anyNA(c(fig_size, fig_res, min_edges)) || fig_size <= 0 || fig_res <= 0) {
  stop("--size, --res and --min-edges take positive integers")
}
if (is.na(min_arc) || min_arc <= 0 || min_arc > 10) {
  stop("--min-arc is a positive number of degrees (a block band; 0.9 is the default)")
}
if (is.na(max_arc) || max_arc < min_arc || max_arc > 30) {
  stop("--max-arc is a number of degrees, at least --min-arc (", min_arc, ")")
}
if (is.na(fill_frac) || fill_frac <= 0 || fill_frac >= 1) {
  stop("--fill is the fraction of the genome the block bands take, in (0, 1)")
}
if (is.na(gap_deg) || gap_deg < 0 || gap_deg > 5) {
  stop("--gap is the gap between chromosome sectors in degrees, in [0, 5]")
}
if (length(args) != 1L) {
  stop("usage: Rscript figure_scripts/fig_string_chord_by_genome.R config/<pheno>.yml",
       " [--chrsizes PATH] [--size N] [--res N] [--weight W] [--min-edges N]",
       " [--outfile NAME] [--palette NAME] [--min-arc D] [--fill F] [--gap D]",
       " [--linked-only] [--genes] [--cex-chr N] [--cex-block N] [--cex-gene N]",
       " [--cex-legend N] [--track-chr N] [--track-block N] [--track-gene N]")
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

# ---- genes without a block cannot be placed on the genome -------------------
no_blk <- igraph::V(g)$gene[is.na(igraph::V(g)$block)]
if (length(no_blk)) {
  message("dropping ", length(no_blk), " gene(s) with no LD block: ",
          paste(sort(no_blk), collapse = ", "))
  g <- igraph::induced_subgraph(g, igraph::V(g)[!is.na(igraph::V(g)$block)])
}
blk  <- igraph::V(g)$block
gene <- igraph::V(g)$gene
m    <- igraph::ecount(g)

# ---- the genome: every chromosome, at its GRCh38 length ---------------------
# The point of this figure. Sector width is chromosome length, so an empty
# chromosome still spends its share of the circle and reads as empty.
if (!file.exists(chrsizes)) {
  stop("chromosome sizes not found: ", chrsizes,
       "\n  expected a two-column table: chr <tab> size (see data/chromosomes/)")
}
cs <- read.delim(chrsizes, sep = "\t", header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("chr", "size") %in% names(cs))) {
  stop("chromosome sizes file must have columns 'chr' and 'size': ", chrsizes)
}
# Sizes arrive quoted and comma-grouped ("248,956,422"); same cleaning as
# fig_karyoploter_trait_ld_blocks.R.
cs$size <- as.numeric(gsub('[",\r]', "", trimws(cs$size)))
cs$chr  <- sub("^chr", "", trimws(cs$chr), ignore.case = TRUE)
if (anyNA(cs$size)) stop("unparseable size(s) in ", chrsizes, " for chr ",
                         paste(cs$chr[is.na(cs$size)], collapse = ", "))

chr_rank <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  ifelse(is.na(n), 100 + match(toupper(x), c("X", "Y", "MT")), n)
}
cs <- cs[order(chr_rank(cs$chr)), ]
genome_bp <- sum(cs$size)

# ---- blocks, and where they sit on it ---------------------------------------
# region_start/region_end are the LD block's own interval, which is what belongs
# on a genome figure — the gene coordinates inside it are detail. A block that
# was merged from several intervals is taken from its outermost bounds.
bt_path <- file.path(cfg$paths$intermediate, s$input$block_table)
if (!file.exists(bt_path)) stop("block table not found: ", bt_path)
bt <- read.delim(bt_path)
gcol <- s$input$block_gene_col
bcol <- s$input$block_block_col
for (col in c("chromosome", "region_start", "region_end")) {
  if (!col %in% names(bt)) {
    stop("block table is missing the '", col, "' column: ", bt_path,
         "\n  this figure needs the LD block intervals, not only gene positions")
  }
}

bl <- do.call(rbind, lapply(split(bt, bt[[bcol]]), function(d) data.frame(
  block = d[[bcol]][1],
  chr   = as.character(d$chromosome[1]),
  lo    = min(d$region_start),
  hi    = max(d$region_end),
  n_genes_table = nrow(d),
  stringsAsFactors = FALSE)))
bl$chr <- sub("^chr", "", bl$chr)
bl$mid <- (bl$lo + bl$hi) / 2

# Blocks the network actually reaches. The rest are drawn (unless --linked-only)
# as a locus with no partner, which is a finding and not a gap in the data.
n_linked_gene <- table(blk)
bl$n_genes <- as.integer(n_linked_gene[bl$block]); bl$n_genes[is.na(bl$n_genes)] <- 0L
bl$linked  <- bl$block %in% blk
if (linked_only) bl <- bl[bl$linked, , drop = FALSE]
bl <- bl[order(chr_rank(bl$chr), bl$mid), ]

missing_chr <- setdiff(bl$chr, cs$chr)
if (length(missing_chr)) {
  stop("block(s) on chromosome(s) absent from ", chrsizes, ": ",
       paste(sort(unique(missing_chr)), collapse = ", "),
       "\n  add the row(s) to that file, or the block cannot be placed")
}
over <- bl$hi > cs$size[match(bl$chr, cs$chr)]
if (any(over)) {
  stop("block(s) past the end of their chromosome: ",
       paste(bl$block[over], collapse = ", "),
       "\n  the block table and ", chrsizes, " are on different assemblies")
}

sectors <- bl$block                       # ribbon bookkeeping is still per block
carried <- cs$chr[cs$chr %in% bl$chr]     # in chromosome order, not lexicographic
empty   <- cs$chr[!cs$chr %in% bl$chr]
message(sprintf("genome: %d chromosomes, %.2f Gb", nrow(cs), genome_bp / 1e9))
message("blocks: ", nrow(bl), " placed (", sum(bl$linked), " with a STRING association",
        if (!linked_only && any(!bl$linked))
          paste0(", ", sum(!bl$linked), " without — drawn as empty loci") else "", ")")
message("chromosomes carrying a block: ", length(carried), " (",
        paste(carried, collapse = " "), ")")
message("chromosomes with none: ",
        if (length(empty)) paste(empty, collapse = " ") else "none")
message("chord input: ", igraph::vcount(g), " genes, ", m, " associations")

# genes of each block, in position order, for the module ring and --genes
g_hit  <- match(paste(gene, blk), paste(bt[[gcol]], bt[[bcol]]))
if (anyNA(g_hit)) {
  stop("no GRCh38 coordinate for: ", paste(gene[is.na(g_hit)], collapse = ", "),
       "\n  every connected gene must appear in ", bt_path)
}
g_pos <- bt$start_position[g_hit]
gene_of <- lapply(bl$block, function(b) {
  d <- data.frame(gene = gene[blk == b], pos = g_pos[blk == b],
                  stringsAsFactors = FALSE)
  d[order(d$pos), ]
})
names(gene_of) <- bl$block

# ---- observed and expected between-block association counts ----------------
# Identical to fig_string_chord_by_block.R; see the note on duplication above.
el  <- igraph::as_edgelist(g, names = FALSE)
bi  <- blk[el[, 1]]; bj <- blk[el[, 2]]
obs <- matrix(0, length(sectors), length(sectors),
              dimnames = list(sectors, sectors))
for (k in seq_len(nrow(el))) {
  a <- bi[k]; b <- bj[k]
  if (a != b && a %in% sectors && b %in% sectors) {
    obs[a, b] <- obs[a, b] + 1; obs[b, a] <- obs[b, a] + 1
  }
}
within  <- sum(bi == bj)
between <- sum(bi != bj)

deg <- igraph::degree(g)
K   <- vapply(sectors, function(b) sum(deg[blk == b]), numeric(1))
n_i <- setNames(as.numeric(bl$n_genes), bl$block)
e_deg <- outer(K, K) / (2 * m)                       # configuration model
dens  <- m / choose(igraph::vcount(g), 2)
e_siz <- outer(n_i, n_i) * dens
dimnames(e_deg) <- dimnames(e_siz) <- dimnames(obs)

safe_ratio <- function(o, e) { r <- o / e; r[!is.finite(r)] <- 0; r }
oe_deg <- safe_ratio(obs, e_deg)
oe_siz <- safe_ratio(obs, e_siz)

W <- switch(weight_by, observed = obs, oe_degree = oe_deg, oe_size = oe_siz)
W[obs < min_edges] <- 0

ut    <- which(upper.tri(obs) & W > 0, arr.ind = TRUE)
pairs <- data.frame(
  block_i  = sectors[ut[, 1]],
  block_j  = sectors[ut[, 2]],
  chr_i    = bl$chr[ut[, 1]],
  chr_j    = bl$chr[ut[, 2]],
  observed = obs[ut],
  exp_degree = round(e_deg[ut], 4),
  oe_degree  = round(oe_deg[ut], 4),
  exp_size   = round(e_siz[ut], 4),
  oe_size    = round(oe_siz[ut], 4),
  width_weight = W[ut],
  stringsAsFactors = FALSE
)
pairs$same_chr <- pairs$chr_i == pairs$chr_j
pairs <- pairs[order(-pairs$observed, -pairs$oe_degree), ]
pairs_path <- file.path(outdir, paste0(outfile, "_pairs.tsv"))
write.table(pairs, pairs_path, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("associations: %d within-block, %d between-block (%.1f%% cross)",
                within, between, 100 * between / m))
message(sprintf("drawing %d block pairs (of %d linked) at --weight %s --min-edges %d",
                nrow(pairs), sum(obs[upper.tri(obs)] > 0), weight_by, min_edges))
if (nrow(pairs)) {
  message(sprintf("of those, %d join blocks on the same chromosome and %d cross chromosomes",
                  sum(pairs$same_chr), sum(!pairs$same_chr)))
}

# ---- block bands: width from weight, floored, then packed -------------------
# One global constant turns ribbon weight into base pairs, so a foot means the
# same thing anywhere on the circle. --fill sets it by fixing the total: all
# block bands together take that fraction of the genome.
Wsum     <- rowSums(W)
floor_bp <- min_arc / 360 * genome_bp
ceil_bp  <- max_arc / 360 * genome_bp

# --fill fixes the total: all bands together take that share of the genome. On a
# sparse network that share is divided among few blocks and nearly every band
# hits --max-arc, which flattens the width encoding to "all the same". So the
# scale is also held to the point where only the top tenth of blocks is capped —
# whichever of the two rules gives the thinner ribbon wins. leptin (21 cross-block
# associations) is set by the second, CLL (214) by the first.
pos_w  <- Wsum[Wsum > 0]
k_fill <- if (sum(Wsum) > 0) fill_frac * genome_bp / sum(Wsum) else 0
k_q90  <- if (length(pos_w)) ceil_bp / stats::quantile(pos_w, 0.9, names = FALSE) else 0
kscale <- min(k_fill, k_q90)
raw    <- kscale * Wsum[bl$block]
want     <- pmin(ceil_bp, pmax(floor_bp, raw))
capped   <- bl$block[raw > ceil_bp]
gap_bp   <- 0.25 * floor_bp     # daylight between adjacent bands, and their labels

# Push overlapping bands apart within a chromosome, keeping them inside it and
# in genomic order, and move them as little as possible while doing it.
#
# "As little as possible" is the whole job, because every base pair of
# displacement is the figure lying about where a block is. Sweeping left to
# right and shoving each band clear of its predecessor satisfies the
# constraints but dumps the entire cost on the rightmost block of a cluster:
# on chromosome 6 that put B23 40 Mb from home while B21, the block that
# actually needed the room, did not move at all.
#
# The least-squares placement is available exactly. Minimising sum (c_i - mid_i)^2
# subject to c_(i+1) - c_i >= (w_i + w_(i+1))/2 + gap becomes, after subtracting
# the cumulative minimum spacing from each target, a plain isotonic regression:
# fit a non-decreasing sequence to the shifted targets, add the spacing back.
# stats::isoreg is that fit. A cluster then opens symmetrically about its own
# centre of mass, and blocks with room to spare do not move at all.
pack_spans <- function(mid, want, size, gap) {
  n <- length(mid)
  need <- sum(want) + gap * max(0, n - 1)
  shrink <- if (need > size) size / need else 1   # the chromosome is the hard bound
  w <- want * shrink
  if (n == 1L) {
    cen <- min(max(mid, w / 2), size - w / 2)
    return(list(lo = cen - w / 2, hi = cen + w / 2, shrink = shrink))
  }

  spacing <- (w[-n] + w[-1]) / 2 + gap            # centre-to-centre minimum
  offs    <- c(0, cumsum(spacing))

  # The ends of the chromosome are two more constraints, and they cannot be
  # applied by sliding the finished solution: on chromosome 19, B48 sits 250 kb
  # from B47 and gets pushed to the telomere, and shifting every band left to
  # keep it inside put B45 ten megabases off the START of the chromosome —
  # outside its own sector, where circlize quietly draws the ribbon anyway.
  #
  # In the transformed coordinates they are bounds on the monotone sequence,
  # and each has to be made monotone itself (a running max from the left, a
  # running min from the right) before it can bind. Clipping the isotonic fit
  # to monotone bounds leaves it monotone and is still the closest feasible
  # placement; `shrink` above is what guarantees a feasible one exists.
  lower <- cummax(w / 2 - offs)
  upper <- rev(cummin(rev(size - w / 2 - offs)))
  cen   <- pmin(pmax(stats::isoreg(mid - offs)$yf, lower), upper) + offs
  list(lo = cen - w / 2, hi = cen + w / 2, shrink = shrink)
}

bl$band_lo <- NA_real_; bl$band_hi <- NA_real_
chr_shrink <- setNames(rep(1, nrow(cs)), cs$chr)
for (c in unique(bl$chr)) {
  idx <- which(bl$chr == c)                     # already in position order
  p   <- pack_spans(bl$mid[idx], want[idx], cs$size[match(c, cs$chr)], gap_bp)
  bl$band_lo[idx] <- p$lo; bl$band_hi[idx] <- p$hi
  chr_shrink[[c]] <- p$shrink
}
bl$band_mid <- (bl$band_lo + bl$band_hi) / 2
bl$shift_mb <- (bl$band_mid - bl$mid) / 1e6

# A band outside its chromosome is a ribbon foot outside its sector, which
# circlize draws anyway with only a note on the console. Catch it here, where
# the cause is still visible, rather than in a figure nobody re-reads.
esc <- which(bl$band_lo < -1 | bl$band_hi > cs$size[match(bl$chr, cs$chr)] + 1)
if (length(esc)) {
  stop("band(s) placed outside their chromosome: ",
       paste(sprintf("%s (chr%s, %.1f-%.1f Mb of %.1f)", bl$block[esc], bl$chr[esc],
                     bl$band_lo[esc] / 1e6, bl$band_hi[esc] / 1e6,
                     cs$size[match(bl$chr[esc], cs$chr)] / 1e6), collapse = ", "),
       "\n  this is a packing bug, not a data problem — lowering --min-arc works around it")
}

moved <- which(abs(bl$shift_mb) > 0.5)
message(sprintf("block bands: %.1f-%.1f deg (%.0f-%.0f Mb), %.2f Mb of arc per association (%s); %d of %d set by weight, %d floored, %d capped",
                min_arc, max_arc, floor_bp / 1e6, ceil_bp / 1e6, kscale / 1e6,
                if (k_fill <= k_q90) sprintf("--fill %.2f", fill_frac) else "held down to cap only the top tenth",
                sum(raw >= floor_bp & raw <= ceil_bp), nrow(bl),
                sum(raw < floor_bp), length(capped)))
if (length(capped)) {
  message("  capped at --max-arc, so their ribbons are compressed relative to the rest: ",
          paste(capped, collapse = " "))
}
if (length(moved)) {
  message(sprintf("  %d band(s) pushed off their true position to avoid a collision (max %.1f Mb): %s",
                  length(moved), max(abs(bl$shift_mb)),
                  paste(bl$block[moved], collapse = " ")))
  message("  ticks on the ideogram mark the true positions; leaders join them to the bands")
}
if (any(chr_shrink < 1)) {
  sh <- names(chr_shrink)[chr_shrink < 1]
  message(sprintf("  chr %s could not hold its bands at full width; they were scaled to %s of it",
                  paste(sh, collapse = " "),
                  paste(sprintf("%.0f%%", 100 * chr_shrink[sh]), collapse = " ")))
}

blocks_path <- file.path(outdir, paste0(outfile, "_blocks.tsv"))
write.table(
  data.frame(block = bl$block, chr = bl$chr,
             region_start = bl$lo, region_end = bl$hi,
             genes_in_block = bl$n_genes_table, genes_in_network = bl$n_genes,
             linked = bl$linked, ribbon_weight = as.numeric(Wsum[bl$block]),
             band_start = round(bl$band_lo), band_end = round(bl$band_hi),
             displaced_mb = round(bl$shift_mb, 3), stringsAsFactors = FALSE),
  blocks_path, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf(
  "\ncaption line: %d of %d associations (%.0f%%) cross a block boundary, spread over %d of the %d chromosomes; %s carry none.\n",
  between, m, 100 * between / m, length(carried), nrow(cs),
  if (length(empty)) paste0("chromosomes ", paste(empty, collapse = ", ")) else "none"))

# ---- colours ---------------------------------------------------------------
pal      <- get_palette(palette)
NEUTRAL  <- pal_col(pal, "neutral", "#898781")
INK      <- pal_col(pal, "ink",     "#0b0b0b")
MUTED    <- pal_col(pal, "muted",   "#4b505a")
COL_PAPER<- pal_col(pal, "paper",   "white")
GRID     <- pal_col(pal, "grid",    "#e8e8e4")
AXIS     <- pal_col(pal, "axis",    "#9aa0ab")
POS      <- pal_col(pal, "pos",     "#7a2738")
NEG      <- pal_col(pal, "neg",     "#00516f")

# Module colours match fig_string_network_by_module.R and the sibling chord
# exactly, so all three can be read against each other: same module, same hue.
keep  <- sort(msum$module_id[!msum$positional_artifact])
SLOT  <- pal_cat(pal, max(3L, length(keep)), c("#2a78d6", "#eb6834", "#1baf7a"))
mid_of<- mod$module_id[match(gene, mod$gene)]
art   <- msum$positional_artifact[match(mid_of, msum$module_id)]
gene_mod_col <- ifelse(is.na(mid_of), GRID,
                       ifelse(art, NEUTRAL, SLOT[match(mid_of, keep)]))
names(gene_mod_col) <- gene

OE_CAP  <- 2
oe_ramp <- circlize::colorRamp2(c(-OE_CAP, 0, OE_CAP), c(NEG, GRID, POS))
ribbon_col <- function(oe, alpha = 0.62) {
  v <- log2(pmax(oe, 1e-6))
  circlize::add_transparency(oe_ramp(pmax(-OE_CAP, pmin(OE_CAP, v))), 1 - alpha)
}

# ---- ribbon feet inside each band ------------------------------------------
# Feet keep the global scale (kscale bp per unit weight) and are centred in the
# band, so a band floored to the minimum arc shows its ribbons narrower than
# itself rather than stretched to fill it. Where a chromosome was scaled down,
# its feet scale with it. Feet are ordered by the partner's position around the
# circle, which keeps ribbons from braiding.
chr_start <- setNames(c(0, cumsum(cs$size)[-nrow(cs)]), cs$chr)   # genome coords
gpos      <- setNames(chr_start[bl$chr] + bl$band_mid, bl$block)  # for ordering
ord       <- setNames(rank(gpos), bl$block)

feet <- list()
for (b in bl$block) {
  ptn <- bl$block[W[b, bl$block] > 0]
  if (!length(ptn)) next
  ptn <- ptn[order(ord[ptn])]
  i   <- match(b, bl$block)
  bw  <- bl$band_hi[i] - bl$band_lo[i]
  w   <- kscale * W[b, ptn] * chr_shrink[[bl$chr[i]]]
  if (sum(w) > bw) w <- w * bw / sum(w)         # never wider than the band
  st  <- bl$band_lo[i] + (bw - sum(w)) / 2 + c(0, cumsum(w)[-length(w)])
  feet[[b]] <- data.frame(partner = ptn, x1 = st, x2 = st + w,
                          stringsAsFactors = FALSE)
}

# gene slot i of n inside a block's band
gx <- function(b, i) {
  j <- match(b, bl$block); n <- nrow(gene_of[[b]])
  bl$band_lo[j] + (i - 0.5) * (bl$band_hi[j] - bl$band_lo[j]) / n
}
ge <- function(b, i) {
  j <- match(b, bl$block); n <- nrow(gene_of[[b]])
  bl$band_lo[j] + c(i - 1, i) * (bl$band_hi[j] - bl$band_lo[j]) / n
}

blocks_of <- split(bl, factor(bl$chr, levels = cs$chr))   # empty chr -> 0 rows

# ---- the drawing ------------------------------------------------------------
REPORTED <- NULL; SKIP_REPORTED <- NULL

draw_chord <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2), bg = COL_PAPER, xpd = NA)
  circos.clear()
  circos.par(start.degree = 90, clock.wise = TRUE,
             gap.degree = gap_deg, cell.padding = c(0, 0, 0, 0),
             track.margin = c(0, 0.004), points.overflow.warning = FALSE,
             canvas.xlim = c(-1.04, 1.04), canvas.ylim = c(-1.26, 1.03))
  circos.initialize(factors = factor(cs$chr, levels = cs$chr),
                    xlim = cbind(0, cs$size))

  PAD <- 0.014
  tb <- if (is.na(trk_block)) max(strwidth(bl$block, cex = cex_block, units = "user")) + PAD else trk_block
  tg <- if (!show_genes) 0 else
        if (is.na(trk_gene)) max(strwidth(gene, cex = cex_gene, units = "user")) + PAD else trk_gene
  if (is.null(REPORTED)) {
    message(sprintf("type: chr cex %.2f in ring %.3f | block cex %.2f in ring %.3f | gene ring %s | legend x%.2f",
                    cex_chr, trk_chr, cex_block, tb,
                    if (show_genes) sprintf("cex %.2f in %.3f", cex_gene, tg) else "off",
                    cex_leg))
    REPORTED <<- TRUE
  }

  # ---- chromosome ideogram, outermost --------------------------------------
  # Alternating ink/muted, the ideogram convention, which separates neighbours
  # without spending a hue. Every chromosome is drawn, which is the point.
  chr_fill <- rep(c(INK, MUTED), length.out = nrow(cs))
  names(chr_fill) <- cs$chr
  skipped <- character(0)
  circos.track(ylim = c(0, 1), track.height = trk_chr, bg.border = NA,
               panel.fun = function(x, y) {
                 ch <- CELL_META$sector.index
                 circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                             col = chr_fill[[ch]], border = NA)
                 # Tangential where the arc can hold the label, radial where it
                 # cannot. Every chromosome gets named either way — an unnamed
                 # sector is exactly the confusion this figure exists to remove.
                 arc <- CELL_META$cell.end.degree - CELL_META$cell.start.degree
                 room <- abs(arc) * pi / 180 * (1 - trk_chr / 2)
                 if (room > strwidth(ch, cex = cex_chr, units = "user") * 1.3) {
                   circos.text(mean(CELL_META$xlim), 0.5, ch, col = COL_PAPER,
                               facing = "bending.inside", niceFacing = TRUE,
                               cex = cex_chr)
                 } else {
                   circos.text(mean(CELL_META$xlim), 0.5, ch, col = COL_PAPER,
                               facing = "clockwise", niceFacing = TRUE,
                               cex = cex_chr)
                   skipped <<- c(skipped, ch)
                 }
               })

  # ---- true position, and the leader to the band it was drawn at -----------
  # The honesty track. A tick at the block's real coordinate, and a straight
  # leader down to the centre of its band. Undisplaced blocks show a plumb line
  # and nothing else; the more a band had to move, the more the leader leans.
  circos.track(ylim = c(0, 1), track.height = 0.030, bg.border = NA,
               panel.fun = function(x, y) {
                 d <- blocks_of[[CELL_META$sector.index]]
                 if (!nrow(d)) return(invisible(NULL))
                 for (i in seq_len(nrow(d))) {
                   col_i <- if (d$linked[i]) INK else NEUTRAL
                   circos.segments(d$mid[i], 1, d$mid[i], 0.62, col = col_i, lwd = 0.7)
                   circos.segments(d$mid[i], 0.62, d$band_mid[i], 0,
                                   col = col_i, lwd = 0.5, straight = TRUE)
                 }
               })

  # ---- block names ---------------------------------------------------------
  circos.track(ylim = c(0, 1), track.height = tb, bg.border = NA,
               panel.fun = function(x, y) {
                 d <- blocks_of[[CELL_META$sector.index]]
                 if (!nrow(d)) return(invisible(NULL))
                 for (i in seq_len(nrow(d))) {
                   circos.text(d$band_mid[i], 0, d$block[i],
                               facing = "clockwise", niceFacing = TRUE,
                               adj = c(0, 0.5), cex = cex_block,
                               col = if (d$linked[i]) INK else NEUTRAL)
                 }
               })

  # ---- gene symbols, on demand --------------------------------------------
  if (show_genes) {
    circos.track(ylim = c(0, 1), track.height = tg, bg.border = NA,
                 panel.fun = function(x, y) {
                   d <- blocks_of[[CELL_META$sector.index]]
                   if (!nrow(d)) return(invisible(NULL))
                   for (i in seq_len(nrow(d))) {
                     gs <- gene_of[[d$block[i]]]$gene
                     for (j in seq_along(gs)) {
                       circos.text(gx(d$block[i], j), 0.04, gs[j],
                                   facing = "clockwise", niceFacing = TRUE,
                                   adj = c(0, 0.5), cex = cex_gene, col = INK)
                     }
                   }
                 })
  }

  # ---- MCL module of each gene in the block -------------------------------
  circos.track(ylim = c(0, 1), track.height = 0.030, bg.border = NA,
               panel.fun = function(x, y) {
                 d <- blocks_of[[CELL_META$sector.index]]
                 if (!nrow(d)) return(invisible(NULL))
                 for (i in seq_len(nrow(d))) {
                   gs <- gene_of[[d$block[i]]]$gene
                   if (!length(gs)) {   # a block with no network gene: hollow
                     circos.rect(d$band_lo[i], 0.25, d$band_hi[i], 0.75,
                                 col = NA, border = NEUTRAL, lwd = 0.5)
                     next
                   }
                   for (j in seq_along(gs)) {
                     e <- ge(d$block[i], j)
                     circos.rect(e[1], 0, e[2], 1,
                                 col = gene_mod_col[[gs[j]]], border = NA)
                   }
                 }
               })

  # ---- the block band itself, where the ribbons land -----------------------
  band_col <- rep(c(MUTED, AXIS), length.out = nrow(bl))
  names(band_col) <- bl$block
  circos.track(ylim = c(0, 1), track.height = 0.022, bg.border = NA,
               panel.fun = function(x, y) {
                 d <- blocks_of[[CELL_META$sector.index]]
                 if (!nrow(d)) return(invisible(NULL))
                 for (i in seq_len(nrow(d))) {
                   circos.rect(d$band_lo[i], 0, d$band_hi[i], 1,
                               col = if (d$linked[i]) band_col[[d$block[i]]] else GRID,
                               border = NA)
                 }
               })

  # ---- ribbons, widest first so the thin ones stay visible on top ----------
  o <- order(-pairs$width_weight)
  for (r in o) {
    a <- pairs$block_i[r]; b <- pairs$block_j[r]
    fa <- feet[[a]][feet[[a]]$partner == b, ]
    fb <- feet[[b]][feet[[b]]$partner == a, ]
    if (!nrow(fa) || !nrow(fb)) next
    ca <- bl$chr[match(a, bl$block)]; cb <- bl$chr[match(b, bl$block)]
    circos.link(ca, c(fa$x1, fa$x2), cb, c(fb$x1, fb$x2),
                col = ribbon_col(pairs$oe_degree[r]), border = NA)
  }

  if (length(skipped) && is.null(SKIP_REPORTED)) {
    message("chromosome band: chr ", paste(unique(skipped), collapse = ", "),
            " labelled radially — too narrow for tangential text at --cex-chr ",
            sprintf("%.2f", cex_chr))
    SKIP_REPORTED <<- TRUE
  }

  # ---- legend band ---------------------------------------------------------
  # Same two keys as the sibling figure, plus the one line this variant owes the
  # reader: the bands are not to scale, and the ticks are where the blocks are.
  hd   <- function(x, s) text(x, -1.005, s, cex = 0.42 * cex_leg, col = INK, adj = c(0, 0.5), font = 2)
  sub  <- function(x, s) text(x, -1.064, s, cex = 0.35 * cex_leg, col = MUTED, adj = c(0, 0.5))
  rule <- function(x0, x1) segments(x0, -1.032, x1, -1.032, col = AXIS, lwd = 1.2)

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

  text(x0, -1.222, sprintf(
    "Chromosomes are to scale (GRCh38). Block bands are not: each is drawn at least %.1f degrees wide, and the tick on the ideogram is its true position.",
    min_arc), cex = 0.30 * cex_leg, col = MUTED, adj = c(0, 0.5))

  circos.clear()
}

png(file.path(outdir, paste0(outfile, ".png")),
    width = fig_size, height = fig_size, res = fig_res)
draw_chord(); dev.off()
message("wrote: ", file.path(outdir, paste0(outfile, ".png")))

# svglite sizes in inches at a fixed 72 px/in and has no `res`; the canvas is
# size/72 inches with the pointsize scaled by res/72, which lands the same pixel
# box and the same proportions as the PNG. (Same treatment as the sibling.)
svg_path <- file.path(outdir, paste0(outfile, ".svg"))
svglite::svglite(svg_path,
                 width = fig_size / 72, height = fig_size / 72,
                 pointsize = 12 * fig_res / 72, bg = COL_PAPER)
draw_chord(); dev.off()

svg <- readLines(svg_path, warn = FALSE)
svg[1:3] <- sub(sprintf("width='%.2fpt' height='%.2fpt'", fig_size, fig_size),
                sprintf("width='%d' height='%d'", fig_size, fig_size),
                svg[1:3], fixed = TRUE)
writeLines(svg, svg_path)
message("wrote: ", svg_path)
message("wrote: ", pairs_path)
message("wrote: ", blocks_path)
