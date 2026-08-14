#!/usr/bin/env Rscript
# ============================================================
# fig_string_network_by_category.R
#
# Trait categories projected onto the STRING network and onto the MCL modules.
#
# READ THIS BEFORE INTERPRETING ANYTHING HERE
# -------------------------------------------
# Trait categories attach to LD BLOCKS, not to genes. Lead SNPs sit in introns
# and intergenic space, and the genes were selected by positional overlap with
# the block region, so no single gene owns a trait. Every gene therefore
# INHERITS the full category profile of its block, and every gene in a block
# carries an identical profile.
#
# Consequence: two block-pure modules from the same block have identical
# category profiles by construction. That is not a finding, it is the data
# model. The heatmap below shows it directly (M1/M6/M8 are all pure B2).
# Only modules spanning several blocks have a profile of their own.
#
# Outputs:
#   string_network_category_facets.png   one panel per category, on the network
#   string_module_category_heatmap.png   module x category profile
#   string_04_module_category_profile.tsv
#
# Usage:
#   Rscript figure_scripts/fig_string_network_by_category.R config/leptin.yml
#   Rscript figure_scripts/fig_string_network_by_category.R config/leptin.yml --top 16
# ============================================================

if (!file.exists("R/utils.R")) stop("run from the repository root")
source("R/utils.R")
source("R/string_utils.R")
source("R/palettes.R")
suppressPackageStartupMessages({ library(igraph); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
pal_opt <- NULL
if ("--palette" %in% args) {
  i <- which(args == "--palette"); pal_opt <- args[i + 1]; args <- args[-c(i, i + 1)]
}
top_n <- 12L
if ("--top" %in% args) {
  i <- which(args == "--top"); top_n <- as.integer(args[i + 1]); args <- args[-c(i, i + 1)]
}
if (length(args) != 1L) stop("usage: ... config/<pheno>.yml [--top N] [--palette NAME]")

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

g    <- readRDS(file.path(paths$stage_dir, "string_02_graph_downstream.rds"))
mod  <- read.delim(file.path(paths$stage_dir, "string_03_modules.tsv"))
msum <- read.delim(file.path(paths$stage_dir, "string_03_module_summary.tsv"))
lt   <- read.delim(file.path(paths$stage_dir, "string_02_layout.tsv"))
lay  <- igraph::norm_coords(
  as.matrix(lt[match(igraph::V(g)$name, lt$string_id), c("x", "y")]), -1, 1, -1, 1)

outdir <- file.path(cfg$paths$figures, "string")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Colours. --palette wins; with no flag get_palette() falls back to the
# FIGURE_PALETTE environment variable, and with neither, to the literals below,
# which are the ones this script used before it knew about palettes. So the
# default rendering is unchanged.
pal      <- get_palette(pal_opt)
SLOT     <- pal_col(pal, "accent",  "#2a78d6")
INK      <- pal_col(pal, "ink",     "#0b0b0b")
COL_PAPER<- pal_col(pal, "paper",   "white")
COL_OFF  <- pal_col(pal, "grid",    "#ebebe8")   # node outside the highlighted set
COL_EDGE <- pal_col(pal, "grid",    "#e8e8e4")
COL_NOTE <- pal_col(pal, "muted",   "#52514e")
COL_BOX  <- pal_col(pal, "axis",    "#d8d7d2")
RAMP_LO  <- pal_col(pal, "lo",      "#f4f4f2")

# ---- block -> category weights ---------------------------------------------
f  <- load_intermediate("06_crosstraits_filtered", cfg)
t7 <- load_intermediate("07_trait_table", cfg)
bc <- f %>% select(Block, MAPPED_TRAIT) %>% distinct() %>%
      left_join(t7 %>% select(MAPPED_TRAIT, Category) %>% distinct(),
                by = "MAPPED_TRAIT", relationship = "many-to-many") %>%
      filter(!is.na(Category)) %>% count(Block, Category, name = "n_traits")

# each block's profile as a proportion of its own trait mass
bprof <- bc %>% group_by(Block) %>% mutate(w = n_traits / sum(n_traits)) %>% ungroup()

cat_rank <- bc %>% count(Category, wt = n_traits, name = "total") %>% arrange(desc(total))
cats     <- cat_rank$Category
show     <- head(cats, top_n)

# ---- 1. categories on the network ------------------------------------------
# A gene is highlighted if ITS BLOCK carries the category. Whole blocks light up
# together -- that is the data model, not a rendering choice.
blk  <- igraph::V(g)$block
hasc <- function(k) {
  b <- bc$Block[bc$Category == k]
  !is.na(blk) & blk %in% b
}

n <- length(show); ncol <- ceiling(sqrt(n)); nrow <- ceiling(n / ncol)
png(file.path(outdir, "string_network_category_facets.png"),
    width = 400 * ncol, height = 400 * nrow + 90, res = 150)
op <- par(mfrow = c(nrow, ncol), mar = c(0.4, 0.4, 2.4, 0.4),
          oma = c(0, 0, 4.0, 0), bg = COL_PAPER)
for (k in show) {
  hit <- hasc(k)
  plot(g, layout = lay, rescale = FALSE, xlim = c(-1, 1), ylim = c(-1, 1),
       vertex.color = ifelse(hit, SLOT, COL_OFF),
       vertex.frame.color = ifelse(hit, COL_PAPER, COL_OFF),
       vertex.size = ifelse(hit, 8, 4), vertex.label = NA,
       edge.color = COL_EDGE, edge.width = 0.8)
  title(sprintf("%s\n%d genes / %d blocks",
                strtrim(k, 38), sum(hit), length(unique(blk[hit]))),
        cex.main = 0.82, col.main = INK, font.main = 1, line = 0.2)
}
if (n < nrow * ncol) for (i in seq_len(nrow * ncol - n)) plot.new()
mtext(sprintf("%s — trait categories on the STRING network (top %d of %d)",
              cfg$phenotype$name, n, length(cats)),
      outer = TRUE, cex = 0.95, col = INK, line = 2.2)
mtext("Categories are inherited from the gene's LD block, so whole blocks light up together — this is not gene-level evidence.",
      outer = TRUE, cex = 0.7, col = COL_NOTE, line = 0.9)
par(op); dev.off()
message("wrote: ", file.path(outdir, "string_network_category_facets.png"))

# ---- 2. module x category profile ------------------------------------------
# Gene-weighted mixture of the member genes' block profiles.
prof <- mod %>% left_join(bprof, by = c("block" = "Block"),
                          relationship = "many-to-many") %>%
  filter(!is.na(Category)) %>%
  count(module_id, Category, wt = w, name = "mass") %>%
  group_by(module_id) %>% mutate(frac = mass / sum(mass)) %>% ungroup()

readr::write_tsv(
  prof %>% left_join(msum %>% select(module_id, size, blocks, positional_artifact),
                     by = "module_id") %>% arrange(module_id, desc(frac)),
  file.path(paths$stage_dir, "string_04_module_category_profile.tsv"), na = "")

ord  <- msum$module_id[order(-msum$size)]
M    <- matrix(0, nrow = length(ord), ncol = length(show),
               dimnames = list(paste0("M", ord), show))
for (i in seq_along(ord)) for (j in seq_along(show)) {
  v <- prof$frac[prof$module_id == ord[i] & prof$Category == show[j]]
  M[i, j] <- if (length(v)) v else 0
}
art <- msum$positional_artifact[match(ord, msum$module_id)]

ramp <- colorRampPalette(c(RAMP_LO, SLOT))(64)
png(file.path(outdir, "string_module_category_heatmap.png"),
    width = 2100, height = 1250, res = 190)
op <- par(mar = c(13.5, 6.5, 4.2, 2.2), bg = COL_PAPER)
image(x = seq_len(ncol(M)), y = seq_len(nrow(M)), z = t(M[nrow(M):1, , drop = FALSE]),
      col = ramp, axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq_len(ncol(M)), labels = colnames(M), las = 2, cex.axis = 0.68,
     col = NA, col.ticks = NA, col.axis = INK)
lab <- sprintf("%s  (n=%d%s)", rownames(M), msum$size[match(ord, msum$module_id)],
               ifelse(art, ", artifact", ""))
axis(2, at = seq_len(nrow(M)), labels = rev(lab), las = 2, cex.axis = 0.72,
     col = NA, col.ticks = NA, col.axis = INK)
grid(nx = ncol(M), ny = nrow(M), col = COL_PAPER, lty = 1, lwd = 1.5)
box(col = COL_BOX)
title(sprintf("%s — trait-category profile of each MCL module", cfg$phenotype$name),
      cex.main = 1.0, col.main = INK, font.main = 1, line = 2.6, adj = 0)
mtext("Share of the module's inherited trait mass. Profiles come from LD blocks, so block-pure modules (M1/M6/M8 = all B2) are identical by construction.",
      side = 3, line = 0.9, adj = 0, cex = 0.66, col = COL_NOTE)
par(op); dev.off()
message("wrote: ", file.path(outdir, "string_module_category_heatmap.png"))

# ---- 3. written description per module -------------------------------------
message("\n", paste(rep("=", 78), collapse = ""))
message("MODULE TRAIT-CATEGORY DESCRIPTIONS  (inherited from LD blocks)")
message(paste(rep("=", 78), collapse = ""))
for (id in ord) {
  r <- msum[msum$module_id == id, ]
  p <- prof %>% filter(module_id == id) %>% arrange(desc(frac)) %>% head(5)
  message(sprintf("\nM%-2d  n=%-2d  blocks: %-28s %s", id, r$size, r$blocks,
                  if (r$positional_artifact) "[POSITIONAL ARTIFACT]" else ""))
  if (!nrow(p)) {
    message("     no cross-traits in any of its blocks")
  } else {
    for (k in seq_len(nrow(p))) {
      message(sprintf("     %5.1f%%  %s", 100 * p$frac[k], p$Category[k]))
    }
  }
}
message("\nprofile table -> ",
        file.path(paths$stage_dir, "string_04_module_category_profile.tsv"))
