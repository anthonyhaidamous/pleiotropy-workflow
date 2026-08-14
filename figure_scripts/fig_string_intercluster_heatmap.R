#!/usr/bin/env Rscript
# =============================================================================
# STRING inter-cluster connectivity heatmap  (phenotype-agnostic)
#
# Builds an N×N matrix of edge counts between k-means protein clusters from a
# STRING network export and renders it as a ggplot2 geom_tile heatmap.
#
# Inputs (STRING "plus50" k-means export, tab-separated, '#'-prefixed header):
#   <indir>/plus50_interactions_short.tsv          edges (node1, node2, scores)
#   <indir>/plus50_kmeans_clusters.tsv             node -> cluster
#   <indir>/plus50_kmeans_cluster_descriptions.tsv cluster names (optional)
#
# Cells are raw edge COUNTS (the STRING "short" file has per-channel sub-scores,
# not a single combined confidence). Diagonal = intra-cluster edges.
#
# Usage:
#   Rscript figure_scripts/fig_string_intercluster_heatmap.R --trait leptin
#   # first six clusters only, log-compressed colour:
#   Rscript figure_scripts/fig_string_intercluster_heatmap.R \
#       --trait chronic_lymphocytic_leukemia --clusters 1-6 --log
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))
if (!requireNamespace("optparse", quietly = TRUE))
  stop("Please install the 'optparse' package.", call. = FALSE)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--trait",      type = "character", default = "chronic_lymphocytic_leukemia",
              help = "Phenotype name; sets default in/out paths [default: %default]"),
  make_option("--indir",      type = "character", default = NULL,
              help = "STRING export dir [default: results/intermediate/<trait>/string]"),
  make_option("--prefix",     type = "character", default = "plus50",
              help = "Filename prefix for the k-means export trio (<prefix>_interactions_short.tsv etc.) [default: %default]"),
  make_option("--outfile",    type = "character", default = NULL,
              help = "Output file; extension picks the format (.png or .svg) [default: figures/<trait>/string_heatmap_intercluster.png]"),
  make_option("--clusters",   type = "character", default = "all",
              help = "Clusters to include: 'all', '1-6', '1,3,5', '1-3,7' [default: %default]"),
  make_option("--log",        action = "store_true", default = FALSE,
              help = "Log-compress the colour scale (log1p) to offset the intra-cluster diagonal"),
  make_option("--highlight",     type = "character", default = "1",
              help = "Cluster number to emphasise, or 'none' [default: %default]"),
  make_option("--names-file",   type = "character", default = NULL,
              help = "Optional TSV (cluster_number <tab> name) overriding axis labels"),
  make_option("--short-labels", action = "store_true", default = FALSE,
              help = "Use 'Cluster N' on axes instead of full names (good for panels)"),
  make_option("--width",        type = "double", default = NULL,
              help = "Figure width in inches (overrides auto-sizing)"),
  make_option("--height",       type = "double", default = NULL,
              help = "Figure height in inches (overrides auto-sizing)"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)))

INDIR   <- if (is.null(opt$indir))   file.path("results/intermediate", opt$trait, "string") else opt$indir
OUTFILE <- if (is.null(opt$outfile)) file.path("figures", opt$trait, "string_heatmap_intercluster.png") else opt$outfile
LOG     <- isTRUE(opt$log)
FONT    <- "Arial"

PREFIX       <- opt$prefix
EDGE_FILE    <- file.path(INDIR, paste0(PREFIX, "_interactions_short.tsv"))
CLUSTER_FILE <- file.path(INDIR, paste0(PREFIX, "_kmeans_clusters.tsv"))
DESC_FILE    <- file.path(INDIR, paste0(PREFIX, "_kmeans_cluster_descriptions.tsv"))
for (f in c(EDGE_FILE, CLUSTER_FILE))
  if (!file.exists(f)) stop("Missing required file: ", f, call. = FALSE)

# ── Read ('comment.char=""' keeps the '#'-prefixed header line) ───────────────
read_string <- function(path)
  read.delim(path, sep = "\t", header = TRUE, quote = "",
             comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)

edges <- read_string(EDGE_FILE)
nodes <- read_string(CLUSTER_FILE)

# node1/node2 = first two cols of edge file; cluster number = col2, gene = col6.
e <- data.frame(a = edges[[1]], b = edges[[2]], stringsAsFactors = FALSE)
node2cluster <- setNames(as.integer(nodes[[2]]), nodes[[6]])
available    <- sort(unique(as.integer(nodes[[2]])))

# ── Cluster selection ────────────────────────────────────────────────────────
parse_clusters <- function(spec, available) {
  spec <- tolower(trimws(spec))
  if (spec %in% c("", "all")) return(available)
  out <- integer(0)
  for (p in strsplit(spec, ",")[[1]]) {
    p <- trimws(p)
    if (grepl("[-:]", p)) { r <- as.integer(strsplit(p, "[-:]")[[1]]); out <- c(out, seq(r[1], r[2])) }
    else                    out <- c(out, as.integer(p))
  }
  sel <- intersect(sort(unique(out)), available)
  if (!length(sel)) stop("No requested clusters exist (available: ",
                          paste(available, collapse = ","), ")", call. = FALSE)
  sel
}
sel <- parse_clusters(opt$clusters, available)
N   <- length(sel)

# ── Cluster names: short-labels > names-file > descriptions file > "Cluster N" ─
trunc <- function(x, n = 38) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
nm <- setNames(paste("Cluster", available), available)
if (!isTRUE(opt[["short-labels"]])) {
  if (file.exists(DESC_FILE)) {
    d <- read_string(DESC_FILE); num <- as.integer(d[[2]]); prim <- trimws(d[[6]])
    for (k in seq_along(num)) {
      lab <- prim[k]
      if (is.na(lab) || lab %in% c("", "-")) lab <- paste("Peripheral", num[k])
      nm[as.character(num[k])] <- lab
    }
  }
  if (!is.null(opt[["names-file"]])) {
    o <- read.delim(opt[["names-file"]], sep = "\t", header = FALSE, quote = "",
                    stringsAsFactors = FALSE, comment.char = "#")
    nm[as.character(as.integer(o[[1]]))] <- o[[2]]
  }
}
labs <- trunc(unname(nm[as.character(sel)]))

# ── Join cluster membership, subset to selected clusters ─────────────────────
e$ca <- node2cluster[e$a]; e$cb <- node2cluster[e$b]
dropped_unmapped <- sum(is.na(e$ca) | is.na(e$cb))
e <- e[!is.na(e$ca) & !is.na(e$cb), ]
key <- paste(pmin(e$a, e$b), pmax(e$a, e$b), sep = "\t")
dups <- sum(duplicated(key)); e <- e[!duplicated(key), ]
e <- e[e$ca %in% sel & e$cb %in% sel, ]

# ── N×N symmetric connectivity matrix (positions follow `sel`) ───────────────
pa <- match(e$ca, sel); pb <- match(e$cb, sel)
lo <- pmin(pa, pb); hi <- pmax(pa, pb)
M  <- matrix(0L, N, N)
for (k in seq_along(lo)) {
  M[lo[k], hi[k]] <- M[lo[k], hi[k]] + 1L
  if (lo[k] != hi[k]) M[hi[k], lo[k]] <- M[hi[k], lo[k]] + 1L
}
cat(sprintf("Clusters: %s | edges kept: %d (%d unmapped, %d dup) | intra: %d | inter: %d\n",
            paste(sel, collapse = ","), nrow(e), dropped_unmapped, dups,
            sum(diag(M)), (sum(M) - sum(diag(M))) / 2))

# ── Long form ────────────────────────────────────────────────────────────────
df <- expand.grid(i = 1:N, j = 1:N)
df$count <- M[cbind(df$i, df$j)]
df$x <- factor(labs[df$i], levels = labs)
df$y <- factor(labs[df$j], levels = rev(labs))      # first selected cluster at top
tv   <- if (LOG) log1p(df$count) else df$count       # transformed value drives text colour
df$white_text <- tv > 0.5 * max(tv)

# ── Highlight ────────────────────────────────────────────────────────────────
hl_pos <- if (tolower(opt$highlight) == "none") NA_integer_ else match(as.integer(opt$highlight), sel)

# Colours, resolved once. Each fallback is the literal this script used
# before it took a --palette, so the default rendering is unchanged.
source("R/palettes.R")
pal        <- get_palette(opt$palette)
hi_col     <- pal_col(pal, "accent",  "#b30000")   # highlighted labels
lab_col    <- pal_col(pal, "muted",   "grey25")    # ordinary labels
grad_lo    <- pal_col(pal, "lo",      "white")     # empty cell
grad_hi    <- pal_col(pal, "hi",      "#8b0000")   # densest cell
cell_line  <- pal_col(pal, "paper",   "white")     # tile separators
txt_on_hi  <- pal_col(pal, "paper",   "white")     # label on a dark cell
txt_on_lo  <- pal_col(pal, "ink",     "grey20")    # label on a light cell
col_paper  <- pal_col(pal, "paper",   "white")
mk_axis <- function(order_labs) {
  is_hl <- !is.na(hl_pos) & order_labs == labs[hl_pos]
  list(col = ifelse(is_hl, hi_col, lab_col), face = ifelse(is_hl, "bold", "plain"))
}
ax <- mk_axis(labs); ay <- mk_axis(rev(labs))

# ── Plot ─────────────────────────────────────────────────────────────────────
fill_scale <- if (LOG) {
  scale_fill_gradient(low = grad_lo, high = grad_hi, trans = "log1p",
                      breaks = c(0, 1, 5, 25, 100, 500), name = "Edges\n(log scale)")
} else {
  scale_fill_gradient(low = grad_lo, high = grad_hi, name = "Edges")
}

p <- ggplot(df, aes(x, y, fill = count)) +
  geom_tile(color = cell_line, linewidth = 0.4) +
  geom_text(aes(label = ifelse(count > 0, count, ""), color = white_text),
            family = FONT, size = 1.8, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = txt_on_hi, `FALSE` = txt_on_lo)) +
  fill_scale +
  scale_x_discrete(position = "top") +
  coord_equal() +
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  theme_minimal(base_size = 6, base_family = FONT) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 0, vjust = 0,
                                    colour = ax$col, face = ax$face, size = 5.5),
    axis.text.y      = element_text(colour = ay$col, face = ay$face,
                                    size = 5.5, hjust = 1),
    legend.title     = element_text(size = 5.5, face = "bold"),
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(0.25, "cm"),
    legend.key.width = unit(0.2,  "cm"),
    legend.position  = "right",
    legend.margin    = margin(0, 0, 0, 1),
    panel.grid       = element_blank(),
    plot.margin      = margin(1, 1, 1, 1)
  )

# ── Save ─────────────────────────────────────────────────────────────────────
dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)
# Device follows the output extension: .svg -> vector (svglite), else PNG (ragg if available).
ext <- tolower(tools::file_ext(OUTFILE))
dev <- if (ext == "svg") {
  if (!requireNamespace("svglite", quietly = TRUE))
    stop("Writing SVG needs the 'svglite' package (install.packages('svglite')).", call. = FALSE)
  # fix_text_size = FALSE -> no textLength/lengthAdjust, so text stays editable.
  function(...) svglite::svglite(..., fix_text_size = FALSE)
} else if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
# Auto-size: compact panel proportions. Short labels need much less margin.
short <- isTRUE(opt[["short-labels"]])
side  <- max(1.6, 0.28 * N + if (short) 0.3 else 0.8)
W <- if (!is.null(opt$width))  opt$width  else side + if (short) 0.9 else 1.4
H <- if (!is.null(opt$height)) opt$height else side + if (short) 0.3 else 0.5
ggsave(OUTFILE, plot = p, device = dev,
       width = W, height = H, units = "in", dpi = 300, bg = col_paper)
cat(sprintf("Saved -> %s  (%.2f × %.2f in)\n", OUTFILE, W, H))
cat("Saved ->", OUTFILE, "\n")
