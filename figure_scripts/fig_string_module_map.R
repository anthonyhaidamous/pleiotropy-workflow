#!/usr/bin/env Rscript
# ============================================================
# fig_string_module_map.R
# The STRING network with its modules ringed and named.
#
# The module figure that answers "what is in there", as opposed to
# fig_string_network_by_module.R, which answers "which gene is where". Nodes
# carry no labels; each module is enclosed by a dashed outline over a
# half-opaque fill, and named with the function its enrichment supports.
#
# EVERY MODULE IS RINGED, INCLUDING THE POSITIONAL ONES
# -----------------------------------------------------
# fig_string_network_by_module.R greys the modules that stage 03 flags as
# positional artifacts, because that figure is about the modules that carry a
# cross-locus claim. This one does the opposite on purpose: a module that is
# 21 of 31 genes from B21 is still a real cluster of the network, and naming it
# "MHC class II antigen presentation" is the clearest possible statement of WHY
# it is flagged. The flag is reported in the console and can be printed in the
# label with --mark-artifacts; it never removes a module from the drawing.
#
# WHICH MODULES GET A NAME
# ------------------------
# The names come from analysis/module_enrichment.R (module_labels.tsv), which
# only accepts a headline term that meets its support rule. A module with no
# qualifying term is drawn and ringed but left unnamed rather than given a weak
# label — several of the small modules are three genes with nothing in common
# but each other, and a GO term stretched over them would be an invention.
# Run module_enrichment.R with --min-blocks 1 to name the positional modules
# too, which is the intended pairing for this figure.
#
# HULLS
# -----
# Each module's outline is the convex hull of discs of radius --expand centred
# on its nodes: rounded corners, guaranteed to contain every member with a
# margin, and degenerating gracefully to a circle for a two-gene module. Hulls
# of large modules can overlap a neighbour's nodes — 11 such intrusions in the
# CLL layout, all in the three biggest modules — which is a property of a
# force-directed layout, not of the clustering. The console reports the count.
#
# Inputs (all under <stage_dir>):
#   string_02_layout.tsv         node coordinates (from fig_string_network_by_block.R)
#   string_02_graph_downstream.rds
#   string_03_modules.tsv, string_03_module_summary.tsv
#   module_enrichment/module_labels.tsv     names, if present
#
# Outputs (in <figures>/string/):
#   string_module_map.svg and .png
#
# Usage:
#   Rscript figure_scripts/fig_string_module_map.R config/cll.yml --palette poster
#
#   --labels FILE     two columns, module_id and label, overriding the
#                     enrichment names. Written as a template on first run so a
#                     curated wording can replace the raw term.
#   --layout MODE     network (default) keeps the stage-02 coordinates, so this
#                     figure is directly comparable with
#                     fig_string_network_by_module.R — same layout, same node
#                     size, same shape. Hulls then overlap where the embedding
#                     interleaves modules (21 genes of one module inside
#                     another's hull, on CLL); that is a property of a
#                     force-directed layout, not of the clustering.
#                     modules re-lays the graph so every module is a separated
#                     disc: clean rings, but an embedding that shares no
#                     landmarks with the sibling figure.
#   --min-genes N     ring modules with at least N genes           [2]
#   --expand N        hull padding, in layout units                [0.25]
#   --gap N           clear space left between two hulls            [0.25]
#   --alpha A         module fill opacity                          [0.5]
#   --mark-artifacts  append "(positional)" to flagged modules' names
#   --no-edges        drop the interaction lines, leaving modules alone
#   --node-size N     node diameter, ggplot units. Left unset it follows
#                     --base-size; set it to match a sibling figure drawn with
#                     igraph (vertex.size 7 on a [-1,1] layout is 4.5 here,
#                     which is the default and matches the sibling figure).
#   --base-size N     type size in points at the requested --width [11]
#   --width, --height, --dpi, --outfile, --palette
# ============================================================

suppressPackageStartupMessages({
  library(igraph); library(ggplot2); library(dplyr)
})

if (!file.exists("R/utils.R")) stop("run from the repository root")
source("R/utils.R")
source("R/string_utils.R")
source("R/palettes.R")

args <- commandArgs(trailingOnly = TRUE)
flag <- function(f) { r <- f %in% args; args <<- args[args != f]; r }
mark_art <- flag("--mark-artifacts")
no_edges <- flag("--no-edges")
take <- function(f, default, cast = as.character) {
  if (!(f %in% args)) return(default)
  i <- which(args == f); v <- cast(args[i + 1]); args <<- args[-c(i, i + 1)]; v
}
layout_mode <- take("--layout",    "network")
if (!layout_mode %in% c("modules", "network")) stop("--layout is modules or network")
labels_file <- take("--labels",    NULL)
min_genes   <- take("--min-genes", 2L,   as.integer)
expand      <- take("--expand",    0.25, as.numeric)
gap_between <- take("--gap",       0.25, as.numeric)
alpha_fill  <- take("--alpha",     0.5,  as.numeric)
base_size   <- take("--base-size", 11,   as.numeric)
node_size   <- take("--node-size", 4.5,  as.numeric)
fig_w       <- take("--width",     9.51, as.numeric)
fig_h       <- take("--height",    8.0,  as.numeric)
fig_dpi     <- take("--dpi",       600L, as.integer)
outfile     <- take("--outfile",   "string_module_map")
palette     <- take("--palette",   NULL)

if (is.na(alpha_fill) || alpha_fill < 0 || alpha_fill > 1) stop("--alpha is in [0, 1]")
if (is.na(expand) || expand <= 0) stop("--expand must be positive")
if (length(args) != 1L) {
  stop("usage: Rscript figure_scripts/fig_string_module_map.R config/<pheno>.yml",
       " [--labels FILE] [--min-genes N] [--expand N] [--alpha A]",
       " [--mark-artifacts] [--no-edges] [--base-size N] [--width N] [--height N]",
       " [--dpi N] [--outfile NAME] [--palette NAME]")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

need <- function(p) { if (!file.exists(p)) stop("not found: ", p); p }
lay  <- read.delim(need(file.path(paths$stage_dir, "string_02_layout.tsv")), stringsAsFactors = FALSE)
mod  <- read.delim(need(file.path(paths$stage_dir, "string_03_modules.tsv")), stringsAsFactors = FALSE)
msum <- read.delim(need(file.path(paths$stage_dir, "string_03_module_summary.tsv")), stringsAsFactors = FALSE)
g    <- readRDS(need(file.path(paths$stage_dir, "string_02_graph_downstream.rds")))

outdir <- file.path(cfg$paths$figures, "string")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

nodes <- merge(lay[, c("gene", "x", "y")], mod[, c("gene", "module_id")], by = "gene")
if (!nrow(nodes)) stop("no gene matched between the layout and the module table")
missing_lay <- setdiff(mod$gene, lay$gene)
if (length(missing_lay)) {
  message("no layout coordinate for ", length(missing_lay), " gene(s), not drawn: ",
          paste(head(missing_lay, 8), collapse = ", "))
}

# ---- layout ----------------------------------------------------------------
# The stage-02 layout is a single force-directed embedding of the whole graph,
# which is right for a gene-level figure and wrong for this one: modules
# interpenetrate, so their hulls overlap — 34 foreign genes inside a hull, on
# the CLL network — and a ringed region stops meaning anything.
#
# `--layout modules` (the default) rebuilds it in two levels instead. Module
# CENTRES come from a force-directed layout of the module graph, so modules that
# interact sit near each other and the picture still carries the network's shape;
# members are then laid out locally inside a disc of radius proportional to
# sqrt(size), and the centres are relaxed apart until no two discs touch. Hulls
# are then separated by construction. `--layout network` keeps the stage-02
# coordinates, which makes this figure directly comparable with
# fig_string_network_by_module.R at the cost of the overlaps.
pack_modules <- function(nodes, g, hull_pad, gap = 0.25, iter = 400L) {
  ms <- sort(unique(nodes$module_id))
  n_of <- vapply(ms, function(m) sum(nodes$module_id == m), integer(1))
  # Two radii, and keeping them distinct is what stops the map sprawling: the
  # members occupy `node_r`, and the HULL reaches node_r + the padding it is
  # drawn with. Relaxation must separate on the second, or every pair is held
  # apart by a margin that is not in the picture — which is what put a disc's
  # width of blank between neighbours.
  node_r <- 0.62 * sqrt(n_of)
  rad    <- node_r + hull_pad

  # module-level graph: an edge wherever two modules interact, weighted by count
  # V(g)$name is the STRING protein id (9606.ENSP...), NOT the gene symbol —
  # matching gene symbols against it fails silently and empties the edge list.
  gname <- igraph::V(g)$gene
  mm <- setNames(nodes$module_id, nodes$gene)
  el <- igraph::as_edgelist(g, names = FALSE)
  a  <- mm[gname[el[, 1]]]; b <- mm[gname[el[, 2]]]
  keep <- !is.na(a) & !is.na(b) & a != b
  cen <- if (any(keep)) {
    pairs <- paste(pmin(a[keep], b[keep]), pmax(a[keep], b[keep]))
    tab   <- table(pairs)
    ends  <- do.call(rbind, strsplit(names(tab), " "))
    mg <- igraph::graph_from_data_frame(
      data.frame(from = ends[, 1], to = ends[, 2], weight = as.numeric(tab)),
      directed = FALSE, vertices = data.frame(name = as.character(ms)))
    set.seed(1)
    igraph::layout_with_fr(mg, weights = igraph::E(mg)$weight)
  } else {
    igraph::layout_in_circle(igraph::make_empty_graph(length(ms)))
  }
  # Start compact and let the relaxation open it only as far as it must.
  cen <- cen / max(1e-9, max(abs(cen))) * (0.55 * sqrt(sum(n_of)))

  # relax the centres until no two discs overlap
  for (it in seq_len(iter)) {
    moved <- FALSE
    for (i in seq_along(ms)) for (j in seq_along(ms)) {
      if (j <= i) next
      d <- cen[j, ] - cen[i, ]
      dist <- sqrt(sum(d^2))
      need <- rad[i] + rad[j] + gap
      if (dist < need) {
        if (dist < 1e-6) { d <- c(cos(i * 2.4), sin(i * 2.4)); dist <- 1e-6 }
        push <- (need - dist) / 2 * d / dist
        cen[i, ] <- cen[i, ] - push; cen[j, ] <- cen[j, ] + push
        moved <- TRUE
      }
    }
    if (!moved) break
  }

  out <- do.call(rbind, lapply(seq_along(ms), function(k) {
    m  <- ms[k]; d <- nodes[nodes$module_id == m, ]
    sub <- igraph::induced_subgraph(g, which(igraph::V(g)$gene %in% d$gene))
    set.seed(k)
    loc <- if (igraph::vcount(sub) > 2 && igraph::ecount(sub) > 0) {
      l <- igraph::layout_with_fr(sub)
      rownames(l) <- igraph::V(sub)$gene
      l[d$gene, , drop = FALSE]
    } else if (nrow(d) > 1) {
      th <- seq(0, 2 * pi, length.out = nrow(d) + 1)[-(nrow(d) + 1)]
      cbind(cos(th), sin(th))
    } else cbind(0, 0)
    span <- max(1e-9, max(abs(loc)))
    loc  <- loc / span * node_r[k]                # the hull pad is added on top
    data.frame(gene = d$gene, module_id = m,
               x = cen[k, 1] + loc[, 1], y = cen[k, 2] + loc[, 2],
               stringsAsFactors = FALSE)
  }))
  out
}

if (identical(layout_mode, "modules")) {
  nodes <- pack_modules(nodes, g, hull_pad = expand, gap = gap_between)
  lay   <- nodes[, c("gene", "x", "y")]           # edges follow the new positions
}

# ---- module names ----------------------------------------------------------
# Raw enrichment terms are written for a table, not a figure: "Peptide Antigen
# Assembly With MHC Protein Complex (GO:0002396)" is nine words and an
# accession. They are shortened here, and the shortened set is written out so a
# curated wording can be supplied with --labels instead.
shorten <- function(x) {
  x <- sub("\\s*\\((GO|R-HSA|WP)[:0-9A-Za-z-]+\\)\\s*$", "", x)   # drop accessions
  x <- sub("\\s+Homo sapiens.*$", "", x)
  # KEGG returns SHOUTED terms; sentence-case them, leave mixed case alone.
  up <- !is.na(x) & x == toupper(x) & nchar(x) > 3
  x[up] <- paste0(substr(x[up], 1, 1), tolower(substring(x[up], 2)))
  trimws(x)
}

lab_path <- file.path(paths$stage_dir, "module_enrichment", "module_labels.tsv")
labs <- data.frame(module_id = msum$module_id, label = NA_character_,
                   stringsAsFactors = FALSE)
if (file.exists(lab_path)) {
  ml <- read.delim(lab_path, stringsAsFactors = FALSE)
  labs$label <- shorten(ml$label[match(labs$module_id, ml$module_id)])
} else {
  message("no module_labels.tsv at ", lab_path,
          "\n  run analysis/module_enrichment.R first; modules will be numbered only")
}
if (!is.null(labels_file)) {
  if (!file.exists(labels_file)) stop("--labels file not found: ", labels_file)
  cur <- read.delim(labels_file, stringsAsFactors = FALSE)
  if (!all(c("module_id", "label") %in% names(cur)))
    stop("--labels needs columns module_id and label")
  hit <- match(labs$module_id, cur$module_id)
  ok  <- !is.na(hit) & nzchar(trimws(cur$label[hit]))
  labs$label[ok] <- trimws(cur$label[hit][ok])
  message("curated labels applied for ", sum(ok), " module(s) from ", labels_file)
}

labs$artifact <- msum$positional_artifact[match(labs$module_id, msum$module_id)]
labs$size     <- msum$size[match(labs$module_id, msum$module_id)]
labs$n_blocks <- msum$n_blocks[match(labs$module_id, msum$module_id)]

# The template is written every run so it tracks a re-run of the enrichment;
# it is never read back automatically, only through --labels.
tmpl <- file.path(outdir, paste0(outfile, "_labels_template.tsv"))
write.table(labs[, c("module_id", "size", "n_blocks", "artifact", "label")],
            tmpl, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ---- hulls -----------------------------------------------------------------
# Convex hull of discs centred on the members: rounded, always containing every
# node with `expand` to spare, and a plain circle when a module has one gene.
disc_hull <- function(x, y, r, k = 36L) {
  th <- seq(0, 2 * pi, length.out = k + 1L)[-(k + 1L)]
  px <- as.vector(outer(x, r * cos(th), "+"))
  py <- as.vector(outer(y, r * sin(th), "+"))
  h  <- grDevices::chull(px, py)
  data.frame(x = px[h], y = py[h])
}
in_poly <- function(px, py, hx, hy) {
  n <- length(hx); j <- n; inside <- FALSE
  for (i in seq_len(n)) {
    if (((hy[i] > py) != (hy[j] > py)) &&
        (px < (hx[j] - hx[i]) * (py - hy[i]) / (hy[j] - hy[i]) + hx[i])) inside <- !inside
    j <- i
  }
  inside
}

drawn <- sort(unique(nodes$module_id[ave(nodes$module_id, nodes$module_id, FUN = length) >= min_genes]))
hulls <- do.call(rbind, lapply(drawn, function(m) {
  d <- nodes[nodes$module_id == m, ]
  cbind(disc_hull(d$x, d$y, expand), module_id = m)
}))

intr <- 0L
for (m in drawn) {
  h <- hulls[hulls$module_id == m, ]; o <- nodes[nodes$module_id != m, ]
  intr <- intr + sum(mapply(function(a, b) in_poly(a, b, h$x, h$y), o$x, o$y))
}
message(sprintf("%d modules drawn (of %d); %d gene(s) of another module fall inside a hull",
                length(drawn), nrow(msum), intr))

# ---- text placed above each hull -------------------------------------------
# Labels are anchored RADIALLY OUTWARD from the centre of the whole map, not
# above each hull. With modules packed as discs, "outward" is by construction
# the direction with nothing in it; "above" walks a label straight onto the
# neighbour that happens to sit there, which is what put M2's name across M6.
gcx <- mean(range(nodes$x)); gcy <- mean(range(nodes$y))
wrap_at <- function(x, n = 24L) paste(strwrap(x, width = n), collapse = "\n")

cap <- do.call(rbind, lapply(drawn, function(m) {
  h <- hulls[hulls$module_id == m, ]; r <- labs[labs$module_id == m, ]
  cx <- mean(range(h$x)); cy <- mean(range(h$y))
  d  <- c(cx - gcx, cy - gcy)
  nd <- sqrt(sum(d^2)); if (nd < 1e-6) { d <- c(0, 1); nd <- 1 }
  d  <- d / nd
  reach <- max(sqrt((h$x - cx)^2 + (h$y - cy)^2))          # hull radius that way
  txt <- sprintf("M%d", m)
  if (!is.na(r$label) && nzchar(r$label)) txt <- paste0(txt, " · ", wrap_at(r$label))
  if (mark_art && isTRUE(r$artifact)) txt <- paste0(txt, "\n(positional)")
  data.frame(module_id = m,
             x = cx + d[1] * (reach + expand * 1.6),
             y = cy + d[2] * (reach + expand * 1.6),
             label = txt, stringsAsFactors = FALSE)
}))
# A label anchored outward can still land on a neighbour — always for the
# module nearest the centre of the map, where "outward" has no meaning (M1, in
# the CLL network, whose anchor fell inside M7). Push any anchor that lands in
# a foreign hull further along its own direction until it is clear.
for (i in seq_len(nrow(cap))) {
  m  <- cap$module_id[i]
  cx <- mean(range(hulls$x[hulls$module_id == m]))
  cy <- mean(range(hulls$y[hulls$module_id == m]))
  d  <- c(cap$x[i] - cx, cap$y[i] - cy)
  nd <- sqrt(sum(d^2)); if (nd < 1e-6) next
  d  <- d / nd
  for (step in 0:24) {
    px <- cap$x[i] + d[1] * step * expand * 0.5
    py <- cap$y[i] + d[2] * step * expand * 0.5
    clash <- any(vapply(setdiff(drawn, m), function(o) {
      h <- hulls[hulls$module_id == o, ]; in_poly(px, py, h$x, h$y)
    }, logical(1)))
    if (!clash) { cap$x[i] <- px; cap$y[i] <- py; break }
  }
}

named <- sum(!is.na(labs$label[labs$module_id %in% drawn]) &
             nzchar(labs$label[labs$module_id %in% drawn]))
message(sprintf("%d of the %d drawn modules carry a name; the rest are numbered only",
                named, length(drawn)))
message("label template written to ", tmpl, " — edit and pass with --labels to curate")

# ---- colours ---------------------------------------------------------------
pal   <- get_palette(palette)
PAPER <- pal_col(pal, "paper", "white")
INK   <- pal_col(pal, "ink",   "#0b0b0b")
MUTED <- pal_col(pal, "muted", "#4b505a")
EDGE  <- pal_col(pal, "axis",  "#9aa0ab")   # grid is too faint under a 50% fill
# Fallback ramp, used only when the palette in force defines no `cat`.
RAMP <- c("#4477AA", "#EE7733", "#117733", "#EE3377", "#88CCEE", "#882255",
          "#DDCC77", "#009988", "#AA4499", "#999933", "#EE6677", "#332288",
          "#994F00", "#CCBB44", "#66CCEE", "#44AA99", "#771155", "#BBBBBB",
          "#0077BB", "#225522")

# Module colours come from the palette's own categorical ramp, and the order is
# chosen so this figure AGREES with fig_string_network_by_module.R: that script
# colours the non-artifact modules with pal_cat(pal, length(keep)) in ascending
# module order and greys the rest, so taking the same slots for the same modules
# makes M1 the same colour in both. The artifact modules — which this figure
# draws in colour rather than greying — take the ramp's remaining entries,
# interpolated if there are more modules than slots.
catv <- if (!is.null(pal$cat)) as.character(pal$cat) else RAMP
keep_mods <- sort(msum$module_id[!msum$positional_artifact])
keep_mods <- intersect(keep_mods, drawn)
rest_mods <- setdiff(drawn, keep_mods)

n_keep <- min(length(keep_mods), length(catv))
cols_keep <- if (length(keep_mods)) {
  grDevices::colorRampPalette(catv[seq_len(max(1L, n_keep))])(length(keep_mods))
} else character()
tail_cols <- if (length(catv) > n_keep) catv[(n_keep + 1):length(catv)] else catv
cols_rest <- if (length(rest_mods)) {
  grDevices::colorRampPalette(tail_cols)(length(rest_mods))
} else character()
cols <- setNames(c(cols_keep, cols_rest), c(keep_mods, rest_mods))
cols <- cols[as.character(drawn)]
names(cols) <- drawn

# ---- edges -----------------------------------------------------------------
el   <- igraph::as_edgelist(g, names = FALSE)
gsym <- igraph::V(g)$gene
e1   <- gsym[el[, 1]]; e2 <- gsym[el[, 2]]
ed <- data.frame(
  x    = lay$x[match(e1, lay$gene)], y    = lay$y[match(e1, lay$gene)],
  xend = lay$x[match(e2, lay$gene)], yend = lay$y[match(e2, lay$gene)])
ed <- ed[complete.cases(ed), ]

# ---- draw ------------------------------------------------------------------
nodes$fill <- cols[as.character(nodes$module_id)]
nodes$fill[is.na(nodes$fill)] <- MUTED
hulls$module_id <- factor(hulls$module_id, levels = drawn)

# Order matters: the hull FILL goes down first, the edges over it, then the
# dashed rim and the nodes. Drawn the other way round a half-opaque fill washes
# out every intra-module edge, which is most of them.
p <- ggplot() +
  geom_polygon(data = hulls, aes(x = x, y = y, group = module_id, fill = module_id),
               alpha = alpha_fill, colour = NA)
if (!no_edges) {
  p <- p + geom_segment(data = ed, aes(x = x, y = y, xend = xend, yend = yend),
                        colour = EDGE, linewidth = 0.22, alpha = 0.75)
}
p <- p +
  geom_polygon(data = hulls, aes(x = x, y = y, group = module_id, colour = module_id),
               fill = NA, linetype = "dashed", linewidth = 0.45) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols, guide = "none") +
  geom_point(data = nodes, aes(x = x, y = y), colour = nodes$fill,
             size = if (is.na(node_size)) base_size * 0.16 else node_size,
             alpha = 0.95) +
  coord_equal(clip = "off") +
  theme_void(base_size = base_size) +
  theme(plot.background  = element_rect(fill = PAPER, colour = NA),
        panel.background = element_rect(fill = PAPER, colour = NA),
        plot.margin      = margin(6, 6, 6, 6))

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(
    data = cap, aes(x = x, y = y, label = label),
    size = base_size * 0.30, colour = INK, seed = 1, fontface = "bold",
    direction = "both", max.overlaps = Inf,
    box.padding = 0.30, point.padding = 0.1, min.segment.length = 0.4,
    segment.colour = MUTED, segment.size = 0.25, lineheight = 0.95)
} else {
  message("NOTE: ggrepel not installed; labels placed without collision avoidance.")
  p <- p + geom_text(data = cap, aes(x = x, y = y + expand, label = label),
                     size = base_size * 0.30, colour = INK, fontface = "bold")
}

svg_path <- file.path(outdir, paste0(outfile, ".svg"))
if (!requireNamespace("svglite", quietly = TRUE)) {
  stop("svglite is required for the SVG output; install it or use --outfile with png only")
}
ggsave(svg_path, p, width = fig_w, height = fig_h, device = svglite::svglite, bg = PAPER)
message("SVG saved -> ", svg_path)

png_path <- file.path(outdir, paste0(outfile, ".png"))
ggsave(png_path, p, width = fig_w, height = fig_h, dpi = fig_dpi, bg = PAPER)
message("PNG saved -> ", png_path)
