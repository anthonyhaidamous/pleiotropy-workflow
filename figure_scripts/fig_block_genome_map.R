#!/usr/bin/env Rscript
# ============================================================
# fig_block_genome_map.R
# Where the LD blocks are, and nothing else.
#
# Twenty-two chromosomes drawn to scale as plain rounded bars, one per row, with
# every block marked at its own position. No cytobands, no centromeres, no
# ideogram staining: the question is which loci and where, and a Giemsa pattern
# answers a different one while costing most of the ink.
#
# WHY NOT karyoploteR
# -------------------
# fig_karyoploter_trait_ld_blocks.R is the right tool when the base-pair detail
# matters — it draws real ideograms, per-chromosome zooms, and LD structure
# inside a locus. On a poster read from two metres none of that survives: the
# cytobands become grey mush, the chromosome outlines carry more visual weight
# than the marks, and 49 blocks on 22 ideograms is a picture of an ideogram, not
# of the result. This script keeps the one thing that has to be legible — the
# position of each block against the whole genome — and drops the rest.
#
# WHAT A MARK ENCODES
# -------------------
# Position is the block's midpoint. Blocks are 0.03-1.2 Mb against a 249 Mb
# chromosome, so drawn to scale a block is a third of a pixel: the mark is a
# dot, deliberately not a scaled interval, and the console reports the true span
# of each so nothing is implied about width.
#
# With --crosstraits the dot AREA is the block's cross-trait count, which puts
# the pleiotropy on the map without a second figure: the eye picks out B17, B23
# and B10 as the large dots and reads their positions in the same glance. Blocks
# with no cross-trait are drawn hollow rather than dropped — an LD block that
# carries none is a result, and this is the figure that can say so.
#
# Inputs:
#   --blocks    04_blocks.{tsv,rds}          Block, Genomic_position (chr:start-end)
#   --crosstraits 06_crosstraits_filtered.{tsv,rds}   optional, sizes the dots
#   --chrsizes  chr <tab> size               default data/chromosomes/chrom_sizes.txt
#
# Outputs (in --outdir): <outfile>.png, and .svg with --svg
#
# Usage:
#   Rscript figure_scripts/fig_block_genome_map.R \
#       --blocks      results/intermediate/chronic_lymphocytic_leukemia/04_blocks.rds \
#       --crosstraits results/intermediate/chronic_lymphocytic_leukemia/06_crosstraits_filtered.tsv \
#       --outdir      figures/chronic_lymphocytic_leukemia --palette poster
#
#   --cols N        chromosome columns: 1 (a tall ladder) or 2 (squarer) [1]
#   --label-top N   label the N most pleiotropic blocks; 0 labels none,
#                   -1 labels every block                                 [12]
#   --no-size       one dot size for every block, ignoring --crosstraits
#   --base-size N   type size in points at the requested --width          [11]
#   --width, --height, --dpi, --svg, --outfile, --palette
# ============================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2); library(gtools)
})

if (!file.exists("R/figure_utils.R")) stop("run from the repository root")
source("R/figure_utils.R")
source("R/palettes.R")

option_list <- list(
  make_option(c("-b", "--blocks"),      type = "character", help = "04_blocks file [required]"),
  make_option(c("-c", "--crosstraits"), type = "character", default = NULL,
              help = "06_crosstraits_filtered file; sizes the dots [optional]"),
  make_option(c("--chrsizes"),          type = "character",
              default = "data/chromosomes/chrom_sizes.txt"),
  make_option(c("-o", "--outdir"),      type = "character", default = "."),
  make_option(c("--outfile"),           type = "character", default = "block_genome_map"),
  make_option(c("--trait-col"),         type = "character", default = "MAPPED_TRAIT"),
  make_option(c("--cols"),              type = "integer",   default = 1L),
  make_option(c("--label-top"),         type = "integer",   default = 12L),
  make_option(c("--no-size"),           action = "store_true", default = FALSE),
  make_option(c("--base-size"),         type = "double",    default = 11),
  make_option(c("--width"),             type = "double",    default = 9.51),
  make_option(c("--height"),            type = "double",    default = 6.5),
  make_option(c("--dpi"),               type = "integer",   default = 600),
  make_option(c("--svg"),               action = "store_true", default = FALSE),
  make_option(c("--palette"),           type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "LD block positions on the genome"))

if (is.null(opt$blocks)) { cat("ERROR: --blocks is required\n", file = stderr()); quit(status = 2) }
for (f in c(opt$blocks, opt$chrsizes, opt$crosstraits)) {
  if (!is.null(f) && !file.exists(f)) {
    cat("ERROR: file not found: ", f, "\n", sep = "", file = stderr()); quit(status = 2)
  }
}
if (!opt$cols %in% c(1L, 2L)) { cat("ERROR: --cols is 1 or 2\n", file = stderr()); quit(status = 2) }
safe_dir_create(opt$outdir)

# ---- the genome ------------------------------------------------------------
cs <- read.delim(opt$chrsizes, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cs$size <- as.numeric(gsub('[",\r]', "", trimws(cs$size)))
cs$chr  <- sub("^chr", "", trimws(cs$chr), ignore.case = TRUE)
chr_rank <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  ifelse(is.na(n), 100 + match(toupper(x), c("X", "Y", "MT")), n)
}
cs <- cs[order(chr_rank(cs$chr)), ]

# ---- the blocks ------------------------------------------------------------
blocks <- read_any(opt$blocks, as_data_frame = TRUE)
for (nm in c("Block", "Genomic_position")) if (!nm %in% names(blocks)) {
  cat("ERROR: --blocks has no '", nm, "' column\n", sep = "", file = stderr()); quit(status = 2)
}
# "chr6:31834261-33006171" -> chr, start, end. Parsed rather than assumed: a
# block that fails to parse is named, not silently placed at position zero.
m <- regmatches(blocks$Genomic_position,
                regexec("^(?:chr)?([0-9XYMT]+):([0-9]+)-([0-9]+)$",
                        trimws(blocks$Genomic_position)))
ok <- lengths(m) == 4L
if (!all(ok)) {
  cat("ERROR: cannot parse Genomic_position for: ",
      paste(blocks$Block[!ok], collapse = ", "), "\n", sep = "", file = stderr())
  quit(status = 2)
}
bl <- data.frame(
  block = blocks$Block,
  chr   = vapply(m, `[`, character(1), 2),
  start = as.numeric(vapply(m, `[`, character(1), 3)),
  end   = as.numeric(vapply(m, `[`, character(1), 4)),
  stringsAsFactors = FALSE
)
bl$mid    <- (bl$start + bl$end) / 2
bl$span_kb <- (bl$end - bl$start) / 1000

missing_chr <- setdiff(bl$chr, cs$chr)
if (length(missing_chr)) {
  cat("ERROR: block(s) on chromosome(s) absent from ", opt$chrsizes, ": ",
      paste(sort(unique(missing_chr)), collapse = ", "), "\n", sep = "", file = stderr())
  quit(status = 2)
}

# ---- cross-traits, if offered ----------------------------------------------
bl$traits <- NA_integer_
if (!is.null(opt$crosstraits)) {
  ct <- read_any(opt$crosstraits, as_data_frame = TRUE)
  tc <- opt[["trait-col"]]
  if (!all(c("Block", tc) %in% names(ct))) {
    cat("ERROR: --crosstraits needs 'Block' and '", tc, "'\n", sep = "", file = stderr())
    quit(status = 2)
  }
  ct[[tc]] <- trimws(ct[[tc]])
  per <- ct %>% group_by(Block) %>% summarise(n = n_distinct(.data[[tc]]), .groups = "drop")
  bl$traits <- per$n[match(bl$block, per$Block)]
  bl$traits[is.na(bl$traits)] <- 0L
}
size_on <- !opt[["no-size"]] && !all(is.na(bl$traits))

message(sprintf("%d blocks on %d of %d chromosomes; span %.0f-%.0f kb (drawn as points, not to scale)",
                nrow(bl), length(unique(bl$chr)), nrow(cs),
                min(bl$span_kb), max(bl$span_kb)))
empty_chr <- cs$chr[!cs$chr %in% bl$chr]
message("chromosomes with no block: ",
        if (length(empty_chr)) paste(empty_chr, collapse = " ") else "none")
if (size_on) {
  nz <- sum(bl$traits == 0)
  message(sprintf("cross-traits per block: %d-%d; %d block(s) with none, drawn hollow",
                  min(bl$traits), max(bl$traits), nz))
}

# ---- layout ----------------------------------------------------------------
# Chromosome 1 at the top. With two columns the split is by rank, so the left
# column holds the long chromosomes and the right the short ones — the columns
# then have roughly equal ink rather than equal counts.
cs$rank <- seq_len(nrow(cs))
cs$col  <- if (opt$cols == 1L) 1L else ifelse(cs$rank <= ceiling(nrow(cs) / 2), 1L, 2L)

# y is a FACTOR of the chromosome name, not a number. A numeric y with
# facet_wrap gives both columns the same eleven positions, and the scale then
# draws every label in both panels: "1" and "12" land on the same row and
# overprint each other. A discrete scale with scales = "free_y" drops the
# levels a panel does not use, so each column labels only its own chromosomes.
cs$y <- factor(cs$chr, levels = rev(cs$chr))
bl <- left_join(bl, cs[, c("chr", "col", "y", "size")], by = "chr")

pal   <- get_palette(opt$palette)
PAPER <- pal_col(pal, "paper",   "white")
INK   <- pal_col(pal, "ink",     "#0b0b0b")
MUTED <- pal_col(pal, "muted",   "#4b505a")
GRIDC <- pal_col(pal, "grid",    "#e8e8e4")
ACC   <- pal_col(pal, "accent",  "#7a2738")

MB <- 1e6
p <- ggplot() +
  # the chromosome: one rounded segment, which is the whole ideogram
  geom_segment(data = cs,
               aes(x = 0, xend = size / MB, y = y, yend = y),
               colour = GRIDC, linewidth = opt[["base-size"]] * 0.30,
               lineend = "round") +
  facet_wrap(~ col, ncol = opt$cols, scales = "free_y")

if (size_on) {
  p <- p +
    geom_point(data = bl[bl$traits > 0, ],
               aes(x = mid / MB, y = y, size = traits),
               colour = ACC, alpha = 0.85) +
    geom_point(data = bl[bl$traits == 0, ],
               aes(x = mid / MB, y = y),
               shape = 21, fill = PAPER, colour = MUTED, stroke = 0.5,
               size = opt[["base-size"]] * 0.16) +
    scale_size_area(max_size = opt[["base-size"]] * 0.62, name = "Cross-traits",
                    breaks = pretty(c(1, max(bl$traits, na.rm = TRUE)), 4))
} else {
  p <- p + geom_point(data = bl, aes(x = mid / MB, y = y),
                      colour = ACC, size = opt[["base-size"]] * 0.20)
}

n_lab <- opt[["label-top"]]
if (n_lab != 0L) {
  keyed <- if (size_on) bl[order(-bl$traits), ] else bl
  keyed <- if (n_lab > 0L) head(keyed, n_lab) else keyed
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = keyed, aes(x = mid / MB, y = y, label = block),
      size = opt[["base-size"]] * 0.26, colour = INK, seed = 1,
      # direction = "x" pins every label to one band above its own chromosome,
      # so a label can slide sideways past a neighbour but can never drift into
      # the row above. max.overlaps = Inf because a dropped label is worse than
      # a tight one: a reader counts the dots and asks which block is missing.
      nudge_y = 0.42, direction = "x", max.overlaps = Inf,
      min.segment.length = 0.2, segment.colour = MUTED, segment.size = 0.25,
      box.padding = 0.16, point.padding = 0.1)
  } else {
    message("NOTE: ggrepel not installed; labels placed without collision avoidance.")
    p <- p + geom_text(data = keyed, aes(x = mid / MB, y = y + 0.34, label = block),
                       size = opt[["base-size"]] * 0.26, colour = INK)
  }
}

p <- p +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
  scale_x_continuous(name = "Position (Mb)",
                     expand = expansion(mult = c(0.035, 0.02))) +
  theme_minimal(base_size = opt[["base-size"]]) +
  theme(
    plot.background   = element_rect(fill = PAPER, colour = NA),
    panel.background  = element_rect(fill = PAPER, colour = NA),
    panel.grid        = element_blank(),
    panel.grid.major.x = element_line(colour = GRIDC, linewidth = 0.3),
    axis.text.y       = element_text(size = rel(0.9), colour = INK, hjust = 1),
    axis.text.x       = element_text(size = rel(0.8), colour = MUTED),
    axis.title.x      = element_text(size = rel(0.85), colour = MUTED),
    axis.title.y      = element_blank(),
    strip.text        = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = rel(0.85), colour = INK),
    legend.text       = element_text(size = rel(0.8), colour = INK),
    panel.spacing.x   = unit(opt[["base-size"]] * 1.6, "pt")
  )

png_path <- file.path(opt$outdir, paste0(opt$outfile, ".png"))
ggsave(png_path, p, width = opt$width, height = opt$height, dpi = opt$dpi, bg = PAPER)
message("PNG saved -> ", png_path)

if (opt$svg) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("NOTE: svglite not installed; skipping the SVG.")
  } else {
    svg_path <- file.path(opt$outdir, paste0(opt$outfile, ".svg"))
    ggsave(svg_path, p, width = opt$width, height = opt$height,
           device = svglite::svglite, bg = PAPER)
    message("SVG saved -> ", svg_path)
  }
}
