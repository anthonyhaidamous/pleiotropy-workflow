#!/usr/bin/env Rscript
# =============================================================================
# Karyoplotter — trait LD blocks (phenotype-agnostic)
# Produces:
#   - whole-genome ideogram (PDF + SVG) with LD blocks (and, optionally, a
#     focal gene) highlighted
#   - per-chromosome ideograms (PDF)
#   - per-block Manhattan-style LD plots (PDF + PNG + SVG)
#   - combined PDFs per chromosome and one master combined PDF
#
# Inputs (pipeline outputs):
#   --blocks    stage-4 blocks table, 04_blocks.rds OR 04_blocks.tsv
#               (needs Block, Genomic_position, Top_lead_SNP, lSNPs, pSNPs, SNPs)
#   --ldlink    stage-2 LD proxy table, 02_proxy_table.tsv OR .rds
#               (needs RS_Number, Coord, R2)
#   --chrsizes  chromosome sizes file (tab-separated: chr <tab> size; size may
#               be quoted / comma-grouped, e.g. data/chromosomes/chrom_sizes.txt)
#
# All genomic coordinates must share one build (the pipeline uses GRCh38).
#
# Usage:
#   Rscript figure_scripts/fig_karyoploter_trait_ld_blocks.R \
#       --blocks   results/intermediate/leptin/04_blocks.rds \
#       --ldlink   results/intermediate/leptin/02_proxy_table.tsv \
#       --chrsizes data/chromosomes/chrom_sizes.txt \
#       --outdir   figures/leptin/blocks_karyoploteR \
#       --gene     chr7:128241278-128257629:LEP
#
#   # ideograms only (skip per-block Manhattan plots):
#   Rscript figure_scripts/fig_karyoploter_trait_ld_blocks.R ... --no-manhattan
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(karyoploteR)
  library(pdftools)
  library(gtools)
})

option_list <- list(
  make_option(c("--blocks"), type = "character",
              help = "Stage-4 blocks table (04_blocks.rds or 04_blocks.tsv) [required]"),
  make_option(c("--ldlink"), type = "character",
              help = "Stage-2 LD proxy table (02_proxy_table.tsv or .rds) [required]"),
  make_option(c("--chrsizes"), type = "character",
              help = "Chromosome sizes file (2-col TSV: chr, size) [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--gene"), type = "character", default = NULL,
              help = "Optional focal gene to highlight, as 'chr:start-end' or 'chr:start-end:LABEL' (e.g. chr7:128241278-128257629:LEP)"),
  make_option(c("--flank"), type = "double", default = 1e6,
              help = "Flanking region in bp for Manhattan plots [default: %default]"),
  make_option(c("--no-manhattan"), action = "store_true", default = FALSE,
              help = "Skip per-block Manhattan plots (ideograms only)"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Karyoplotter: trait LD block ideograms and Manhattan plots")
opt <- parse_args(parser)

source("R/palettes.R")

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal        <- get_palette(opt$palette)
COL_LEAD   <- pal_col(pal, "pos",     "#B5191A")   # lead-SNP marker
COL_PROXY  <- pal_col(pal, "accent2", "#E8941A")   # proxy span
COL_INK    <- pal_col(pal, "ink",     "black")
COL_SUBTLE <- pal_col(pal, "grid",    "grey")

for (arg in c("blocks", "ldlink", "chrsizes", "outdir")) {
  if (is.null(opt[[arg]])) {
    cat("ERROR: --", arg, " is required\n\n", sep = "", file = stderr())
    print_help(parser)
    quit(status = 2)
  }
}
for (arg in c("blocks", "ldlink", "chrsizes")) {
  if (!file.exists(opt[[arg]])) {
    cat("ERROR: File not found:", opt[[arg]], "\n", file = stderr())
    quit(status = 2)
  }
}

OUTDIR  <- opt$outdir
DIR_PDF <- file.path(OUTDIR, "pdf")
DIR_PNG <- file.path(OUTDIR, "png")
DIR_SVG <- file.path(OUTDIR, "svg")
for (d in c(OUTDIR, DIR_PDF, DIR_PNG, DIR_SVG))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# =============================================================================
# 0. INPUT READERS (accept .rds or tab-separated, flatten list-columns)
# =============================================================================

read_table_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    df <- as.data.frame(readRDS(path), check.names = FALSE, stringsAsFactors = FALSE)
    list_cols <- vapply(df, is.list, logical(1))
    df[list_cols] <- lapply(df[list_cols], function(col)
      vapply(col, function(x) paste(x, collapse = ", "), character(1)))
    df
  } else {
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
               quote = "", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
  }
}

require_cols <- function(df, cols, what) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    cat("ERROR: ", what, " is missing expected column(s): ",
        paste(missing, collapse = ", "), "\n", sep = "", file = stderr())
    quit(status = 2)
  }
}

# Optional focal gene -------------------------------------------------------
gene_gr <- NULL
gene_label <- NULL
if (!is.null(opt$gene)) {
  parts <- trimws(strsplit(opt$gene, ":")[[1]])
  if (length(parts) < 2 || !grepl("^\\d+-\\d+$", parts[2])) {
    cat("ERROR: --gene must look like 'chr:start-end' or 'chr:start-end:LABEL'\n",
        file = stderr())
    quit(status = 2)
  }
  gene_chr   <- parts[1]
  se         <- as.numeric(strsplit(parts[2], "-")[[1]])
  gene_label <- if (length(parts) >= 3) parts[3] else "gene"
  gene_gr    <- toGRanges(data.frame(chr = gene_chr, start = se[1], end = se[2]))
}

# =============================================================================
# 1. READ & PREPROCESS BLOCKS
# =============================================================================

blocks_raw_df <- read_table_any(opt$blocks)
require_cols(blocks_raw_df,
            c("Block", "Genomic_position", "Top_lead_SNP", "lSNPs", "pSNPs", "SNPs"),
            "blocks table")

blocks <- blocks_raw_df %>%
  select(Block, Genomic_position) %>%
  separate(Genomic_position, into = c("chr", "pos"), sep = ":") %>%
  separate(pos, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start), end = as.numeric(end))

blocks_gr    <- makeGRangesFromDataFrame(blocks)
blocks_label <- blocks$Block

# =============================================================================
# 2. CHROMOSOME SIZES
# =============================================================================

chr_data <- read.delim(opt$chrsizes, sep = "\t", header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
names(chr_data)[1:2] <- c("chr", "size")
# Sizes may be quoted and comma-grouped, with trailing CR from CRLF files.
chr_data$size <- as.numeric(gsub('[",\r]', "", trimws(chr_data$size)))
chr_data <- data.frame(chr = chr_data$chr, start = 1, end = chr_data$size,
                       stringsAsFactors = FALSE)
chr_data <- chr_data %>% filter(chr %in% blocks$chr)
genome_gr <- makeGRangesFromDataFrame(chr_data,
                                      seqnames.field = "chr",
                                      start.field    = "start",
                                      end.field      = "end")

# =============================================================================
# SHARED DRAW HELPER — ideogram
# =============================================================================

draw_ideogram <- function(kp) {
  kpRect(kp, chr = blocks$chr, x0 = blocks$start, x1 = blocks$end,
         y0 = 0, y1 = 1, data.panel = "ideogram", col = COL_LEAD, border = NA)
  if (!is.null(gene_gr)) {
    kpRect(kp, chr = as.character(seqnames(gene_gr)),
           x0 = start(gene_gr), x1 = end(gene_gr),
           y0 = 0, y1 = 1, data.panel = "ideogram", col = COL_PROXY, border = NA)
  }
  kpPlotMarkers(kp, data = blocks_gr, labels = blocks_label,
                r0 = 0, r1 = 0.25, marker.parts = c(0.1, 0.05, 0.01),
                label.margin = 1, cex = 1,
                adjust.label.position = TRUE, ignore.chromosome.ends = TRUE)
  if (!is.null(gene_gr)) {
    kpPlotMarkers(kp, data = gene_gr, labels = gene_label,
                  r0 = 0.25, r1 = 0.55, marker.parts = c(0.15, 0.05, 0.01),
                  label.margin = 0.8, cex = 0.55,
                  label.color = COL_INK, line.color = COL_INK,
                  text.orientation = "horizontal",
                  adjust.label.position = FALSE, ignore.chromosome.ends = TRUE)
  }
}

# =============================================================================
# 3. WHOLE-GENOME IDEOGRAM
# =============================================================================

pp <- getDefaultPlotParams(plot.type = 1)
pp$data1outmargin <- 100

pdf(file.path(DIR_PDF, "ideogram_whole_genome.pdf"), height = 14, width = 16)
kp <- plotKaryotype(genome = genome_gr, plot.type = 1,
                    chromosomes = chr_data$chr, plot.params = pp)
draw_ideogram(kp)
dev.off()
message("Saved: ideogram_whole_genome.pdf")

svg(file.path(DIR_SVG, "ideogram_whole_genome.svg"), height = 14, width = 16)
kp <- plotKaryotype(genome = genome_gr, plot.type = 1,
                    chromosomes = chr_data$chr, plot.params = pp)
draw_ideogram(kp)
dev.off()
message("Saved: ideogram_whole_genome.svg")

# =============================================================================
# 4. PER-CHROMOSOME IDEOGRAMS
# =============================================================================

for (chr in unique(blocks$chr)) {
  sel              <- blocks$chr == chr
  blocks_chr       <- blocks[sel, ]
  blocks_gr_chr    <- blocks_gr[seqnames(blocks_gr) == chr]
  blocks_label_chr <- blocks_label[sel]
  genome_gr_chr    <- genome_gr[seqnames(genome_gr) == chr]

  pdf(file.path(DIR_PDF, paste0("ideogram_", chr, ".pdf")), height = 2, width = 16)
  pp_chr <- getDefaultPlotParams(plot.type = 1)
  pp_chr$data1outmargin <- 100
  kp <- plotKaryotype(genome = genome_gr_chr, plot.type = 1,
                      chromosomes = chr, plot.params = pp_chr)
  kpRect(kp, chr = blocks_chr$chr, x0 = blocks_chr$start, x1 = blocks_chr$end,
         y0 = 0, y1 = 1, data.panel = "ideogram", col = COL_INK, border = NA)
  kpPlotMarkers(kp, data = blocks_gr_chr, labels = blocks_label_chr,
                r0 = 0, r1 = 0.25, marker.parts = c(0.1, 0.05, 0.01),
                label.margin = 1, cex = 1,
                adjust.label.position = TRUE, ignore.chromosome.ends = TRUE)
  dev.off()
}
message("Saved per-chromosome ideograms")

if (opt[["no-manhattan"]]) {
  message("--no-manhattan set; skipping block Manhattan plots.")
  quit(status = 0)
}

# =============================================================================
# 5. LEAD SNP LIST
# =============================================================================

# Top_lead_SNP is "rsID, pos, pval, ORbeta" — keep the rsID only.
main_leads <- trimws(sub(",.*", "", blocks_raw_df$Top_lead_SNP))

# =============================================================================
# 6. LD PROXIES → GRANGES
# =============================================================================

ldlink_raw <- read_table_any(opt$ldlink)
require_cols(ldlink_raw, c("RS_Number", "Coord", "R2"), "LD proxy table")
ldlink <- ldlink_raw[, c("RS_Number", "Coord", "R2")]
names(ldlink) <- c("SNP", "Coord", "R2")
ldlink <- ldlink %>%
  separate(Coord, into = c("seqnames", "ranges"), sep = ":") %>%
  mutate(ranges = as.numeric(ranges), seqnames = as.character(seqnames), R2 = as.numeric(R2))

snps <- GRanges(seqnames = ldlink$seqnames,
                ranges   = IRanges(start = ldlink$ranges, width = 1),
                SNP      = ldlink$SNP,
                R2       = ldlink$R2)

# =============================================================================
# 7. BLOCKS WITH SNP MEMBERSHIP
# =============================================================================

blocks_with_SNPs <- blocks_raw_df %>%
  select(Block, Genomic_position, lSNPs, pSNPs, SNPs) %>%
  separate(Genomic_position, into = c("chr", "pos"), sep = ":") %>%
  separate(pos, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start), end = as.numeric(end))

# =============================================================================
# 8. BLOCK PLOT HELPER
# =============================================================================

draw_block_plot <- function(kp, bchr, bstart, bend, zoom_reg,
                             outside_snps, inside_snps, others, main_snp, block_name) {
  kpAddMainTitle(kp, paste0(block_name, " — chr", bchr), cex = 2)
  kpRect(kp, chr = bchr, x0 = start(zoom_reg), x1 = end(zoom_reg),
         y0 = 0, y1 = 1, col = COL_SUBTLE, border = NA, data.panel = "ideogram")
  kpAddBaseNumbers(kp, tick.dist = (end(zoom_reg) - start(zoom_reg)) / 2, add.units = TRUE)
  kpAxis(kp, ymax = 1.1, tick.pos = c(0.1, 0.9, 1))
  kpAddLabels(kp, labels = "Correlation Coefficient (R²)", srt = 90, pos = 1, cex = 1.3,
              r0 = 0, r1 = 1, label.margin = 0.02)
  kpAbline(kp, h = c(0.9, 0.1), col = "darkgrey", lty = c(2, 3), ymax = 1.1)
  kpAbline(kp, v = c(bstart, bend), col = COL_SUBTLE, lty = 2)
  kpRect(kp, chr = bchr, x0 = bstart, x1 = bend,
         y0 = 0, y1 = 1, col = COL_INK, border = NA, data.panel = "ideogram")
  if (length(outside_snps) > 0)
    kpPoints(kp, data = outside_snps, pch = 19, cex = 1.2,
             col = "lightgrey", y = outside_snps$R2, ymax = 1.1)
  if (length(inside_snps) > 0) {
    cols <- ifelse(inside_snps$R2 >= 0.9, "red2",
                   ifelse(inside_snps$R2 > 0.1, "steelblue", "lightskyblue"))
    kpPoints(kp, data = inside_snps, pch = 19, cex = 1.2,
             col = cols, y = inside_snps$R2, ymax = 1.1)
  }
  if (length(others) > 0)
    kpPoints(kp, data = others, pch = 19, cex = 1.5, col = "red4",
             y = others$R2, ymax = 1.1)
  if (!is.null(main_snp)) {
    kpPoints(kp, data = main_snp, pch = 19, cex = 1.5, col = COL_INK,
             y = main_snp$R2, ymax = 1.1)
    kpText(kp, data = main_snp, labels = main_snp$SNP,
           y = main_snp$R2, ymax = 1.1, pos = 3, col = COL_INK)
  }
}

# =============================================================================
# 9. GENERATE MANHATTAN PLOTS
# =============================================================================

mk_ids <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  ids <- trimws(unlist(strsplit(x, ","), use.names = FALSE))
  ids[ids != ""]
}

logfile <- file.path(OUTDIR, "logfile.log")
if (file.exists(logfile)) file.remove(logfile)
cat("Manhattan plot generation log\nStarted at:", format(Sys.time()), "\n\n",
    file = logfile, append = TRUE)

pb <- utils::txtProgressBar(min = 0, max = nrow(blocks), style = 3)
start_time <- Sys.time()

for (i in seq_len(nrow(blocks))) {
  utils::setTxtProgressBar(pb, i)
  block_name <- blocks$Block[i]
  bchr       <- blocks$chr[i]
  bstart     <- blocks$start[i]
  bend       <- blocks$end[i]

  cat(sprintf("--- Block %d (%s) ---\nchr: %s   start-end: %d-%d\n",
              i, block_name, bchr, bstart, bend), file = logfile, append = TRUE)

  row        <- blocks_with_SNPs[i, ]
  lead_ids   <- mk_ids(row$lSNPs)
  proxy_ids  <- mk_ids(row$pSNPs)
  other_ids  <- mk_ids(row$SNPs)
  inside_ids <- unique(c(lead_ids, proxy_ids, other_ids))

  snps_chr     <- snps[seqnames(snps) == bchr]
  inside_snps  <- snps_chr[snps_chr$SNP %in% inside_ids]
  outside_snps <- snps_chr[!snps_chr$SNP %in% inside_ids]
  top_snps     <- inside_snps[inside_snps$R2 == 1]
  cand         <- top_snps[top_snps$SNP %in% main_leads]

  if (length(cand) > 0) {
    center   <- (bstart + bend) / 2
    dists    <- abs(start(cand) - center)
    main_snp <- cand[which.min(dists)]
  } else {
    main_snp <- NULL
  }

  lead_top <- top_snps[top_snps$SNP %in% lead_ids]
  others   <- if (!is.null(main_snp)) lead_top[lead_top$SNP != main_snp$SNP] else lead_top

  zoom_reg  <- toGRanges(data.frame(chr = bchr, start = bstart - opt$flank, end = bend + opt$flank))
  pp_block  <- getDefaultPlotParams(plot.type = 4)
  pp_block$ideogram.height <- 0.2
  base_name <- paste0("manhattan_", block_name, "_chr", bchr)

  pdf(file.path(DIR_PDF, paste0(base_name, ".pdf")), width = 12, height = 8)
  kp <- plotKaryotype(chromosomes = bchr, zoom = zoom_reg, cytobands = NULL, plot.params = pp_block)
  draw_block_plot(kp, bchr, bstart, bend, zoom_reg, outside_snps, inside_snps, others, main_snp, block_name)
  dev.off()

  png(file.path(DIR_PNG, paste0(base_name, ".png")), width = 12, height = 8, units = "in", res = 300)
  kp <- plotKaryotype(chromosomes = bchr, zoom = zoom_reg, cytobands = NULL, plot.params = pp_block)
  draw_block_plot(kp, bchr, bstart, bend, zoom_reg, outside_snps, inside_snps, others, main_snp, block_name)
  dev.off()

  svg(file.path(DIR_SVG, paste0(base_name, ".svg")), width = 12, height = 8)
  kp <- plotKaryotype(chromosomes = bchr, zoom = zoom_reg, cytobands = NULL, plot.params = pp_block)
  draw_block_plot(kp, bchr, bstart, bend, zoom_reg, outside_snps, inside_snps, others, main_snp, block_name)
  dev.off()
}

close(pb)
elapsed <- difftime(Sys.time(), start_time, units = "secs")
cat(sprintf("Finished at: %s\nElapsed: %.2f seconds\n", format(Sys.time()), as.numeric(elapsed)),
    file = logfile, append = TRUE)
message("Done. Log saved to: ", logfile)

# =============================================================================
# 10. COMBINE PDFs
# =============================================================================

chroms <- unique(blocks$chr)

for (chr in chroms) {
  chr_pdfs <- list.files(DIR_PDF, pattern = paste0("manhattan_B.*_chr", chr, "\\.pdf$"),
                         full.names = TRUE)
  if (length(chr_pdfs) > 0) {
    chr_pdfs <- mixedsort(chr_pdfs)
    pdf_combine(chr_pdfs, output = file.path(DIR_PDF, paste0("combined_", chr, ".pdf")))
  }
}

for (chr in chroms) {
  ideogram_pdf <- file.path(DIR_PDF, paste0("ideogram_", chr, ".pdf"))
  combined_pdf <- file.path(DIR_PDF, paste0("combined_", chr, ".pdf"))
  merged_pdf   <- file.path(DIR_PDF, paste0("combined_ideo_", chr, ".pdf"))
  if (file.exists(ideogram_pdf) && file.exists(combined_pdf))
    pdf_combine(c(ideogram_pdf, combined_pdf), output = merged_pdf)
}

all_pdfs <- list.files(DIR_PDF, pattern = "manhattan_B.*\\.pdf$", full.names = TRUE)
all_pdfs <- mixedsort(all_pdfs)
pdf_combine(all_pdfs, output = file.path(DIR_PDF, "combined_all_blocks.pdf"))
message("All done. Combined PDF: ", file.path(DIR_PDF, "combined_all_blocks.pdf"))
