#!/usr/bin/env Rscript
# ============================================================
# 10_extract_genes.R
#
# Builds a protein-coding gene list for each LD block from two sources:
#   List A  biomaRt query on the block's Genomic_position (± flank)
#   List B  gene symbols already in the 04_blocks Genes column
# Both lists are merged and deduplicated per block.
#
# Inputs:
#   --blocks  04_blocks.{tsv,rds}   (required)
#
# Outputs (written to --outdir):
#   protein_coding_genes_by_block.tsv   one row per gene × block
#   protein_coding_genes_collapsed.tsv  one row per block, genes comma-joined
#   genes_list_complete.txt             flat unique gene list (one per line)
#                                       → feed this into Enrichr
#
# Usage:
#   Rscript pipeline/10_extract_genes.R \
#       --blocks results/intermediate/leptin/04_blocks.rds \
#       --outdir results/intermediate/leptin/enrichr
#
#   # wider flanks, specific Ensembl version, useast mirror:
#   Rscript pipeline/10_extract_genes.R \
#       --blocks results/intermediate/leptin/04_blocks.rds \
#       --outdir results/intermediate/leptin/enrichr \
#       --flank 25000 --ensembl_version 112 --mirror useast
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt", ask = FALSE, update = FALSE)

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

source("R/figure_utils.R")   # safe_dir_create, read_any
source("R/utils.R")          # write_stage_manifest

# ============================================================
# Helpers
# ============================================================

# safe_dir_create() and read_any() live in R/figure_utils.R.

as_int <- function(x) suppressWarnings(as.integer(as.character(x)))

stop_if_missing <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    cat("ERROR: Missing column(s):\n - ", paste(miss, collapse = "\n - "),
        "\nAvailable columns:\n - ", paste(names(df), collapse = "\n - "), "\n",
        sep = "", file = stderr())
    quit(status = 2)
  }
}

parse_genomic_pos <- function(x) {
  x <- gsub("\\s+", "", x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  m <- regexec("^([^:]+):([0-9]+)-([0-9]+)$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) != 4)
    return(list(chr = NA_character_, start = NA_integer_, end = NA_integer_))
  list(chr = r[2], start = as_int(r[3]), end = as_int(r[4]))
}

keep_gene <- function(gstart, gend, rstart, rend, mode = c("overlap", "contained")) {
  mode <- match.arg(mode)
  if (mode == "overlap") !(gend < rstart | gstart > rend)
  else                   gstart >= rstart & gend <= rend
}

getBM_retry <- function(..., retries = 3, sleep_sec = 1) {
  last_err <- NULL
  for (k in seq_len(retries)) {
    out <- tryCatch(getBM(...), error = function(e) { last_err <<- e; NULL })
    if (!is.null(out)) return(out)
    Sys.sleep(sleep_sec)
  }
  stop(last_err)
}

# Connect to Ensembl with automatic mirror fallback.
#
# biomaRt limitation: `version` and `mirror` are mutually exclusive — when a
# version is set, biomaRt silently ignores the mirror and hits the archive host.
# Valid current-release mirrors: www (default/""), useast, asia.
connect_ensembl <- function(version = NULL, mirror = "", retries = 3) {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 20)  # fail fast per attempt instead of the default 60s

  if (!is.null(version)) {
    message("  Archive version — mirrors not supported; retrying main site only.")
    last_err <- NULL
    for (attempt in seq_len(retries)) {
      message("  Attempt ", attempt, "/", retries, " (v", version, ") ...")
      result <- tryCatch(
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = version),
        error = function(e) { last_err <<- e; NULL }
      )
      if (!is.null(result)) { message("  Connected."); return(result) }
      if (attempt < retries) Sys.sleep(2)
    }
    stop("Archive v", version, " unreachable. Omit --ensembl_version to use current ",
         "release (mirror fallback supported), or use --method ensdb.\n",
         "Last error: ", conditionMessage(last_err))
  }

  # Validate mirror; valid values for biomaRt: www (=""), useast, asia
  if (nzchar(mirror) && !mirror %in% c("www", "useast", "asia"))
    stop("Invalid --mirror '", mirror, "'. Choose from: www, useast, asia")
  mirrors <- if (nzchar(mirror)) mirror else c("", "useast", "asia")

  last_err <- NULL
  for (m in mirrors) {
    label <- if (nzchar(m)) m else "www"
    for (attempt in seq_len(retries)) {
      message("  Trying ", label, " (attempt ", attempt, "/", retries, ") ...")
      result <- tryCatch(
        if (nzchar(m))
          useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = m)
        else
          useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl"),
        error = function(e) { last_err <<- e; NULL }
      )
      if (!is.null(result)) { message("  Connected via ", label); return(result) }
      if (attempt < retries) Sys.sleep(2)
    }
    message("  Mirror '", label, "' unreachable.")
  }
  stop("All Ensembl mirrors unreachable. Try --method ensdb for offline annotation.\n",
       "Last error: ", conditionMessage(last_err))
}

# Query genes from a local EnsDb via AnnotationHub (downloaded + cached on first use).
# Returns a data.frame with the same columns expected by the downstream merge.
query_genes_ensdb <- function(blocks, mode = "overlap", ensembl_version = NULL) {
  for (p in c("AnnotationHub", "ensembldb", "GenomicRanges")) {
    if (!requireNamespace(p, quietly = TRUE))
      BiocManager::install(p, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages({
    library(AnnotationHub)
    library(ensembldb)
    library(GenomicRanges)
  })

  message("  Loading EnsDb from AnnotationHub (downloaded + cached on first use) ...")
  ah  <- AnnotationHub(ask = FALSE)
  recs <- query(ah, c("EnsDb", "Homo sapiens", "GRCh38"))
  if (length(recs) == 0) stop("No EnsDb GRCh38 records found in AnnotationHub.")

  # Pick exact version if requested, otherwise most recent
  ahid <- if (!is.null(ensembl_version)) {
    hits <- names(recs)[grepl(paste0("v", ensembl_version, "\\."), recs$title, ignore.case = TRUE)]
    if (length(hits)) hits[1] else { message("  v", ensembl_version, " not found; using latest."); names(recs)[length(recs)] }
  } else {
    names(recs)[length(recs)]
  }
  edb <- ah[[ahid]]
  message("  Using: ", recs[ahid]$title)

  regions_gr <- GRanges(
    seqnames = blocks$chromosome,
    ranges   = IRanges(start = blocks$query_start, end = blocks$query_end)
  )

  message("  Querying ", nrow(blocks), " block regions ...")
  genes_gr <- genes(edb,
                    filter  = GRangesFilter(regions_gr),
                    columns = c("gene_id", "gene_name", "gene_biotype"))

  if (length(genes_gr) == 0) {
    message("  No genes returned by EnsDb query.")
    return(data.frame(chromosome_name = character(0), start_position = integer(0),
                      end_position = integer(0), ensembl_gene_id = character(0),
                      external_gene_name = character(0), gene_biotype = character(0),
                      stringsAsFactors = FALSE))
  }

  df <- as.data.frame(genes_gr)
  df$chromosome_name    <- as.character(df$seqnames)
  df$start_position     <- as.integer(df$start)
  df$end_position       <- as.integer(df$end)
  df$ensembl_gene_id    <- as.character(df$gene_id)
  df$external_gene_name <- ifelse(is.na(df$gene_name) | df$gene_name == "",
                                   df$gene_id, as.character(df$gene_name))
  df$gene_biotype       <- as.character(df$gene_biotype)
  df[, c("chromosome_name", "start_position", "end_position",
         "ensembl_gene_id", "external_gene_name", "gene_biotype")]
}


# ------------------------------------------------------------
# Biotype lookup for List B symbols
#
# List A comes from a region query and is biotype-filtered at source. List B is
# whatever sits in the catalog's MAPPED_GENE column, which is NOT restricted to
# protein-coding genes -- for leptin it contributed 12 pseudogenes, lncRNAs and
# a miRNA (LINC-PINT, MIR129-1, RSL24D1P4, SEPTIN14P21, ...) to a file named
# "protein_coding". They then appear as STRING "unmapped" genes, inflating the
# reported submission from 93 to 105 and the unmapped count from 5 to 17.
#
# This resolves each List B symbol's biotype so the union can be filtered.
# Symbols that cannot be resolved at all are RETAINED and reported, so a
# lookup failure never silently deletes a real gene.
# ------------------------------------------------------------
protein_coding_symbols <- function(symbols, mart) {
  symbols <- unique(symbols[nzchar(symbols) & !is.na(symbols)])
  if (length(symbols) == 0) return(character(0))
  bt <- tryCatch(
    getBM_retry(attributes = c("external_gene_name", "gene_biotype"),
                filters    = "external_gene_name",
                values     = symbols,
                mart       = mart),
    error = function(e) NULL
  )
  if (is.null(bt) || nrow(bt) == 0) {
    message("  WARNING: biotype lookup for List B returned nothing; ",
            "keeping all ", length(symbols), " symbol(s) unfiltered.")
    return(symbols)
  }
  coding     <- unique(bt$external_gene_name[bt$gene_biotype == "protein_coding"])
  resolved   <- unique(bt$external_gene_name)
  unresolved <- setdiff(symbols, resolved)
  dropped    <- setdiff(intersect(symbols, resolved), coding)

  if (length(dropped) > 0) {
    message("  List B biotype filter: dropped ", length(dropped),
            " non-protein-coding symbol(s):")
    message("    ", paste(sort(dropped), collapse = ", "))
  }
  if (length(unresolved) > 0) {
    message("  List B biotype filter: ", length(unresolved),
            " symbol(s) could not be resolved and are KEPT: ",
            paste(sort(unresolved), collapse = ", "))
  }
  sort(unique(c(intersect(symbols, coding), unresolved)))
}

# Parses one Genes cell — handles:
#   RDS list-column element : character vector, e.g. c("RSL24D1P4 - SEPTIN14P21", "GCKR")
#   TSV string              : e.g. "RSL24D1P4 - SEPTIN14P21, GCKR"
# Splits on comma and " - " (flanking-gene notation), deduplicates.
parse_block_genes <- function(x) {
  if (is.list(x)) x <- unlist(x)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0) return(character(0))
  parts <- unlist(strsplit(x, "\\s*,\\s*|\\s+-\\s+"))
  parts <- trimws(parts)
  unique(parts[nzchar(parts)])
}

# ============================================================
# CLI
# ============================================================

option_list <- list(
  make_option(c("-b", "--blocks"), type = "character",
              help = "04_blocks file (.tsv or .rds) [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option(c("--flank"), type = "integer", default = 10000L,
              help = "Flank in bp added to both sides of each block [default: %default]"),
  make_option(c("--mode"), type = "character", default = "overlap",
              help = "Gene selection: overlap|contained [default: %default]"),
  make_option(c("--ensembl_version"), type = "integer", default = NULL,
              help = "Ensembl archive version (e.g. 112). Omit to use the current release [default: current]"),
  make_option(c("--mirror"), type = "character", default = "",
              help = "Optional Ensembl mirror: useast|uswest|asia [default: main]"),
  make_option(c("--retries"), type = "integer", default = 3L,
              help = "Number of retries for Ensembl queries [default: %default]"),
  make_option(c("--method"), type = "character", default = "biomart",
              help = "Annotation source: biomart (live API) or ensdb (AnnotationHub local cache) [default: %default]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Extract protein-coding genes per LD block (biomaRt + Genes column)")
opt <- parse_args(parser)

if (is.null(opt$blocks)) {
  cat("ERROR: --blocks is required\n\n", file = stderr())
  print_help(parser); quit(status = 2)
}
if (!file.exists(opt$blocks)) {
  cat("ERROR: File not found:", opt$blocks, "\n", file = stderr()); quit(status = 2)
}
opt$mode <- match.arg(opt$mode, c("overlap", "contained"))
safe_dir_create(opt$outdir)

message("blocks          = ", opt$blocks)
message("outdir          = ", opt$outdir)
message("flank (bp)      = ", opt$flank)
message("mode            = ", opt$mode)
message("ensembl_version = ", if (is.null(opt$ensembl_version)) "current" else opt$ensembl_version)
message("mirror          = ", if (nzchar(opt$mirror)) opt$mirror else "(default)")
message("retries         = ", opt$retries)
message("method          = ", opt$method)

opt$method <- match.arg(opt$method, c("biomart", "ensdb"))

# ============================================================
# Read 04_blocks
# ============================================================

blocks_in <- read_any(opt$blocks)

# Tolerate "Genomic Position" (space) as well as "Genomic_position"
if (!("Genomic_position" %in% names(blocks_in)) &&
    "Genomic Position" %in% names(blocks_in)) {
  names(blocks_in)[names(blocks_in) == "Genomic Position"] <- "Genomic_position"
}

stop_if_missing(blocks_in, c("Block", "Genomic_position", "Genes"))

# Parse coordinate strings
parsed <- lapply(blocks_in$Genomic_position, parse_genomic_pos)
blocks <- data.frame(
  Block      = blocks_in$Block,
  chromosome = vapply(parsed, `[[`, character(1), "chr"),
  start      = vapply(parsed, `[[`, integer(1),   "start"),
  end        = vapply(parsed, `[[`, integer(1),   "end"),
  stringsAsFactors = FALSE
)

bad <- which(is.na(blocks$chromosome) | is.na(blocks$start) | is.na(blocks$end))
if (length(bad) > 0) {
  cat("ERROR: Failed to parse Genomic_position for block(s):\n", file = stderr())
  print(head(blocks_in[bad, c("Block", "Genomic_position")], 10))
  quit(status = 2)
}

flank <- as_int(opt$flank)
blocks$query_start <- pmax(1L, blocks$start - flank)
blocks$query_end   <- blocks$end + flank

# List B: genes from the Genes column (handles list column or string)
block_gene_map <- data.frame(Block = blocks_in$Block, stringsAsFactors = FALSE)
block_gene_map$genes_from_block <- lapply(blocks_in$Genes, parse_block_genes)

# ============================================================
# List A: query annotation source for all blocks
# ============================================================

if (opt$method == "biomart") {
  message("\nConnecting to Ensembl (",
          if (is.null(opt$ensembl_version)) "current release" else paste0("v", opt$ensembl_version),
          ") ...")
  ensembl <- connect_ensembl(version = opt$ensembl_version,
                             mirror  = opt$mirror,
                             retries = opt$retries)

  regions <- paste0(blocks$chromosome, ":", blocks$query_start, ":", blocks$query_end)
  message("Querying ", nrow(blocks), " blocks in one request ...")
  message("  regions: ", paste(regions, collapse = "; "))

  raw_all <- getBM_retry(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name",
                   "start_position", "end_position", "gene_biotype"),
    filters    = "chromosomal_region",
    values     = regions,
    mart       = ensembl,
    retries    = opt$retries,
    sleep_sec  = 2
  )
  message("  Received ", nrow(raw_all), " gene records across all blocks")

  raw_all$ensembl_gene_id    <- as.character(raw_all$ensembl_gene_id)
  raw_all$external_gene_name <- as.character(raw_all$external_gene_name)
  raw_all$chromosome_name    <- as.character(raw_all$chromosome_name)
  raw_all$gene_biotype       <- as.character(raw_all$gene_biotype)
  raw_all$start_position     <- as_int(raw_all$start_position)
  raw_all$end_position       <- as_int(raw_all$end_position)

} else {
  message("\nUsing EnsDb (AnnotationHub) ...")
  raw_all <- query_genes_ensdb(blocks,
                               mode             = opt$mode,
                               ensembl_version  = opt$ensembl_version)
  message("  Received ", nrow(raw_all), " gene records across all blocks")
}

empty_frame <- data.frame(
  Block = character(0), chromosome = character(0),
  region_start = integer(0), region_end = integer(0),
  ensembl_gene_id = character(0), external_gene_name = character(0),
  chromosome_name = character(0), start_position = integer(0),
  end_position = integer(0), gene_biotype = character(0),
  stringsAsFactors = FALSE
)

# Assign each gene to the block(s) whose region it overlaps
res_list <- lapply(seq_len(nrow(blocks)), function(i) {
  b      <- blocks[i, ]
  rstart <- b$query_start
  rend   <- b$query_end

  sub <- raw_all[
    raw_all$chromosome_name == b$chromosome &
    !is.na(raw_all$start_position) & !is.na(raw_all$end_position),
  , drop = FALSE]

  if (nrow(sub) > 0) {
    keep <- mapply(keep_gene, sub$start_position, sub$end_position,
                   MoreArgs = list(rstart = rstart, rend = rend, mode = opt$mode))
    sub <- sub[keep, , drop = FALSE]
  }

  if (nrow(sub) == 0) return(empty_frame)

  sub$Block        <- b$Block
  sub$chromosome   <- b$chromosome
  sub$region_start <- rstart
  sub$region_end   <- rend
  sub
})

# Long table: protein-coding only, one row per gene×block
genes_long <- bind_rows(res_list) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  mutate(
    external_gene_name = ifelse(
      is.na(external_gene_name) | external_gene_name == "",
      ensembl_gene_id, external_gene_name
    )
  ) %>%
  distinct(Block, ensembl_gene_id, .keep_all = TRUE) %>%
  arrange(Block, chromosome_name, start_position, external_gene_name)

out_long <- file.path(opt$outdir, "protein_coding_genes_by_block.tsv")
write.table(genes_long, out_long, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("\nWrote: ", out_long)

# ============================================================
# Merge List A (biomaRt) + List B (Genes column)
# ============================================================

biomart_collapsed <- genes_long %>%
  group_by(Block) %>%
  summarise(biomart_genes = list(unique(external_gene_name)), .groups = "drop")

# Biotype-filter List B before the union. List A is already protein-coding-only
# (filtered at source); without this step the merged output is labelled
# "protein_coding" while carrying whatever biotype the catalog happened to name.
listB_all <- unique(unlist(block_gene_map$genes_from_block))
listB_all <- listB_all[!is.na(listB_all) & nzchar(listB_all)]
already_coding <- unique(genes_long$external_gene_name)
listB_to_check <- setdiff(listB_all, already_coding)

message("\nList B (catalog Genes column): ", length(listB_all), " distinct symbol(s), ",
        length(listB_to_check), " not already in List A.")
listB_keep <- if (opt$method == "biomart") {
  c(intersect(listB_all, already_coding),
    protein_coding_symbols(listB_to_check, ensembl))
} else {
  # The EnsDb object is scoped inside query_genes_ensdb(), so there is no
  # handle to look symbols up against here. Retain everything rather than
  # guess, and say so — see .context/todo.md.
  message("  NOTE: --method ensdb does not biotype-filter List B; all ",
          length(listB_to_check), " extra symbol(s) retained.")
  listB_all
}

collapsed <- blocks %>%
  left_join(biomart_collapsed, by = "Block") %>%
  left_join(block_gene_map,   by = "Block") %>%
  mutate(
    biomart_genes    = lapply(biomart_genes,    function(x)
      if (is.null(x) || all(is.na(x))) character(0) else x),
    genes_from_block = lapply(genes_from_block, function(x)
      if (is.null(x) || all(is.na(x))) character(0) else intersect(x, listB_keep)),
    merged_genes = Map(function(a, b) unique(c(a, b)), biomart_genes, genes_from_block),
    n_protein_coding_genes = vapply(merged_genes, length, integer(1)),
    protein_coding_genes   = vapply(merged_genes,
      function(x) paste(sort(unique(x)), collapse = ", "), character(1))
  ) %>%
  dplyr::select(
    Block,
    chromosome,
    region_start = query_start,
    region_end   = query_end,
    n_protein_coding_genes,
    protein_coding_genes
  ) %>%
  arrange(as.integer(sub("^\\D+", "", Block)), Block)

out_collapsed <- file.path(opt$outdir, "protein_coding_genes_collapsed.tsv")
write.table(collapsed, out_collapsed, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote: ", out_collapsed)

# ============================================================
# Flat gene list (for Enrichr)
# ============================================================

all_genes <- collapsed %>%
  dplyr::filter(nzchar(protein_coding_genes)) %>%
  separate_rows(protein_coding_genes, sep = ",\\s*") %>%
  pull(protein_coding_genes) %>%
  trimws() %>%
  unique() %>%
  sort()

out_genes <- file.path(opt$outdir, "genes_list_complete.txt")
writeLines(all_genes, out_genes)

# ============================================================
# STRING submission list
#
# A bare Ensembl gene ID is a real gene and counts toward the gene total, but
# it is NOT a valid STRING input: STRING is queried by symbol, so an ENSG
# reaches the mapper, fails, and lands in the "unmapped" bucket. Passing it
# silently makes the submission count look larger than the number of genes
# that could ever have mapped, and hides a data-quality problem (a gene whose
# symbol was never resolved) inside a number that looks like a STRING result.
#
# Same principle as the exclusion columns in stage 6: the gene is retained and
# counted, the exclusion is by rule, and it is named in the log and in a file
# that reaches the supplement rather than being dropped invisibly.
#
# Affects both phenotypes: 1 ID in leptin, 8 in CLL.
# ============================================================

# A submittable identifier is a gene SYMBOL. Bare Ensembl accessions are the
# only non-symbol form present in either phenotype today; the pattern is
# anchored so that symbols legitimately beginning with "ENSG"-like text are
# unaffected, and it is deliberately narrow — an unknown identifier form is
# submitted (and shows up as STRING-unmapped) rather than silently dropped
# here, so a new form is visible rather than absorbed.
is_bare_ensembl <- function(x) grepl("^ENS[A-Z]*G[0-9]{6,}(\\.[0-9]+)?$", x)

string_excluded  <- all_genes[is_bare_ensembl(all_genes)]
genes_submission <- setdiff(all_genes, string_excluded)

out_submission <- file.path(opt$outdir, "genes_string_submission.txt")
writeLines(genes_submission, out_submission)

out_string_excl <- file.path(opt$outdir, "genes_excluded_from_string.tsv")
write.table(
  data.frame(
    gene_id = string_excluded,
    reason  = rep("bare_ensembl_id_not_a_string_symbol", length(string_excluded)),
    retained_in_gene_table = rep(TRUE, length(string_excluded)),
    stringsAsFactors = FALSE
  ),
  out_string_excl, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

message("\nSTRING submission list:")
message("  genes retained in gene table : ", length(all_genes))
message("  submitted to STRING          : ", length(genes_submission))
if (length(string_excluded) > 0) {
  message("  excluded from submission (", length(string_excluded),
          "), retained in the gene table and written to the unmapped list:")
  for (g in string_excluded) message("    - ", g)
} else {
  message("  excluded from submission     : 0")
}
message("  wrote: ", out_submission)
message("  wrote: ", out_string_excl)

# Stage 10 is CLI-driven rather than config-driven, so the manifest gets a
# minimal cfg shim carrying the phenotype slug implied by --outdir.
.cfg10 <- list(
  phenotype = list(name = basename(opt$outdir), slug = basename(opt$outdir)),
  paths     = list(intermediate = opt$outdir)
)
attr(.cfg10, "config_path") <- paste("(CLI)", paste(commandArgs(TRUE), collapse = " "))
attr(.cfg10, "config_md5")  <- NA_character_
write_stage_manifest(
  stage   = 10,
  cfg     = .cfg10,
  inputs  = opt$blocks,
  outputs = c(out_long, out_collapsed, out_genes, out_submission, out_string_excl),
  counts  = list(
    n_blocks                    = nrow(blocks),
    n_listA_protein_coding      = length(unique(genes_long$external_gene_name)),
    n_listB_symbols             = length(listB_all),
    n_listB_checked             = length(listB_to_check),
    n_listB_kept                = length(listB_keep),
    n_genes_final               = length(all_genes),
    n_genes_string_submitted    = length(genes_submission),
    n_genes_string_excluded     = length(string_excluded),
    string_excluded_ids         = as.list(string_excluded),
    flank                       = as_int(opt$flank),
    method                      = opt$method
  ),
  notes = paste("List B biotype-filtered against Ensembl; unresolvable symbols",
                "retained. Bare Ensembl IDs counted in the gene table but",
                "excluded from the STRING submission."))
message("Wrote: ", out_genes)

message("\nDone.  Blocks: ", nrow(blocks),
        " | Unique protein-coding genes across all blocks: ", length(all_genes))
