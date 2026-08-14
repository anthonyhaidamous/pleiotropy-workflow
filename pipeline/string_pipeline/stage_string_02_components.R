#!/usr/bin/env Rscript
# =============================================================================
# stage_string_02_components.R
#
# Builds the igraph object from the frozen edge list and produces the full
# gene accounting. Every gene that leaves the analysis appears by name in the
# output — nothing is dropped invisibly.
#
# The graph is built over ALL mapped proteins, including zero-degree ones.
# Building from the edge list alone would drop them before they could be
# counted, which is exactly the silent filter this stage exists to prevent.
#
# Which node set flows downstream (the largest connected component, or every
# gene with degree >= 1) is a config parameter, not a decision baked into code.
#
# Inputs:
#   config/<phenotype>.yml
#   <edges_dir>/<slug>_edges_v<ver>_s<thr>.tsv      + its .provenance.json
#   <edges_dir>/<slug>_gene_mapping_v<ver>.tsv
#   string.input.block_table                        (optional, for LD blocks)
#
# Outputs (to results/intermediate/<slug>/<stage_subdir>/):
#   string_02_graph_full.rds          igraph over all mapped proteins
#   string_02_graph_downstream.rds    the selected node set
#   string_02_accounting.tsv          one row per submitted gene
#   string_02_accounting_summary.tsv  the six Q1b quantities
#
# Usage:
#   Rscript pipeline/string_pipeline/stage_string_02_components.R config/leptin.yml
#
# Run from the repository root.
# =============================================================================

if (!file.exists("R/utils.R")) {
  stop("run this script from the repository root (R/utils.R not found)")
}
source("R/utils.R")
source("R/string_utils.R")
suppressPackageStartupMessages(library(igraph))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript pipeline/string_pipeline/stage_string_02_components.R config/<pheno>.yml")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

stage_banner("string-02", "Graph construction and gene accounting")

node_set <- s$components$downstream_node_set
if (is.null(node_set)) {
  stop("config must define string.components.downstream_node_set ('lcc' or 'degree_ge1')")
}

message("Edge file:        ", paths$edges)
message("Downstream set:   ", node_set)

# ---- inputs -----------------------------------------------------------------
edges <- read_frozen_edges(paths$edges, paths$provenance, verify = TRUE)
message("  edge file verified against its provenance checksum")

if (!file.exists(paths$mapping)) {
  stop("gene mapping table not found: ", paths$mapping, "\n  Re-run stage 01.")
}
mapping <- readr::read_tsv(paths$mapping, quote = "", progress = FALSE,
                           col_types = readr::cols(.default = readr::col_character(),
                                                   match_rank = readr::col_integer(),
                                                   n_candidate_ids = readr::col_integer()))

block_map <- NULL
if (!is.null(s$input$block_table)) {
  bt <- s$input$block_table
  if (!grepl("^(/|\\./|\\.\\./|[A-Za-z]:)", bt)) bt <- file.path(cfg$paths$intermediate, bt)
  block_map <- load_block_map(bt,
                              gene_col  = s$input$block_gene_col  %||% "external_gene_name",
                              block_col = s$input$block_block_col %||% "Block",
                              policy    = s$input$multi_block_policy %||% "error")
  message("  block table:    ", bt, "  (", length(unique(block_map)), " blocks)")
}

# ---- graph ------------------------------------------------------------------
g <- build_string_graph(edges, mapping, block_map)
message("\nGraph: ", igraph::vcount(g), " vertices, ", igraph::ecount(g), " edges")

acc  <- component_accounting(g, mapping, block_map)
summ <- summarise_accounting(acc)

# ---- accounting report ------------------------------------------------------
message("\n", paste(rep("=", 70), collapse = ""))
message("GENE ACCOUNTING")
message(paste(rep("=", 70), collapse = ""))
for (i in seq_len(nrow(summ))) {
  message(sprintf("  %-12s %4d   %s", summ$quantity[i], summ$n[i], summ$meaning[i]))
}

named <- function(label, syms) {
  if (!length(syms)) return(invisible(NULL))
  message("\n", label, " (", length(syms), "):")
  for (x in syms) message("    - ", x)
}
named("UNMAPPED — no STRING protein",        acc$input_symbol[acc$status == "unmapped"])
named("AMBIGUOUS — >1 protein",              acc$input_symbol[acc$status == "ambiguous"])
named("DUPLICATE — same protein as another", acc$input_symbol[acc$status == "duplicate_protein"])
named("ZERO-DEGREE — mapped, no association at this threshold",
      acc$input_symbol[acc$status == "singleton"])

isl <- acc[acc$status == "island", ]
if (nrow(isl)) {
  message("\nNON-LCC COMPONENTS (", nrow(isl), " genes):")
  for (cid in sort(unique(isl$component_id))) {
    m <- isl[isl$component_id == cid, ]
    message("    component ", cid, " (size ", m$component_size[1], "): ",
            paste(m$input_symbol, collapse = ", "))
  }
}

# ---- LD block composition ---------------------------------------------------
if (!is.null(block_map)) {
  inplay <- acc[acc$status %in% c("lcc", "island"), ]
  message("\nLD BLOCK COMPOSITION of the ", nrow(inplay), " connected genes:")
  tb <- sort(table(inplay$block), decreasing = TRUE)
  for (nm in names(tb)) {
    message(sprintf("    %-6s %3d  (%4.1f%%)", nm, tb[[nm]], 100 * tb[[nm]] / nrow(inplay)))
  }
}

# ---- downstream selection ---------------------------------------------------
gd <- select_downstream(g, node_set)
message("\nDownstream graph ('", node_set, "'): ", igraph::vcount(gd),
        " vertices, ", igraph::ecount(gd), " edges")

# ---- write ------------------------------------------------------------------
dir.create(paths$stage_dir, recursive = TRUE, showWarnings = FALSE)
p_full <- file.path(paths$stage_dir, "string_02_graph_full.rds")
p_down <- file.path(paths$stage_dir, "string_02_graph_downstream.rds")
p_acc  <- file.path(paths$stage_dir, "string_02_accounting.tsv")
p_sum  <- file.path(paths$stage_dir, "string_02_accounting_summary.tsv")

saveRDS(g,  p_full)
saveRDS(gd, p_down)
readr::write_tsv(acc,  p_acc,  na = "")
readr::write_tsv(summ, p_sum,  na = "")

message("\n  ", p_full)
message("  ", p_down)
message("  ", p_acc)
message("  ", p_sum)
message("\nEvery submitted gene is accounted for by name in ", basename(p_acc), ".")
