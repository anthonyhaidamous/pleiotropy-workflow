#!/usr/bin/env Rscript
# ============================================================================
# stage_string_03_modules.R
#
# Markov clustering of the STRING network into modules, and the LD block
# composition of each module.
#
# MCL is used because it is deterministic (no RNG, so a run is reproducible
# without seeding), assumes no number of clusters, and partitions by density
# rather than by path. This calls van Dongen's reference `mcl` binary, not the
# CRAN MCL package.
#
# The positional-artifact rule comes from the config and is fixed in advance:
# any module with more than string.modules.max_single_block_fraction of its
# members from one LD block is flagged. Setting that threshold after seeing the
# numbers is how you talk yourself into keeping a module you like.
#
# Inputs:
#   config/<phenotype>.yml
#   <stage_dir>/string_02_graph_downstream.rds
#
# Outputs (to <stage_dir>/):
#   string_03_modules.tsv            one row per gene: module, block
#   string_03_module_summary.tsv     one row per module: size, blocks, flag
#   string_03_block_by_module.tsv    block x module contingency table
#   string_03_mcl_input.abc          exact input handed to mcl
#   string_03_mcl_raw.txt            exact mcl output
#   string_03_provenance.json        mcl version, inflation, input checksums
#
# Usage:
#   Rscript pipeline/string_pipeline/stage_string_03_modules.R config/leptin.yml
# ============================================================================

if (!file.exists("R/utils.R")) stop("run this script from the repository root")
source("R/utils.R")
source("R/string_utils.R")
suppressPackageStartupMessages(library(igraph))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript pipeline/string_pipeline/stage_string_03_modules.R config/<pheno>.yml")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

stage_banner("string-03", "MCL module detection and block composition")

# ---- config, no silent defaults for anything that changes the answer -------
if (is.null(s$mcl$inflation)) stop("config must define string.mcl.inflation")
if (is.null(s$modules$max_single_block_fraction)) {
  stop("config must define string.modules.max_single_block_fraction")
}
infl    <- as.numeric(s$mcl$inflation)
maxfrac <- as.numeric(s$modules$max_single_block_fraction)
binary  <- s$mcl$binary %||% "mcl"
wscale  <- as.numeric(s$mcl$weight_scale %||% 1000)

mcl_path <- Sys.which(binary)
if (!nzchar(mcl_path)) {
  stop("mcl binary not found on PATH: '", binary, "'\n",
       "  Install van Dongen's reference implementation:  sudo apt-get install -y mcl")
}
mcl_ver <- tryCatch(system2(binary, "--version", stdout = TRUE, stderr = TRUE)[1],
                    error = function(e) NA_character_)

g_path <- file.path(paths$stage_dir, "string_02_graph_downstream.rds")
if (!file.exists(g_path)) {
  stop("graph not found: ", g_path, "\n  Run stage_string_02_components.R first.")
}
g <- readRDS(g_path)

message("mcl:              ", mcl_path, "  (", mcl_ver, ")")
message("inflation:        ", infl)
message("graph:            ", igraph::vcount(g), " genes, ", igraph::ecount(g), " associations")
message("artifact rule:    >", round(100 * maxfrac), "% from one LD block")

# ---- write the mcl input ---------------------------------------------------
# .abc format: one edge per line, "labelA labelB weight". Gene symbols are used
# as labels so the raw mcl output is readable on its own.
el <- igraph::as_data_frame(g, what = "edges")
el$from_g <- igraph::V(g)$gene[match(el$from, igraph::V(g)$name)]
el$to_g   <- igraph::V(g)$gene[match(el$to,   igraph::V(g)$name)]
el$weight <- el$combined_score / wscale

abc_path <- file.path(paths$stage_dir, "string_03_mcl_input.abc")
raw_path <- file.path(paths$stage_dir, "string_03_mcl_raw.txt")
write.table(el[, c("from_g", "to_g", "weight")], abc_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- run mcl ---------------------------------------------------------------
mcl_args <- c(shQuote(abc_path), "--abc", "-I", infl,
              "-o", shQuote(raw_path), unlist(s$mcl$extra_args))
message("\n$ ", binary, " ", paste(mcl_args, collapse = " "))
st <- system2(binary, mcl_args, stdout = TRUE, stderr = TRUE)
if (!is.null(attr(st, "status")) && attr(st, "status") != 0) {
  stop("mcl failed:\n", paste(st, collapse = "\n"))
}
if (!file.exists(raw_path)) stop("mcl produced no output file: ", raw_path)

# ---- parse -----------------------------------------------------------------
lines <- readLines(raw_path, warn = FALSE)
lines <- lines[nzchar(trimws(lines))]
clusters <- lapply(strsplit(lines, "\t"), function(z) trimws(z[nzchar(trimws(z))]))
clusters <- clusters[order(-lengths(clusters))]      # largest first, stable ids

modules <- tibble::tibble(
  gene      = unlist(clusters),
  module_id = rep(seq_along(clusters), lengths(clusters))
)
message("\nmcl returned ", length(clusters), " modules covering ",
        nrow(modules), " genes")

placed <- igraph::V(g)$gene
if (!setequal(modules$gene, placed)) {
  miss <- setdiff(placed, modules$gene)
  stop(length(miss), " gene(s) absent from the mcl output, e.g. ",
       paste(utils::head(miss, 5), collapse = ", "))
}

# ---- attach blocks ---------------------------------------------------------
blk <- setNames(igraph::V(g)$block, igraph::V(g)$gene)
modules$block <- unname(blk[modules$gene])
modules$degree <- as.integer(igraph::degree(g))[match(modules$gene, igraph::V(g)$gene)]
modules <- modules[order(modules$module_id, modules$block, modules$gene), ]

# ---- per-module summary ----------------------------------------------------
summarise_module <- function(m) {
  tb <- sort(table(m$block[!is.na(m$block)]), decreasing = TRUE)
  dom <- if (length(tb)) names(tb)[1] else NA_character_
  frac <- if (length(tb)) as.numeric(tb[1]) / nrow(m) else NA_real_
  tibble::tibble(
    module_id      = m$module_id[1],
    size           = nrow(m),
    n_blocks       = length(tb),
    dominant_block = dom,
    dominant_n     = if (length(tb)) as.integer(tb[1]) else NA_integer_,
    dominant_frac  = frac,
    blocks         = paste(sprintf("%s(%d)", names(tb), as.integer(tb)), collapse = " "),
    genes          = paste(m$gene, collapse = ", ")
  )
}
msum <- do.call(rbind, lapply(split(modules, modules$module_id), summarise_module))
msum$positional_artifact <- !is.na(msum$dominant_frac) & msum$dominant_frac > maxfrac
msum <- msum[order(-msum$size), ]

# ---- block x module contingency -------------------------------------------
ct <- as.data.frame.matrix(table(block = modules$block, module = modules$module_id))
ct <- cbind(block = rownames(ct), ct, total = rowSums(ct))
ct <- ct[order(-ct$total), ]

# ---- report ----------------------------------------------------------------
message("\n", paste(rep("=", 74), collapse = ""))
message("MODULES  (inflation ", infl, ")")
message(paste(rep("=", 74), collapse = ""))
for (i in seq_len(nrow(msum))) {
  flag <- if (msum$positional_artifact[i]) "  <-- POSITIONAL ARTIFACT" else ""
  message(sprintf("\nM%-2d  n=%-3d  %d block(s)  %s%s",
                  msum$module_id[i], msum$size[i], msum$n_blocks[i],
                  msum$blocks[i], flag))
  for (k in strwrap(msum$genes[i], width = 68)) message("     ", k)
}

n_art <- sum(msum$positional_artifact)
message("\n", paste(rep("-", 74), collapse = ""))
message(sprintf("%d modules | %d flagged as positional artifact (>%d%% one block) | %d genes in them",
                nrow(msum), n_art, round(100 * maxfrac),
                sum(msum$size[msum$positional_artifact])))
message(sprintf("%d modules survive for downstream hypotheses, covering %d genes",
                nrow(msum) - n_art, sum(msum$size[!msum$positional_artifact])))

message("\nBLOCK x MODULE")
print(ct, row.names = FALSE)

# ---- write -----------------------------------------------------------------
p_mod <- file.path(paths$stage_dir, "string_03_modules.tsv")
p_sum <- file.path(paths$stage_dir, "string_03_module_summary.tsv")
p_ct  <- file.path(paths$stage_dir, "string_03_block_by_module.tsv")
readr::write_tsv(modules, p_mod, na = "")
readr::write_tsv(msum,    p_sum, na = "")
readr::write_tsv(ct,      p_ct,  na = "")

prov <- list(
  script          = "pipeline/string_pipeline/stage_string_03_modules.R",
  mcl_binary      = unname(mcl_path),
  mcl_version     = mcl_ver,
  inflation       = infl,
  weight_scale    = wscale,
  extra_args      = as.list(unlist(s$mcl$extra_args)),
  node_set        = s$components$downstream_node_set,
  n_genes         = igraph::vcount(g),
  n_edges         = igraph::ecount(g),
  n_modules       = nrow(msum),
  max_single_block_fraction = maxfrac,
  n_positional_artifact     = n_art,
  edge_file       = paths$edges,
  edge_md5        = file_md5(paths$edges),
  mcl_input_md5   = file_md5(abc_path),
  mcl_output_md5  = file_md5(raw_path),
  run_on          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version       = R.version.string
)
jsonlite::write_json(prov, file.path(paths$stage_dir, "string_03_provenance.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")

message("\n  ", p_mod, "\n  ", p_sum, "\n  ", p_ct)
