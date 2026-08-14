#!/usr/bin/env Rscript
# ============================================================================
# stage_string_04_stability.R
#
# Is the module definition an artifact of inflation 1.5?
#
# Stage 03 clusters the frozen network at one inflation value and applies one
# positional-artifact threshold. Both were carried over from an earlier manual
# run and neither has been validated. Every downstream claim — above all M3,
# the single cross-block module — rests on that one partition. This stage asks
# whether it survives being poked.
#
# Three questions, one pass:
#
#   1. INFLATION. Re-cluster the identical graph across a grid of inflation
#      values. How many modules, how many survive the artifact rule, and how
#      does the partition score, at each?
#
#   2. PERSISTENCE. For each module stage 03 found, do its members stay
#      together as inflation moves? Reported as the fraction of within-module
#      gene PAIRS still co-clustered — a module that sheds one gene scores
#      differently from one that shatters. "Intact" is also reported with the
#      size of the host cluster, because at low inflation everything merges and
#      a module can be intact only by being swallowed.
#
#   3. THE 60% RULE. Which modules change verdict at 55% / 60% / 65%? M4 clears
#      the rule at 57.1% — one gene either way. A threshold doing arbitrary
#      work should be visible, not discovered by a reviewer.
#
# This stage READS the frozen network and re-analyses it. It defines no new
# dataset, writes nothing stage 03 owns, and is safe to re-run.
#
# Two integrity checks guard against sweeping something other than what stage
# 03 clustered:
#   - the mcl input built here must be byte-identical to string_03_mcl_input.abc
#   - the sweep's partition AT the reference inflation must reproduce
#     string_03_modules.tsv exactly (mcl carries no RNG, so it must)
# Either failing means stage 03's outputs are stale. The stage stops.
#
# Inputs:
#   config/<phenotype>.yml            (string.stability block required)
#   <stage_dir>/string_02_graph_downstream.rds
#   <stage_dir>/string_03_modules.tsv
#   <stage_dir>/string_03_module_summary.tsv
#
# Outputs (to <stage_dir>/):
#   string_04_stability_sweep.tsv         one row per inflation value
#   string_04_stability_persistence.tsv   one row per (stage-03 module x inflation)
#   string_04_stability_modules.tsv       one row per stage-03 module, summarised
#   string_04_stability_artifact_rule.tsv one row per (module x artifact fraction)
#   string_04_stability_gene_pairs.tsv    consensus co-clustering, one row per pair
#   string_04_stability_provenance.json   grid, mcl version, checksums
#   string_04_sweep/mcl_I<value>.txt      exact mcl output at each inflation
#
# Usage:
#   Rscript pipeline/string_pipeline/stage_string_04_stability.R config/leptin.yml
# ============================================================================

if (!file.exists("R/utils.R")) stop("run this script from the repository root")
source("R/utils.R")
source("R/string_utils.R")
suppressPackageStartupMessages(library(igraph))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript pipeline/string_pipeline/stage_string_04_stability.R config/<pheno>.yml")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)

stage_banner("string-04", "Module stability: inflation sweep and artifact-rule sensitivity")

# ---- config ----------------------------------------------------------------
# Same rule as stage 03: nothing that changes a reported number gets a silent
# default. The grid defines the coverage of the robustness claim, so it is
# declared in the config and travels with the phenotype.
if (is.null(s$mcl$inflation)) stop("config must define string.mcl.inflation")
if (is.null(s$modules$max_single_block_fraction)) {
  stop("config must define string.modules.max_single_block_fraction")
}
if (is.null(s$stability$inflation_grid)) {
  stop("config must define string.stability.inflation_grid\n",
       "  e.g.  inflation_grid: [1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0]")
}
if (is.null(s$stability$artifact_fractions)) {
  stop("config must define string.stability.artifact_fractions\n",
       "  e.g.  artifact_fractions: [0.55, 0.60, 0.65]")
}
if (is.null(s$stability$core_consensus_levels)) {
  stop("config must define string.stability.core_consensus_levels\n",
       "  e.g.  core_consensus_levels: [0.5, 0.6, 0.7, 0.8]")
}

infl_ref <- as.numeric(s$mcl$inflation)
maxfrac  <- as.numeric(s$modules$max_single_block_fraction)
binary   <- s$mcl$binary %||% "mcl"
wscale   <- as.numeric(s$mcl$weight_scale %||% 1000)

# The reference inflation is always swept, whether or not it was listed: the
# persistence and artifact-rule sections are both anchored to it.
grid  <- sort(unique(round(c(as.numeric(unlist(s$stability$inflation_grid)), infl_ref), 3)))
fracs <- sort(unique(round(as.numeric(unlist(s$stability$artifact_fractions)), 4)))
if (!maxfrac %in% fracs) fracs <- sort(c(fracs, maxfrac))
levels_core <- sort(unique(round(as.numeric(unlist(s$stability$core_consensus_levels)), 4)))
sweep_subdir <- s$stability$sweep_subdir %||% "string_04_sweep"

mcl_path <- Sys.which(binary)
if (!nzchar(mcl_path)) {
  stop("mcl binary not found on PATH: '", binary, "'\n",
       "  Install van Dongen's reference implementation:  sudo apt-get install -y mcl")
}
mcl_ver <- tryCatch(system2(binary, "--version", stdout = TRUE, stderr = TRUE)[1],
                    error = function(e) NA_character_)

# ---- inputs ----------------------------------------------------------------
g_path   <- file.path(paths$stage_dir, "string_02_graph_downstream.rds")
mod_path <- file.path(paths$stage_dir, "string_03_modules.tsv")
sum_path <- file.path(paths$stage_dir, "string_03_module_summary.tsv")
abc_03   <- file.path(paths$stage_dir, "string_03_mcl_input.abc")
for (p in c(g_path, mod_path, sum_path)) {
  if (!file.exists(p)) {
    stop("input not found: ", p, "\n  Run stage_string_03_modules.R first.")
  }
}
g     <- readRDS(g_path)
m03   <- read.delim(mod_path, stringsAsFactors = FALSE)
msum03 <- read.delim(sum_path, stringsAsFactors = FALSE)

message("mcl:              ", mcl_path, "  (", mcl_ver, ")")
message("graph:            ", igraph::vcount(g), " genes, ", igraph::ecount(g), " associations")
message("reference I:      ", infl_ref, "   (stage 03's partition)")
message("inflation grid:   ", paste(format(grid, nsmall = 1), collapse = ", "))
message("artifact rule:    ", paste0(round(100 * fracs), "%", collapse = " / "),
        "   (config: >", round(100 * maxfrac), "%)")

sweep_dir <- file.path(paths$stage_dir, sweep_subdir)
dir.create(sweep_dir, recursive = TRUE, showWarnings = FALSE)

# ---- the mcl input, rebuilt exactly as stage 03 builds it ------------------
el <- igraph::as_data_frame(g, what = "edges")
el$from_g <- igraph::V(g)$gene[match(el$from, igraph::V(g)$name)]
el$to_g   <- igraph::V(g)$gene[match(el$to,   igraph::V(g)$name)]
el$weight <- el$combined_score / wscale

abc_path <- file.path(sweep_dir, "string_04_mcl_input.abc")
write.table(el[, c("from_g", "to_g", "weight")], abc_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

if (file.exists(abc_03) && !identical(file_md5(abc_path), file_md5(abc_03))) {
  stop("the mcl input rebuilt here differs from stage 03's:\n",
       "    ", abc_03, "  (", file_md5(abc_03), ")\n",
       "    ", abc_path, "  (", file_md5(abc_path), ")\n",
       "  The graph or the weight scale has changed since stage 03 ran, so a\n",
       "  sweep would not be describing stage 03's partition. Re-run stage 03.")
}
message("  mcl input verified identical to stage 03's")

# ---- helpers ---------------------------------------------------------------
parse_mcl <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  cl <- lapply(strsplit(ln, "\t"), function(z) trimws(z[nzchar(trimws(z))]))
  cl[order(-lengths(cl))]                 # largest first — stable module ids
}

run_mcl_at <- function(infl) {
  out <- file.path(sweep_dir, sprintf("mcl_I%.2f.txt", infl))
  a <- c(shQuote(abc_path), "--abc", "-I", format(infl, trim = TRUE),
         "-o", shQuote(out), unlist(s$mcl$extra_args))
  st <- system2(binary, a, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) {
    stop("mcl failed at I=", infl, ":\n", paste(st, collapse = "\n"))
  }
  if (!file.exists(out)) stop("mcl produced no output at I=", infl)
  list(file = out, clusters = parse_mcl(out))
}

membership_of <- function(clusters) {
  stats::setNames(rep(seq_along(clusters), lengths(clusters)), unlist(clusters))
}

# Block composition. Mirrors summarise_module() in stage 03 — including the
# denominator, which is the module's full size, not the number of genes that
# carry a block. Asserted against stage 03's own output below.
block_profile <- function(genes, blk) {
  b  <- blk[genes]
  tb <- sort(table(b[!is.na(b)]), decreasing = TRUE)
  list(
    n_blocks       = length(tb),
    dominant_block = if (length(tb)) names(tb)[1] else NA_character_,
    dominant_n     = if (length(tb)) as.integer(tb[1]) else NA_integer_,
    dominant_frac  = if (length(tb)) as.numeric(tb[1]) / length(genes) else NA_real_,
    blocks         = paste(sprintf("%s(%d)", names(tb), as.integer(tb)), collapse = " ")
  )
}

# Fraction of the module's gene PAIRS still in one cluster. A module that loses
# a single gene and one that fragments completely both stop being "intact";
# only the pair fraction separates them.
pair_cocluster <- function(genes, memb) {
  if (length(genes) < 2L) return(NA_real_)
  m <- memb[genes]
  n <- length(m)
  (sum(outer(m, m, "==")) - n) / (n * (n - 1))
}

# A partition as a comparable object: member sets, each sorted, themselves
# sorted. unname() because split() names its groups and mcl's output does not,
# and identical() would compare the names rather than the partition.
partition_key <- function(clusters) {
  unname(sort(vapply(clusters, function(z) paste(sort(z), collapse = "|"),
                     character(1))))
}

blk <- stats::setNames(igraph::V(g)$block, igraph::V(g)$gene)
w   <- igraph::E(g)$combined_score / wscale

# ---- sweep -----------------------------------------------------------------
message("\nrunning mcl at ", length(grid), " inflation values ...")
runs <- stats::setNames(lapply(grid, run_mcl_at), sprintf("%.2f", grid))
membs <- lapply(runs, function(r) membership_of(r$clusters))

for (i in seq_along(grid)) {
  placed <- names(membs[[i]])
  if (!setequal(placed, igraph::V(g)$gene)) {
    stop("mcl at I=", grid[i], " did not place every gene (", length(placed),
         " of ", igraph::vcount(g), ")")
  }
}

# ---- integrity: the reference run must reproduce stage 03 ------------------
ref_key <- sprintf("%.2f", infl_ref)
ref_clusters <- runs[[ref_key]]$clusters
m03_sets <- split(m03$gene, m03$module_id)
if (!identical(partition_key(ref_clusters), partition_key(m03_sets))) {
  stop("the sweep at I=", infl_ref, " does not reproduce string_03_modules.tsv.\n",
       "  mcl is deterministic, so this means stage 03's outputs are stale\n",
       "  relative to the current graph or config. Re-run stage_string_03_modules.R.")
}
message("  reference partition reproduces stage 03 exactly (",
        length(ref_clusters), " modules)")

# Stage 03's module ids, keyed by gene set, so this stage's labels are its labels.
m03_id_by_key <- stats::setNames(
  names(m03_sets),
  vapply(m03_sets, function(z) paste(sort(z), collapse = "|"), character(1))
)
ref_ids <- as.integer(m03_id_by_key[
  vapply(ref_clusters, function(z) paste(sort(z), collapse = "|"), character(1))
])

# ---- per-inflation summary -------------------------------------------------
sweep_rows <- lapply(seq_along(grid), function(i) {
  cl   <- runs[[i]]$clusters
  memb <- membs[[i]]
  prof <- lapply(cl, block_profile, blk = blk)
  dfr  <- vapply(prof, function(p) p$dominant_frac, numeric(1))
  sizes <- lengths(cl)

  row <- tibble::tibble(
    inflation      = grid[i],
    n_modules      = length(cl),
    n_singletons   = sum(sizes == 1L),
    largest_module = max(sizes),
    median_size    = stats::median(sizes),
    n_genes        = sum(sizes),
    modularity     = igraph::modularity(g, memb[igraph::V(g)$gene], weights = w),
    mcl_output     = basename(runs[[i]]$file)
  )
  # Survivors at each candidate artifact threshold, not just the configured one.
  for (f in fracs) {
    art <- !is.na(dfr) & dfr > f
    row[[sprintf("n_artifact_f%02d", round(100 * f))]]     <- sum(art)
    row[[sprintf("n_survive_f%02d", round(100 * f))]]      <- sum(!art)
    row[[sprintf("genes_survive_f%02d", round(100 * f))]]  <- sum(sizes[!art])
  }
  row
})
sweep <- do.call(rbind, sweep_rows)

# ---- reference modules: block profile, and the artifact rule across fracs --
ref_prof <- lapply(ref_clusters, block_profile, blk = blk)
ref_tbl <- tibble::tibble(
  module_id      = ref_ids,
  size           = lengths(ref_clusters),
  n_blocks       = vapply(ref_prof, function(p) p$n_blocks, integer(1)),
  dominant_block = vapply(ref_prof, function(p) p$dominant_block, character(1)),
  dominant_n     = vapply(ref_prof, function(p) p$dominant_n, integer(1)),
  dominant_frac  = vapply(ref_prof, function(p) p$dominant_frac, numeric(1)),
  blocks         = vapply(ref_prof, function(p) p$blocks, character(1)),
  genes          = vapply(ref_clusters, function(z) paste(z, collapse = ", "), character(1))
)
ref_tbl <- ref_tbl[order(ref_tbl$module_id), ]

# Assert this stage computes block composition the way stage 03 does. If the
# two ever drift, every verdict below would be quietly measured on a different
# rule than the one stage 03 published.
chk <- merge(ref_tbl[, c("module_id", "dominant_frac")],
             msum03[, c("module_id", "dominant_frac")],
             by = "module_id", suffixes = c("_04", "_03"))
if (nrow(chk) != nrow(ref_tbl) ||
    !isTRUE(all.equal(chk$dominant_frac_04, chk$dominant_frac_03, tolerance = 1e-9))) {
  stop("block composition computed here disagrees with string_03_module_summary.tsv")
}

art_rows <- do.call(rbind, lapply(fracs, function(f) {
  tibble::tibble(
    module_id     = ref_tbl$module_id,
    size          = ref_tbl$size,
    dominant_block = ref_tbl$dominant_block,
    dominant_frac = ref_tbl$dominant_frac,
    artifact_fraction = f,
    is_configured = isTRUE(all.equal(f, maxfrac)),
    verdict       = ifelse(!is.na(ref_tbl$dominant_frac) & ref_tbl$dominant_frac > f,
                           "artifact", "survives")
  )
}))
art_rows <- art_rows[order(art_rows$module_id, art_rows$artifact_fraction), ]
flips <- vapply(split(art_rows$verdict, art_rows$module_id),
                function(v) length(unique(v)) > 1L, logical(1))
art_rows$flips <- unname(flips[as.character(art_rows$module_id)])

# ---- persistence of each reference module across the grid ------------------
pers_rows <- list()
for (k in seq_along(ref_clusters)) {
  genes <- ref_clusters[[k]]
  for (i in seq_along(grid)) {
    memb <- membs[[i]]
    tb   <- sort(table(memb[genes]), decreasing = TRUE)
    host <- as.integer(names(tb)[1])
    host_size <- sum(memb == host)
    pers_rows[[length(pers_rows) + 1L]] <- tibble::tibble(
      module_id    = ref_ids[k],
      size         = length(genes),
      inflation    = grid[i],
      n_parts      = length(tb),
      largest_part = as.integer(tb[1]),
      intact       = length(tb) == 1L,
      pair_cocluster = pair_cocluster(genes, memb),
      host_size    = host_size,           # cluster holding the largest part
      host_purity  = as.numeric(tb[1]) / host_size,
      absorbed     = host_size > length(genes)
    )
  }
}
pers <- do.call(rbind, pers_rows)

off <- pers[abs(pers$inflation - infl_ref) > 1e-9, ]
agg <- function(f) vapply(split(off, off$module_id), f, numeric(1))
ref_tbl$cohesion    <- unname(agg(function(d) mean(d$pair_cocluster, na.rm = TRUE))[
                                as.character(ref_tbl$module_id)])
ref_tbl$intact_frac <- unname(agg(function(d) mean(d$intact))[
                                as.character(ref_tbl$module_id)])
ref_tbl$absorbed_frac <- unname(agg(function(d) mean(d$absorbed))[
                                as.character(ref_tbl$module_id)])
ref_tbl$artifact_configured <- !is.na(ref_tbl$dominant_frac) & ref_tbl$dominant_frac > maxfrac
ref_tbl$rule_flips <- unname(flips[as.character(ref_tbl$module_id)])

# ---- consensus co-clustering, every gene pair ------------------------------
# An inflation that returns a single cluster covering the whole graph asserts
# nothing about which genes belong together — every pair scores 1 trivially. It
# is excluded from the consensus rather than allowed to lift every pair by
# 1/|grid|. Excluded by a rule that fires and is named, never dropped silently.
informative <- vapply(runs, function(r) length(r$clusters) > 1L, logical(1))
degenerate  <- grid[!informative]
if (all(!informative)) stop("every inflation in the grid returned a single cluster")
if (length(degenerate)) {
  message("\n  note: I = ", paste(sprintf("%.2f", degenerate), collapse = ", "),
          " returned a single cluster covering all ", igraph::vcount(g),
          " genes.\n  Excluded from the consensus and cores below (it separates nothing).",
          "\n  It is still reported in the sweep and persistence tables.")
}
membs_co <- membs[informative]

genes_all <- sort(igraph::V(g)$gene)
pair_idx  <- utils::combn(length(genes_all), 2L)
co <- rowMeans(vapply(membs_co, function(m) {
  as.numeric(m[genes_all][pair_idx[1, ]] == m[genes_all][pair_idx[2, ]])
}, numeric(ncol(pair_idx))))
ref_mem <- membs[[ref_key]]
pairs <- tibble::tibble(
  gene1 = genes_all[pair_idx[1, ]],
  gene2 = genes_all[pair_idx[2, ]],
  block1 = unname(blk[genes_all[pair_idx[1, ]]]),
  block2 = unname(blk[genes_all[pair_idx[2, ]]]),
  same_block = !is.na(blk[genes_all[pair_idx[1, ]]]) &
               blk[genes_all[pair_idx[1, ]]] == blk[genes_all[pair_idx[2, ]]],
  ref_together = ref_mem[genes_all[pair_idx[1, ]]] == ref_mem[genes_all[pair_idx[2, ]]],
  consensus = co
)
pairs <- pairs[pairs$consensus > 0 | pairs$ref_together, ]
pairs <- pairs[order(-pairs$consensus, pairs$gene1, pairs$gene2), ]

# ---- stable cores ----------------------------------------------------------
# A module's pair-cohesion says how much of it holds together, not WHICH part.
# For M3 the distinction is the whole result: a module can look cross-block at
# one inflation because it glues several within-block sub-units together.
#
# So: within each stage-03 module, keep the gene pairs whose consensus clears a
# level, and take the connected components. Each component is a set of genes
# that travels together across the grid. A core spanning >1 LD block is the
# only thing in this analysis that positional co-location cannot explain.
#
# The level is not a hypothesis test and is not picked — every level in the
# config is reported, so its influence is visible rather than chosen.
pair_key <- paste(pairs$gene1, pairs$gene2, sep = "\r")
cons_of <- function(a, b) {
  k <- ifelse(a < b, paste(a, b, sep = "\r"), paste(b, a, sep = "\r"))
  out <- pairs$consensus[match(k, pair_key)]
  ifelse(is.na(out), 0, out)
}

core_rows <- list()
for (k in seq_along(ref_clusters)) {
  gs <- sort(ref_clusters[[k]])
  for (lv in levels_core) {
    if (length(gs) == 1L) {
      comps <- list(gs)
    } else {
      idx <- utils::combn(length(gs), 2L)
      keep <- cons_of(gs[idx[1, ]], gs[idx[2, ]]) >= lv
      sub <- igraph::make_empty_graph(n = length(gs), directed = FALSE)
      igraph::V(sub)$name <- gs
      if (any(keep)) sub <- igraph::add_edges(sub, as.vector(idx[, keep, drop = FALSE]))
      cm <- igraph::components(sub)
      comps <- split(gs, cm$membership)
      comps <- comps[order(-lengths(comps))]
    }
    for (ci in seq_along(comps)) {
      cg <- comps[[ci]]
      bp <- block_profile(cg, blk)
      minc <- if (length(cg) > 1L) {
        ii <- utils::combn(length(cg), 2L)
        min(cons_of(cg[ii[1, ]], cg[ii[2, ]]))
      } else NA_real_
      core_rows[[length(core_rows) + 1L]] <- tibble::tibble(
        module_id      = ref_ids[k],
        consensus_level = lv,
        core_id        = ci,
        size           = length(cg),
        n_blocks       = bp$n_blocks,
        dominant_block = bp$dominant_block,
        dominant_frac  = bp$dominant_frac,
        blocks         = bp$blocks,
        cross_block    = bp$n_blocks > 1L,
        min_consensus  = minc,
        genes          = paste(cg, collapse = ", ")
      )
    }
  }
}
cores <- do.call(rbind, core_rows)

# ---- report ----------------------------------------------------------------
bar <- function(n) paste(rep("=", n), collapse = "")
fcols <- sprintf("n_survive_f%02d", round(100 * fracs))

message("\n", bar(74))
message("INFLATION SWEEP")
message(bar(74))
hdr <- sprintf("%6s %9s %9s %9s %11s   %s", "I", "modules", "largest", "singl.",
               "modularity", paste(sprintf("surv@%d%%", round(100 * fracs)), collapse = " "))
message(hdr)
for (i in seq_len(nrow(sweep))) {
  message(sprintf("%6.2f %9d %9d %9d %11.3f   %s",
                  sweep$inflation[i], sweep$n_modules[i], sweep$largest_module[i],
                  sweep$n_singletons[i], sweep$modularity[i],
                  paste(sprintf("%7d", unlist(sweep[i, fcols])), collapse = " ")))
}
message("\n  modularity is descriptive, not a selection criterion — it has a known")
message("  resolution limit and will favour whatever granularity the graph's size")
message("  suits. Read it alongside the module counts, not instead of them.")

message("\n", bar(74))
message("PERSISTENCE OF STAGE-03 MODULES   (reference I = ", infl_ref, ")")
message(bar(74))
message("fraction of each module's gene PAIRS still co-clustered\n")
message(sprintf("%-6s %4s  %s  %9s %7s", "module", "n",
                paste(sprintf("%5.2f", grid), collapse = " "), "cohesion", "intact"))
for (k in seq_len(nrow(ref_tbl))) {
  id <- ref_tbl$module_id[k]
  d  <- pers[pers$module_id == id, ]
  d  <- d[match(grid, d$inflation), ]
  cells <- vapply(d$pair_cocluster, function(x)
    if (is.na(x)) "    ." else sprintf("%5.2f", x), character(1))
  tag <- if (ref_tbl$artifact_configured[k]) " (artifact)" else ""
  message(sprintf("%-6s %4d  %s  %9.2f %7.0f%%%s",
                  paste0("M", id), ref_tbl$size[k], paste(cells, collapse = " "),
                  ref_tbl$cohesion[k], 100 * ref_tbl$intact_frac[k], tag))
}
message("\n  '.' = single-gene module, no pairs to score.")
message("  cohesion and intact% are means over the grid EXCLUDING the reference,")
message("  where both are 1 by construction.")

absorbed_any <- ref_tbl$module_id[!is.na(ref_tbl$absorbed_frac) & ref_tbl$absorbed_frac > 0]
if (length(absorbed_any)) {
  message("\n  absorbed at some inflation (intact only by being swallowed by a")
  message("  larger cluster — check host_size in the persistence table):")
  for (id in absorbed_any) {
    d <- pers[pers$module_id == id & pers$absorbed, ]
    message(sprintf("    M%-3s at I = %s", id,
                    paste(sprintf("%.2f (host %d)", d$inflation, d$host_size),
                          collapse = ", ")))
  }
}

message("\n", bar(74))
message("ARTIFACT-RULE SENSITIVITY   (reference I = ", infl_ref, ")")
message(bar(74))
message(sprintf("%-6s %4s %8s %8s   %s", "module", "n", "block", "frac",
                paste(sprintf("%8s", paste0(round(100 * fracs), "%")), collapse = " ")))
for (k in seq_len(nrow(ref_tbl))) {
  id <- ref_tbl$module_id[k]
  d  <- art_rows[art_rows$module_id == id, ]
  d  <- d[match(fracs, d$artifact_fraction), ]
  cells <- ifelse(d$verdict == "survives", "keep", "drop")
  cells[d$is_configured] <- toupper(cells[d$is_configured])   # the configured rule
  flag <- if (isTRUE(ref_tbl$rule_flips[k])) "   <-- FLIPS" else ""
  message(sprintf("%-6s %4d %8s %7.1f%%   %s%s",
                  paste0("M", id), ref_tbl$size[k],
                  ref_tbl$dominant_block[k], 100 * ref_tbl$dominant_frac[k],
                  paste(sprintf("%8s", cells), collapse = " "), flag))
}
message("\n  UPPERCASE = the verdict at the configured rule (>",
        round(100 * maxfrac), "%).")
n_flip <- sum(ref_tbl$rule_flips, na.rm = TRUE)
if (n_flip == 0L) {
  message("  No module changes verdict across ",
          paste0(round(100 * fracs), "%", collapse = "/"),
          " — the threshold is not doing arbitrary work.")
} else {
  message("  ", n_flip, " module(s) change verdict within the tested range. Their",
          "\n  status is a property of the threshold, not of the network.")
}

message("\n", bar(74))
message("STABLE CORES   (consensus over ", sum(informative), " informative inflation values)")
message(bar(74))
message("sub-units of each module that travel together across the grid.")
message("a module can look cross-block at one inflation by gluing several")
message("within-block sub-units together — this is where that shows up.\n")
survivors <- ref_tbl$module_id[!ref_tbl$artifact_configured]
for (id in survivors) {
  rt <- ref_tbl[ref_tbl$module_id == id, ]
  message(sprintf("M%s  n=%d  %s", id, rt$size, rt$blocks))
  for (lv in levels_core) {
    d <- cores[cores$module_id == id & cores$consensus_level == lv, ]
    d <- d[order(-d$size), ]
    shed <- sum(d$size == 1L)
    for (j in seq_len(nrow(d))) {
      if (d$size[j] == 1L) next
      lab <- if (j == 1L) sprintf("  >=%.2f", lv) else "        "
      message(sprintf("%s  %2d genes  %d block%s  %-26s %s%s",
                      lab, d$size[j], d$n_blocks[j],
                      if (d$n_blocks[j] == 1L) " " else "s", d$blocks[j],
                      d$genes[j],
                      if (d$cross_block[j]) "   <-- cross-block" else ""))
    }
    if (all(d$size == 1L)) {
      message(sprintf("  >=%.2f   no pair holds together — the module is not a unit above this level", lv))
    } else if (shed > 0L) {
      message(sprintf("        (%d gene%s shed)", shed, if (shed == 1L) "" else "s"))
    }
  }
  message("")
}
xb <- cores[cores$cross_block & cores$size >= 2L, ]
if (nrow(xb)) {
  top <- xb[xb$consensus_level == max(xb$consensus_level), ]
  message("  cross-block cores holding at the strictest level tested (>=",
          sprintf("%.2f", max(xb$consensus_level)), "):")
  if (nrow(top)) {
    # Scanned across every stage-03 module, not only the survivors: a
    # positional module can still contain a cross-block pair, and that pair is
    # evidence the artifact rule discards wholesale at module level.
    for (j in seq_len(nrow(top))) {
      art <- ref_tbl$artifact_configured[match(top$module_id[j], ref_tbl$module_id)]
      message(sprintf("    M%-3s %-52s [%s]%s", top$module_id[j], top$genes[j],
                      top$blocks[j],
                      if (isTRUE(art)) "  (from a flagged module)" else ""))
    }
  } else {
    message("    none")
  }
} else {
  message("  no cross-block core survives any tested consensus level. Every")
  message("  sub-unit that holds together across the grid is within one LD block.")
}

# ---- write -----------------------------------------------------------------
p_sweep <- file.path(paths$stage_dir, "string_04_stability_sweep.tsv")
p_pers  <- file.path(paths$stage_dir, "string_04_stability_persistence.tsv")
p_mods  <- file.path(paths$stage_dir, "string_04_stability_modules.tsv")
p_art   <- file.path(paths$stage_dir, "string_04_stability_artifact_rule.tsv")
p_pair  <- file.path(paths$stage_dir, "string_04_stability_gene_pairs.tsv")
p_core  <- file.path(paths$stage_dir, "string_04_stability_cores.tsv")
readr::write_tsv(cores,    p_core,  na = "")
readr::write_tsv(sweep,    p_sweep, na = "")
readr::write_tsv(pers,     p_pers,  na = "")
readr::write_tsv(ref_tbl,  p_mods,  na = "")
readr::write_tsv(art_rows, p_art,   na = "")
readr::write_tsv(pairs,    p_pair,  na = "")

prov <- list(
  script            = "pipeline/string_pipeline/stage_string_04_stability.R",
  mcl_binary        = unname(mcl_path),
  mcl_version       = mcl_ver,
  reference_inflation = infl_ref,
  inflation_grid    = as.list(grid),
  artifact_fractions = as.list(fracs),
  core_consensus_levels = as.list(levels_core),
  consensus_inflations  = as.list(grid[informative]),
  degenerate_inflations = as.list(degenerate),
  configured_max_single_block_fraction = maxfrac,
  weight_scale      = wscale,
  extra_args        = as.list(unlist(s$mcl$extra_args)),
  node_set          = s$components$downstream_node_set,
  n_genes           = igraph::vcount(g),
  n_edges           = igraph::ecount(g),
  edge_file         = paths$edges,
  edge_md5          = file_md5(paths$edges),
  mcl_input_md5     = file_md5(abc_path),
  stage_03_mcl_input_md5 = if (file.exists(abc_03)) file_md5(abc_03) else NULL,
  mcl_output_md5    = as.list(stats::setNames(
                        vapply(runs, function(r) file_md5(r$file), character(1)),
                        sprintf("I%.2f", grid))),
  n_modules_by_inflation = as.list(stats::setNames(
                        as.integer(sweep$n_modules), sprintf("I%.2f", grid))),
  n_rule_flips      = as.integer(n_flip),
  run_on            = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version         = R.version.string,
  package_versions  = package_pins(c("igraph", "readr", "tibble", "jsonlite", "yaml"))
)
jsonlite::write_json(prov, file.path(paths$stage_dir, "string_04_stability_provenance.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")

message("\n  ", p_sweep, "\n  ", p_pers, "\n  ", p_mods,
        "\n  ", p_art, "\n  ", p_core, "\n  ", p_pair,
        "\n  ", sweep_dir, "/  (", length(grid), " raw mcl outputs)")
