#!/usr/bin/env Rscript
# =============================================================================
# stage_string_01_fetch_edges.R
#
# Freezes the STRING network induced on this phenotype's protein-coding gene
# list, at a declared release and a declared score threshold.
#
# The edge file this writes is the reproducibility anchor for everything
# downstream. Once written it is not silently regenerated: an
# existing edge file causes the script to refuse unless --force is given, and
# --force moves the old file aside rather than destroying it - pretty neat, huh?
#
# What it does
#   1. Reads the gene list and maps every symbol to a STRING protein id via
#      preferred_name then aliases. Ambiguous and unmapped symbols are reported
#      by name, avoiding them being silently dropped or silently resolved.
#   2. Scans the version-pinned bulk links file in chunks, keeping rows whose
#      BOTH endpoints are in the mapped set and whose combined_score meets the
#      YAML threshold. All per-channel subscores are retained, so the
#      literature-bias check (textmining-only edges) is possible later without
#      re-deriving the network.
#   3. De-duplicates STRING's reciprocal A-B / B-A storage.
#   4. Writes the edge file, a gene-mapping audit table, and a provenance
#      sidecar carrying the md5 of the edge file.
#
# md5 is a 32-char fingerprint of the edge file. Same bytes → same hash; a
# single changed byte induced a completely different hash. It is not a random
# seed, meaning it does not control or regenerate anything: it lets anyone 
# re-hash the frozen network later and prove it is bit-for-bit the one recorded
# in the provenance side-car, so every downstream result is traceable to this
# exact same file.
# 
# Inputs:
#   config/<phenotype>.yml
#   string.bulk.dir/*.txt.gz          (run stage_string_00_download_bulk.R first)
#   string.input.gene_list            (default: stage 10's genes_list_complete.txt)
#
# Outputs (to string.output.edges_dir):
#   <slug>_edges_v<ver>_s<thr>.tsv                  the frozen network
#   <slug>_edges_v<ver>_s<thr>.provenance.json      how it was made + md5
#   <slug>_gene_mapping_v<ver>.tsv                  every submitted symbol
#
# Usage:
#   Rscript pipeline/string_pipeline/stage_string_01_fetch_edges.R config/leptin.yml
#   Rscript pipeline/string_pipeline/stage_string_01_fetch_edges.R config/leptin.yml --force
#
# Run from the repository root.
# =============================================================================

if (!file.exists("R/utils.R")) {
  stop("run this script from the repository root (R/utils.R not found)")
}
source("R/utils.R")
source("R/string_utils.R")

PINNED <- c("yaml", "readr", "dplyr", "tibble", "jsonlite")

# ---- arguments --------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
args  <- args[args != "--force"]
if (length(args) != 1L) {
  stop("usage: Rscript pipeline/string_pipeline/stage_string_01_fetch_edges.R config/<pheno>.yml [--force]")
}

cfg   <- load_config(args[1])
s     <- string_cfg(cfg)
paths <- string_paths(cfg, s)
spec  <- string_bulk_spec(s)

stage_banner("string-01", "Freeze the induced STRING edge list")

message("Phenotype:        ", cfg$phenotype$name, "  (slug: ", cfg$phenotype$slug, ")")
message("STRING version:   ", s$version)
message("Species taxon:    ", s$species_taxon)
message("Network type:     ", s$network_type)
message("Score threshold:  ", s$score_threshold, "   (combined_score >= this)")
message("Edge file:        ", paths$edges)

# ---- refuse early -----------------------------------------------------------
# Checked before any expensive work: there is no point scanning 13M rows to
# then discover we are not allowed to write the result.
if (file.exists(paths$edges) && !force) {
  stop("the frozen edge file already exists:\n    ", paths$edges, "\n",
       "It is the reproducibility anchor for every downstream stage and will\n",
       "not be regenerated silently. Re-run with --force only if you intend to\n",
       "redefine the network (the existing file is moved aside, not deleted).")
}

# ---- bulk files present? ----------------------------------------------------
absent <- spec$path[!file.exists(spec$path)]
if (length(absent)) {
  stop("missing STRING bulk file(s):\n    ", paste(absent, collapse = "\n    "),
       "\n  Run first:  Rscript pipeline/string_pipeline/stage_string_00_download_bulk.R ",
       args[1])
}
if (!file.exists(paths$manifest)) {
  stop("missing bulk manifest: ", paths$manifest,
       "\n  Run stage_string_00_download_bulk.R — provenance cites its checksums.")
}
manifest <- jsonlite::fromJSON(paths$manifest, simplifyVector = FALSE)

# The manifest describes a specific release/species/network. If the config has
# moved on since the download, the bulk files on disk are the wrong ones.
if (!identical(as.character(manifest$string_version), s$version) ||
    !identical(as.integer(manifest$species_taxon), s$species_taxon) ||
    !identical(as.character(manifest$network_type), s$network_type)) {
  stop("bulk_manifest.json describes v", manifest$string_version, " / taxon ",
       manifest$species_taxon, " / ", manifest$network_type,
       ", but the config asks for v", s$version, " / taxon ", s$species_taxon,
       " / ", s$network_type, ".\n  Re-run stage_string_00_download_bulk.R.")
}

# ---- 1. gene list -----------------------------------------------------------
gene_list_path <- resolve_gene_list_path(cfg, s)
message("\n[1/4] Gene list: ", gene_list_path)
genes <- read_gene_list(gene_list_path)
message("  ", length(genes), " unique symbols submitted")

# ---- 2. symbol -> STRING protein id ----------------------------------------
message("\n[2/4] Mapping symbols to STRING protein ids")
info    <- load_string_info(spec$path[spec$key == "info"])
aliases <- load_string_aliases(spec$path[spec$key == "aliases"])
message("  reference: ", nrow(info), " proteins, ", nrow(aliases), " aliases")

mp     <- map_symbols_to_string(genes, info, aliases)
mapped <- mp$mapped

message("  mapped:    ", nrow(mapped), " / ", length(genes))

if (length(mp$unmapped)) {
  message("  UNMAPPED (", length(mp$unmapped), "):")
  for (g in mp$unmapped) message("    - ", g)
}
if (nrow(mp$ambiguous)) {
  amb <- unique(mp$ambiguous$input_symbol)
  message("  AMBIGUOUS (", length(amb), ") — resolve to >1 protein, excluded:")
  for (g in amb) {
    ids <- mp$ambiguous$string_id[mp$ambiguous$input_symbol == g]
    message("    - ", g, " -> ", paste(ids, collapse = ", "))
  }
}
if (nrow(mp$duplicates)) {
  message("  DUPLICATE PROTEIN (", nrow(mp$duplicates),
          ") — collide with an already-mapped symbol, excluded:")
  for (i in seq_len(nrow(mp$duplicates))) {
    message("    - ", mp$duplicates$input_symbol[i], " = ",
            mp$duplicates$duplicate_of[i], "  (", mp$duplicates$string_id[i], ")")
  }
}
n_fuzzy <- sum(mapped$match_rank >= 3L)
if (n_fuzzy > 0L) {
  message("  note: ", n_fuzzy, " symbol(s) matched only case-insensitively — ",
          "see match_rank in the mapping table")
}

if (nrow(mapped) == 0L) {
  stop("no symbols mapped to STRING protein ids — check the species taxon ",
       "and the gene list format")
}

dir.create(paths$edges_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(mp$table, paths$mapping, na = "")
message("  mapping table -> ", paths$mapping)

# ---- 3. induced edges -------------------------------------------------------
links_path <- spec$path[spec$key == "links"]
message("\n[3/4] Scanning ", basename(links_path))
edges_raw <- read_links_induced(links_path,
                                keep_ids   = mapped$string_id,
                                threshold  = s$score_threshold,
                                chunk_size = as.integer(s$read$chunk_size))
rows_scanned <- attr(edges_raw, "rows_scanned")
message("  ", format(rows_scanned, big.mark = ","), " rows scanned, ",
        nrow(edges_raw), " retained at score >= ", s$score_threshold)

edges <- dedup_reciprocal(edges_raw)
# Captured now: dplyr verbs below drop non-standard attributes.
rows_before_dedup <- attr(edges, "rows_before_dedup")
message("  ", nrow(edges), " unique associations after reciprocal de-duplication")

# Attach the submitted symbols, so the frozen file is readable without a join.
sym_of <- setNames(mapped$input_symbol, mapped$string_id)
edges  <- tibble::add_column(edges,
                             gene1 = unname(sym_of[edges$protein1]),
                             gene2 = unname(sym_of[edges$protein2]),
                             .after = "protein2")
edges <- dplyr::arrange(edges, dplyr::desc(combined_score), gene1, gene2)

# ---- 4. freeze --------------------------------------------------------------
message("\n[4/4] Writing the frozen network")
write_frozen_edges(edges, paths$edges, force = force)
edge_md5 <- file_md5(paths$edges)
message("  edges -> ", paths$edges)
message("  md5     ", edge_md5)

bulk_records <- lapply(manifest$files, function(f) {
  list(key = f$key, filename = f$filename, url = f$url,
       bytes = f$bytes, md5 = f$md5, retrieved_on = f$retrieved_on)
})

provenance <- list(
  edge_file            = paths$edges,
  edge_md5             = edge_md5,
  n_edges              = nrow(edges),

  string_version       = s$version,
  network_type         = s$network_type,
  species_taxon        = s$species_taxon,
  score_threshold      = s$score_threshold,
  source               = "bulk download",
  url_bases_tried      = as.list(string_url_bases(s)),
  bulk_files           = bulk_records,
  live_version_at_download = manifest$live_version_at_download,

  phenotype            = cfg$phenotype$name,
  phenotype_slug       = cfg$phenotype$slug,
  gene_list_file       = gene_list_path,
  gene_list_md5        = file_md5(gene_list_path),
  n_genes_submitted    = length(genes),
  n_genes_mapped       = nrow(mapped),
  n_genes_unmapped     = length(mp$unmapped),
  unmapped_symbols     = as.list(mp$unmapped),
  n_genes_ambiguous    = length(unique(mp$ambiguous$input_symbol)),
  ambiguous_symbols    = as.list(unique(mp$ambiguous$input_symbol)),
  n_genes_duplicate    = nrow(mp$duplicates),
  duplicate_symbols    = as.list(mp$duplicates$input_symbol),
  mapping_file         = paths$mapping,

  rows_scanned         = rows_scanned,
  rows_before_dedup    = rows_before_dedup,

  config_file          = args[1],
  retrieved_on         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version            = R.version.string,
  package_versions     = package_pins(PINNED),
  script               = "pipeline/string_pipeline/stage_string_01_fetch_edges.R"
)
jsonlite::write_json(provenance, paths$provenance, auto_unbox = TRUE,
                     pretty = TRUE, null = "null", na = "null")
message("  provenance -> ", paths$provenance)

# ---- summary ----------------------------------------------------------------
# Node-level accounting (degree, components, who leaves the analysis) is
# stage 02's job. This is only what stage 01 can state on its own.
in_edges <- length(unique(c(edges$protein1, edges$protein2)))
message("\n", paste(rep("-", 70), collapse = ""))
message("submitted symbols       ", length(genes))
message("mapped to STRING        ", nrow(mapped))
message("unmapped                ", length(mp$unmapped))
message("ambiguous (excluded)    ", length(unique(mp$ambiguous$input_symbol)))
message("duplicate (excluded)    ", nrow(mp$duplicates))
message("associations frozen     ", nrow(edges))
message("proteins with >=1 edge  ", in_edges, "  (full accounting is stage 02)")
message(paste(rep("-", 70), collapse = ""))
message("Frozen. Do not regenerate this file without a reason you can write down.")
