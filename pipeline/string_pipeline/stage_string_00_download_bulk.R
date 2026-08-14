#!/usr/bin/env Rscript
# =============================================================================
# stage_string_00_download_bulk.R
#
# Retrieves the version-pinned STRING bulk files and records a checksummed
# manifest. Bulk download rather than the REST API: no rate limits, no version
# drift between calls, and the source file itself becomes an archivable
# artifact.
#
# The files are keyed by species + release, not by phenotype, so this runs
# once per (species, version, network_type) and is shared by every phenotype
# pinned to the same release. It is idempotent: a file already present whose
# md5 matches the manifest is left alone.
#
# Inputs:
#   config/<phenotype>.yml   (string: block — version, species_taxon,
#                             network_type, bulk.dir, bulk.source_url_base)
#
# Outputs (to string.bulk.dir):
#   <taxon>.protein.links.full.v<ver>.txt.gz     per-channel subscores + combined
#   <taxon>.protein.info.v<ver>.txt.gz           preferred display names
#   <taxon>.protein.aliases.v<ver>.txt.gz        identifier mapping
#   bulk_manifest.json                           url, bytes, md5, retrieval date
#
# Usage:
#   Rscript pipeline/string_pipeline/stage_string_00_download_bulk.R config/leptin.yml
#   Rscript pipeline/string_pipeline/stage_string_00_download_bulk.R config/leptin.yml --force
#
#   --force  re-download files that are already present and verified
#
# Run from the repository root.
# =============================================================================

if (!file.exists("R/utils.R")) {
  stop("run this script from the repository root (R/utils.R not found)")
}
source("R/utils.R")
source("R/string_utils.R")

# ---- arguments --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
args  <- args[args != "--force"]
if (length(args) != 1L) {
  stop("usage: Rscript pipeline/string_pipeline/stage_string_00_download_bulk.R config/<pheno>.yml [--force]")
}

cfg <- load_config(args[1])
s   <- string_cfg(cfg)

stage_banner("string-00", "Download version-pinned STRING bulk files")

message("STRING version pin: ", s$version)
message("Species taxon:      ", s$species_taxon)
message("Network type:       ", s$network_type)
message("Destination:        ", s$bulk$dir)

# ---- live release check (informational) -------------------------------------
# A mismatch is not an error. The pin is deliberate: it is what makes the
# analysis reproducible after STRING moves on. But an unannounced release is
# worth knowing about, so it is reported and recorded.
live <- NA_character_
if (isTRUE(s$bulk$verify_live_version)) {
  api  <- if (!is.null(s$bulk$version_api)) s$bulk$version_api else "https://string-db.org/api/json/version"
  live <- string_live_version(api)
  if (is.na(live)) {
    message("\nLive version check: unreachable (offline?) — continuing with the pin.")
  } else if (!identical(live, s$version)) {
    message("\nLive version check: string-db.org now serves ", live,
            ", this analysis is pinned to ", s$version, ".")
    message("  That is expected and intended. Nothing to do.")
  } else {
    message("\nLive version check: ", live, " — matches the pin.")
  }
}

# ---- resolve the file set ---------------------------------------------------
spec <- string_bulk_spec(s)
dir.create(s$bulk$dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(s$bulk$dir, "bulk_manifest.json")
prev <- if (file.exists(manifest_path)) {
  jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
} else NULL

prev_field <- function(filename, field) {
  if (is.null(prev) || is.null(prev$files)) return(NA_character_)
  for (f in prev$files) {
    if (identical(f$filename, filename)) {
      v <- f[[field]]
      return(if (is.null(v)) NA_character_ else as.character(v))
    }
  }
  NA_character_
}

# ---- fetch ------------------------------------------------------------------
records <- vector("list", nrow(spec))

for (i in seq_len(nrow(spec))) {
  key  <- spec$key[i]
  path <- spec$path[i]
  url  <- spec$url[i]
  message("\n[", key, "] ", spec$filename[i])

  known  <- prev_field(spec$filename[i], "md5")
  action <- "downloaded"

  if (file.exists(path) && !force) {
    have <- file_md5(path)
    if (!is.na(known) && identical(have, known)) {
      message("  present, md5 matches manifest — skipping")
      action <- "cached"
    } else if (is.na(known)) {
      message("  present but not in any manifest — adopting and checksumming it")
      action <- "adopted"
    } else {
      # Present, but a different file than the manifest describes. Do not
      # silently accept: the manifest is what downstream provenance cites.
      stop("file present but its md5 does not match bulk_manifest.json:\n    ", path,
           "\n  manifest: ", known, "\n  on disk:  ", have,
           "\n  Delete the file, or re-run with --force to replace it.")
    }
  } else {
    got <- download_bulk_file(string_urls_for(s, spec$dataset[i], spec$filename[i]),
                              path, timeout_s = s$bulk$download_timeout,
                              quiet = FALSE)
    url <- got$url          # the host that actually answered, for the manifest
  }

  sz  <- file.size(path)
  md5 <- file_md5(path)
  message("  ", format(sz, big.mark = ","), " bytes   md5 ", md5)

  # Retrieval date: the moment of download for a fresh file; for a cached one
  # the date the manifest already recorded, so re-running never rewrites
  # history. A file adopted without a manifest falls back to its mtime.
  retrieved <- if (action == "downloaded") {
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  } else {
    r <- prev_field(spec$filename[i], "retrieved_on")
    if (is.na(r)) format(file.mtime(path), "%Y-%m-%dT%H:%M:%S%z") else r
  }

  records[[i]] <- list(
    key          = key,
    dataset      = spec$dataset[i],
    filename     = spec$filename[i],
    url          = url,
    path         = path,
    bytes        = sz,
    md5          = md5,
    action       = action,
    retrieved_on = retrieved
  )
}

# ---- manifest ---------------------------------------------------------------
manifest <- list(
  string_version = s$version,
  species_taxon  = s$species_taxon,
  network_type   = s$network_type,
  url_bases_tried = as.list(string_url_bases(s)),
  live_version_at_download = live,
  written_on     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version      = R.version.string,
  files          = records
)
jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE,
                     null = "null", na = "null")

message("\n", paste(rep("-", 70), collapse = ""))
message("manifest -> ", manifest_path)
message("Bulk files ready. Next: pipeline/string_pipeline/stage_string_01_fetch_edges.R")
