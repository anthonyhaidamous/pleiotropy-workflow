# =============================================================================
# R/string_utils.R
# Shared helpers for the STRING module-detection stage group
# (pipeline/string_pipeline/stage_string_*.R). Source AFTER R/utils.R:
#     source("R/utils.R")
#     source("R/string_utils.R")
#
# Scope: version-pinned bulk retrieval, identifier mapping, induced-edge
# construction, provenance. No clustering, no statistics, no plotting —
# those belong to later stages.
#
# Conventions (inherited from R/utils.R)
# --------------------------------------
#   - Functions receive config as an argument, never via globals.
#   - No setwd(); every path is resolved relative to the repository root.
#   - Nothing that changes the analysis has a default in code. Anything the
#     analysis depends on comes from YAML and aborts loudly if absent.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
})


# -----------------------------------------------------------------------------
# Config: validate the `string:` block
#
# score_threshold has NO default. STRING's combined score is the definition of
# the dataset, not a display setting — a silently defaulted threshold would
# silently redefine the network, so its absence is a hard error.
# -----------------------------------------------------------------------------
string_cfg <- function(cfg) {
  s <- cfg$string
  if (is.null(s)) {
    stop("config has no 'string:' block — the STRING stage group cannot run.\n",
         "  See config/leptin.yml for the expected schema.")
  }

  required <- c("version", "species_taxon", "network_type", "score_threshold")
  absent <- required[vapply(required, function(k) is.null(s[[k]]), logical(1))]
  if (length(absent)) {
    stop("config string: block is missing required key(s): ",
         paste0("string.", absent, collapse = ", "))
  }

  if (!s$network_type %in% c("functional", "physical")) {
    stop("string.network_type must be 'functional' or 'physical', got: ",
         s$network_type)
  }

  thr <- suppressWarnings(as.integer(s$score_threshold))
  if (is.na(thr) || thr < 0 || thr > 1000) {
    stop("string.score_threshold must be an integer in [0, 1000], got: ",
         s$score_threshold)
  }
  s$score_threshold <- thr
  s$species_taxon   <- as.integer(s$species_taxon)
  s$version         <- as.character(s$version)

  # Structural defaults only — these change where files land, never what the
  # analysis says.
  if (is.null(s$bulk$dir))              s$bulk$dir <- "data/string/bulk"
  if (is.null(s$bulk$download_timeout)) s$bulk$download_timeout <- 3600
  if (is.null(s$output$edges_dir))     s$output$edges_dir <- "data/string"
  if (is.null(s$output$stage_subdir))  s$output$stage_subdir <- "string_modules"
  if (is.null(s$read$chunk_size))      s$read$chunk_size <- 1e6

  s
}


# STRING version as it appears in filenames: "12.0" -> "12-0".
# Dots are replaced so the version never reads as a file extension.
string_version_tag <- function(version) {
  gsub("[^A-Za-z0-9]+", "-", as.character(version))
}


# -----------------------------------------------------------------------------
# Download hosts, in the order they are tried.
#
# The default is DERIVED from the version pin: string-db.org's own version API
# advertises a per-release stable address (version-12-0.string-db.org), which
# keeps serving v12.0 after the project moves to v13. Pointing a pinned
# analysis at a host that follows the current release would quietly break the
# pin the first time STRING ships a new version.
#
# Which host served a file does not change the file — the manifest records the
# md5 either way — so falling back on transport failure is safe. It records
# which host actually answered.
# -----------------------------------------------------------------------------
string_url_bases <- function(s) {
  primary <- s$bulk$source_url_base
  if (is.null(primary) || !nzchar(primary)) {
    primary <- sprintf("https://version-%s.string-db.org/download",
                       string_version_tag(s$version))
  }
  bases <- c(primary, unlist(s$bulk$fallback_url_bases, use.names = FALSE))
  unique(sub("/+$", "", bases[nzchar(bases)]))
}


# -----------------------------------------------------------------------------
# Bulk source specification
#
# Filenames and dataset directories are DERIVED from version + taxon +
# network_type rather than typed into YAML. A hand-entered filename can
# disagree with the declared network_type; a derived one cannot.
# -----------------------------------------------------------------------------
string_bulk_spec <- function(s) {
  v   <- s$version
  tax <- s$species_taxon

  links_dataset <- if (identical(s$network_type, "physical")) {
    sprintf("protein.physical.links.full.v%s", v)
  } else {
    sprintf("protein.links.full.v%s", v)
  }

  spec <- tibble::tibble(
    key     = c("links", "info", "aliases"),
    dataset = c(links_dataset,
                sprintf("protein.info.v%s", v),
                sprintf("protein.aliases.v%s", v))
  )
  spec$filename <- sprintf("%d.%s.txt.gz", tax, spec$dataset)
  spec$path     <- file.path(s$bulk$dir, spec$filename)
  spec$url      <- sprintf("%s/%s/%s", string_url_bases(s)[1],
                           spec$dataset, spec$filename)
  spec
}

# All candidate URLs for one dataset, in try order.
string_urls_for <- function(s, dataset, filename) {
  sprintf("%s/%s/%s", string_url_bases(s), dataset, filename)
}


# -----------------------------------------------------------------------------
# Output paths for the frozen network bundle
#
# The three files below travel together and describe one network at one
# version at one threshold. They are written once and never regenerated
# silently.
# -----------------------------------------------------------------------------
string_paths <- function(cfg, s) {
  tag  <- string_version_tag(s$version)
  slug <- cfg$phenotype$slug
  stem <- sprintf("%s_edges_v%s_s%d", slug, tag, s$score_threshold)

  list(
    tag        = tag,
    stem       = stem,
    edges_dir  = s$output$edges_dir,
    edges      = file.path(s$output$edges_dir, paste0(stem, ".tsv")),
    provenance = file.path(s$output$edges_dir, paste0(stem, ".provenance.json")),
    mapping    = file.path(s$output$edges_dir,
                           sprintf("%s_gene_mapping_v%s.tsv", slug, tag)),
    manifest   = file.path(s$bulk$dir, "bulk_manifest.json"),
    stage_dir  = file.path(cfg$paths$intermediate, s$output$stage_subdir)
  )
}


# Gene list path: absolute / explicitly-relative paths are used as given,
# a bare filename is resolved under the phenotype's intermediate directory.
resolve_gene_list_path <- function(cfg, s) {
  p <- s$input$gene_list
  if (is.null(p) || !nzchar(p)) {
    stop("config must define string.input.gene_list")
  }
  if (grepl("^(/|\\./|\\.\\./|[A-Za-z]:)", p)) p else file.path(cfg$paths$intermediate, p)
}


# -----------------------------------------------------------------------------
# Gene list
#
# Strips CR explicitly. Gene lists in this repo are LF, but files that have
# been through a Windows editor are CRLF, and readLines() would then yield
# "ABHD1\r" — which fails to match any STRING alias and would be reported as
# a genuinely unmappable gene. A wrong answer that looks like a real result.
# -----------------------------------------------------------------------------
read_gene_list <- function(path) {
  if (!file.exists(path)) stop("gene list not found: ", path)
  x <- readLines(path, warn = FALSE)
  x <- gsub("\r", "", x, fixed = TRUE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  n_dup <- sum(duplicated(x))
  if (n_dup > 0L) {
    message("  note: ", n_dup, " duplicate symbol(s) in the gene list, collapsed")
  }
  unique(x)
}


# -----------------------------------------------------------------------------
# Checksums and provenance
# tools::md5sum() is base R — one fewer dependency between the reproducibility
# anchor and the claim that it has not changed.
# -----------------------------------------------------------------------------
file_md5 <- function(path) unname(tools::md5sum(path))

package_pins <- function(pkgs) {
  as.list(vapply(pkgs, function(p) as.character(utils::packageVersion(p)),
                 character(1)))
}


# -----------------------------------------------------------------------------
# Download
#
# The gzip magic-number check is not paranoia: a mistyped dataset directory
# returns an HTML error page with HTTP 200 at some mirrors. Saved as .txt.gz
# it would fail much later, inside a chunked read, with an unrelated message.
# -----------------------------------------------------------------------------
download_bulk_file <- function(urls, dest, timeout_s = 3600, quiet = FALSE) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".part")

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = timeout_s)

  failures <- character(0)

  for (url in urls) {
    if (file.exists(tmp)) unlink(tmp)
    message("  GET ", url)

    err <- tryCatch({
      if (requireNamespace("curl", quietly = TRUE)) {
        curl::curl_download(url, tmp, quiet = quiet, mode = "wb")
      } else {
        utils::download.file(url, tmp, mode = "wb", quiet = quiet)
      }
      NULL
    }, error = function(e) conditionMessage(e))

    if (!is.null(err)) {
      message("    failed: ", err)
      failures <- c(failures, paste0(url, "\n      ", err))
      next
    }

    con <- file(tmp, "rb"); magic <- readBin(con, "raw", 2L); close(con)
    if (!identical(magic, as.raw(c(0x1f, 0x8b)))) {
      unlink(tmp)
      message("    failed: response is not gzip (error page?)")
      failures <- c(failures, paste0(url, "\n      not gzip-compressed"))
      next
    }

    if (!file.rename(tmp, dest)) {
      file.copy(tmp, dest, overwrite = TRUE); unlink(tmp)
    }
    return(list(path = dest, url = url))
  }

  stop("could not retrieve ", basename(dest), " from any configured host:\n    ",
       paste(failures, collapse = "\n    "))
}


# Live release check. A mismatch against the pin is INFORMATION, not an error:
# the pin is deliberate. Network failure is also not an error — an offline
# rerun of a pinned analysis must still work.
string_live_version <- function(api_url = "https://string-db.org/api/json/version",
                                timeout_s = 20) {
  tryCatch({
    old <- getOption("timeout"); on.exit(options(timeout = old), add = TRUE)
    options(timeout = timeout_s)
    txt <- paste(readLines(api_url, warn = FALSE), collapse = "")
    as.character(jsonlite::fromJSON(txt)$string_version[1])
  }, error = function(e) NA_character_)
}


# -----------------------------------------------------------------------------
# Bulk table readers
#
# quote = "" throughout: the annotation column in protein.info contains
# apostrophes and quote marks that would otherwise swallow field boundaries.
# Columns are renamed by position because the header's first field is
# "#string_protein_id" and the leading "#" has moved between releases.
# -----------------------------------------------------------------------------
load_string_info <- function(path) {
  if (!file.exists(path)) stop("STRING info file not found: ", path)
  x <- readr::read_tsv(path, quote = "", progress = FALSE,
                       col_types = readr::cols(.default = readr::col_character()))
  if (ncol(x) < 2L) stop("unexpected column count in ", path)
  names(x)[1:2] <- c("string_id", "preferred_name")
  x[, c("string_id", "preferred_name")]
}

load_string_aliases <- function(path) {
  if (!file.exists(path)) stop("STRING aliases file not found: ", path)
  x <- readr::read_tsv(path, quote = "", progress = FALSE,
                       col_types = readr::cols(.default = readr::col_character()))
  if (ncol(x) < 3L) stop("unexpected column count in ", path)
  names(x)[1:3] <- c("string_id", "alias", "source")
  x[, c("string_id", "alias", "source")]
}


# -----------------------------------------------------------------------------
# Symbol -> STRING protein id
#
# Resolution is ranked, and the best available rank wins:
#   1  exact match on preferred_name (protein.info)
#   2  exact match on any alias      (protein.aliases)
#   3  case-insensitive preferred_name
#   4  case-insensitive alias
#
# A symbol that resolves to MORE THAN ONE protein at its best rank is reported
# as ambiguous and dropped — it is not silently resolved by picking the first.
# Which alias source produced a match is recorded so ranks 3-4 in particular
# can be eyeballed.
#
# Returns: list(mapped, unmapped, ambiguous, table)
#   table  one row per submitted symbol, in submission order — the audit trail.
# -----------------------------------------------------------------------------
map_symbols_to_string <- function(symbols, info, aliases) {
  sym <- unique(symbols)
  q   <- tibble::tibble(input_symbol = sym, .key = sym)
  qi  <- tibble::tibble(input_symbol = sym, .key = toupper(sym))

  hit <- function(query, lookup, key_col, src_col, rank) {
    l <- tibble::tibble(.key      = lookup[[key_col]],
                        string_id = lookup$string_id,
                        match_source = if (is.null(src_col)) "preferred_name"
                                       else paste0("alias:", lookup[[src_col]]))
    out <- dplyr::inner_join(query, l, by = ".key", relationship = "many-to-many")
    out$match_rank <- rank
    out[, c("input_symbol", "string_id", "match_source", "match_rank")]
  }

  info_ci    <- dplyr::mutate(info,    .u = toupper(preferred_name))
  aliases_ci <- dplyr::mutate(aliases, .u = toupper(alias))

  cand <- dplyr::bind_rows(
    hit(q,  info,       "preferred_name", NULL,     1L),
    hit(q,  aliases,    "alias",          "source", 2L),
    hit(qi, info_ci,    ".u",             NULL,     3L),
    hit(qi, aliases_ci, ".u",             "source", 4L)
  )

  if (nrow(cand) == 0L) {
    best <- cand[0, ]
  } else {
    cand <- dplyr::distinct(cand, input_symbol, string_id, match_rank, .keep_all = TRUE)
    best <- cand |>
      dplyr::group_by(input_symbol) |>
      dplyr::filter(match_rank == min(match_rank)) |>
      dplyr::ungroup()
  }

  counts <- best |>
    dplyr::group_by(input_symbol) |>
    dplyr::summarise(n_ids = dplyr::n_distinct(string_id), .groups = "drop")

  ok_syms  <- counts$input_symbol[counts$n_ids == 1L]
  amb_syms <- counts$input_symbol[counts$n_ids  > 1L]

  resolved <- best |>
    dplyr::filter(input_symbol %in% ok_syms) |>
    dplyr::arrange(input_symbol, match_source) |>   # deterministic source label
    dplyr::group_by(input_symbol) |>
    dplyr::slice(1L) |>                             # string_id identical within group
    dplyr::ungroup() |>
    dplyr::left_join(info, by = "string_id")

  # Two submitted symbols can resolve to the SAME protein — a current symbol
  # and a retired alias of it, say. Left alone that protein would enter the
  # network once but be counted twice in the gene accounting, and any
  # symbol-keyed lookup would silently keep whichever came first. The second
  # symbol is excluded and named, with the symbol it collides with.
  dup_idx    <- duplicated(resolved$string_id)
  duplicates <- resolved[dup_idx, , drop = FALSE]
  if (nrow(duplicates) > 0L) {
    keeper <- resolved$input_symbol[!dup_idx]
    keyed  <- resolved$string_id[!dup_idx]
    duplicates$duplicate_of <- keeper[match(duplicates$string_id, keyed)]
  } else {
    duplicates$duplicate_of <- character(0)
  }
  mapped <- resolved[!dup_idx, , drop = FALSE]

  ambiguous <- best |>
    dplyr::filter(input_symbol %in% amb_syms) |>
    dplyr::left_join(info, by = "string_id") |>
    dplyr::arrange(input_symbol, string_id)

  unmapped <- setdiff(sym, counts$input_symbol)

  # Audit trail: every submitted symbol, in submission order, with the reason
  # it is or is not in the network.
  tbl <- tibble::tibble(input_symbol = sym) |>
    dplyr::left_join(
      dplyr::select(resolved, input_symbol, string_id, preferred_name,
                    match_source, match_rank),
      by = "input_symbol") |>
    dplyr::mutate(
      status = dplyr::case_when(
        input_symbol %in% amb_syms                ~ "ambiguous",
        input_symbol %in% duplicates$input_symbol ~ "duplicate_protein",
        is.na(string_id)                          ~ "unmapped",
        TRUE                                      ~ "mapped"),
      n_candidate_ids = counts$n_ids[match(input_symbol, counts$input_symbol)],
      n_candidate_ids = ifelse(is.na(n_candidate_ids), 0L, n_candidate_ids),
      note = dplyr::case_when(
        status == "duplicate_protein" ~ paste0(
          "same STRING protein as ",
          duplicates$duplicate_of[match(input_symbol, duplicates$input_symbol)]),
        status == "ambiguous" ~ "resolves to >1 protein at its best match rank",
        TRUE ~ NA_character_)
    )

  list(mapped = mapped, unmapped = unmapped, ambiguous = ambiguous,
       duplicates = duplicates, table = tbl)
}


# -----------------------------------------------------------------------------
# Induced edges from the bulk links file
#
# Read in chunks and filter each chunk down to the submitted proteins before
# binding. The unfiltered human links.full file is ~13M rows x 16 columns;
# reading it whole costs several GB. Chunking bounds memory regardless of how
# large the source is, and works unchanged on a pre-filtered download.
# -----------------------------------------------------------------------------
read_links_induced <- function(path, keep_ids, threshold, chunk_size = 1e6,
                               quiet = FALSE) {
  if (!file.exists(path)) stop("STRING links file not found: ", path)
  keep_ids <- unique(keep_ids)
  n_seen <- 0

  cb <- function(chunk, pos) {
    n_seen <<- n_seen + nrow(chunk)
    if (!quiet) {
      message("    scanned ", format(n_seen, big.mark = ","), " rows")
    }
    chunk[chunk$protein1 %in% keep_ids &
          chunk$protein2 %in% keep_ids &
          chunk$combined_score >= threshold, , drop = FALSE]
  }

  out <- readr::read_delim_chunked(
    path,
    callback   = readr::DataFrameCallback$new(cb),
    chunk_size = chunk_size,
    delim      = " ",
    quote      = "",
    col_types  = readr::cols(protein1 = readr::col_character(),
                             protein2 = readr::col_character(),
                             .default = readr::col_integer()),
    progress   = FALSE
  )
  if (is.null(out)) out <- tibble::tibble()
  attr(out, "rows_scanned") <- n_seen
  out
}


# -----------------------------------------------------------------------------
# Reciprocal de-duplication
#
# STRING stores each association twice (A-B and B-A). The second row carries
# no information. Endpoints are put in canonical order, then duplicates
# dropped. Reciprocal rows are asserted to agree on combined_score first: if
# they disagree the file is not the symmetric storage we assume, and every
# downstream degree and density would be wrong.
# -----------------------------------------------------------------------------
dedup_reciprocal <- function(edges) {
  if (nrow(edges) == 0L) return(edges)

  a <- pmin(edges$protein1, edges$protein2)
  b <- pmax(edges$protein1, edges$protein2)

  n_self <- sum(a == b)
  if (n_self > 0L) {
    message("  dropped ", n_self, " self-loop row(s)")
    keep <- a != b
    edges <- edges[keep, , drop = FALSE]; a <- a[keep]; b <- b[keep]
  }

  edges$protein1 <- a
  edges$protein2 <- b
  key <- paste(a, b, sep = "\t")

  disagree <- tapply(edges$combined_score, key, function(z) length(unique(z)) > 1L)
  if (any(disagree)) {
    bad <- names(disagree)[disagree]
    stop("reciprocal rows disagree on combined_score for ", length(bad),
         " pair(s), e.g. ", sub("\t", " / ", bad[1]),
         "\nThe links file is not symmetric as assumed; aborting.")
  }

  out <- edges[!duplicated(key), , drop = FALSE]
  attr(out, "rows_before_dedup") <- nrow(edges)
  out
}


# -----------------------------------------------------------------------------
# Freeze
#
# An existing edge file is never destroyed. --force moves it aside with a
# timestamp rather than overwriting, so a network that something downstream
# was already built against remains recoverable.
# -----------------------------------------------------------------------------
write_frozen_edges <- function(edges, path, force = FALSE) {
  if (file.exists(path)) {
    if (!force) {
      stop("refusing to overwrite the frozen edge file:\n    ", path, "\n",
           "This file is the reproducibility anchor for every downstream stage.\n",
           "Re-run with --force only if you intend to redefine the network.")
    }
    superseded <- sprintf("%s.superseded_%s", path,
                          format(Sys.time(), "%Y%m%d_%H%M%S"))
    file.rename(path, superseded)
    message("  existing edge file moved aside -> ", superseded)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(edges, path, na = "")
  path
}


# =============================================================================
# Stage 02 helpers — graph construction and component accounting
# =============================================================================

# -----------------------------------------------------------------------------
# Read the frozen edge list, verifying it against its provenance sidecar.
# A checksum mismatch means the anchor moved under us and every number computed
# downstream of it is invalid, so it aborts rather than warns.
# -----------------------------------------------------------------------------
read_frozen_edges <- function(edge_path, provenance_path, verify = TRUE) {
  if (!file.exists(edge_path)) {
    stop("frozen edge file not found: ", edge_path,
         "\n  Run pipeline/string_pipeline/stage_string_01_fetch_edges.R first.")
  }
  if (verify) {
    if (!file.exists(provenance_path)) {
      stop("provenance sidecar not found: ", provenance_path,
           "\n  The edge file cannot be verified. Re-run stage 01.")
    }
    prov <- jsonlite::fromJSON(provenance_path, simplifyVector = FALSE)
    have <- file_md5(edge_path)
    if (!identical(have, prov$edge_md5)) {
      stop("frozen edge file does not match its provenance record:\n",
           "    file:      ", edge_path, "\n",
           "    recorded:  ", prov$edge_md5, "\n",
           "    on disk:   ", have, "\n",
           "  The reproducibility anchor has changed. Nothing downstream of it\n",
           "  can be trusted until this is resolved.")
    }
  }
  readr::read_tsv(edge_path, quote = "", progress = FALSE,
                  col_types = readr::cols(
                    protein1 = readr::col_character(),
                    protein2 = readr::col_character(),
                    gene1    = readr::col_character(),
                    gene2    = readr::col_character(),
                    .default = readr::col_integer()))
}


# LD block of origin per gene symbol, from stage 10's by-block table.
# Returns a named character vector (symbol -> block), or NULL if not configured.
load_block_map <- function(path, gene_col = "external_gene_name",
                           block_col = "Block",
                           policy = "error",
                           gene_start_col   = "start_position",
                           gene_end_col     = "end_position",
                           region_start_col = "region_start",
                           region_end_col   = "region_end") {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) stop("block table not found: ", path)
  b <- readr::read_tsv(path, quote = "", progress = FALSE,
                       col_types = readr::cols(.default = readr::col_character()))
  for (cl in c(gene_col, block_col)) {
    if (!cl %in% names(b)) {
      stop("block table is missing column '", cl, "': ", path)
    }
  }
  dup <- unique(b[[gene_col]][duplicated(b[[gene_col]])])
  if (length(dup)) {
    # A gene in two blocks makes "block of origin" ill-defined and would
    # double-count in any block x module table. This is not hypothetical: a
    # gene long enough can contain two independent association signals (CLL's
    # LPP spans 738 kb and holds B11 and B12 outright), and the +/-10 kb
    # extraction padding can make two adjacent blocks' windows overlap.
    #
    # Leptin's 15 blocks are far enough apart that it never arises. CLL's 49
    # are not. There is no defensible default here, so the policy is declared
    # in the config and whichever way it resolves is named on the console.
    if (identical(policy, "error")) {
      stop(length(dup), " gene(s) appear in more than one LD block, e.g. ",
           paste(utils::head(dup, 3), collapse = ", "),
           "\n  Block of origin must be unique per gene for the overlay to mean anything.",
           "\n  Set string.input.multi_block_policy to 'largest_overlap' or 'drop'",
           "\n  if this phenotype legitimately has genes spanning two blocks.")
    }
    if (!policy %in% c("largest_overlap", "drop")) {
      stop("unknown string.input.multi_block_policy: '", policy,
           "'  (expected 'error', 'largest_overlap' or 'drop')")
    }

    message("  ", length(dup), " gene(s) span more than one LD block — policy '",
            policy, "':")
    drop_rows <- rep(FALSE, nrow(b))
    for (gname in dup) {
      i <- which(b[[gene_col]] == gname)
      if (identical(policy, "drop")) {
        drop_rows[i] <- TRUE
        message(sprintf("    %-12s %s  -> dropped from the block overlay",
                        gname, paste(b[[block_col]][i], collapse = ", ")))
        next
      }
      # largest_overlap: base-pair overlap between the gene body and each
      # block's extraction window. Ties go to the lowest block number, so the
      # result does not depend on row order in the table.
      gs <- as.numeric(b[[gene_start_col]][i]);   ge <- as.numeric(b[[gene_end_col]][i])
      rs <- as.numeric(b[[region_start_col]][i]); re <- as.numeric(b[[region_end_col]][i])
      ov <- pmax(0, pmin(ge, re) - pmax(gs, rs))
      bn <- suppressWarnings(as.integer(sub("^[^0-9]*", "", b[[block_col]][i])))
      win <- order(-ov, bn)[1]
      drop_rows[i[-win]] <- TRUE
      message(sprintf("    %-12s %s  -> %s", gname,
                      paste(sprintf("%s (%.0f kb)", b[[block_col]][i], ov / 1000),
                            collapse = "  "),
                      b[[block_col]][i[win]]))
    }
    b <- b[!drop_rows, , drop = FALSE]
  }
  setNames(b[[block_col]], b[[gene_col]])
}


# -----------------------------------------------------------------------------
# Build the graph.
#
# Vertices are ALL mapped proteins, including those with no edges. Building
# from the edge list alone would drop zero-degree genes before they could be
# counted — precisely the silent filter the accounting exists to prevent.
# -----------------------------------------------------------------------------
build_string_graph <- function(edges, mapping, block_map = NULL) {
  mapped <- mapping[mapping$status == "mapped", ]
  v <- tibble::tibble(name         = mapped$string_id,
                      gene         = mapped$input_symbol,
                      preferred    = mapped$preferred_name)
  v$block <- if (is.null(block_map)) NA_character_ else unname(block_map[v$gene])

  el <- as.matrix(edges[, c("protein1", "protein2")])
  g  <- igraph::graph_from_data_frame(
    d = data.frame(from = el[, 1], to = el[, 2],
                   combined_score = edges$combined_score,
                   stringsAsFactors = FALSE),
    directed = FALSE, vertices = as.data.frame(v))
  g
}


# -----------------------------------------------------------------------------
# Per-gene accounting: one row per SUBMITTED symbol, with the reason it is or
# is not in the analysis. Nothing leaves invisibly.
# -----------------------------------------------------------------------------
component_accounting <- function(g, mapping, block_map = NULL) {
  comp <- igraph::components(g)
  lcc  <- which.max(comp$csize)

  vt <- tibble::tibble(
    string_id      = igraph::V(g)$name,
    input_symbol   = igraph::V(g)$gene,
    degree         = as.integer(igraph::degree(g)),
    component_id   = as.integer(comp$membership),
    component_size = as.integer(comp$csize[comp$membership])
  )
  vt$in_lcc <- vt$component_id == lcc

  acc <- mapping |>
    dplyr::select(input_symbol, string_id, preferred_name, status, match_rank, note) |>
    dplyr::left_join(dplyr::select(vt, -input_symbol), by = "string_id") |>
    dplyr::mutate(
      block  = if (is.null(block_map)) NA_character_ else unname(block_map[input_symbol]),
      status = dplyr::case_when(
        status != "mapped"          ~ status,      # unmapped / ambiguous / duplicate
        degree == 0L                ~ "singleton", # mapped but no association
        in_lcc                      ~ "lcc",
        TRUE                        ~ "island"     # in a smaller component
      )
    ) |>
    dplyr::arrange(dplyr::desc(!is.na(component_size)), component_size, input_symbol)

  attr(acc, "lcc_component") <- lcc
  acc
}


# The Q1b table: the six quantities, in the order the method notes fix.
summarise_accounting <- function(acc) {
  tibble::tibble(
    quantity = c("n_in", "n_mapped", "n_deg_ge1", "n_lcc", "n_singleton", "n_island"),
    n = c(
      nrow(acc),
      sum(!acc$status %in% c("unmapped", "ambiguous", "duplicate_protein")),
      sum(acc$status %in% c("lcc", "island")),
      sum(acc$status == "lcc"),
      sum(acc$status == "singleton"),
      sum(acc$status == "island")
    ),
    meaning = c(
      "gene symbols submitted",
      "resolved to a unique STRING protein",
      "mapped and carrying at least one association",
      "in the largest connected component",
      "mapped but with no association at this threshold",
      "in a connected component other than the largest"
    )
  )
}


# Which node set flows downstream. No default at the call site: the value comes
# from string.components.downstream_node_set in the config.
select_downstream <- function(g, which) {
  if (!which %in% c("lcc", "degree_ge1")) {
    stop("string.components.downstream_node_set must be 'lcc' or 'degree_ge1', got: ",
         which)
  }
  if (identical(which, "lcc")) {
    comp <- igraph::components(g)
    keep <- igraph::V(g)[comp$membership == which.max(comp$csize)]
  } else {
    keep <- igraph::V(g)[igraph::degree(g) > 0]
  }
  igraph::induced_subgraph(g, keep)
}
