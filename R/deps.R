# =============================================================================
# R/deps.R
# The project's package manifest, and the two questions every new machine asks
# of it: what is missing here, and how do I get it.
#
#     source("R/deps.R")     # from the project root
#
# Scope and rules
# ---------------
#   - Base R only, and no side effects on source. This file has to run on a
#     machine where nothing is installed yet, which rules out optparse, yaml,
#     and even BiocManager (install_deps() installs that one itself when a
#     Bioconductor package is actually wanted).
#   - PROJECT_DEPS is the single place a dependency is declared. tools/setup_env.R
#     is a thin CLI over it; the README package list is prose that mirrors it.
#     Add a library() call to a script, add the package here.
#   - `required = FALSE` means the scripts degrade gracefully without it and say
#     so on the console (ggchicklet -> square bars, ragg -> the base png device,
#     curl -> download.file). Never make a figure silently different: if the
#     absence changes a number, the package is required.
#
# Two machines, one repo
# ----------------------
# On the Linux box the pleio container already carries every package in its own
# site library; there, renv is skipped (PLEIO_NO_RENV=1, see .Rprofile) and this
# file only ever reports "all present". On the laptop renv owns a project
# library under renv/library/, and setup_env.R fills it. Both paths install into
# whatever .libPaths()[1] is at the time, so neither needs to know about the
# other.
# =============================================================================


# -----------------------------------------------------------------------------
# The manifest
#
# source:  "CRAN"  installed with install.packages()
#          "Bioc"  installed with BiocManager::install() / renv's bioc:: prefix
#          "<user>/<repo>"  a GitHub remote (ggchicklet is not on CRAN)
# -----------------------------------------------------------------------------
PROJECT_DEPS <- rbind(
  # -- core: attached by R/utils.R and R/string_utils.R, so every stage needs them
  data.frame(package = c("yaml", "dplyr", "tidyr", "stringr", "jsonlite",
                         "gtools", "readr", "tibble"),
             source = "CRAN", group = "core", required = TRUE),

  # -- pipeline stages
  data.frame(package = "LDlinkR", source = "CRAN", group = "stage 2 (LDlink)",
             required = TRUE),
  data.frame(package = c("biomaRt", "GenomicRanges", "ensembldb", "AnnotationHub"),
             source = "Bioc", group = "stage 10 (gene extraction)", required = TRUE),
  data.frame(package = "igraph", source = "CRAN", group = "STRING stages",
             required = TRUE),

  # -- figures, tables and analysis scripts
  data.frame(package = c("optparse", "ggplot2", "ggrepel", "scales", "patchwork",
                         "svglite", "httr"),
             source = "CRAN", group = "figures", required = TRUE),
  data.frame(package = "circlize", source = "CRAN",
             group = "figures (chord diagram)", required = TRUE),
  data.frame(package = c("pdftools"), source = "CRAN",
             group = "figures (karyoplot)", required = TRUE),
  data.frame(package = "karyoploteR", source = "Bioc",
             group = "figures (karyoplot)", required = TRUE),

  # -- nice to have; the scripts check for these themselves and carry on without
  data.frame(package = "ragg", source = "CRAN", group = "optional", required = FALSE),
  data.frame(package = "curl", source = "CRAN", group = "optional", required = FALSE),
  data.frame(package = "ggchicklet", source = "hrbrmstr/ggchicklet",
             group = "optional", required = FALSE),

  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# Which of them this machine actually has
#
# `find.package` rather than `require`: asking whether a package is INSTALLED
# must not attach it, load its namespace, or run .onLoad. Several of the
# Bioconductor packages print banners and mask base generics on attach, and a
# status report has no business doing that.
# -----------------------------------------------------------------------------
deps_status <- function(deps = PROJECT_DEPS) {
  deps$installed <- vapply(deps$package, function(p) {
    length(find.package(p, quiet = TRUE)) > 0L
  }, logical(1))
  deps
}

missing_deps <- function(deps = PROJECT_DEPS, required_only = TRUE) {
  st <- deps_status(deps)
  if (required_only) st <- st[st$required, , drop = FALSE]
  st[!st$installed, , drop = FALSE]
}


# -----------------------------------------------------------------------------
# Where packages come from
#
# On Linux, plain CRAN means compiling every package from source, and that is
# what makes a fresh laptop painful: karyoploteR reaches biovizBase -> Hmisc ->
# htmlTable -> htmlwidgets -> rmarkdown -> bslib -> sass -> fs, and `fs` stops
# on a missing libuv1-dev the same way pdftools stops on libpoppler-cpp-dev. A
# missing system header, three packages deep, for a chord diagram.
#
# Posit Package Manager serves prebuilt binaries for the exact distribution from
# the ordinary src/contrib path, so the whole chain becomes downloads and no -dev
# package is needed (the binary links against the runtime .so, which a desktop
# already has). It hands out binaries only to a client that identifies its R
# version and platform, which R does through the HTTP user agent — hence
# deps_use_repos() below, and hence this being useless without it.
#
# On Windows and macOS CRAN already ships binaries, so this returns plain CRAN.
# Bioconductor has no PPM equivalent here; those few still build from source,
# which works because their compiled dependencies are already in place.
# -----------------------------------------------------------------------------
deps_repos <- function() {
  cran <- c(CRAN = "https://cloud.r-project.org")
  if (Sys.info()[["sysname"]] != "Linux") return(cran)
  codename <- linux_codename()
  if (is.na(codename)) return(cran)
  c(PPM = sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename),
    cran)
}

# Point this session at those repos, keeping any others (the Bioconductor ones
# renv restores from the lockfile) and replacing only the CRAN mirror.
deps_use_repos <- function() {
  repos <- deps_repos()
  if ("PPM" %in% names(repos)) {
    options(HTTPUserAgent = sprintf(
      "R/%s R (%s)", getRversion(),
      paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))
  }
  keep <- getOption("repos", character())
  nm <- names(keep); if (is.null(nm)) nm <- rep("", length(keep))
  keep <- keep[!nm %in% c("CRAN", "PPM", "P3M", "RSPM") &
                 !grepl("cran|packagemanager", keep, ignore.case = TRUE)]
  options(repos = c(repos, keep))
  invisible(getOption("repos"))
}

linux_codename <- function() {
  if (!file.exists("/etc/os-release")) return(NA_character_)
  kv <- tryCatch(readLines("/etc/os-release", warn = FALSE),
                 error = function(e) character())
  code <- sub("^VERSION_CODENAME=", "", grep("^VERSION_CODENAME=", kv, value = TRUE))
  if (!length(code)) return(NA_character_)
  gsub('"', "", code[1])
}


# -----------------------------------------------------------------------------
# Install
#
# Split by whether the lockfile has an opinion, not by machine:
#
#   in renv.lock       renv::restore(), which pins the version that produced the
#                      locked numbers. This is the reproducibility path and it
#                      ignores the binary repos above on purpose — the recorded
#                      version wins over a fast install.
#   everything else    install.packages() / BiocManager into .libPaths()[1],
#                      which is the renv project library when renv is active and
#                      the container's site library when it is not. Both cases
#                      want the binaries, so both go through deps_use_repos().
#
# renv::install() would also serve the second group, but the route below is the
# one verified to come back with binaries here: renv left getOption("repos") on
# cloud.r-project.org and compiled the whole karyoploteR chain from source.
# -----------------------------------------------------------------------------
renv_active <- function() {
  isTRUE(requireNamespace("renv", quietly = TRUE)) &&
    nzchar(Sys.getenv("RENV_PROJECT"))
}

# Packages recorded in the lockfile, i.e. the ones restore() can pin.
locked_packages <- function(lockfile = "renv.lock") {
  if (!file.exists(lockfile)) return(character())
  txt <- readLines(lockfile, warn = FALSE)
  # Package entries are the "Package": "name" fields under Packages; reading the
  # JSON would mean depending on jsonlite, which may itself be missing here.
  m <- regmatches(txt, regexec('^\\s*"Package":\\s*"([^"]+)"', txt))
  unique(vapply(m[lengths(m) == 2L], `[`, character(1), 2L))
}

install_deps <- function(pkgs, deps = PROJECT_DEPS, quiet = FALSE) {
  if (!length(pkgs)) return(invisible(character()))
  say <- function(...) if (!quiet) message(...)

  locked <- if (renv_active()) intersect(pkgs, locked_packages()) else character()
  if (length(locked)) {
    say("- renv::restore (pinned by renv.lock): ", paste(locked, collapse = ", "))
    renv::restore(packages = locked, prompt = FALSE)
  }

  fresh <- setdiff(pkgs, locked)
  if (!length(fresh)) return(invisible(pkgs))

  deps_use_repos()
  rows <- deps[match(fresh, deps$package), , drop = FALSE]
  cran <- rows$package[rows$source == "CRAN"]
  bioc <- rows$package[rows$source == "Bioc"]
  gh   <- rows[grepl("/", rows$source, fixed = TRUE), , drop = FALSE]

  if (length(cran)) {
    say("- install.packages: ", paste(cran, collapse = ", "))
    install.packages(cran)
  }
  if (nrow(gh)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    say("- remotes::install_github: ", paste(gh$source, collapse = ", "))
    remotes::install_github(gh$source)
  }
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    say("- BiocManager::install: ", paste(bioc, collapse = ", "))
    BiocManager::install(bioc, ask = FALSE, update = FALSE)
  }
  invisible(pkgs)
}


# -----------------------------------------------------------------------------
# System dependencies: not R packages, and either one stops a fresh clone cold
#
#   mcl           the STRING module stage shells out to it (config string.mcl).
#                 Nothing installs this from R; it has to come from apt.
#   poppler-cpp   pdftools links against it. Two different requirements wear the
#                 same name, and conflating them is what sends people to apt for
#                 no reason: a SOURCE build needs the headers (libpoppler-cpp-dev),
#                 a binary needs only the shared library, which any desktop with
#                 a PDF viewer already has. We report the runtime, and mention
#                 the headers only when a source build is what is on offer.
# -----------------------------------------------------------------------------
system_deps_status <- function() {
  linux <- Sys.info()[["sysname"]] == "Linux"
  mcl   <- Sys.which("mcl")

  headers <- if (linux && nzchar(Sys.which("pkg-config"))) {
    system2("pkg-config", c("--exists", "poppler-cpp"), stdout = NULL, stderr = NULL) == 0L
  } else NA
  runtime <- if (linux && nzchar(Sys.which("ldconfig"))) {
    any(grepl("libpoppler-cpp",
              suppressWarnings(system2("ldconfig", "-p", stdout = TRUE, stderr = NULL))))
  } else NA

  list(mcl         = if (nzchar(mcl)) mcl else NA_character_,
       poppler_hdr = headers,
       poppler_lib = runtime,
       # On Windows and macOS CRAN ships a binary pdftools with poppler bundled,
       # so neither check applies and neither is reported.
       poppler_relevant = linux)
}
