#!/usr/bin/env Rscript
# =============================================================================
# tools/setup_env.R
# One command that makes a machine able to run this repo, and one command that
# tells you whether it already can.
#
#   Rscript tools/setup_env.R            install every missing REQUIRED package
#   Rscript tools/setup_env.R --check    report only; exit 1 if something required
#                                        is missing (this is what run_pipeline.sh
#                                        pre-flights with)
#   Rscript tools/setup_env.R --all      also install the optional three
#                                        (ragg, curl, ggchicklet)
#   Rscript tools/setup_env.R --dry-run  print what would be installed
#
# Why this exists
# ---------------
# The pipeline runs on two machines that resolve packages in completely
# different ways: the Linux box runs it inside the pleio container, where the
# image already carries every package; the laptop has a bare R and an renv
# project library. Before this script the second machine only discovered a gap
# by hitting it mid-run — `Error in library(circlize) : there is no package
# called 'circlize'`, forty minutes into a stage, with no statement anywhere of
# what the full set even was.
#
# So: R/deps.R holds the manifest, this reports and fills it, and .Rprofile
# points at this script whenever a library() call fails. Nothing here knows
# which of the two machines it is on — it installs into .libPaths()[1], which
# renv or the container has already pointed at the right place.
#
# Base R only, deliberately: it must run on a machine where the first package
# has not been installed yet, so no optparse for the flags and no yaml for the
# config.
# =============================================================================

if (!file.exists("R/deps.R")) stop("run from the repository root: Rscript tools/setup_env.R")
source("R/deps.R")

args      <- commandArgs(trailingOnly = TRUE)
known     <- c("--check", "--all", "--dry-run", "--quiet", "--help", "-h")
unknown   <- setdiff(args, known)
if (length(unknown)) {
  stop("unknown flag(s): ", paste(unknown, collapse = " "),
       "\nusage: Rscript tools/setup_env.R [--check] [--all] [--dry-run] [--quiet]")
}
if (any(c("--help", "-h") %in% args)) {
  writeLines(c(
    "usage: Rscript tools/setup_env.R [--check] [--all] [--dry-run] [--quiet]",
    "  --check    report only, exit 1 if a required package is missing",
    "  --all      include the optional packages (ragg, curl, ggchicklet)",
    "  --dry-run  print what would be installed, install nothing",
    "  --quiet    drop the inventory, print only what is wrong (for pre-flights)"))
  quit(save = "no", status = 0L)
}
check_only  <- "--check"   %in% args
want_all    <- "--all"     %in% args
dry_run     <- "--dry-run" %in% args
quiet       <- "--quiet"   %in% args

# --quiet keeps the verdict and drops the inventory: run_pipeline.sh pre-flights
# with it, and a clean environment should add nothing to the log.
say  <- function(...) if (!quiet) cat(...)
rule <- function() say(strrep("-", 72), "\n", sep = "")
pad  <- function(x, n) formatC(x, width = -n, flag = " ")


# -----------------------------------------------------------------------------
# Where this session would install to, and from
# -----------------------------------------------------------------------------
repos  <- deps_repos()
in_renv <- renv_active()

rule()
say("Environment\n")
rule()
say("  R           ", R.version.string, "\n")
say("  platform    ", R.version$platform, "\n")
if (Sys.info()[["sysname"]] == "Linux" && file.exists("/etc/os-release")) {
  pretty <- grep("^PRETTY_NAME=", readLines("/etc/os-release", warn = FALSE), value = TRUE)
  if (length(pretty)) say("  distro      ", gsub('"', "", sub("^PRETTY_NAME=", "", pretty[1])), "\n")
}
say("  renv        ", if (in_renv) "active" else "not active",
    if (!in_renv && identical(Sys.getenv("PLEIO_NO_RENV"), "1"))
      " (PLEIO_NO_RENV=1)" else "", "\n", sep = "")
say("  library     ", .libPaths()[1], "\n")
say("  repository  ", repos[[1]], "\n")
say("              ", if (names(repos)[1] == "PPM") "binaries, not source builds" else
                      "CRAN", "\n", sep = "")
if (in_renv) say("              packages in renv.lock are restored at the pinned version\n")
say("\n")


# -----------------------------------------------------------------------------
# Package status
# -----------------------------------------------------------------------------
st <- deps_status()
st <- st[order(!st$required, st$group, st$package), , drop = FALSE]

rule()
say("R packages\n")
rule()
group <- ""
for (i in seq_len(nrow(st))) {
  if (!identical(st$group[i], group)) {
    group <- st$group[i]
    say("  ", group, "\n", sep = "")
  }
  mark <- if (st$installed[i]) "ok     " else if (st$required[i]) "MISSING" else "absent "
  say("    ", mark, "  ", pad(st$package[i], 16),
      if (st$source[i] == "CRAN") "" else paste0("[", st$source[i], "]"), "\n", sep = "")
}
say("\n")

missing_req <- st$package[!st$installed & st$required]
missing_opt <- st$package[!st$installed & !st$required]


# -----------------------------------------------------------------------------
# Everything that is not an R package
# -----------------------------------------------------------------------------
sysd     <- system_deps_status()
binaries <- names(repos)[1] == "PPM" || Sys.info()[["sysname"]] != "Linux"

rule()
say("System\n")
rule()
say("  ", if (is.na(sysd$mcl)) "MISSING" else "ok     ", "   mcl",
    if (is.na(sysd$mcl)) "            not on PATH; the STRING module stage needs it"
    else paste0("            ", sysd$mcl), "\n", sep = "")
if (isTRUE(sysd$poppler_relevant)) {
  # pdftools is the only package here with a system library behind it, so say
  # which of the two poppler requirements is actually in play.
  ok <- if (binaries) isTRUE(sysd$poppler_lib) else isTRUE(sysd$poppler_hdr)
  say("  ", if (ok) "ok     " else "MISSING", "   poppler-cpp    ",
      if (binaries) "runtime library (pdftools installs as a binary)"
      else "headers (pdftools would be built from source)", "\n", sep = "")
}
say("  ", if (nzchar(Sys.getenv("LDLINK_TOKEN"))) "ok     " else "unset  ", "   LDLINK_TOKEN",
    if (nzchar(Sys.getenv("LDLINK_TOKEN"))) "" else "   stage 2 only; export before running it",
    "\n", sep = "")
say("\n")

hint_system <- function() {
  poppler_gap <- isTRUE(sysd$poppler_relevant) &&
    !isTRUE(if (binaries) sysd$poppler_lib else sysd$poppler_hdr)
  if (is.na(sysd$mcl) || poppler_gap) {
    cat("These come from the OS, not from R:\n")
    if (is.na(sysd$mcl)) cat("    sudo apt install mcl\n")
    if (poppler_gap)     cat("    sudo apt install ",
                             if (binaries) "libpoppler-cpp0v5" else "libpoppler-cpp-dev",
                             "   # pdftools links against it\n", sep = "")
    cat("\n")
  }
}


# -----------------------------------------------------------------------------
# Report, or act
# -----------------------------------------------------------------------------
if (check_only) {
  hint_system()
  if (length(missing_req)) {
    cat("MISSING ", length(missing_req), " required package(s): ",
        paste(missing_req, collapse = ", "), "\n", sep = "")
    cat("Install them with:  Rscript tools/setup_env.R\n")
    quit(save = "no", status = 1L)
  }
  say("All required packages present.\n")
  if (length(missing_opt)) {
    say("Optional, absent (scripts fall back and say so): ",
        paste(missing_opt, collapse = ", "), "\n", sep = "")
  }
  quit(save = "no", status = 0L)
}

to_install <- if (want_all) c(missing_req, missing_opt) else missing_req
hint_system()

if (!length(to_install)) {
  say("Nothing to install: every required package is present",
      if (want_all) " and so are the optional ones" else "", ".\n", sep = "")
  if (!want_all && length(missing_opt)) {
    say("Optional, absent: ", paste(missing_opt, collapse = ", "),
        "   (--all installs these too)\n", sep = "")
  }
  quit(save = "no", status = 0L)
}

rule()
cat("Installing ", length(to_install), " package(s) into\n  ", .libPaths()[1], "\n", sep = "")
rule()
if (dry_run) {
  cat("--dry-run: would install ", paste(to_install, collapse = ", "), "\n", sep = "")
  quit(save = "no", status = 0L)
}

install_deps(to_install)

# Report against a fresh look at the library rather than assuming success: a
# Bioconductor build can fail on a missing system header and still leave the
# rest of the batch installed.
after <- deps_status()
still <- after$package[!after$installed & after$package %in% to_install]
cat("\n")
rule()
if (length(still)) {
  cat("FAILED to install: ", paste(still, collapse = ", "), "\n", sep = "")
  cat("Scroll up for the compiler or download error. A missing system library is\n")
  cat("the usual cause on Linux; see the System section above.\n")
  quit(save = "no", status = 1L)
}
cat("Installed: ", paste(to_install, collapse = ", "), "\n", sep = "")
if (in_renv) {
  fresh <- setdiff(to_install, locked_packages())
  if (length(fresh)) {
    cat("\n", length(fresh), " of these are not in renv.lock (",
        paste(fresh, collapse = ", "), ").\n", sep = "")
    cat("Record the versions so the other machine gets the same ones:\n")
    cat("    Rscript -e 'renv::snapshot()'\n")
  }
}
