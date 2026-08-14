# =============================================================================
# .Rprofile — project startup, sourced by every R session started from the
# repository root (including every Rscript the pipeline runs).
#
# Two jobs, both of them about running the same repo on two machines that
# resolve packages in completely different ways.
#
# 1. renv, unless told otherwise.
#    Inside the pleio container every package is already in the image's site
#    library, and activating renv there would put an empty project library in
#    front of it and hide the lot. Set PLEIO_NO_RENV=1 in the container (or for
#    any one-off session) and R uses the library it was given:
#        PLEIO_NO_RENV=1 Rscript figure_scripts/fig_string_chord_by_block.R config/cll.yml
#    Everywhere else renv activates as before, so the laptop keeps the pinned
#    library under renv/library/ and the lockfile still governs versions.
#
# 2. A missing package should print its own fix.
#    `Error in library(circlize) : there is no package called 'circlize'` says
#    nothing about where the packages for this project come from. The handler
#    below appends that, and only for that error — everything else keeps R's
#    ordinary behaviour and its exit status.
# =============================================================================

if (!identical(Sys.getenv("PLEIO_NO_RENV"), "1") && file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

local({

  # Pull the package name out of R's message. R quotes it with directional
  # marks in a UTF-8 locale and with ASCII quotes otherwise, so accept both.
  missing_package <- function(msg) {
    m <- regmatches(msg, regexec(
      "there is no package called [‘'\"]?([^’'\"\n]+)", msg))[[1]]
    if (length(m) == 2L) trimws(m[2]) else NA_character_
  }

  hint <- function() {
    pkg <- tryCatch(missing_package(geterrmessage()), error = function(e) NA_character_)
    if (is.na(pkg)) return(invisible(NULL))
    cat(sep = "", file = stderr(),
        "\n", strrep("-", 72), "\n",
        "  '", pkg, "' is not installed in ", .libPaths()[1], "\n\n",
        "  What this machine is missing, in full:\n",
        "      Rscript tools/setup_env.R --check\n",
        "  Install it:\n",
        "      Rscript tools/setup_env.R\n",
        "  In the pleio container, use the container's own library instead:\n",
        "      PLEIO_NO_RENV=1 Rscript <script> <config>\n",
        strrep("-", 72), "\n")
  }

  options(error = function() {
    try(hint(), silent = TRUE)
    # The handler replaces R's default action, so restore the half of it that
    # matters: a non-interactive session must still stop, and still exit
    # non-zero, or run_pipeline.sh's `set -e` would sail past a failed stage.
    if (!interactive()) quit(save = "no", status = 1L)
  })

})
