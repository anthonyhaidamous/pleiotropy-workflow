# =============================================================================
# R/palettes.R
# Named colour palettes for the figure scripts, so the same figure can be
# rendered for different destinations — a journal article, a conference poster,
# a slide deck — without editing the script.
#
#     source("R/palettes.R")         # from the project root
#
# Scope and rules
# ---------------
#   - Base R only. Loads NO packages and has NO side effects, matching
#     R/figure_utils.R, so it is safe to source from any script. The one
#     exception is reading a custom palette from YAML, which requires the
#     `yaml` package lazily, only when a file path is actually passed.
#   - Opt-in by construction. With no palette selected every accessor returns
#     the fallback the calling script passes in — that fallback is the literal
#     the script already used, so existing output is unchanged byte for byte.
#     A figure only changes colour when someone asks for a palette.
#   - Partial palettes are legal. A palette (especially a user's own) may
#     define any subset of the roles; unset roles fall back per above.
#
# Using it from a figure script
# -----------------------------
# Add the option, resolve once, then read colours through the accessors:
#
#     make_option("--palette", type = "character", default = NULL,
#                 help = palette_help())
#
#     pal     <- get_palette(opt$palette)
#     col_ink <- pal_col(pal, "ink",    "#2c3e50")   # <- the script's own literal
#     col_pt  <- pal_col(pal, "accent", "steelblue")
#     cols    <- pal_cat(pal, n_groups, tableau10)
#
# Selecting a palette
# -------------------
#     --palette poster                 a registered palette
#     --palette my_colours.yml         your own, as role: "#hex" pairs
#     FIGURE_PALETTE=poster Rscript …  applies to any script that has no
#                                      explicit --palette (handy for rendering
#                                      a whole figure set for one destination)
#
# Roles
# -----
#   paper    plot and panel background
#   ink      primary text, axis labels, dark marks
#   muted    secondary text, captions, annotations
#   grid     panel gridlines
#   axis     axis lines, ticks, panel borders
#   accent   the primary mark colour for a single-series figure
#   accent2  a second mark colour
#   pos/neg  the two poles of a signed scale (enriched/depleted,
#            concordant/discordant, up/down)
#   neutral  non-significant, zero, "other"
#   lo/mid/hi  endpoints and midpoint of a continuous scale
#   cat      categorical vector, for modules / categories / blocks
# =============================================================================


PALETTE_ROLES <- c("paper", "ink", "muted", "grid", "axis", "accent", "accent2",
                   "pos", "neg", "neutral", "lo", "mid", "hi", "cat")


# -----------------------------------------------------------------------------
# The registry
#
# Colour choices are perceptual, not eyeballed. Within a palette, `pos` and
# `neg` are matched in CIE L* and chroma so neither pole reads as louder than
# the other, and each `cat` vector walks a constant-L*, constant-chroma circle
# so no category dominates. See the note on each palette.
# -----------------------------------------------------------------------------
.FIGURE_PALETTES <- list(

  # The JOBIM 2026 poster (poster/index.html). Cream paper, near-black navy
  # ink, crimson accent. The poster itself is monochrome-plus-crimson and has
  # no second pole, so `neg` is derived: the crimson's exact L* (29.2) and
  # chroma (38.6) rotated to hue 225 deg.
  poster = list(
    paper   = "#F7F6F3",
    ink     = "#1B2330",
    muted   = "#4B505A",
    grid    = "#D8D4CB",
    axis    = "#9AA0AB",
    accent  = "#7A2738",
    accent2 = "#00516F",
    pos     = "#7A2738",
    neg     = "#00516F",
    neutral = "#A39D8D",
    lo      = "#F7F6F3",
    mid     = "#EFECE5",
    hi      = "#7A2738",
    cat     = c("#A85B65", "#A2614B", "#906B3A", "#767539", "#577C47", "#2E8061",
                "#00817F", "#007F99", "#1C79A9", "#5C70A9", "#86659B", "#9F5C82")
  ),

  # Print journal. White paper, near-black ink, and the Okabe-Ito qualitative
  # set, which stays distinguishable under all common colour-vision
  # deficiencies and survives greyscale conversion.
  article = list(
    paper   = "#FFFFFF",
    ink     = "#1A1A1A",
    muted   = "#4D4D4D",
    grid    = "#E6E6E6",
    axis    = "#808080",
    accent  = "#0072B2",
    accent2 = "#D55E00",
    pos     = "#D55E00",
    neg     = "#0072B2",
    neutral = "#999999",
    lo      = "#F7FBFF",
    mid     = "#F0F0F0",
    hi      = "#08306B",
    cat     = c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#E69F00", "#56B4E9", "#F0E442", "#000000")
  ),

  # Slides. Same structure, pushed brighter and higher-contrast so it survives
  # a projector and a lit room.
  ppt = list(
    paper   = "#FFFFFF",
    ink     = "#12232E",
    muted   = "#41576B",
    grid    = "#DFE6EC",
    axis    = "#8FA3B3",
    accent  = "#0F6FC6",
    accent2 = "#F2760B",
    pos     = "#E8453C",
    neg     = "#1F7AE0",
    neutral = "#9AA7B0",
    lo      = "#EAF2FB",
    mid     = "#F2F2F2",
    hi      = "#0F4C81",
    cat     = c("#DD6568", "#C5793A", "#988D2A", "#5B9B48", "#00A17D",
                "#00A1B5", "#009ADD", "#598AE5", "#B173C8", "#D96395")
  )
)


# -----------------------------------------------------------------------------
# Names of the registered palettes
# -----------------------------------------------------------------------------
list_palettes <- function() names(.FIGURE_PALETTES)


# -----------------------------------------------------------------------------
# Help text for the --palette option, so every script documents the same thing
# -----------------------------------------------------------------------------
palette_help <- function() {
  paste0("Colour palette: ", paste(list_palettes(), collapse = " | "),
         ", or a path to a YAML file of role: \"#hex\" pairs. ",
         "Defaults to the script's own colours, or $FIGURE_PALETTE if set.")
}


# -----------------------------------------------------------------------------
# Resolve a palette selection into a list of roles, or NULL
#
# NULL means "no palette selected" — every accessor then yields the calling
# script's own fallback, which is how default output stays unchanged.
#
# `name` may be a registered name, a path to a YAML file, or NULL/""/"default"
# for none. When it is NULL the FIGURE_PALETTE environment variable is
# consulted, so a whole figure set can be rendered for one destination without
# touching any command line.
# -----------------------------------------------------------------------------
get_palette <- function(name = NULL) {
  if (is.null(name) || !nzchar(name)) name <- Sys.getenv("FIGURE_PALETTE", "")
  if (!nzchar(name) || identical(name, "default") || identical(name, "none")) {
    return(NULL)
  }

  if (name %in% names(.FIGURE_PALETTES)) return(.FIGURE_PALETTES[[name]])

  if (grepl("\\.(ya?ml)$", name, ignore.case = TRUE)) {
    if (!file.exists(name)) {
      stop("Palette file not found: ", name, call. = FALSE)
    }
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Reading a palette from YAML needs the 'yaml' package.", call. = FALSE)
    }
    pal <- yaml::yaml.load_file(name)
    if (!is.list(pal)) {
      stop("Palette file must be a mapping of role: value — ", name, call. = FALSE)
    }
    unknown <- setdiff(names(pal), PALETTE_ROLES)
    if (length(unknown)) {
      warning("Ignoring unknown palette role(s) in ", name, ": ",
              paste(unknown, collapse = ", "), call. = FALSE)
      pal <- pal[intersect(names(pal), PALETTE_ROLES)]
    }
    # cat may arrive as a YAML list; flatten to a character vector
    if (!is.null(pal$cat)) pal$cat <- as.character(unlist(pal$cat, use.names = FALSE))
    return(pal)
  }

  stop("Unknown palette '", name, "'. Available: ",
       paste(list_palettes(), collapse = ", "),
       ", or a path to a .yml file.", call. = FALSE)
}


# -----------------------------------------------------------------------------
# One colour, by role
#
# `fallback` is what the script used before it knew about palettes, and is what
# it still gets when no palette is active or the active palette leaves the role
# unset. It is required — a role with no fallback would make the script's
# default rendering depend on the registry.
# -----------------------------------------------------------------------------
pal_col <- function(pal, role, fallback) {
  if (is.null(pal) || is.null(pal[[role]])) return(fallback)
  as.character(pal[[role]])[1]
}


# -----------------------------------------------------------------------------
# `n` categorical colours
#
# Returns the palette's `cat` vector trimmed to n. When n exceeds what the
# palette carries, the vector is interpolated to length — the in-between shades
# are less distinct than the anchors, which is a fair warning that the figure
# is asking colour to carry more categories than it can.
# -----------------------------------------------------------------------------
pal_cat <- function(pal, n, fallback) {
  cols <- if (is.null(pal) || is.null(pal$cat)) fallback else as.character(pal$cat)
  if (length(cols) == 0L) return(fallback)
  if (n <= length(cols)) return(cols[seq_len(n)])
  grDevices::colorRampPalette(cols)(n)
}
