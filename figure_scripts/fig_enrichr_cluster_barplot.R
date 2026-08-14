#!/usr/bin/env Rscript
# =============================================================================
# Enrichr pathway enrichment — faceted bar chart (phenotype-agnostic)
# Reads Enrichr result files from a directory and produces a faceted OR bar
# chart (one facet per database), labelled with significance stars.
#
# Output : <outdir>/<outfile>.{png,svg}
#          <outdir>/enrichr_table_top_terms.txt   (data behind the plot)
#
# Expected files inside --enrichr-dir follow the pipeline's prefix-based naming
# "<prefix>_<library>.txt", selected by --prefix (the cluster id that precedes
# the library name, e.g. "cluster1"):
#   <prefix>_GO_Biological_Process_*   <prefix>_GO_Molecular_Function_*
#   <prefix>_GO_Cellular_Component_*   <prefix>_BioPlanet_*
#   <prefix>_BioCarta_*                <prefix>_Reactome_*   <prefix>_KEGG_*
#
# Nothing is phenotype-specific: the term-shortening map is supplied via
# --label-map (optional) instead of being hardcoded. The database panel and
# per-database term counts live in the editable db_specs list below.
#
# Usage:
#   Rscript fig_enrichr_cluster_barplot.R \
#       --enrichr-dir results/intermediate/chronic_lymphocytic_leukemia/enrichr \
#       --prefix cluster1 \
#       --outdir figures/chronic_lymphocytic_leukemia
#
#   # with a curated label map:
#   Rscript fig_enrichr_cluster_barplot.R \
#       --enrichr-dir results/intermediate/leptin/enrichr \
#       --prefix cluster1 \
#       --outdir figures/leptin \
#       --label-map data/enrichr/leptin_cluster_label_map.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

option_list <- list(
  make_option(c("-d", "--enrichr-dir"), type = "character",
              help = "Directory containing Enrichr result files [required]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory [required]"),
  make_option(c("--outfile"), type = "character", default = NULL,
              help = "Output base filename (no extension) [default: enrichr_<prefix>_barplot]"),
  make_option(c("--prefix"), type = "character",
              help = "Cluster prefix preceding the library name in each filename, e.g. 'cluster1' (matches <prefix>_<library>.txt) [required]"),
  make_option(c("--top-n"), type = "integer", default = 3L,
              help = "Top N significant terms per database to plot [default: %default]"),
  make_option(c("--exclude-dbs"), type = "character",
              default = "GO Cellular Component,BioCarta 2016",
              help = "Comma-separated database display names to omit from the plot (kept in the table) [default: %default]"),
  make_option(c("--label-map"), type = "character", default = NULL,
              help = "Optional TSV (columns: term, label) to shorten term labels; use literal \\n in 'label' for line breaks"),
  make_option(c("--width"), type = "double", default = 16.0,
              help = "Figure width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 4.5,
              help = "Figure height in inches [default: %default]"),
  make_option(c("--dpi"), type = "integer", default = 300L,
              help = "PNG resolution [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Enrichr faceted OR bar chart")
opt <- parse_args(parser)

enrichr_dir <- opt[["enrichr-dir"]]
if (is.null(enrichr_dir)) {
  cat("ERROR: --enrichr-dir is required\n\n", file = stderr()); print_help(parser); quit(status = 2)
}
if (!dir.exists(enrichr_dir)) {
  cat("ERROR: Directory not found:", enrichr_dir, "\n", file = stderr()); quit(status = 2)
}
if (is.null(opt$outdir)) {
  cat("ERROR: --outdir is required\n\n", file = stderr()); print_help(parser); quit(status = 2)
}
if (is.null(opt$prefix)) {
  cat("ERROR: --prefix is required\n\n", file = stderr()); print_help(parser); quit(status = 2)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

source("R/palettes.R")

# Colours, resolved once. Each fallback is the literal this script used before
# it took a --palette, so the default rendering is unchanged.
pal        <- get_palette(opt$palette)
col_bar    <- pal_col(pal, "accent",  "darkgrey")  # OR bars
col_ink    <- pal_col(pal, "ink",     "black")     # outlines, axis text, value labels
col_muted  <- pal_col(pal, "muted",   "grey30")    # secondary annotation
col_grid   <- pal_col(pal, "grid",    "grey85")
col_strip  <- pal_col(pal, "ink",     "grey20")    # facet strip fill
col_striptx<- pal_col(pal, "paper",   "white")     # facet strip text, on that fill
col_paper  <- pal_col(pal, "paper",   "white")

cluster_prefix <- opt$prefix
if (is.null(opt$outfile)) opt$outfile <- paste0("enrichr_", cluster_prefix, "_barplot")
top_n       <- opt[["top-n"]]
exclude_dbs <- trimws(strsplit(opt[["exclude-dbs"]], ",")[[1]])
exclude_dbs <- exclude_dbs[nzchar(exclude_dbs)]

# ── optional term-label overrides ─────────────────────────────
label_map <- character(0)
if (!is.null(opt[["label-map"]])) {
  lm_path <- opt[["label-map"]]
  if (!file.exists(lm_path)) {
    cat("ERROR: --label-map file not found:", lm_path, "\n", file = stderr()); quit(status = 2)
  }
  lm <- read.delim(lm_path, sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("term", "label") %in% names(lm))) {
    cat("ERROR: --label-map must have columns 'term' and 'label'\n", file = stderr()); quit(status = 2)
  }
  # Literal "\n" in the file becomes a real newline for multi-line labels.
  label_map <- setNames(gsub("\\\\n", "\n", lm$label), lm$term)
  message("Loaded ", length(label_map), " label override(s) from ", lm_path)
}

# =============================================================================
# DATABASE PANEL — edit here to change which DBs appear / how many terms each
#   pattern : library token matched after "<prefix>_" inside --enrichr-dir
#   name    : display name used in the facet strip and the table
#   n       : number of top terms (by Combined Score) to keep for this DB
# =============================================================================
db_specs <- list(
  list(pattern = "GO_Biological_Process", name = "GO Biological Process", n = top_n),
  list(pattern = "GO_Molecular_Function", name = "GO Molecular Function", n = top_n),
  list(pattern = "GO_Cellular_Component", name = "GO Cellular Component", n = min(top_n, 2L)),
  list(pattern = "BioPlanet",             name = "BioPlanet 2019",        n = top_n),
  list(pattern = "BioCarta",              name = "BioCarta 2016",         n = 1L),
  list(pattern = "Reactome",              name = "Reactome 2024",         n = top_n),
  list(pattern = "KEGG",                  name = "KEGG 2026",             n = top_n)
)

# =============================================================================
# LOAD FILES
# =============================================================================

# Find the Enrichr file for `db_pattern` belonging to cluster `prefix`, i.e. a
# basename "<prefix>_<library>.txt" whose library begins with db_pattern.
# Anchoring on "^<prefix>_<db_pattern>" stops "cluster1" from also grabbing
# "cluster1_plus30_..." (the stepping files share the same cluster stem).
find_file <- function(prefix, db_pattern) {
  all_files <- list.files(enrichr_dir, full.names = TRUE)
  bn   <- basename(all_files)
  rx   <- paste0("^", prefix, "_", db_pattern)
  hits <- all_files[grepl(rx, bn) & endsWith(bn, ".txt")]
  if (length(hits) == 0) {
    cat("ERROR: no file matching '", rx, "' (in ", enrichr_dir, ")\n",
        sep = "", file = stderr())
    quit(status = 2)
  }
  hits[1]
}

top_terms <- function(df, db_name, n) {
  df %>%
    filter(Adjusted.P.value < 0.05) %>%
    arrange(desc(Combined.Score)) %>%
    slice_head(n = n) %>%
    mutate(
      Database = db_name,
      Term     = gsub("\\s*\\(GO:\\d+\\)", "", Term),
      FDR      = Adjusted.P.value,
      OR       = round(Odds.Ratio, 2),
      CS       = round(Combined.Score, 1),
      N_sig    = sum(df$Adjusted.P.value < 0.05)
    ) %>%
    select(Database, Term, OR, FDR, CS, N_sig)
}

df <- bind_rows(lapply(db_specs, function(spec) {
  f <- find_file(cluster_prefix, spec$pattern)
  top_terms(read.delim(f), spec$name, spec$n)
}))

df$Term <- str_to_sentence(df$Term)

table_path <- file.path(opt$outdir, "enrichr_table_top_terms.txt")
write.table(df, table_path, row.names = FALSE, sep = "\t")
message("Table saved -> ", table_path)

# Exclude requested databases from the plot (they stay in the table above)
df <- df[!(df$Database %in% exclude_dbs), ]
if (nrow(df) == 0) {
  cat("ERROR: no terms left to plot after exclusions / filtering.\n", file = stderr())
  quit(status = 2)
}

df <- df %>%
  mutate(stars = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    TRUE        ~ ""
  ))

# Facet order follows db_specs, minus the excluded databases
db_order <- setdiff(vapply(db_specs, `[[`, character(1), "name"), exclude_dbs)
df$Database <- factor(df$Database, levels = db_order)

df <- df %>% arrange(Database, desc(OR))

# Display labels: --label-map override (raw term -> label) else wrap at width 20
make_label <- function(t) {
  if (t %in% names(label_map)) label_map[[t]]
  else paste(strwrap(t, width = 20), collapse = "\n")
}
df$TermLabel <- vapply(df$Term, make_label, character(1))
df$TermLabel <- factor(df$TermLabel, levels = unique(df$TermLabel))

or_max <- max(df$OR)
df <- df %>%
  mutate(or_label_y   = OR + or_max * 0.06,
         star_label_y = OR + or_max * 0.15)

# =============================================================================
# PLOT
# =============================================================================

p <- ggplot(df, aes(x = TermLabel, y = OR)) +
  geom_col(fill = col_bar, width = 0.7, colour = col_ink, linewidth = 0.4) +
  geom_text(aes(y = or_label_y, label = round(OR, 0)),
            size = 5, colour = col_muted, vjust = 0) +
  geom_text(aes(y = star_label_y, label = stars),
            size = 6, colour = col_ink, vjust = 0, fontface = "bold") +
  facet_wrap(~ Database, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Odds Ratio (OR)") +
  theme_classic(base_size = 12) +
  theme(
    strip.background    = element_rect(fill = col_strip, colour = col_ink, linewidth = 0.5),
    strip.text          = element_text(size = 12, face = "bold", colour = col_striptx),
    axis.text.x         = element_text(size = 12, angle = 40, hjust = 1, colour = col_ink),
    axis.text.y         = element_text(size = 12, colour = col_ink),
    axis.title.y        = element_text(size = 12, colour = col_ink, margin = margin(r = 6)),
    axis.ticks          = element_line(colour = col_ink, linewidth = 0.4),
    axis.ticks.length   = unit(2, "pt"),
    panel.grid.major.y  = element_line(colour = col_grid, linewidth = 0.3),
    panel.grid.minor    = element_blank(),
    panel.spacing       = unit(0.4, "lines"),
    panel.border        = element_rect(colour = col_ink, fill = NA, linewidth = 0.5),
    plot.margin         = margin(10, 15, 10, 10)
  )

# =============================================================================
# EXPORT
# =============================================================================

base <- file.path(opt$outdir, opt$outfile)

ggsave(paste0(base, ".png"), plot = p,
       width = opt$width, height = opt$height, dpi = opt$dpi, bg = col_paper)
message("PNG saved -> ", paste0(base, ".png"))

ggsave(paste0(base, ".svg"), plot = p, device = "svg",
       width = opt$width, height = opt$height, fix_text_size = FALSE)
message("SVG saved -> ", paste0(base, ".svg"))
