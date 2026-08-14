#!/usr/bin/env Rscript
# =============================================================================
# Table 1 — GWAS studies for a phenotype of interest
# Pipeline-adapted: sources R/utils.R and applies the same study-level filters
# as pipeline/06_filter (ethnicity, CNV, sample size, GxE), instead of the old
# inline ancestry map / ad-hoc sample-size parsing.
#
# Inputs : a GWAS Catalog studies TSV (data/studies/*.tsv) + the studied
#          trait term(s) (MAPPED_TRAIT values that define the phenotype).
#
# Filters (mirroring pipeline/06_filter, in the same order):
#   (1) Ethnicity   — ancestry encoded via R/utils.R encode_ancestries_for_rows.
#                     Annotate-only: there is no LD block to be concordant with
#                     at the study level, so no rows are dropped.
#   (2) Sample size — initial + replication >= --min-sample-size (default 1000,
#                     matching cfg$filters$min_sample_size).
#   (3) CNV         — skipped with a warning: the studies TSV has no CNV column
#                     (CNV is an association-level field in 06_filter).
#   (4) GxE         — two-tier regex (default_gxe_patterns) on DISEASE/TRAIT,
#                     MAPPED_TRAIT and STUDY. auto_exclude rows are dropped;
#                     multi-mapped rows are kept provisionally (as in 06's
#                     first run).
#
# Style  : Nature/Cell 3-line booktabs, black & white, no row shading
# Output : <outdir>/table1_studies.tsv   (data)
#          <outdir>/table1_studies.svg   (publication figure)
#
# Run from the project root, e.g.:
#   Rscript figure_scripts/table1_gwas_studies.R \
#       -i data/studies/gwas-catalog-v1.0.3.1-studies-r2026-02-16.tsv \
#       --traits "leptin measurement,BMI-adjusted leptin measurement" \
#       --trait-labels "Leptin,BMI-adj. leptin" \
#       -o figures/leptin
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(svglite)
  library(grid)
})

# Locate and source the shared pipeline helpers (encode_ancestries_for_rows,
# sum_sample_sizes, default_gxe_patterns) whether run from the project root or
# from the figure_scripts/ directory.
get_script_dir <- function() {
  ca <- commandArgs(FALSE)
  m  <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  getwd()
}
utils_candidates <- c("R/utils.R",
                      file.path(get_script_dir(), "..", "..", "R", "utils.R"))
utils_path <- utils_candidates[file.exists(utils_candidates)][1]
if (is.na(utils_path)) {
  cat("ERROR: cannot locate R/utils.R (run from the project root).\n",
      file = stderr())
  quit(status = 2)
}
source(utils_path)

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Path to GWAS Catalog studies TSV [required]"),
  make_option(c("--traits"), type = "character",
              help = "Comma-separated MAPPED_TRAIT term(s) defining the studied trait [required]"),
  make_option(c("--trait-labels"), type = "character", default = NULL,
              help = "Optional comma-separated display labels, parallel to --traits"),
  make_option(c("--min-sample-size"), type = "integer", default = 1000L,
              help = "Sample-size filter: initial + replication >= this [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory required"),
  make_option(c("--svg-width"), type = "double", default = 24.0,
              help = "SVG canvas width in cm [default: %default]"),
  make_option(c("--row-height"), type = "double", default = 0.52,
              help = "Row height in cm; raise to 0.65 if text clips [default: %default]"),
  make_option(c("--palette"), type = "character", default = NULL,
              help = "Colour palette: article | poster | ppt, or a path to a YAML file [default: the script's own colours]")
)

parser <- OptionParser(option_list = option_list,
                       description = "Table 1: GWAS studies summary (pipeline-filtered)")
opt <- parse_args(parser)

if (is.null(opt$input)) {
  cat("ERROR: --input is required\n\n", file = stderr())
  print_help(parser)
  quit(status = 2)
}
if (is.null(opt$traits)) {
  cat("ERROR: --traits is required (the MAPPED_TRAIT term(s) of the studied trait)\n\n",
      file = stderr())
  print_help(parser)
  quit(status = 2)
}
if (is.null(opt$outdir)) {
  cat("ERROR: --outdir is required\n\n", file = stderr())
  print_help(parser)
  quit(status = 2)
}
if (!file.exists(opt$input)) {
  cat("ERROR: File not found:", opt$input, "\n", file = stderr())
  quit(status = 2)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

trait_terms     <- trimws(strsplit(opt$traits, ",")[[1]])
min_sample_size <- opt[["min-sample-size"]]

trait_labels_map <- NULL
if (!is.null(opt[["trait-labels"]])) {
  labs <- trimws(strsplit(opt[["trait-labels"]], ",")[[1]])
  if (length(labs) != length(trait_terms)) {
    cat("ERROR: --trait-labels must have the same number of comma-separated",
        "entries as --traits\n", file = stderr())
    quit(status = 2)
  }
  trait_labels_map <- setNames(labs, trait_terms)
}
message("Studied trait term(s): ", paste(trait_terms, collapse = " | "))

# =============================================================================
# 1. DATA PIPELINE
# =============================================================================

studies <- read.delim(opt$input, sep = "\t", quote = "",
                      check.names = FALSE, stringsAsFactors = FALSE)
studies <- subset(studies, MAPPED_TRAIT %in% trait_terms)
if (nrow(studies) == 0) {
  cat("ERROR: no studies match the given trait term(s) in", opt$input, "\n",
      file = stderr())
  quit(status = 2)
}

n_strata_total <- nrow(studies)
n_pubmed_total <- n_distinct(studies$`PUBMED ID`)

# Annotations needed for both filtering and display. Sample sizes and ancestry
# are computed with the canonical pipeline helpers (R/utils.R), so they stay in
# lockstep with pipeline/01 and pipeline/06 — in particular this drops the old
# inline ancestry map that still misclassified Thai/Filipino/Zulu.
studies$YEAR          <- sapply(strsplit(studies$DATE, "-"), `[[`, 1)
studies$N_discovery   <- sum_sample_sizes(studies$`INITIAL SAMPLE SIZE`)
studies$N_replication <- sum_sample_sizes(studies$`REPLICATION SAMPLE SIZE`)
studies$N_total       <- studies$N_discovery + studies$N_replication
studies$Ancestry      <- encode_ancestries_for_rows(studies)

# -----------------------------------------------------------------------------
# Filtering — mirrors pipeline/06_filter, same order, with a per-filter log.
# -----------------------------------------------------------------------------
initial_n <- nrow(studies)
log_rows <- list(
  list(filter = "0_input", n_before = NA_integer_, n_after = initial_n,
       n_dropped = 0L)
)
message("Starting from ", initial_n, " study strata.")

# (1) Ethnicity — annotate only. At the study level there is no LD block for the
#     trait's own studies to be concordant with, so ancestry is encoded for
#     display but no rows are dropped.
message("(1) Ethnicity: ancestry encoded (annotate-only, no rows dropped).")

# (2) Sample size: initial + replication >= min_sample_size
n_before <- nrow(studies)
studies  <- studies[studies$N_total >= min_sample_size, ]
n_after  <- nrow(studies)
message("(2) Sample size ≥ ", min_sample_size, ": kept ", n_after, " / ",
        n_before, " (dropped ", n_before - n_after, ").")
log_rows[[length(log_rows) + 1]] <- list(
  filter = "2_sample_size", n_before = n_before, n_after = n_after,
  n_dropped = n_before - n_after
)

# (3) CNV: not available at the study level (CNV is an association-level field
#     in 06_filter), so this filter is skipped.
warning("CNV filter skipped: the studies TSV has no CNV column ",
        "(CNV is an association-level field in 06_filter).")
message("(3) CNV: skipped (not available at study level).")

# (4) GxE / multi-trait: same two-tier regex as 06_filter.
if (nrow(studies) > 0) {
  patterns <- default_gxe_patterns()

  na0 <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  disease_trait <- na0(studies$`DISEASE/TRAIT`)
  mapped_trait  <- na0(studies$MAPPED_TRAIT)
  study_txt     <- na0(studies$STUDY)

  flag_matrix <- vapply(patterns, function(p) {
    grepl(p, disease_trait, ignore.case = TRUE, perl = TRUE) |
    grepl(p, mapped_trait,  ignore.case = TRUE, perl = TRUE) |
    grepl(p, study_txt,     ignore.case = TRUE, perl = TRUE)
  }, logical(nrow(studies)))

  flag_auto <- if (is.matrix(flag_matrix)) apply(flag_matrix, 1, any) else any(flag_matrix)
  flag_manual <- grepl(",", mapped_trait) & !flag_auto

  studies$gxe_decision <- ifelse(flag_auto, "auto_exclude",
                          ifelse(flag_manual, "manual_check", "keep"))

  message("(4) GxE / multi-trait classification:")
  message("       auto_exclude : ", sum(studies$gxe_decision == "auto_exclude"))
  message("       manual_check : ", sum(studies$gxe_decision == "manual_check"),
          " (kept provisionally, as in 06's first run)")
  message("       keep         : ", sum(studies$gxe_decision == "keep"))

  # Like 06's first run, multi-mapped rows are provisionally kept.
  studies$gxe_decision[studies$gxe_decision == "manual_check"] <- "keep"

  n_before <- nrow(studies)
  studies  <- studies[studies$gxe_decision == "keep", ]
  n_after  <- nrow(studies)
  message("(4) GxE / multi-trait: kept ", n_after, " / ", n_before,
          " (dropped ", n_before - n_after, ").")
  log_rows[[length(log_rows) + 1]] <- list(
    filter = "4_gxe", n_before = n_before, n_after = n_after,
    n_dropped = n_before - n_after
  )
}

# -----------------------------------------------------------------------------
# Display columns + phenotype labels
# -----------------------------------------------------------------------------
studies <- studies %>%
  mutate(`Array type` = str_extract(`GENOTYPING TECHNOLOGY`, "^[^\\[]+") %>% str_trim()) %>%
  rename(Associations = `ASSOCIATION COUNT`)

studies$Phenotype <- studies$MAPPED_TRAIT
if (!is.null(trait_labels_map)) {
  studies$Phenotype <- ifelse(studies$Phenotype %in% names(trait_labels_map),
                              trait_labels_map[studies$Phenotype],
                              studies$Phenotype)
}

cat(sprintf(
  "Strata: %d -> %d kept after filters | PubMed IDs (pre-filter): %d | 0-hit strata among kept: %d\n",
  n_strata_total, nrow(studies), n_pubmed_total, sum(studies$Associations == 0)
))

pheno_levels <- unique(studies$Phenotype)
studies_final <- studies %>%
  mutate(
    `Author (Year)` = paste0(`FIRST AUTHOR`, " et al. (", YEAR, ")"),
    Phenotype = factor(Phenotype, levels = pheno_levels)
  ) %>%
  filter(Associations > 0) %>%
  group_by(`PUBMED ID`) %>%
  summarize(
    `Author (Year)`       = first(`Author (Year)`),
    Phenotype             = paste(levels(Phenotype)[levels(Phenotype) %in% unique(Phenotype)],
                                  collapse = ", "),
    `N (total)`           = format(max(N_total), big.mark = ",", trim = TRUE),
    Ancestry              = paste(unique(unlist(strsplit(Ancestry, ",\\s*"))), collapse = ", "),
    `Genotyping platform` = paste(unique(`Array type`), collapse = ", "),
    `GWAS hits`           = sum(Associations),
    .groups = "drop"
  ) %>%
  select(`Author (Year)`, Phenotype, `GWAS hits`, `N (total)`, Ancestry, `Genotyping platform`)

# -----------------------------------------------------------------------------
# Filtering summary (06-style)
# -----------------------------------------------------------------------------
bar <- paste(rep("─", 60), collapse = "")
message("\n", bar)
message("Study-level filtering summary")
message(bar)
message(sprintf("%-25s  %10s  %10s  %10s", "filter", "in", "out", "dropped"))
for (r in log_rows[-1]) {
  message(sprintf("%-25s  %10d  %10d  %10d",
                  r$filter, r$n_before, r$n_after, r$n_dropped))
}
message(bar)
message(sprintf("Initial -> final : %d -> %d (%d dropped)",
                initial_n, nrow(studies), initial_n - nrow(studies)))
message(bar)

# =============================================================================
# 2. TSV EXPORT
# =============================================================================

tsv_path <- file.path(opt$outdir, "table1_studies.tsv")
write.table(studies_final, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("TSV saved -> ", tsv_path)

# =============================================================================
# 3. SVG FIGURE
# =============================================================================

FONT      <- "Arial"
FS_HEADER <- 8
FS_BODY   <- 8
LWD_THICK <- 1.5
LWD_THIN  <- 0.75

# Colours, resolved once. Each fallback is the literal this script used
# before it took a --palette, so the default rendering is unchanged.
source("R/palettes.R")
pal       <- get_palette(opt$palette)
COL_LINE  <- pal_col(pal, "axis", "black")
COL_TEXT  <- pal_col(pal, "ink",  "black")

svg_width     <- opt[["svg-width"]]
row_h         <- opt[["row-height"]]
header_h      <- 0.62
margin_top    <- 0.12
margin_bottom <- 0.20
pad_l         <- 0.28
pad_r         <- 0.28

n_rows <- nrow(studies_final)
n_cols <- ncol(studies_final)
W      <- svg_width
H      <- margin_top + header_h + n_rows * row_h + margin_bottom

# Column proportions (must sum to 1):
#   Author (Year)       27 chars
#   Phenotype           23 chars
#   GWAS hits            9 chars [right]
#   N (total)            8 chars [right]
#   Ancestry            18 chars
#   Genotyping platform 28 chars
col_w <- c(0.24, 0.20, 0.08, 0.09, 0.16, 0.23)
stopifnot("col_w must sum to 1" = abs(sum(col_w) - 1) < 1e-9)

right_cols <- c(3L, 4L)

x_left       <- cumsum(c(0, col_w[-n_cols]))
x_right      <- x_left + col_w
y_top_rule   <- 1 - margin_top / H
y_midrule    <- y_top_rule - header_h / H
y_header_mid <- (y_top_rule + y_midrule) / 2
y_bot_rule   <- margin_bottom / H

svg_path <- file.path(opt$outdir, "table1_studies.svg")
svglite(svg_path, width = W / 2.54, height = H / 2.54,
        system_fonts = list(sans = FONT), fix_text_size = FALSE)
grid.newpage()

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_top_rule,y_top_rule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

for (j in seq_len(n_cols)) {
  is_r  <- j %in% right_cols
  x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
  grid.text(label = names(studies_final)[j],
            x = unit(x_txt,"npc"), y = unit(y_header_mid,"npc"),
            hjust = if (is_r) 1 else 0, vjust = 0.5,
            gp = gpar(fontfamily = FONT, fontface = "bold",
                      fontsize = FS_HEADER, col = COL_TEXT))
}

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_midrule,y_midrule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THIN))

for (i in seq_len(n_rows)) {
  y_mid <- y_midrule - ((i - 0.5) * row_h) / H
  for (j in seq_len(n_cols)) {
    is_r  <- j %in% right_cols
    x_txt <- if (is_r) (x_right[j] - pad_r / W) else (x_left[j] + pad_l / W)
    grid.text(label = as.character(studies_final[[j]][i]),
              x = unit(x_txt,"npc"), y = unit(y_mid,"npc"),
              hjust = if (is_r) 1 else 0, vjust = 0.5,
              gp = gpar(fontfamily = FONT, fontface = "plain",
                        fontsize = FS_BODY, col = COL_TEXT))
  }
}

grid.lines(x = unit(c(0,1),"npc"), y = unit(c(y_bot_rule,y_bot_rule),"npc"),
           gp = gpar(col = COL_LINE, lwd = LWD_THICK))

dev.off()
message("SVG saved -> ", svg_path)
