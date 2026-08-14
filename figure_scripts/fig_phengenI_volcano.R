#!/usr/bin/env Rscript
# ============================================================
# PhenGenI Association 2021 — Volcano Plot (CLI)
# Phenotype-agnostic: panel labels are derived from each input file's name
# (or supplied explicitly); nothing trait-specific is hardcoded.
#
# Usage:
#   Rscript figure_scripts/fig_phengenI_volcano.R <input_sets> <out_dir> <top_n> <pheno_int> [labels]
#
# Arguments:
#   input_sets   semicolon-separated sets; within each set, files are
#                comma-separated. Each set is assembled as its own panel
#                figure with shared axis limits internal to that set.
#   out_dir      output directory for pdf/png/svg
#   top_n        number of top terms to label (by Combined Score)
#   pheno_int    phenotype of interest to force-label (use "NULL" to disable)
#   labels       OPTIONAL panel labels, same nested structure as input_sets
#                (';' between sets, ',' within). Use "NULL" or omit to derive
#                each label from its filename; a name starting with "network"
#                becomes "Whole network", otherwise the prefix before "_step"
#                is used.
#
# Examples:
#   Single set:
#     Rscript figure_scripts/fig_phengenI_volcano.R \
#       "network_step0.txt,network_step1.txt" figures 20 Leptin
#   Two sets with explicit labels (whole network + a cluster):
#     Rscript figure_scripts/fig_phengenI_volcano.R \
#       "network_step0.txt,network_step1.txt;cluster2_step0.txt,cluster2_step1.txt" \
#       figures 20 Leptin \
#       "Whole network,Whole network;Cluster 2,Cluster 2"
# ============================================================

# ── ARGUMENTS ─────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript figure_scripts/fig_phengenI_volcano.R <input_sets> <out_dir> <top_n> <pheno_int> [labels]\n",
    "  input_sets: semicolon-separated sets, comma-separated files within each set\n",
    "  pheno_int:  phenotype term to force-label, or NULL to disable\n",
    "  labels:     optional panel labels mirroring input_sets, or NULL to auto-derive",
    call. = FALSE
  )
}

# Parse sets: list of character vectors
IN_SETS   <- lapply(strsplit(args[1], ";")[[1]], function(s) strsplit(s, ",")[[1]])
OUT_DIR   <- args[2]
TOP_N     <- as.integer(args[3])
PHENO_INT <- if (args[4] == "NULL") NULL else args[4]

# Optional explicit panel labels, mirroring IN_SETS' nested structure.
LABELS <- if (length(args) >= 5 && args[5] != "NULL") {
  lapply(strsplit(args[5], ";")[[1]], function(s) strsplit(s, ",")[[1]])
} else NULL

# ── PACKAGES ──────────────────────────────────────────────────
pkgs <- c("ggplot2", "ggrepel", "scales", "patchwork", "svglite")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
  library(p, character.only = TRUE)
}))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── CONSTANTS ─────────────────────────────────────────────────
PADJ_THRESH     <- 0.05
padj_thresh_log <- -log10(PADJ_THRESH)

# Colours. This script takes positional arguments rather than flags, so it has
# no --palette; it honours the FIGURE_PALETTE environment variable instead.
# Each fallback is the literal used before, so the default rendering is
# unchanged.
source("R/palettes.R")
pal      <- get_palette(NULL)
sig_lvls <- pal_cat(pal, 3, c("#d62728", "#ff7f0e", "#1f77b4"))

sig_colours <- c(
  "FDR < 0.001" = sig_lvls[1],
  "FDR < 0.01"  = sig_lvls[2],
  "FDR < 0.05"  = sig_lvls[3],
  "ns"          = pal_col(pal, "neutral", "#aec7e8")
)
col_ink   <- pal_col(pal, "ink",   "black")
col_ref   <- pal_col(pal, "axis",  "grey50")
col_grid  <- pal_col(pal, "grid",  "grey93")
# Three shades of secondary text. A palette collapses them onto its single
# `muted`; with no palette each keeps the distinct grey it always had.
col_seg   <- pal_col(pal, "muted", "grey30")   # ggrepel leader lines
col_note  <- pal_col(pal, "muted", "grey45")   # italic in-panel annotation
col_sub   <- pal_col(pal, "muted", "grey40")   # subtitle

# ── HELPERS ───────────────────────────────────────────────────
parse_suffix <- function(filepath) {
  base_in <- tools::file_path_sans_ext(basename(filepath))
  suffix  <- sub("^PhenGenI_Association_2021_table_", "", base_in)
  suffix  <- gsub("_step_", "_step", suffix)
  suffix
}

# Natural (human) ordering: zero-pad every digit run so that embedded
# numbers sort numerically ("plus30" < "plus120") rather than lexically.
natural_order <- function(x) {
  m <- gregexpr("[0-9]+", x)
  padded <- regmatches(x, m)
  regmatches(x, m) <- lapply(padded, function(nums)
    formatC(as.numeric(nums), width = 18, flag = "0", format = "d"))
  order(x)
}

# Derive a human panel label from a file's suffix when none is supplied.
# Strips the step/plus marker (and any trailing assay cruft like
# "_PhenGenI_Association_2021"), so "cluster1_plus30_..." -> "Cluster 1" and
# "network_step0" -> "Whole network".
derive_label <- function(suffix) {
  prefix <- sub("[_-]?step[_-]?[0-9]+$", "", suffix, ignore.case = TRUE)  # _stepN
  prefix <- sub("[_-]plus[0-9]+.*$",     "", prefix, ignore.case = TRUE)  # _plusN + trailing cruft
  prefix <- gsub("([A-Za-z])([0-9])", "\\1 \\2", prefix)                  # "cluster1" -> "cluster 1"
  prefix <- trimws(gsub("[_-]+", " ", prefix))
  if (prefix == "") return("Set")
  if (grepl("network", prefix, ignore.case = TRUE)) return("Whole network")
  paste0(toupper(substring(prefix, 1, 1)), substring(prefix, 2))
}

# Describe the panel's step/depth for the title: "plusN" -> "+N genes",
# "stepN" -> "step N"; falls back to a trailing number if neither marker exists.
parse_step_desc <- function(suffix) {
  m <- regmatches(suffix, regexec("[_-]plus([0-9]+)", suffix, ignore.case = TRUE))[[1]]
  if (length(m) == 2) return(paste0("+", m[2], " genes"))
  m <- regmatches(suffix, regexec("step[_-]?([0-9]+)", suffix, ignore.case = TRUE))[[1]]
  if (length(m) == 2) return(paste0("step ", m[2]))
  num <- regmatches(suffix, regexpr("[0-9]+$", suffix))
  if (length(num) && nzchar(num)) paste0("step ", num) else ""
}

build_df <- function(filepath) {
  tbl <- read.delim(filepath, sep = "\t", header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
  tbl$neg_log10_pval    <- -log10(tbl[["P-value"]])
  tbl$neg_log10_adjpval <- -log10(tbl[["Adjusted P-value"]])
  tbl$overlap_num       <- as.integer(sub("/.*", "", tbl$Overlap))
  tbl$overlap_den       <- as.integer(sub(".*/", "", tbl$Overlap))
  tbl$overlap_frac      <- tbl$overlap_num / tbl$overlap_den
  tbl$sig <- ifelse(tbl[["Adjusted P-value"]] < 0.001, "FDR < 0.001",
                    ifelse(tbl[["Adjusted P-value"]] < 0.01,  "FDR < 0.01",
                           ifelse(tbl[["Adjusted P-value"]] < 0.05,  "FDR < 0.05", "ns")))
  tbl$sig <- factor(tbl$sig, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns"))
  tbl
}

save_plot <- function(plot_obj, base_path, width, height) {
  ggsave(paste0(base_path, ".pdf"),
         plot = plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(base_path, ".png"),
         plot = plot_obj, width = width, height = height, dpi = 300)
  ggsave(paste0(base_path, ".svg"),
         plot = plot_obj, width = width, height = height, device = "svg",
         fix_text_size = FALSE)
}

# ── PLOT FUNCTION ─────────────────────────────────────────────
make_plot <- function(df, suffix, panel_label, shared_x = NULL, shared_y = NULL, label_terms = TRUE) {
  step_desc  <- parse_step_desc(suffix)
  step_label <- panel_label

  sig_terms <- df[df[["Adjusted P-value"]] < PADJ_THRESH, ]
  top_terms <- sig_terms[order(-sig_terms[["Combined Score"]]), ][seq_len(min(TOP_N, nrow(sig_terms))), ]
  
  if (!is.null(PHENO_INT)) {
    pheno_row <- df[df$Term == PHENO_INT, ]
    if (nrow(pheno_row) == 0) warning("PHENO_INT '", PHENO_INT, "' not found in: ", suffix)
    top_terms <- unique(rbind(top_terms, pheno_row))
  }
  
  p <- ggplot(df, aes(
    x     = .data[["Odds Ratio"]],
    y     = neg_log10_adjpval,
    size  = overlap_num,
    color = sig,
    fill  = sig
  )) +
    
    geom_hline(yintercept = padj_thresh_log, linetype = "dashed",
               color = col_ref, linewidth = 0.5) +
    
    geom_point(shape = 21, alpha = 0.85, stroke = 0.8) +
    
    {
      if (!is.null(PHENO_INT)) {
        pheno_df <- df[df$Term == PHENO_INT, ]
        if (nrow(pheno_df) > 0) {
          geom_point(
            data        = pheno_df,
            shape       = 21,
            stroke      = 2.2,
            colour      = col_ink,
            fill        = NA,
            alpha       = 1,
            size        = 5,
            show.legend = FALSE,
            inherit.aes = FALSE,
            mapping     = aes(x = .data[["Odds Ratio"]], y = neg_log10_adjpval)
          )
        }
      }
    } +
    
    (if (label_terms) {
      geom_text_repel(
        data               = top_terms,
        aes(label          = Term),
        size               = 2.6,
        fontface           = "plain",
        color              = col_ink,
        box.padding        = unit(0.9, "lines"),
        point.padding      = unit(0.6, "lines"),
        segment.color      = col_seg,
        segment.size       = 0.35,
        segment.curvature  = -0.1,
        min.segment.length = 0,
        max.overlaps       = Inf,
        force              = 2.5,
        force_pull         = 0.5,
        max.iter           = 20000,
        max.time           = 2,
        seed               = 42
      )
    } else NULL) +
    
    annotate("text",
             x     = if (!is.null(shared_x)) shared_x[2] else max(df[["Odds Ratio"]], na.rm = TRUE),
             y     = padj_thresh_log + 0.12,
             label = paste0("p[adj] == ", PADJ_THRESH),
             parse = TRUE,
             hjust = 1, vjust = 0,
             size  = 2.8, color = col_note, fontface = "italic") +
    
    scale_color_manual(values = sig_colours, name = "Significance (BH-FDR)") +
    scale_fill_manual( values = sig_colours, name = "Significance (BH-FDR)") +
    scale_size_continuous(
      name   = "Genes in overlap",
      range  = c(1.5, 8),
      breaks = c(2, 5, 10, 15)
    ) +
    scale_x_continuous(
      name   = "Odds Ratio",
      breaks = pretty_breaks(n = 8),
      expand = expansion(mult = c(0.02, 0.05)),
      limits = shared_x
    ) +
    scale_y_continuous(
      name   = expression(-log[10]~"adjusted p-value"),
      breaks = pretty_breaks(n = 8),
      expand = expansion(mult = c(0.02, 0.08)),
      limits = shared_y
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle    = element_text(size = 9, color = col_sub, hjust = 0,
                                      margin = margin(b = 8)),
      axis.title       = element_text(size = 10),
      axis.text        = element_text(size = 9),
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.55, "cm"),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = col_grid),
      plot.margin      = margin(12, 14, 10, 12)
    ) +
    
    labs(
      title    = if (nzchar(step_desc)) paste0(step_label, "  |  ", step_desc) else step_label,
      subtitle = paste0("n = ", nrow(df), " phenotype terms")
    )
  
  p
}

# ── PROCESS EACH SET ──────────────────────────────────────────
for (set_idx in seq_along(IN_SETS)) {
  IN_FILES <- IN_SETS[[set_idx]]

  # Order files within the set by their embedded step/plus number so panels
  # read 30 → 60 → 90 → 120, regardless of the caller's (often lexical) order.
  ord      <- natural_order(IN_FILES)
  IN_FILES <- IN_FILES[ord]
  # Keep any explicit labels aligned with their (now reordered) files.
  if (!is.null(LABELS) && length(LABELS) >= set_idx &&
      length(LABELS[[set_idx]]) == length(ord))
    LABELS[[set_idx]] <- LABELS[[set_idx]][ord]

  message("\n── Processing set ", set_idx, " (", length(IN_FILES), " file(s)) ──")
  
  all_dfs <- lapply(IN_FILES, build_df)
  
  # Axis limits shared within this set only
  set_xlim <- range(sapply(all_dfs, function(d) range(d[["Odds Ratio"]],   na.rm = TRUE)))
  set_ylim <- range(sapply(all_dfs, function(d) range(d$neg_log10_adjpval, na.rm = TRUE)))
  set_xlim[2] <- set_xlim[2] * 1.05
  set_ylim[2] <- set_ylim[2] * 1.08
  
  plots        <- vector("list", length(IN_FILES))
  plots_no_ann <- vector("list", length(IN_FILES))
  
  for (i in seq_along(IN_FILES)) {
    df        <- all_dfs[[i]]
    suffix    <- parse_suffix(IN_FILES[[i]])
    base_name <- paste0("volcano_phengenI_set", set_idx, "_", suffix)

    # Panel label: explicit (if supplied) else derived from the filename
    panel_label <- if (!is.null(LABELS) && length(LABELS) >= set_idx &&
                       length(LABELS[[set_idx]]) >= i) {
      trimws(LABELS[[set_idx]][i])
    } else {
      derive_label(suffix)
    }

    # Individual plot with its own axis limits
    p_ind <- make_plot(df, suffix, panel_label, shared_x = NULL, shared_y = NULL)
    save_plot(p_ind, file.path(OUT_DIR, base_name), width = 10, height = 7)
    message("Saved individual: ", base_name, ".{pdf,png,svg}")

    # Shared-axis versions for panel assembly
    plots[[i]]        <- make_plot(df, suffix, panel_label,
                                   shared_x = set_xlim, shared_y = set_ylim)
    plots_no_ann[[i]] <- make_plot(df, suffix, panel_label,
                                   shared_x = set_xlim, shared_y = set_ylim,
                                   label_terms = FALSE)
  }
  
  # Panel assembly for this set
  if (length(IN_FILES) > 1) {
    ncols <- min(length(plots), 2)
    nrows <- ceiling(length(plots) / ncols)
    
    panel_name <- paste0("volcano_phengenI_set", set_idx, "_panel_",
                         paste(sapply(IN_FILES, parse_suffix), collapse = "_vs_"))
    
    panel_plot        <- wrap_plots(plots,        ncol = ncols, nrow = nrows)
    panel_plot_no_ann <- wrap_plots(plots_no_ann, ncol = ncols, nrow = nrows)
    
    save_plot(panel_plot,
              file.path(OUT_DIR, panel_name),
              width = 10 * ncols, height = 7 * nrows)
    
    save_plot(panel_plot_no_ann,
              file.path(OUT_DIR, paste0(panel_name, "_no_ann")),
              width = 10 * ncols, height = 7 * nrows)
    
    message("Saved panel: ", panel_name, ".{pdf,png,svg}")
    
    # Inkscape-ready individual SVGs
    inkscape_dir <- file.path(OUT_DIR, paste0("inkscape_panels_set", set_idx))
    dir.create(inkscape_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (i in seq_along(IN_FILES)) {
      suffix    <- parse_suffix(IN_FILES[[i]])
      base_name <- paste0("panel_", suffix)
      
      save_plot(plots[[i]],
                file.path(inkscape_dir, base_name),
                width = 10, height = 7)
      
      save_plot(plots_no_ann[[i]],
                file.path(inkscape_dir, paste0(base_name, "_no_ann")),
                width = 10, height = 7)
    }
    
    message("Saved Inkscape-ready panels to: ", inkscape_dir)
  }
}