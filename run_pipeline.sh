#!/usr/bin/env bash
# =============================================================================
# Run the full GWAS pleiotropy pipeline (stages 1-10) for one phenotype.
#
# Stages 1-9 take the YAML config; stage 10 takes CLI flags, so its paths are
# derived from the config's phenotype.name slug + paths.intermediate.
#
# MANUAL CHECKPOINTS
#   Stages 1, 5, 6, 7 pause for analyst review (and stage 4 reads block merges
#   from the YAML). They run unattended ONLY when the matching *_reviewed.tsv
#   already exists — true for an already-reviewed phenotype (leptin, CLL).
#
#   By default this script PRE-FLIGHTS the checkpoint state and refuses to run
#   if a review is pending (a stage wrote *_review.tsv but no *_reviewed.tsv yet),
#   or if there is no review directory at all (a genuine first run). It tells you
#   exactly which file to complete. A checkpoint that produced no review file
#   simply never needed review and is not a blocker. Pass --first-run to proceed
#   regardless (stages emit *_review.tsv for you to fill in, then re-run).
#
# LOGGING
#   By default the full console output (stdout+stderr) is mirrored to a
#   timestamped file logs/<slug>_<YYYYmmdd_HHMMSS>.log. Disable with --no-log.
#   Dry runs are not logged.
#
# Usage:
#   export LDLINK_TOKEN=<your_token>           # needed by stage 2
#   ./run_pipeline.sh config/cll.yml           # re-run a reviewed phenotype (logs to logs/)
#   ./run_pipeline.sh config/leptin.yml --first-run
#   ./run_pipeline.sh config/cll.yml --dry-run # print commands only
#   ./run_pipeline.sh config/cll.yml --no-log  # don't write a log file
#   ./run_pipeline.sh config/cll.yml --skip 2  # skip stage 2 (proxy files already present)
#   ./run_pipeline.sh config/cll.yml --skip 2,3 --skip 4  # skip multiple stages
# =============================================================================
set -euo pipefail

CONFIG=""; DRY_RUN=0; FIRST_RUN=0; NO_LOG=0; SKIP_STAGES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)   DRY_RUN=1 ;;
    --first-run) FIRST_RUN=1 ;;
    --no-log)    NO_LOG=1 ;;
    --skip)
      shift
      if [[ $# -eq 0 || "$1" == -* ]]; then
        echo "ERROR: --skip requires a stage number (e.g. --skip 2 or --skip 2,3)" >&2; exit 2
      fi
      IFS=',' read -ra _nums <<< "$1"
      SKIP_STAGES+=("${_nums[@]}")
      ;;
    --skip=*)
      IFS=',' read -ra _nums <<< "${1#--skip=}"
      SKIP_STAGES+=("${_nums[@]}")
      ;;
    -*)  echo "ERROR: unknown flag: $1" >&2; exit 2 ;;
    *)   CONFIG="$1" ;;
  esac
  shift
done
if [[ -z "$CONFIG" ]]; then
  echo "usage: $0 config/<trait>.yml [--first-run] [--dry-run] [--skip N[,M...]]" >&2; exit 2
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: config not found: $CONFIG" >&2; exit 2
fi

# Derive the per-phenotype intermediate + review dirs (same slug rule as R/utils.R).
read -r INTERMEDIATE REVIEW_DIR < <(Rscript -e '
  cfg <- yaml::read_yaml(commandArgs(TRUE)[1])
  s <- tolower(cfg$phenotype$name)
  s <- gsub("[^a-z0-9]+", "_", s); s <- gsub("^_+|_+$", "", s)
  ib <- if (!is.null(cfg$paths$intermediate)) cfg$paths$intermediate else "results/intermediate"
  rb <- if (!is.null(cfg$paths$review))       cfg$paths$review       else "results/review"
  cat(file.path(ib, s), file.path(rb, s), "\n")' "$CONFIG")

# ── Logging: mirror everything below to a timestamped file in logs/ ──────────
# Simple capture of the console output (stdout+stderr) — nothing more. Skipped
# for --dry-run, or disable with --no-log.
LOG=""
if (( ! DRY_RUN && ! NO_LOG )); then
  mkdir -p logs
  LOG="logs/$(basename "$INTERMEDIATE")_$(date +%Y%m%d_%H%M%S).log"
  exec > >(tee -a "$LOG") 2>&1
fi

echo "Config:        $CONFIG"
echo "Intermediate:  $INTERMEDIATE"
echo "Review dir:    $REVIEW_DIR"
[[ -n "$LOG" ]] && echo "Log:           $LOG"
(( ${#SKIP_STAGES[@]} )) && echo "Skip stages:   ${SKIP_STAGES[*]}"

# ── Pre-flight: manual-checkpoint state ──────────────────────────────────────
# A review is PENDING when a stage wrote <label>_review.tsv but no matching
# <label>_reviewed.tsv exists yet. A label with NEITHER file simply never needed
# review (e.g. CLL had no unmapped ancestries, so stage 5 wrote nothing) — that
# is not a blocker. No review dir at all = a genuine first run for this phenotype.
CHECKPOINTS=("1:01_pvalue" "5:05_unmapped_ancestry" "6:06_gxe" "7:07_trait_table")

if [[ ! -d "$REVIEW_DIR" ]]; then
  echo
  if (( FIRST_RUN )); then
    echo "first-run mode: no review dir yet — stages 1/5/6/7 will write *_review.tsv as needed."
  else
    echo "HALTED: no review directory yet ($REVIEW_DIR)." >&2
    echo "This looks like a first run for this phenotype. Re-run with --first-run to begin;" >&2
    echo "stages 1/5/6/7 will emit *_review.tsv for you to complete, then re-run normally." >&2
    exit 3
  fi
else
  pending=()
  for entry in "${CHECKPOINTS[@]}"; do
    stage="${entry%%:*}"; label="${entry#*:}"
    if [[ -f "$REVIEW_DIR/${label}_review.tsv" && ! -f "$REVIEW_DIR/${label}_reviewed.tsv" ]]; then
      pending+=("$stage:$label")
    fi
  done
  if (( ${#pending[@]} )); then
    echo
    if (( FIRST_RUN )); then
      echo "first-run mode: proceeding with ${#pending[@]} pending review(s)."
    else
      echo "HALTED: pending manual review(s) in $REVIEW_DIR:" >&2
      for m in "${pending[@]}"; do
        st="${m%%:*}"; lb="${m#*:}"
        echo "    stage $st — edit ${lb}_review.tsv, save as ${lb}_reviewed.tsv" >&2
      done
      echo >&2
      echo "Complete the file(s) above (save each as *_reviewed.tsv), then re-run." >&2
      echo "Or pass --first-run to proceed without them." >&2
      exit 3
    fi
  fi
fi

[[ -n "${LDLINK_TOKEN:-}" ]] || echo "WARNING: LDLINK_TOKEN is not set — stage 2 will fail." >&2

_is_skipped() {
  local n="$1"
  for s in "${SKIP_STAGES[@]+"${SKIP_STAGES[@]}"}"; do
    [[ "$s" == "$n" ]] && return 0
  done
  return 1
}

run() {  # run <stage_num> <cmd> [args...]
  local stage="$1"; shift
  echo; echo "=============================================================="
  if _is_skipped "$stage"; then
    echo ">> [SKIP $stage] $*"
    echo "=============================================================="
    return 0
  fi
  echo ">> [stage $stage] $*"
  echo "=============================================================="
  if (( DRY_RUN )); then return 0; fi
  "$@"
}

run 1  Rscript pipeline/01_load_snps.R   "$CONFIG"
run 2  Rscript pipeline/02_ldlink.R      "$CONFIG"
run 3  Rscript pipeline/03_ld_proxies.R  "$CONFIG"
run 4  Rscript pipeline/04_blocks.R      "$CONFIG"
run 5  Rscript pipeline/05_crosstraits.R "$CONFIG"
run 6  Rscript pipeline/06_filter.R      "$CONFIG"
run 7  Rscript pipeline/07_categorize.R  "$CONFIG"
run 8  Rscript pipeline/08_enrichment.R  "$CONFIG"
run 9  Rscript pipeline/09_harmonize.R   "$CONFIG"
run 10 Rscript pipeline/10_extract_genes.R --blocks "$INTERMEDIATE/04_blocks.rds" --outdir "$INTERMEDIATE"

echo; echo "Pipeline complete -> $INTERMEDIATE"
[[ -n "$LOG" ]] && echo "Log saved -> $LOG"
