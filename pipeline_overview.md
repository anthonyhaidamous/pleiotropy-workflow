# Pipeline overview

The 10 stages in `pipeline/` run sequentially from the project root, each taking the phenotype's YAML
config as its argument (stage 10 takes CLI flags instead). Stages hand off data as RDS files under
`results/intermediate/<phenotype>/`, with TSV mirrors for inspection. Manual-checkpoint stages additionally
write a `_review.tsv` to `results/review/<phenotype>/` that the analyst edits, saves as `_reviewed.tsv`, and
the stage re-reads on the next run.

Most stages obtain their packages by sourcing `R/utils.R`, which loads `yaml`, `dplyr`, `tidyr`, `stringr`,
`gtools`, and `jsonlite`. Stage-specific packages are noted per stage below.

Cross-cutting rules that several stages share — trait-name normalisation, and the explicit-exclusion
convention — are in [`.context/conventions.md`](../.context/conventions.md), and the rulings
behind them in [`.context/decisions.md`](../.context/decisions.md).

## Figure render order

**F2 must be rendered last**, after stages 1–10 have completed and the block table is final.
`figure_scripts/fig_karyoploter_trait_ld_blocks.R` draws directly from the canonical block table, so
rendering it earlier bakes in a block set that a later stage-4 re-run can still change. Every other
figure can be rendered in any order once its input stage has run.

Requires `karyoploteR` (Bioconductor). The repo has **no lockfile of any kind** (see `AUDIT.md` §4.1),
so this is declared in `README.md`'s install block rather than pinned; reproducing F2 from a clean
environment depends on that install line being run.

---

## Stage 1 — `01_load_snps.R`

**What it does:** Loads and filters the GWAS Catalog studies + associations for the phenotype, producing the
lead-SNP (lSNP) set with encoded ancestry plus the full significant-association table.

- **Input:** GWAS Catalog studies TSV and associations TSV (paths from the YAML config); ancestry file.
- **Output:** `01_lSNPs.rds` (lead SNPs with ancestry), `01_gwassos_sigP.rds` (full filtered associations).
- **Manual checkpoint:** **Yes.** Writes `01_lSNPs_review.tsv` (ancestry encoding failures, annotate-only) and
  `01_pvalue_review.tsv` (lSNPs at the likely-rounded p-value threshold). Re-run reads `01_pvalue_reviewed.tsv`
  and applies the analyst's drop/correct decisions.
- **Packages:** via `R/utils.R`.

## Stage 2 — `02_ldlink.R`

**What it does:** Retrieves LDproxy results for each lead SNP via the LDlink API.

- **Input:** `01_lSNPs.rds`; `LDLINK_TOKEN` environment variable.
- **Output:** `02_proxy_table.rds` (+ `02_proxy_table.tsv` mirror), `02_lSNPs.rds` (lSNPs after dropping any
  SNP that failed in every queried population).
- **Manual checkpoint:** No — writes `02_failures_review.tsv` as a report of failed SNP–population pairs, but
  does not gate on a `_reviewed.tsv`.
- **Packages:** `LDlinkR`, plus `R/utils.R`.

## Stage 3 — `03_ld_proxies.R`

**What it does:** Splits LDproxy results into the full LD set used for block boundaries (`leadSNPs_LDs`, R² >
`ld_block_r2`) and the high-LD proxy set aggregated per lead SNP (`pSNPs`, R² ≥ `ld_proxy_r2`).

- **Input:** `02_proxy_table.rds`, `02_lSNPs.rds`.
- **Output:** `03_leadSNPs_LDs.rds`, `03_pSNPs_per_lSNP.rds`.
- **Manual checkpoint:** No.
- **Packages:** via `R/utils.R`.

## Stage 4 — `04_blocks.R`

**What it does:** Aggregates lead SNPs into LD-defined blocks via connected components (LD-fusion), then
optionally merges block pairs listed in `cfg$block_merges`.

- **Input:** `02_lSNPs.rds`, `03_leadSNPs_LDs.rds`, `03_pSNPs_per_lSNP.rds`.
- **Output:** `04_blocks.rds` (rich tibble with list-columns) + `04_blocks.tsv` (flattened mirror).
- **Manual checkpoint:** **Yes.** Writes `04_block_overlaps_review.tsv` (overlapping block pairs). The analyst
  decides which to merge and records them under the YAML `block_merges` key, then re-runs.
- **Packages:** via `R/utils.R`.

## Stage 5 — `05_crosstraits.R`

**What it does:** For every SNP in every block (lSNPs and pSNPs), retrieves every other significant GWAS
Catalog association — the raw pleiotropic set, before quality filtering.

- **Input:** `04_blocks.rds`, `01_gwassos_sigP.rds`.
- **Output:** `05_crosstraits_raw.rds` (+ TSV mirror).
- **Manual checkpoint:** **Yes.** Writes `05_unmapped_ancestry_review.tsv` for cross-trait studies that can't
  be ancestry-mapped; re-run reads `05_unmapped_ancestry_reviewed.tsv` (after the analyst updates
  `cfg$ancestry_mapping` as needed).
- **Packages:** via `R/utils.R`.

## Stage 6 — `06_filter.R`

**What it does:** Applies four sequential filters to the raw cross-traits: (1) ancestry concordance with the
block, (2) total sample size ≥ `cfg$filters$min_sample_size`, (3) CNV exclusion (`CNV == "N"`), (4) GxE /
multi-trait two-tier regex (auto-exclude + manual review).

- **Input:** `05_crosstraits_raw.rds`, `04_blocks.rds`.
- **Output:** `06_crosstraits_filtered.rds` (+ TSV mirror), `06_filtering_log.tsv`.
- **Manual checkpoint:** **Yes.** Filter (4) writes `06_gxe_review.tsv` (rows flagged `manual_check`); re-run
  reads `06_gxe_reviewed.tsv` with the `gxe_decision` column set to `keep` or `auto_exclude`.
- **Packages:** via `R/utils.R`.

## Stage 7 — `07_categorize.R`

**What it does:** Builds a unique-trait frequency table from the filtered cross-traits, matches it against an
optional trait reference of categories / hierarchized labels, and flags unmatched traits as `new`.

- **Input:** `06_crosstraits_filtered.rds`; optional trait reference TSV (`cfg$paths$trait_ref`).
- **Output:** `07_trait_table.rds` (+ TSV mirror).
- **Manual checkpoint:** **Yes.** Writes `07_trait_table_review.tsv` with `Category` and `Hierarchized Trait`
  set to `new` for unmatched traits; re-run reads `07_trait_table_reviewed.tsv`. If no reference is configured,
  every trait is `new` on first run (expected for a brand-new phenotype).
- **Packages:** via `R/utils.R`.

## Stage 8 — `08_enrichment.R`

**What it does:** For each trait and category, compares the within-phenotype frequency against the catalog-wide
background; computes fold enrichment, log2-normalized enrichment, Fisher's exact p-value, and BH-FDR adjusted
p-value. Adds columns only — no filtering or reshaping.

- **Input:** `07_trait_table.rds`, `01_gwassos_sigP.rds`.
- **Output:** `08_trait_enrichment.rds`, `08_category_enrichment.rds` (+ TSV mirrors).
- **Manual checkpoint:** No.
- **Packages:** via `R/utils.R`.

## Stage 9 — `09_harmonize.R`

**What it does:** Harmonizes effect direction and allele scale across studies: OR → beta where the CI text
lacks increase/decrease wording; applies sign from "decrease" wording; flips beta/CI to a common reference
allele (CA2 convention) using the LDproxy allele pairs.

- **Input:** `06_crosstraits_filtered.rds`, `04_blocks.rds`, `03_leadSNPs_LDs.rds`, `07_trait_table.rds`,
  `08_trait_enrichment.rds`, `01_lSNPs.rds`.
- **Output:** `09_harmonized_trait_table.rds`, `09_block_directions.rds` (+ TSV mirrors). This is the core
  pipeline's end product.
- **Manual checkpoint:** No.
- **Packages:** via `R/utils.R`.

## Stage 10 — `10_extract_genes.R`

**What it does:** Builds a protein-coding gene list per LD block from two sources — a biomaRt query on each
block's genomic position (± flank) and the gene symbols already in the `04_blocks` Genes column — then merges
and deduplicates.

- **Input:** `04_blocks.{tsv,rds}` (via the `--blocks` flag; this stage takes CLI flags, not the YAML).
- **Output (to `--outdir`):** `protein_coding_genes_by_block.tsv` (one row per gene × block),
  `protein_coding_genes_collapsed.tsv` (one row per block, genes comma-joined),
  `genes_list_complete.txt` (flat unique gene list — the gene total, e.g. leptin's 93),
  `genes_string_submission.txt` (the subset that is valid STRING input — this is what
  `string.input.gene_list` points at), `genes_excluded_from_string.tsv` (retained genes held back
  from the submission, with reason, for the supplement).
- **Manual checkpoint:** No.
- **Note:** a bare Ensembl ID is a real gene and counts toward the gene total, but is not a valid
  STRING input — STRING is queried by symbol. Such genes are retained and counted, excluded from the
  submission by rule, named in the run log, and written to `genes_excluded_from_string.tsv`
  (1 gene in leptin, 8 in CLL). See `.context/conventions.md` §3.
- **Packages:** `biomaRt`, `GenomicRanges`, `ensembldb`, `AnnotationHub`, `dplyr`, `tidyr`, `optparse`
  (Bioconductor; installs `BiocManager` if absent).
