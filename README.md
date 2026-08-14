# GWAS Pleiotropy Pipeline

A modular, phenotype-agnostic pipeline that takes raw GWAS Catalog inputs and produces a harmonized
cross-trait table for any phenotype of interest. Built around a leptin pleiotropy analysis but designed
to run on any trait by editing a single YAML config.

Developed and validated on **leptin** (single trait) and **chronic lymphocytic leukemia (CLL)** (complex disease).

## What it does

Given a phenotype (defined by one or more `MAPPED_TRAIT` terms in the GWAS Catalog), the pipeline:

1. Identifies lead SNPs significantly associated with the phenotype.
2. Pulls their LD proxies from LDlink.
3. Defines LD blocks at the locus level, merging blocks that overlap genomically.
4. Retrieves every other significant association for SNPs in those blocks (cross-trait associations).
5. Filters cross-traits for ancestry concordance, sample size, CNV exclusion, and gene-by-environment / multi-trait artifacts.
6. Categorizes each cross-trait into curated trait categories.
7. Computes enrichment of each category and trait against the catalog-wide background.
8. Harmonizes effect sizes (OR → beta, risk allele → reference allele) across studies.
9. Extracts the protein-coding genes in each LD block (biomaRt query + the catalog's own gene symbols).

The core pipeline's end product is `09_harmonized_trait_table.rds`: a single table with every harmonized
cross-trait association annotated by block, category, ancestry, sample size, and enrichment statistics.
Stage 10 adds a per-block gene list - later used in the STRING section of the study. See `pipeline/string_pipeline/STRING_README.md` for full details.

## Repository structure

| Directory | Contents |
|-----------|----------|
| `pipeline/` | 10 numbered core analysis stages, run sequentially, including the STRING section in `pipeline/string_pipeline/` |
| `figure_scripts/` | Figure and table generation scripts |
| `figures/` | Figure/table outputs per phenotype |
| `R/` | Shared utility functions (`utils.R`, `figure_utils.R`, `string_utils.R`, `palettes.R`) and R-specific environment setup (`deps.R`) |
| `config/` | Phenotype-specific YAML config files |
| `data/` | Input data (GWAS Catalog studies/associations TSVs, chromosome sizes, trait_ref) |
| `results/` | Pipeline outputs: `intermediate/<phenotype>/` RDS handoffs + TSV mirrors, and `review/` checkpoint files |
| `logs/` | Per-phenotype run logs |

See `pipeline_overview.md` for full per-stage I/O and package details.

## Data

`data/` is not tracked in git (apart from the chromosome sizes table) — the GWAS Catalog
files are large and versioned by release date, so you download the ones you want and point
the config at them. A clone therefore ships with `data/associations/`, `data/studies/`, and
`data/trait_ref/` empty.

Download both files from the [GWAS Catalog downloads page](https://www.ebi.ac.uk/gwas/docs/file-downloads)
and place them as follows:

| File | Put it at | Config key |
|------|-----------|------------|
| All studies | `data/studies/<studies-file>.tsv` | `paths.studies` |
| All associations | `data/associations/<associations-file>.tsv` | `paths.associations` |

The shipped configs name the exact releases they were run on, e.g.

```yaml
studies:      "data/studies/gwas-catalog-v1.0.3.1-studies-r2026-02-16.tsv"
associations: "data/associations/gwas-catalog-download-associations-alt-full.tsv"
```

Edit those two lines to match whatever release you downloaded — results are release-dependent,
so record which one you used.

- **Trait references** (`data/trait_ref/<phenotype>/`): the curated trait → category mapping
  produced by stage 7's review loop. Not shipped; for a new phenotype leave `trait_ref` unset
  and stage 7 flags every trait as `new` for you to categorize (see *Adding a new phenotype*).
- **Chromosome sizes** (`data/chromosomes/chrom_sizes.txt`): GRCh38 chromosome lengths, used
  for block mapping on the genome. This one *is* tracked — no download needed.

## Requirements

- R ≥ 4.1 (developed and run on **4.5.0**; 4.1 is not verified)
- About thirty R packages across CRAN and Bioconductor. `R/deps.R` is the
  authoritative list — package, where it comes from, and which stage wants it —
  so there is no second list here to fall out of date with it.
- **System dependencies (not R packages, and neither one installs from R):**
  - `mcl` on `PATH` — the STRING module stage (`config/<pheno>.yml` →
    `string.mcl`).
  - `poppler-cpp`, behind `pdftools`. A *binary* pdftools needs only the runtime
    library (`libpoppler-cpp0v5`, already present on most desktops); a *source*
    build needs the headers (`libpoppler-cpp-dev`).
- An [LDlink API token](https://ldlink.nih.gov/?tab=apiaccess) exported as `LDLINK_TOKEN` (stage 2)
  - Enter your email address in the registration form
  - LDlink emails you a personal access token (a string of letters/numbers)
  - Use that token in your API calls

### Setting a machine up

```bash
Rscript tools/setup_env.R --check    # what is present, what is missing, where it would install
Rscript tools/setup_env.R            # install everything missing
```

`setup_env.R` reads the manifest in `R/deps.R` and installs into whatever
library the session is already using, so the same command works in both places
this pipeline runs:

- **On a bare machine** renv activates and owns `renv/library/`. Packages
  recorded in `renv.lock` are restored at the pinned version that produced the
  locked numbers; anything newer than the lockfile is installed and then wants a
  `renv::snapshot()`, which the script tells you.

```bash
Rscript tools/setup_env.R --all
```

## Running the pipeline

### One command (recommended)

`run_pipeline.sh` runs all 10 stages in order, stops on the first error, and pre-flights the manual
checkpoints so it won't silently skip a pending review:

```bash
export LDLINK_TOKEN=<your_token>              # required by stage 2

./run_pipeline.sh config/cll.yml --first-run  # first pass on any phenotype (see checkpoints below)
./run_pipeline.sh config/cll.yml --dry-run    # print the commands without running them
./run_pipeline.sh config/cll.yml              # re-run once your reviews are saved
```

If a stage has a review pending (it wrote `*_review.tsv` but you haven't saved `*_reviewed.tsv`), or the
phenotype has never been run, the script halts and tells you exactly which file to complete. Pass
`--first-run` to start a new phenotype anyway; the stages will emit their review files for you to fill in.

**A fresh clone always starts at `--first-run`, including for the shipped `leptin.yml` and `cll.yml`
configs.** The manual-review decisions and trait references from the original leptin and CLL analyses
are deliberately not distributed: they are analyst judgements about specific traits, and yours will
differ from ours. What this repository publishes is the method — the stages, the rules they apply, and
the review interface — not one particular set of rulings. Expect to work through the checkpoints below
on your first pass over any phenotype.

### Stage by stage

Equivalently, run each stage yourself from the project root (this is what you do on a new phenotype's first
pass, completing each review file as it appears):

```bash
export LDLINK_TOKEN=<your_token>          # required by stage 2

Rscript pipeline/01_load_snps.R   config/<phenotype>.yml
Rscript pipeline/02_ldlink.R      config/<phenotype>.yml
Rscript pipeline/03_ld_proxies.R  config/<phenotype>.yml
Rscript pipeline/04_blocks.R      config/<phenotype>.yml
Rscript pipeline/05_crosstraits.R config/<phenotype>.yml
Rscript pipeline/06_filter.R      config/<phenotype>.yml
Rscript pipeline/07_categorize.R  config/<phenotype>.yml
Rscript pipeline/08_enrichment.R  config/<phenotype>.yml
Rscript pipeline/09_harmonize.R   config/<phenotype>.yml
Rscript pipeline/10_extract_genes.R --blocks results/intermediate/<phenotype>/04_blocks.rds \
                                    --outdir results/intermediate/<phenotype>
```

Stage 10 takes CLI flags (not the YAML) — see `--help`.

### Stages and their outputs

| Stage | Script                  | Input                       | Output                              | Manual? |
|-------|-------------------------|-----------------------------|-------------------------------------|---------|
| 1     | `01_load_snps.R`        | studies + associations TSVs | `01_lSNPs.rds`, `01_gwassos_sigP.rds` | yes — p-value boundary |
| 2     | `02_ldlink.R`           | `01_lSNPs.rds` + LDLINK_TOKEN | `02_proxy_table.rds`, `02_lSNPs.rds` | no (writes failures report) |
| 3     | `03_ld_proxies.R`       | `02_proxy_table.rds`        | `03_leadSNPs_LDs.rds`, `03_pSNPs_per_lSNP.rds` | no |
| 4     | `04_blocks.R`           | Stage 2/3 outputs           | `04_blocks.rds`, `04_blocks.tsv`    | yes — block merges |
| 5     | `05_crosstraits.R`      | blocks + associations       | `05_crosstraits_raw.rds`            | yes — ancestry overrides |
| 6     | `06_filter.R`           | raw cross-traits + blocks   | `06_crosstraits_filtered.rds`, `06_filtering_log.tsv` | yes — GxE manual review |
| 7     | `07_categorize.R`       | filtered cross-traits + trait reference | `07_trait_table.rds`     | yes — new trait categorization |
| 8     | `08_enrichment.R`       | trait table + full catalog  | `08_trait_enrichment.rds`, `08_category_enrichment.rds` | no |
| 9     | `09_harmonize.R`        | filtered cross-traits + blocks + LDproxy + trait table + enrichment | `09_harmonized_trait_table.rds`, `09_block_directions.rds` | no |
| 10    | `10_extract_genes.R`    | `04_blocks` + biomaRt        | `protein_coding_genes_by_block.tsv`, `protein_coding_genes_collapsed.tsv`, `genes_list_complete.txt` | no |

### Manual checkpoints

Several stages run end-to-end on first invocation but produce review files you can edit and reload to
refine the analysis. The convention is the same for all of them:

- The script writes `<label>_review.tsv` to `results/review/<phenotype>/`.
- You edit it, save your changes as `<label>_reviewed.tsv` (same folder) - be sure to save the file as a .tsv file.
- On the next run, the script picks up your edits via a stable join key (`STUDY ACCESSION`, `trait_id`, etc.) that survives spreadsheet save quirks.

| Stage | Review file                          | What you decide                              |
|-------|--------------------------------------|----------------------------------------------|
| 1     | `01_pvalue_review.tsv`               | lead SNPs sitting at the (likely rounded) p-value threshold: drop or correct |
| 4     | `04_block_overlaps_review.tsv`       | which overlapping block pairs to merge (then add to the YAML's `block_merges`) |
| 5     | `05_unmapped_ancestry_review.tsv`    | superpop assignment for studies whose sample-size text the parser couldn't decode |
| 6     | `06_gxe_review.tsv`                  | for each row flagged `manual_check`: keep or auto_exclude |
| 7     | `07_trait_table_review.tsv`          | category and hierarchized label for each trait flagged `new` |

## Adding a new phenotype

```bash
cp config/leptin.yml config/<phenotype>.yml
# Edit phenotype.name, the trait terms, and self_trait_categories.
# Run the same stages pointing at the new config.
```

The trait reference (`data/trait_ref/...`) is leptin-specific; for a new phenotype, leave `trait_ref`
unset in the YAML and Stage 7 will flag every trait as `new` for review on first run, after which
you can categorize and sub-categorize based on your analysis needs.

## Design choices worth knowing about

**RDS for stage handoffs, TSV for review files.** Stage outputs use R's native binary format
(`saveRDS` / `readRDS`) because they preserve list-columns, factors, and types cleanly. Manual-review
files use TSV so they open in any text editor or spreadsheet.

**Ancestry encoding via cohort-proportion rule.** For every cross-trait row, the parser reads
`INITIAL SAMPLE SIZE` and `REPLICATION SAMPLE SIZE`, parses each into (count, ancestry-descriptor)
cohorts, maps descriptors to 1000G superpopulations (EUR/AFR/EAS/SAS/AMR/NR), and applies the rule:
if one superpop is ≥90% of a field's total, emit just that one; otherwise emit all. Then unions discovery
and replication results. The descriptor → superpop mapping is hardcoded in `R/utils.R`.

**No `setwd()` calls in any script.** All paths come from the YAML config. The pipeline assumes scripts
run from the project root.

**Defensive about catalog data.** The sample-size parser handles `99,900` (thousands separators),
`up to`, `approximately`, `at least`, and similar softeners. Catalog rows with no allele info from LDproxy
are dropped at Stage 9.

## Provenance

Built by refactoring an interactive analysis script developed for a leptin pleiotropy study (Université de
Lorraine bioinformatics internship work). The original linear script ran to ~2,000 lines in one R file; the refactor
splits it into ten stages, externalizes phenotype-specific parameters into YAML, replaces buggy parsers, and
standardizes the manual-review interface across stages.

## License

Released under the [MIT License](LICENSE).