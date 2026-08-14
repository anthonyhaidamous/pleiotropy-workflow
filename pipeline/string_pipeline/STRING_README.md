# STRING module detection — what this is and how it works

This folder takes the genes that came out of your LD blocks, looks them up in the
STRING database, builds a network out of them, and cuts that network into modules.
Then it asks whether those modules are real biology or just an echo of where the
genes sit on the chromosome.

Everything is driven by the `string:` block in `config/<phenotype>.yml`. Nothing is
hardcoded, nothing is phenotype-specific. Run it on CLL tomorrow by pointing it at
`config/cll.yml`.

---

## 1. The idea, in plain terms

You have 15 LD blocks around leptin-associated lead SNPs. Stage 10 of the main
pipeline pulled every protein-coding gene sitting in those regions — 105 genes.

Those 105 genes were selected **because of where they are**, not because they do
anything together. They're neighbours of a GWAS hit. That's the whole problem this
stage group exists to deal with.

The question is: *do these genes actually work together, or are they just
neighbours?*

To answer it:

```
LD blocks  →  genes in those blocks  →  look them up in STRING  →  network
                                                                     ↓
                                             cut the network into modules (MCL)
                                                                     ↓
                                    check: is each module one block, or many?
```

If a module is 100% one LD block, it isn't telling you about biology. It's telling
you those genes are next to each other on a chromosome, which you already knew. We
call that a **positional artifact**.

If a module pulls genes from four or five *different* blocks, that's interesting —
separate loci that converge on the same function.

---

## 2. Concepts you need (and only these)

### STRING

A database that scores how likely two proteins are to be **functionally
associated**. Not "physically touching" — functionally associated. Two proteins get
a high score if they're in the same pathway, co-expressed, mentioned together in
papers, conserved together across species, and so on.

STRING combines **seven channels** of evidence:

| channel | what it means |
|---|---|
| `neighborhood` | the genes sit together in bacterial genomes |
| `fusion` | the two genes are fused into one gene in some other species |
| `cooccurence` | the genes appear/disappear together across species trees |
| `homology` | the proteins are similar in sequence (i.e. **paralogs**) |
| `coexpression` | their mRNA levels rise and fall together |
| `experiments` | physical assays (Y2H, pulldowns, …) |
| `database` | someone curated them into the same pathway |
| `textmining` | they're named together in the literature |

(`*_transferred` versions mean the evidence came from another species and was
mapped across.)

The **combined score** merges all channels into one number from 0 to 1000. It's a
probability estimate ×1000: **400 means "about 40% confident this association is
real."**

### Why the score threshold matters so much

`score_threshold: 400` is not a display setting. It *defines your dataset*. At 400
you get 111 edges; at 700 you'd get far fewer and a different network. Every
downstream number — module count, block composition, everything — moves when it
moves.

That's why the code **refuses to run** if the threshold isn't in the config. There
is no default. A silently defaulted threshold would silently redefine your results.

### Reciprocal edges

STRING stores every association twice: once as `A→B` and once as `B→A`. The second
row carries no new information — it's just how the file is laid out. If you don't
remove the duplicates you double every edge count and every density calculation.
Stage 01 removes them (222 rows → 111 real associations).

### Modules, not pathways

A **module** is a densely-connected clump of nodes in your graph. That's it. It's a
property of the network you built.

A **pathway** is a curated biological object someone published. Do not call your
modules pathways — you'll get taken apart for it in a viva. Same for "complex."
Stick to **module**.

### MCL (Markov Clustering)

The algorithm that finds the modules. The intuition:

> Imagine a random walker on the network. If you start walking from a node, you'll
> tend to stay inside a densely-connected clump for a while before escaping to
> another one. MCL simulates that, then repeatedly exaggerates the difference
> between likely and unlikely paths until the clumps fall apart from each other.

Two operations, repeated until it stops changing:

- **expansion** — square the matrix (simulate walking further)
- **inflation** — raise every value to a power, then renormalise (exaggerate the
  strong paths, kill the weak ones)

Why MCL and not k-means or hierarchical clustering:

1. **It's deterministic.** No random seed. Same input → same output, forever. That
   matters enormously for reproducibility.
2. **You don't tell it how many clusters to find.** k-means makes you guess k.
3. **It partitions by density**, which is what a "module" actually means.

The **inflation** parameter (`-I`) is the one knob. Higher inflation = more,
smaller modules. Lower = fewer, bigger ones. We use **1.5** because that's what
your earlier manual run used. This value has not been validated — see §8.

We use van Dongen's original `mcl` binary (`sudo apt-get install -y mcl`), not the
CRAN `MCL` package.

### The positional artifact rule

**Any module with more than 60% of its members from a single LD block is flagged as
a positional artifact and excluded from downstream hypothesis testing.**

The 60% was fixed *before* looking at any results — that's deliberate and it
matters. If you pick the threshold after seeing the numbers, you'll unconsciously
pick the one that saves the module you like. Setting it in advance is what makes it
a real test instead of a rationalisation.

It lives in the config as `modules.max_single_block_fraction`.

---

## 3. The five stages

Run them in order, from the **repository root** (never `cd` into this folder — every
path is relative to the repo root).

### Stage 00 — download the STRING data

```bash
Rscript pipeline/string_pipeline/stage_string_00_download_bulk.R config/leptin.yml
```

Downloads three files (~180 MB total) into `data/string/bulk/`:

| file | what's in it |
|---|---|
| `9606.protein.links.full.v12.0.txt.gz` | every association in human, with all channel scores |
| `9606.protein.info.v12.0.txt.gz` | STRING protein ID → display name |
| `9606.protein.aliases.v12.0.txt.gz` | every other name each protein goes by |

`9606` is the NCBI taxonomy ID for human. `v12.0` is the pinned STRING release.

**Why bulk download instead of the web API:** no rate limits, no risk of STRING
updating mid-analysis and giving you different answers on different days, and the
source file itself is something you can archive.

It writes `bulk_manifest.json` recording each file's URL, size, MD5 checksum and
download date. Re-running skips files that are already there and match their
checksum — it won't re-download 160 MB.

**Version pinning:** the download URL is *derived* from the pinned version:
`https://version-12-0.string-db.org/download/...`. That address keeps serving v12.0
forever, even after STRING releases v13. If we pointed at the generic host, your
analysis would silently change the day they update.

### Stage 01 — freeze the network

```bash
Rscript pipeline/string_pipeline/stage_string_01_fetch_edges.R config/leptin.yml
```

This is the most important stage. It:

1. Reads your gene list (`genes_list_complete.txt`, 105 genes)
2. Maps each gene symbol to a STRING protein ID
3. Scans all 13.7 million human associations, keeping only those where **both**
   ends are in your gene list **and** the score ≥ 400
4. Removes the reciprocal duplicates
5. Writes the result and never touches it again

**Output — `data/string/leptin_edges_v12-0_s400.tsv`.** The filename encodes
everything that defines it: phenotype, STRING version, score threshold.

This file is the **reproducibility anchor**. Every module, every statistic, every
figure downstream is derived from it. So:

> **If it already exists, the script refuses to run.**
> Re-run with `--force` only if you genuinely mean to redefine the network — and
> even then it renames the old file to `*.superseded_<timestamp>` rather than
> deleting it.

It also writes a **provenance sidecar** (`*.provenance.json`) recording the STRING
version, the source URLs, the MD5 of every bulk file, the MD5 of your gene list, the
MD5 of the edge file itself, package versions, and the date. If someone asks in two
years how you made this network, that file is the answer.

And a **mapping table** (`leptin_gene_mapping_v12-0.tsv`) — one row per submitted
gene, saying whether it mapped and how.

**Gene mapping results for leptin:** 88 of 105 mapped. The 17 that didn't are all
pseudogenes, lncRNAs and miRNAs (`LINC-PINT`, `MIR129-1`, `OR9A1P`, …). That's
correct behaviour — STRING is a *protein* database, and those aren't proteins.

### Stage 02 — build the graph and account for every gene

```bash
Rscript pipeline/string_pipeline/stage_string_02_components.R config/leptin.yml
```

Turns the edge list into an igraph object and produces the full accounting.

The important design point: the graph includes **all 88 mapped proteins, even the
ones with no edges**. If you built the graph from the edge list alone, the
zero-degree genes would vanish before you could count them. The whole point of this
stage is that nothing disappears silently.

**The accounting (leptin):**

```
n_in         105   gene symbols submitted
n_mapped      88   resolved to a unique STRING protein
n_deg_ge1     65   mapped and carrying at least one association
n_lcc         61   in the largest connected component
n_singleton   23   mapped but no association at score >= 400
n_island       4   in a smaller separate component
```

Some vocabulary:

- **degree** = how many edges a node has. Degree 0 = connected to nothing.
- **connected component** = an island of nodes all reachable from each other.
- **LCC** = largest connected component, the biggest island. Here 61 genes.
- The other two islands are `CCNL1–DRD1` and `SSBP1–WEE2`, 2 genes each.

Every single one of those 23 zero-degree genes and 4 island genes is **named
individually** in `string_02_accounting.tsv`. Nothing leaves the analysis
anonymously.

Whether the next stage uses the LCC or all degree≥1 genes is the config setting
`components.downstream_node_set` (currently `lcc`, so 61 genes).

This stage also attaches each gene's **LD block** to the graph as a node attribute,
which is what makes the block figures and the block × module table possible.

### Stage 03 — find the modules

```bash
Rscript pipeline/string_pipeline/stage_string_03_modules.R config/leptin.yml
```

Writes the graph out in `.abc` format (`geneA  geneB  weight`), calls `mcl`, reads
the clusters back, works out each module's block composition, and applies the 60%
rule.

Edge weights are `combined_score / 1000`, so a 900-score edge weighs 0.9.

It saves the *exact* input handed to mcl (`string_03_mcl_input.abc`) and the *exact*
output (`string_03_mcl_raw.txt`), plus a provenance file with the mcl version and
the checksums of both. You can re-run the clustering by hand from those two files.

### Stage 04 — is any of that stable?

```bash
Rscript pipeline/string_pipeline/stage_string_04_stability.R config/leptin.yml
```

Stage 03 clusters at one inflation and applies one artifact threshold, and both
numbers were carried over from an earlier manual run. Stage 04 re-clusters the
**identical** graph across a grid of inflation values and reports what happens.
It defines no dataset, writes nothing stage 03 owns, and is safe to re-run.

Four questions:

**The sweep.** How many modules, how many survive the artifact rule, and what
does the partition score, at each inflation. Weighted modularity is reported but
is *descriptive only* — it has a resolution limit and will favour whatever
granularity the graph's size happens to suit. Don't use it to pick an inflation.

**Persistence.** For each stage-03 module, the fraction of its gene *pairs* still
co-clustered at each inflation. Pairs rather than a yes/no, because a module that
sheds one gene and a module that shatters are not the same event. Reported
alongside `host_size`: at low inflation everything merges, and a module can be
"intact" only by being swallowed. `absorbed` marks exactly that.

**The 60% rule.** Every module's verdict at 55 / 60 / 65%, with any that change
verdict marked `<-- FLIPS`. A threshold doing arbitrary work should be visible
in the output, not discovered by a reviewer.

**Stable cores.** The one that matters. Within each module, gene pairs whose
consensus across the grid clears a level are linked and the connected components
taken. Each component is a set of genes that *travels together*. A module can
look cross-block at one inflation purely by gluing several within-block sub-units
together — this is where that shows up, and for M3 it is the entire result (§7).

Two integrity checks run first, and either failing stops the stage: the mcl input
built here must be byte-identical to `string_03_mcl_input.abc`, and the sweep's
partition **at** the reference inflation must reproduce `string_03_modules.tsv`
exactly. mcl carries no RNG, so it must; if it doesn't, stage 03's outputs are
stale and everything below would be describing a different partition.

An inflation that returns a single cluster covering the whole graph is excluded
from the consensus — it separates nothing, and would otherwise lift every pair by
1/|grid|. It is named on the console and in the provenance, and still appears in
the sweep and persistence tables. Excluded by a rule that fires and is logged.

---

## 4. Making the figures

```bash
Rscript figure_scripts/fig_string_network_by_block.R  config/leptin.yml
Rscript figure_scripts/fig_string_network_by_module.R config/leptin.yml
```

Run the block one **first** — it computes the layout (where each node sits on the
page) and freezes it to `string_02_layout.tsv`. The module script reuses those
coordinates, so the two pictures line up node for node and you can flip between
them.

**Why freeze the layout instead of using `set.seed()`:** graph layout is random, and
`set.seed()` only guarantees the same result within the same igraph version. igraph
changes its RNG use between releases. Writing the actual x/y coordinates to disk is
the only freeze that survives. Use `--relayout` if you want a different arrangement.

Four figures in `figures/leptin/string/`:

- `string_network_by_block.png` — the network, genes coloured by LD block
- `string_network_block_facets.png` — one small panel per block, same layout
- `string_network_modules.png` — surviving modules coloured, artifacts greyed out
- `string_network_module_facets.png` — one panel per module

**Why only 3 colours?** On a network, any two nodes can end up next to each other,
so every pair of colours has to be distinguishable — including for colourblind
readers. Only about three categorical hues clear that bar reliably. Twelve blocks in
twelve colours would be unreadable. That's what the facet figures are for: they name
every block and module individually, one per panel.

---

## 5. Every output file

### `data/string/` — the frozen network (shared, phenotype in filename)

| file | what it is |
|---|---|
| `leptin_edges_v12-0_s400.tsv` | **the anchor.** 111 rows, one per association, all channel scores kept |
| `leptin_edges_v12-0_s400.provenance.json` | how it was made, with checksums |
| `leptin_gene_mapping_v12-0.tsv` | all 105 genes, mapped or not, and why |
| `bulk/` | the three downloaded STRING files + `bulk_manifest.json` |

### `results/intermediate/leptin/string_modules/` — per-phenotype results

| file | what it is |
|---|---|
| `string_02_accounting.tsv` | **one row per submitted gene**, with status, degree, component, block |
| `string_02_accounting_summary.tsv` | the six headline counts |
| `string_02_graph_full.rds` | igraph object, all 88 mapped genes |
| `string_02_graph_downstream.rds` | igraph object, the 61 LCC genes |
| `string_02_layout.tsv` | frozen x/y coordinates for the figures |
| `string_03_modules.tsv` | **one row per gene**: which module, which block |
| `string_03_module_summary.tsv` | **one row per module**: size, blocks, artifact flag, member list |
| `string_03_block_by_module.tsv` | the block × module contingency table |
| `string_03_mcl_input.abc` | exactly what was fed to mcl |
| `string_03_mcl_raw.txt` | exactly what mcl returned |
| `string_03_provenance.json` | mcl version, inflation, checksums |
| `string_04_stability_sweep.tsv` | **one row per inflation**: modules, survivors at each artifact threshold, modularity |
| `string_04_stability_persistence.tsv` | one row per (module × inflation): parts, pair cohesion, host size, absorbed |
| `string_04_stability_modules.tsv` | one row per stage-03 module: cohesion, intact fraction, rule verdict, whether it flips |
| `string_04_stability_artifact_rule.tsv` | one row per (module × artifact threshold): keep/drop |
| `string_04_stability_cores.tsv` | **one row per (module × consensus level × core)**: the sub-units that hold together, and whether each is cross-block |
| `string_04_stability_gene_pairs.tsv` | consensus co-clustering for every gene pair, with both blocks |
| `string_04_stability_provenance.json` | grid, levels, mcl version, per-inflation output checksums |
| `string_04_sweep/` | the exact mcl output at each inflation, plus the input |

### Columns in `string_02_accounting.tsv`

| column | meaning |
|---|---|
| `input_symbol` | the gene as it appeared in your list |
| `string_id` | the STRING protein it resolved to |
| `status` | `lcc` / `island` / `singleton` / `unmapped` / `ambiguous` / `duplicate_protein` |
| `degree` | number of associations |
| `component_id`, `component_size` | which island it's on and how big |
| `in_lcc` | is it in the largest component |
| `block` | its LD block |
| `match_rank` | 1 = exact name match, 2 = matched via alias, 3–4 = case-insensitive |

`ambiguous` = the symbol matched more than one STRING protein, so we refused to
guess. `duplicate_protein` = two of your symbols resolved to the *same* protein
(e.g. a gene and its old name), so the second was dropped and named. Neither
occurred for leptin, but they're checked.

---

## 6. The config block, key by key

```yaml
string:
  version:       "12.0"        # STRING release. Pinned. Appears in every filename.
  species_taxon: 9606          # human
  network_type:  "functional"  # functional | physical
  score_threshold: 400         # REQUIRED. no default. defines the dataset.

  bulk:
    dir: "data/string/bulk"
    source_url_base: null      # null = derive the version-stable host from the pin
    fallback_url_bases: [...]  # tried in order if the primary is down
    download_timeout: 3600
    verify_live_version: true  # compare pin to live release, warn only

  input:
    gene_list:  "genes_list_complete.txt"           # under results/intermediate/<slug>/
    block_table: "protein_coding_genes_by_block.tsv" # supplies each gene's LD block

  output:
    edges_dir:    "data/string"
    stage_subdir: "string_modules"

  components:
    downstream_node_set: "lcc"   # lcc | degree_ge1

  mcl:
    binary: "mcl"
    inflation: 1.5               # the granularity knob
    weight_scale: 1000           # combined_score / this = edge weight

  modules:
    max_single_block_fraction: 0.60   # >60% one block = positional artifact

  stability:                          # stage 04. no defaults in code.
    inflation_grid:        [1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0]
    artifact_fractions:    [0.55, 0.60, 0.65]
    core_consensus_levels: [0.5, 0.6, 0.7, 0.8]
    sweep_subdir: "string_04_sweep"
```

`functional` vs `physical`: functional includes all seven channels (pathway
membership, co-expression, literature). Physical only includes evidence the proteins
actually touch. We use functional — the question is whether these genes *work
together*, not whether they bind.

---

## 7. What the leptin results actually say

**The network:** 105 genes in, 88 recognised by STRING, 65 with at least one
association, 111 associations between them.

**The blocks are all over the network.** 89 of the 111 edges (80%) connect two genes
from the *same* LD block. If the network were blind to genomic position you'd expect
about 24%. So the network is heavily positional.

But that's mostly one block. **B2 alone contributes 62 of those 89 edges** — it has
28 of the 65 connected genes. Its internal *density* is actually the lowest of any
multi-gene block (0.16 vs 0.20 for B8, 0.22 for B10). B2 doesn't dominate because
its genes are unusually connected; it dominates because there are simply a lot of
them. Strip B2 out and the within-block share drops from 80% to ~64%.

Every one of the 12 connected blocks has at least one edge leaving it. It is *not*
true that blocks don't talk to each other.

**The modules:** MCL at inflation 1.5 gives **9 modules**. Seven fail the 60% rule.

| module | n | blocks | verdict |
|---|---|---|---|
| M1 | 14 | B2 ×14 | artifact (100% B2) |
| M2 | 13 | B8 ×10, +3 others | artifact (77% B8) |
| **M3** | **8** | B5 ×3, B2 ×2, B14, B8, B9 | **survives** (38%) |
| **M4** | **7** | B10 ×4, B13, B4, B8 | **survives** (57%) |
| M5 | 6 | B2 ×4, B10 ×2 | artifact |
| M6 | 4 | B2 ×4 | artifact |
| M7 | 4 | B15 ×4 | artifact |
| M8 | 3 | B2 ×3 | artifact |
| M9 | 2 | B1 ×2 | artifact |

**Two modules survive, covering 15 of 61 genes.**

**M3 = FTO, LEP, KLF14, COBLL1, GRB14, SLC38A11, MPV17, UCN.** Five different LD
blocks, no block above 38%. This is the leptin/adiposity/type-2-diabetes axis
assembled from five separate loci. If anything here is real biology, it's M3.

**M4 = TAS2R3/4/5/38, GRM8, NRXN3, LRP1B.** Four blocks, clears the rule at 57%.
Treat it carefully: TAS2R and GRM are both GPCR families, so this may be sequence
homology (paralogy) dressed up as functional association. Check the `homology`
column on those edges before relying on it.

**Note what MCL did to B2:** it didn't recover B2 as one module. It shattered it
across six (14 / 4 / 4 / 3 / 2 / 1) — and each fragment is still B2-pure. So the
clustering isn't recovering the block as a unit; it's carving a gene-dense
chromosomal region into sub-clumps that are each still positional.

---

## 7b. What the stability sweep says (stage 04)

Run on leptin over I = 1.2 … 4.0. The headline: **the module count is entirely a
function of inflation — 1 module at 1.2, 25 at 4.0 — and M3 as an eight-gene
module exists only at I = 1.4 and 1.5.** Its pair cohesion is 0.42 and it is
intact at 2 of the 9 non-reference values. It is the *least* stable of the nine.

That sounds fatal and isn't, because the cores decomposition says what M3
actually is. M3 is three sub-units that MCL glues together at low inflation:

| core | blocks | consensus |
|---|---|---|
| COBLL1, GRB14, SLC38A11 | B5 ×3 | 1.00 — **one block** |
| **FTO, KLF14, LEP** | B14, B9, B8 | ≥0.67 — **three blocks** |
| MPV17, UCN | B2 ×2 | 1.00 — **one block** |

So M3's cross-block character at I = 1.5 was partly manufactured by the merge.
Two of its three sub-units are within-block, i.e. exactly the positional artifact
the 60% rule exists to catch — the rule passed M3 at 37.5% because the mixture
diluted each block below the threshold, not because the module is positionally
clean.

**What survives is FTO + KLF14 + LEP.** Three genes, three different LD blocks,
co-clustered in ≥67% of informative runs, with FTO–KLF14 at 1.00 across the whole
grid. That is a smaller claim than "an eight-gene leptin/adiposity/T2D axis" but
it is a claim the sweep supports rather than one it happens not to have tested.
COBLL1/GRB14/SLC38A11 remain interesting biology — they are simply not *evidence
of cross-locus convergence*, because they are one locus.

**M4 flips.** It is dropped at 55%, kept at 60% and 65%. Its status is a property
of the threshold, not of the network, and it should be reported that way. Its
stable core is GRM8 + TAS2R3/4/5/38 (cohesion holds to ≥0.80) — which is the
paralogy concern from §7 surviving intact, now with GRM8 attached to the bitter-
receptor cluster. The channel-filtering check is still the right next test.

**Two cross-block pairs sit inside *flagged* modules**: AGK–SFXN1 (B10/B7) and
CAD–IMPDH1 (B2/B8), both at consensus ≥0.80 in M2. The artifact rule discards
these wholesale because it operates on whole modules. Worth knowing before
quoting "the rule removes 75% of the network" as if nothing were lost.

**Modularity peaks at I = 1.4 (0.629), barely above 1.5 (0.628).** Do not read
that as validating 1.5. Modularity has a resolution limit and prefers whatever
granularity a 61-node graph suits; it is reported to be looked at, not obeyed.

---

## 8. Things that will bite you

**The 60% rule kills 75% of your network.** That's the honest result, not a bug. It
was predicted by the 80% within-block edge figure before clustering ever ran.

**Inflation 1.5 is unvalidated as a choice** — stage 04 has now swept it, so the
consequences are documented (§7b), but nothing in the sweep argues *for* 1.5. It
remains the value the earlier manual run used. Quote the stable cores, which hold
across the grid, rather than the I = 1.5 module boundaries, which do not.

**Traits cannot be mapped to genes.** Only to blocks. The genes were picked because
they *overlap* the block region (± 10 kb) — purely positional. Lead SNPs sit in
introns and intergenic space, so no single gene owns the trait. `MAPPED_GENE` in the
GWAS Catalog looks like a gene-level link but it only names the lead SNP's nearest
neighbour (GCKR for all 1,883 B2 rows) — it's a block label wearing a gene's name.
Any trait analysis has to be at block level.

**Only 8 of 15 blocks have cross-traits at all**, and B2 holds 824 of them. Anything
trait-related will be underpowered and B2-driven.

**`data/` is in `.gitignore`.** The frozen edge file — the reproducibility anchor —
will not be committed. For a public release you need to add `!data/string/` to
`.gitignore`, or the repo ships without the one file everything depends on.

**Line endings.** A gene list saved from a Windows editor has `\r\n`. Read naively,
`ABHD1` becomes `ABHD1\r`, fails to match anything in STRING, and gets reported as a
genuinely unmappable gene — a wrong answer that looks like a real result. The code
strips `\r` explicitly. Don't remove that.

---

## 9. What isn't built yet

- ~~**Inflation sweep**~~ — **done**, stage 04. See §7b.
- **Channel robustness** — re-run modules on the network with `textmining`
  excluded, and with `homology` excluded. Both are filters on the frozen TSV's
  per-channel columns, not a rebuild. This is now the top item: it is the
  outstanding test on M4's GRM8 + TAS2R core, and the one check that could
  strengthen or kill FTO–KLF14–LEP.
- **Null model** — is 80% within-block actually more than a degree-preserving random
  graph would give? The crude 24% baseline ignores degree.
- **Score-threshold sensitivity** — 23 of 88 mapped genes (26%) are singletons at
  400. A reviewer will ask about 150 and 700. Note this one *does* rebuild the
  frozen edge file, so it needs versioned side-anchors rather than a re-freeze.
- **Enrichment** — functional annotation of the stable cores (not the I = 1.5
  modules), offline against `protein.enrichment.terms` rather than a drifting
  web service
- **Trait/category overlay** at block level
- **FSI augmentation** — adding interactors to each module separately
- **Direction layer / OmniPath** — signed edges
- **H1–H7 hypothesis tests**

The gate is no longer the sweep — it has run. It is now the channel filtering,
because M4's surviving core is a GPCR paralog family and the sweep cannot tell
sequence homology from functional association. Anything built on the module
definitions should be built on the stable cores (§7b), not on the I = 1.5
boundaries.

---

## 10. Running the whole thing from scratch

```bash
# once per STRING version — downloads ~180 MB
Rscript pipeline/string_pipeline/stage_string_00_download_bulk.R config/leptin.yml

# freeze the network (refuses if the edge file already exists)
Rscript pipeline/string_pipeline/stage_string_01_fetch_edges.R  config/leptin.yml

# graph + accounting
Rscript pipeline/string_pipeline/stage_string_02_components.R   config/leptin.yml

# modules (needs: sudo apt-get install -y mcl)
Rscript pipeline/string_pipeline/stage_string_03_modules.R      config/leptin.yml

# stability: inflation sweep, artifact-rule sensitivity, stable cores
# re-analyses the same graph; safe to re-run, touches nothing stage 03 owns
Rscript pipeline/string_pipeline/stage_string_04_stability.R    config/leptin.yml

# figures — block script first, it freezes the layout
Rscript figure_scripts/fig_string_network_by_block.R  config/leptin.yml
Rscript figure_scripts/fig_string_network_by_module.R config/leptin.yml
```

**Dependencies:** R 4.5.0, `yaml` 2.3.12, `readr` 2.2.0, `dplyr` 1.2.1, `tibble`
3.3.0, `jsonlite` 2.0.0, `curl` 7.0.0, `igraph` 2.3.3, and the `mcl` binary
(22-282). Checksums use base R's `tools::md5sum()`.

There's no `renv.lock` yet, so those versions are *declared*, not enforced. Worth
adding before public release.
