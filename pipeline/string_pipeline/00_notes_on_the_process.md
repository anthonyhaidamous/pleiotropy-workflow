# STRING pipeline notes

**The goal** 

I have a set of genes that came out of the LD blocks of the pleiotropy pipeline - genes that happen to sit near main-trait-associated lead SNPs. The idea is to put them into STRING and get a protein interaction network, then cluster it. The natural next question is: *do these clusters mean something biological? Is module <#> a real functional unit, a pathway that leptin genetics (for example) converges on?*

Before claiming that, we have to rule out the boring explanation. So I asked: when STRING draws an edge between two of my genes, are those two genes from the same LD block or from different ones? If my modules were real biology, I would expect plenty of edges reaching *across* blocks - different locim converging on shared function. That's the interesting picture.

---

## The process

### Stage 00: Download STRING bulk files

The network is queried live. Stage 00 pulls the STRING release as bulk files - the full per-channel link table, the protein into file, and the alias map - pinned to a specific species and version (v12.0, 9606 - homo sapiens), and records the URL, byte count, and md5 of each in a manifest. Nothing downstream ever talks to the REST API. That choice buys three things: no rate limits, no risk of the network shifting under me between one call and the next, and a source file I can archive alongside the eventual paper. The sstage is idempotent - a file already on disk whose checksum matches the manifest is left untouched - and it's keyed by species and release rather than by phenotype, so leptin and CLL and anything else pinned to the same version share one download.

### Stage 01: Freezing the network

Stage 01 takes the phenotype's protein-coding gene list and extracts the subgraph STRING induces on it, at a declared release and a declared combined-score cutoff. Every symbol is resolved to a STRING protein identifier through the preferred-name table first, then the alias table; symbols that map ambiguously or not at all are written out by name rather than quietly dropped, so the gene list that entered the network can be reconciled against the one that left it. The bulk links file is then scanned in chunks, retaining only edges whose two endpoints are both in the mapped set and whose combined score clears the threshold, and STRING's reciprocal A-B/B-A duplication is collapsed. The result is written once. An existing edge file makes the script refuse to run; `--force` moves the old file aside rather than deleting it, and the accompanuing provenance file records the md5 of the edges - a fingerprint that lets anyone re-hash the file later and confirm the network behind a given result is bit-for-bit the one that was declared.

In the process of building the tables, **the per-channel scores are retained**. This makes it possible to later filter for, say, textmining-only edges (those with zero score in all other channels) and so test how robust the modules are to the exclusion of literature-supported connections. So it is worth keeping, and is kept. 

> ### The outputs
> 
> **The TSV.** It's basically the graph without its geometry. Shows the interactions, with their scores per-channel
> 
> **The JSON.** Nothing describes the network itself except counts. Read it as four blocks:
> 
> - **What was produced** - `edge_file`, `edge_md5`, `n_edges`. The output and its fingerprint.
> - **What it was built from** - `bulk_files` (each with its own url, bytes, md5, timestamp), `gene_list_file` + `gene_list_md5`, `config_file`. Every input, hashed.
> - **What parameters were applied** - `string_version`, `network_type`, `species_taxon`, `score_threshold`. 
> - **What ran it** - `r_version`, `package_versions`, `script`, `retrieved_on`.
>
> MD5 is a fingerprint that takes input of **any size** and returns 128 bits - 16 bytes, written as 32 hex characters. The 111-edge TSV (leptin) is a few kilobytes; the links file is 161 MB. Both produce 32 characters. Information is destroyed by contruction, so nothing can be recovered from a has. It is not a seed, not compression, not an encoding. It cannot regenerate, position, or control anything. 
>
> Its only function is comparison. Someone runs `md5sum data/string/leptin_edges_v12-0_s400.tsv`, gets `fbe4ac46...`, checks it against the JSON. Match means the file is byte identical to the one that produced the results. Any difference - one score changed, one edge added, a line ending converted - produces a completely unrelated hash. That's the whole mechanism: cheap, verifiable, no trust required.
>
> An md5 is ***derived*** from an input: it is computed from the butes themselves by a fixed procedure, so it can't be attached to the wrong file - anyone can recompute it from the file and check. The file carries its own idetifier implicitly.

### Stage 02: Building the graph and accounting for every gene

Stage 02 verifies the edge file against the checksum (md5) recorded in its provenance sidecar, then constructs the igraph object. The graph is built over all mapped proteins rather than over the edge list, because a graph induced from edges alone contains only genes that have at least one partner - the 23 genes with no interaction above the threshold would vanish bedore anything counted them, in the case of leptin. Building over the full mapped set and adding edges afterwards keeps them present as isolated vertices, so they can be reported rather than lost.

The stage then resolves the submitted gene list into six disjoint quantities: submitted, mapped to a unique STRING protein, carrying at least one interaction, in the largest connected component, mapped but isolated, and in a smaller component. Every gene is named in the accounting table, including those that leave. Which node set continues downstream - the largest connected component, or every gene of degree ≥ 1 - is set in the configuration file rather than fixed in the code, so the choice is declared and changeable rather than implicit.

### Stage 03: Module detection

The downstream graph is partitioned by Markov clustering, using van Dongen's redference `mcl` binary rather than a reimplementation. MCL was chosen for three properties:

- It is deterministic - the algorithm carries no random component, so a run reproduces exactly without seeding, and the module assignments in the paper are not conditional on a seed anyone has to record
- It infers the number of clusters from the graph rather than requiring one to be specified in advance
- It partitions by local density, which is the property a protein module is supposed to have, rather than by path distance. 

The exact input handed to the binary and the exact output it returned are both written to disk, alongside a provenance record carrying the mcl version, the inflation parameter, and the input checksums.

Each resulting module is then decomposed by LD block. A module drawing more than a declared fraction of its members from a single block is flagged as a positional artifact. That threshold is read from the configuration file and fixed before the modules are seen.

> ### MCL
>
> **MCL computes what random walks would do on average. It never performs one.** The mechanism is matrix algebra. You take the adjacency amtrix and normalise each column to sum to 1, which tunrs edge weights into ***transition probabilities*** - from gene A, the chance of stepping to each neighbour is "$norm(weight_{edge A})$". That matrix *is* a description of a random walk. Then two operations alternate:
> - **Expansion** - square the matrix. Entry (i, j) of $M^2$ is the probability of getting from j to i in exactly two steps, summed over every possible intermediate. Not a sample path: the exact total over all paths at once.
> - **Inflation** - raise every entry to the power $r$, then renormalise. This sharpens contrast amplifying already-likely transitions and suppressing unlikely ones.
>
> Iterate until the matrix stops changing, and it converges to a form where the surviving structure reads off directly as clusters.
>
> So the random walk is *interpretive*. It tells you what the numbers mean, and it's why the algorithm works: a walker starting inside a dense region tends to stay there, because more apths lead back **in** than **out**. Inflation exaggerates that tendency until the regions separate. But the computation is over the full probability distribution, not over sampled trajectories, and no RNG is ever consulted.
>
> In short: you can either do a random walk, or compute one. MCL does the second.

---

## The run

### stage_01_fetch_edges.R

```
======================================================================
Stage string-01 — Freeze the induced STRING edge list
======================================================================
Phenotype:        leptin  (slug: leptin)
STRING version:   12.0
Species taxon:    9606
Network type:     functional
Score threshold:  400   (combined_score >= this)
Edge file:        data/string/leptin_edges_v12-0_s400.tsv
Error: the frozen edge file already exists:
    data/string/leptin_edges_v12-0_s400.tsv
It is the reproducibility anchor for every downstream stage and will
not be regenerated silently. Re-run with --force only if you intend to
redefine the network (the existing file is moved aside, not deleted).
Execution halted
```

Execution was halted due to the existence of a file already.

### stage_02_components.R

```
======================================================================
Stage string-02 — Graph construction and gene accounting
======================================================================
Edge file:        data/string/leptin_edges_v12-0_s400.tsv
Downstream set:   lcc
  edge file verified against its provenance checksum
  block table:    results/intermediate/leptin/protein_coding_genes_by_block.tsv  (13 blocks)

Graph: 88 vertices, 111 edges

======================================================================
GENE ACCOUNTING
======================================================================
  n_in          105   gene symbols submitted
  n_mapped       88   resolved to a unique STRING protein
  n_deg_ge1      65   mapped and carrying at least one association
  n_lcc          61   in the largest connected component
  n_singleton    23   mapped but with no association at this threshold
  n_island        4   in a connected component other than the largest

UNMAPPED — no STRING protein (17)

ZERO-DEGREE — mapped, no association at this threshold (23) - disconnected nodes

NON-LCC COMPONENTS (4 genes) - disconnected modules (2 modules of size 2)

LD BLOCK COMPOSITION of the 65 connected genes:
    B2      28  (43.1%)
    B8      12  (18.5%)
    B10      9  (13.8%)
    B15      4  ( 6.2%)
    B5       3  ( 4.6%)
    B1       2  ( 3.1%)
    B7       2  ( 3.1%)
    B13      1  ( 1.5%)
    B14      1  ( 1.5%)
    B4       1  ( 1.5%)
    B6       1  ( 1.5%)
    B9       1  ( 1.5%)

Downstream graph ('lcc'): 61 vertices, 109 edges
  results/intermediate/leptin/string_modules/string_02_graph_full.rds
  results/intermediate/leptin/string_modules/string_02_graph_downstream.rds
  results/intermediate/leptin/string_modules/string_02_accounting.tsv
  results/intermediate/leptin/string_modules/string_02_accounting_summary.tsv
```

### stage_03_modules.R

```
======================================================================
Stage string-03 — MCL module detection and block composition
======================================================================
mcl:              /usr/bin/mcl  (mcl 22-282)
inflation:        1.5
graph:            61 genes, 109 associations
artifact rule:    >60% from one LD block

$ mcl 'results/intermediate/leptin/string_modules/string_03_mcl_input.abc' --abc -I 1.5 -o 'results/intermediate/leptin/string_modules/string_03_mcl_raw.txt'

mcl returned 9 modules covering 61 genes

==========================================================================
MODULES  (inflation 1.5)
==========================================================================

M1   n=14   1 block(s)  B2(14)  <-- POSITIONAL ARTIFACT
     CCDC121, EIF2B4, FNDC4, GPN1, GTF3C2, IFT172, KRTCAP3, NRBP1,
     PPM1G, SLC4A1AP, SNX17, SUPT7L, ZNF512, ZNF513
M2   n=13   4 block(s)  B8(10) B10(1) B2(1) B7(1)  <-- POSITIONAL ARTIFACT
     AGK, CAD, SFXN1, ARF5, FSCN3, GCC1, IMPDH1, LRRC4, METTL2B, PAX4,
     PRRT4, RBM28, SND1
M3   n=8    5 block(s)  B5(3) B2(2) B14(1) B8(1) B9(1)
     FTO, MPV17, UCN, COBLL1, GRB14, SLC38A11, LEP, KLF14
M4   n=7    4 block(s)  B10(4) B13(1) B4(1) B8(1)
     TAS2R3, TAS2R38, TAS2R4, TAS2R5, NRXN3, LRP1B, GRM8
M5   n=6    2 block(s)  B2(4) B10(2)  <-- POSITIONAL ARTIFACT
     CLEC5A, MGAM, BABAM2, MRPL33, RBKS, TCF23
M6   n=4    1 block(s)  B2(4)  <-- POSITIONAL ARTIFACT
     ABHD1, ATRAID, CGREF1, PREB
M7   n=4    1 block(s)  B15(4)  <-- POSITIONAL ARTIFACT
     ADA2, HDHD5, TMEM121B, XKR3
M8   n=3    1 block(s)  B2(3)  <-- POSITIONAL ARTIFACT
     EMILIN1, GCKR, KHK
M9   n=2    1 block(s)  B1(2)  <-- POSITIONAL ARTIFACT
     CEP170, SDCCAG8
--------------------------------------------------------------------------
9 modules | 7 flagged as positional artifact (>60% one block) | 46 genes in them
2 modules survive for downstream hypotheses, covering 15 genes

BLOCK x MODULE
 block  1  2 3 4 5 6 7 8 9 total
    B2 14  1 2 0 4 4 0 3 0    28
    B8  0 10 1 1 0 0 0 0 0    12
   B10  0  1 0 4 2 0 0 0 0     7
   B15  0  0 0 0 0 0 4 0 0     4
    B5  0  0 3 0 0 0 0 0 0     3
    B1  0  0 0 0 0 0 0 0 2     2
   B13  0  0 0 1 0 0 0 0 0     1
   B14  0  0 1 0 0 0 0 0 0     1
    B4  0  0 0 1 0 0 0 0 0     1
    B7  0  1 0 0 0 0 0 0 0     1
    B9  0  0 1 0 0 0 0 0 0     1

  results/intermediate/leptin/string_modules/string_03_modules.tsv
  results/intermediate/leptin/string_modules/string_03_module_summary.tsv
  results/intermediate/leptin/string_modules/string_03_block_by_module.tsv
```

## The results

### The result is M3

Setting aside the "2 modules survive" framing for a moment. M3 is not just an unflagged module, it's the specific thing your analysis was designed to find and had so far failed to produce: 8 genes drawn from 5 different LD blocks - FTO, MPV17, UCN, COBLL1, GRB14, SLC38A11, LEP, KLF14. 

FTO, GRB14, COBLL1, and KLF14 are canonical adiposity and insulin-resistance loci and LEP itself is in there. That is genuinely cross-locus convergence: independent genomic regions whose genes interact functionally. Position cannot explain it, because the genes aren't co-located. Study bias is a weaker objection than usual too - these are famous loci, but they're famous *separately*, from different GWAS traditions.

One module out of nine is not a disappointing yield here. It's the correct yield, and now you can say why the other seven aren't: they're named, counted, and mechanistically explained.

### M4 is fragile

B10 contributes 4 of 7 genes - 57.1%, just under the 60% cutoff. And the four are TAS2R3, TAS2R38, TAS2R4, and TAS2R5: a tandem-duplicated bitter receptor cluster on chr7. That's the textbook case of the confound your rule exists to catch, and it escaped on a roudning margin. The ramining three (NRXN3, LRP1B, GRM8) are large neuronal genes, which is exactly thr profile where text-mining edges dominate.

Two things follow. 

- I need to rerun with the artifact rule at 55%, 60% and then 65% and report M4's status at each - if it flips at 55%, the threshold is doing arbitrary work.
- Check M4's edges with text-mining excluded; the per-channel subscores are already in the frozen TSV, so this is a filter, not a rebuild. If M4 collapses without text-mining and M3 doesnt, that's a strong sentence against the literature-heavy artifact category, and a viable result eventually

### Still outstanding

Inflation is fixed at 1.5 and the sweep hasn't run. M3's membership is what evertyhign now rests on, so its stability across I = 1.2-4.0 is no longer a housekeeping item - it's a robustness check on the only positive finding. If LEP, FTO, GRB14, COBLL1 and KLF14 stay together from 1.2 to 4.0, that's amuch stronger claim than "at I = 1.5 we observed."

23 of 88 mapped genes are singletons at s = 400 - 26%. That's high enough that a reviewer will ask what happens at 150 or 700. **Worth exploring.**

Note: don't read the "7 out of 9 flagged" as the analysis failing. The direction-permutation test doesn't need modules to survive at all; it runs on the frozen network with everything held fixed. The module result is the motivation section. M3 is a bonus that gives something concrete to point at while making that argument.