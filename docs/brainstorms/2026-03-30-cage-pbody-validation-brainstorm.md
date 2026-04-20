# Brainstorm: Independent Validation of CAGE x P-body Analysis

**Date:** 2026-03-30
**Status:** Reviewed

## What We're Building

An independent reproduction of Jason's 9-step analysis pipeline that integrates FANTOM5 CAGE-seq with Hubstenberger et al. P-body RNA-seq data. The goal is methodology validation for paper inclusion — if our independent analysis produces the same results, we have confidence in the interpretation. If not, we identify and resolve discrepancies.

### Scope

- **In scope:** Reproduce all 9 processing steps, quantitatively compare against Jason's reference outputs at each step, critically evaluate three flagged methodological decisions, reproduce Jason's scatter plots in R
- **Out of scope:** Novel analyses beyond what Jason produced, production pipeline engineering, alternative biological interpretations

## Why This Approach

### Language Choice: Python + R hybrid
- **Python (Jupyter notebooks):** Data downloading, filtering, Ensembl ID matching, TPM calculation, dataset merging, unit conversions, quantitative validation
- **R (Rmd):** Final scatter plot reproduction and any plot modifications, to match Jason's plotting code and allow direct visual comparison
- **Rationale:** Python is better suited for the data processing validation (pandas for transparent step-by-step transforms, clear numeric comparisons). R is necessary for reproducing Jason's exact plot aesthetics and for modifying his plotting scripts.

### Project Structure: Step-numbered Jupyter notebooks with inline validation
```
CAGE_P-body_RNA-seq_analysis/
├── notebooks/
│   ├── 01_download_source_data.ipynb
│   ├── 02_extract_protein_coding_genes.ipynb
│   ├── 03_compute_gene_lengths.ipynb
│   ├── 04_calculate_tpm.ipynb
│   ├── 05_process_fantom5_cage.ipynb
│   ├── 06_inner_join_hub_fantom.ipynb
│   ├── 07_calculate_cap_index.ipynb
│   ├── 08_integrate_depmap_impc.ipynb
│   ├── 09_filter_decapping_candidates.ipynb
│   └── 10_plots.Rmd
├── data/
│   ├── raw/          # Downloaded source files
│   └── processed/    # Intermediate outputs from each step
├── reference/        # Jason's outputs (moved here for clarity)
├── docs/
│   └── brainstorms/
├── method_summary.md
└── CLAUDE.md
```

Each notebook:
1. Documents what it does and why
2. Downloads or loads required input data
3. Performs the processing step
4. Exports intermediate output to `data/processed/`
5. Compares output quantitatively against Jason's reference (gene counts, value correlations, exact matches)
6. Documents findings and any discrepancies

### Validation Strategy: Quantitative diffs at every step
- Gene-by-gene comparison where applicable
- Correlation coefficients and summary statistics
- Explicit flagging of any gene or value that differs
- This is the most rigorous option — appropriate because the purpose is methodology review

**Important constraint:** Jason's intermediate files are not in this repo — only his final merged CSVs (12,544 genes and 374 candidates). For early pipeline steps (protein-coding extraction, gene lengths, CAGE processing), validation must extract reference values from columns in Jason's final CSV (e.g., `gene_length_bp`, `FANTOM_Total_CAGE_TPM`, replicate TPMs). Full end-to-end comparison is only possible at the join step (notebook 06) and beyond.

## Key Decisions

1. **Python for processing, R for plotting** — ensures independent reproduction (different tooling from Jason's R pipeline) while allowing exact plot reproduction in R
2. **Scripts download data automatically** — full reproducibility from scratch; each notebook fetches its own source data
3. **Inline validation** — each notebook validates its own output against Jason's reference before proceeding. Findings are co-located with the code that produced them.
4. **Step-by-step planning** — this brainstorm covers the full project; each notebook will be planned individually via `/workflows:plan`

## Critical Methodology Questions to Investigate

### 1. Cap Index Formula Inconsistency
- **Issue:** There appear to be two distinct metrics conflated under "Cap Index":
  - `Cap_index` column in the CSV: appears to be `FANTOM_Total_CAGE_TPM / Cytosol_TPM` (a ratio of CAGE signal to cytosolic RNA abundance)
  - `delta_log10_Hub_minus_FANTOM` column: `log_Hub - log_FANTOM` (difference between two CAGE measures, not involving Cytosol TPM at all)
- **Question:** These are not mathematically equivalent — they use different denominators and measure different things. Which one is the biologically meaningful "Cap Index" used for interpretation? Does the workflow doc's formula (step 11: `Total CAGE TPM / Cytosol TPM`) match what's in the `Cap_index` column?
- **Validation:** Recompute both formulas from the raw columns (`FANTOM_Total_CAGE_TPM`, `Cytosol_TPM`, `log_Hub`, `log_FANTOM`) and verify which matches each CSV column exactly.

### 2. TPM Calculation Order (Before vs After Join)
- **Issue:** Jason flagged this as a procedural correction — TPM should be calculated on the full ~14,690 protein-coding gene set, not the 12,544 post-join set, because the scaling factor denominator changes.
- **Question:** How much does this actually matter? What's the quantitative impact on TPM values and downstream Log2FC?
- **Validation:** Calculate TPM both ways (before and after join) and compare. Report the magnitude of the difference and whether it changes any biological conclusions.

### 3. CAGE Peak Annotation Order
- **Issue:** Jason annotated CAGE peaks with Ensembl IDs *before* summing peaks per gene, correcting an earlier approach that summed by genomic position first.
- **Question:** How many genes are affected? Does the MT pseudogene exclusion change gene counts or Cap Index values meaningfully?
- **Validation:** Process CAGE data both ways (annotate-then-sum vs sum-then-annotate) and compare gene counts and CAGE TPM values.

## Notebook-by-Notebook Outline

### 01: Download Source Data
- Fetch: Hubstenberger mmc3.xlsx, GENCODE v19 GTF, FANTOM5 CAGE HEK293, DepMap files, IMPC lethality data
- Verify file integrity (checksums if available)
- No validation against Jason (data sourcing step)

### 02: Extract Protein-Coding Genes
- Load mmc3.xlsx, strip Ensembl version suffixes
- Map to gene biotype, filter protein-coding
- **Validate:** Compare gene count (~14,692); cross-check gene list against the 12,544 Ensembl IDs in Jason's final CSV (our set should be a superset)

### 03: Compute Gene Lengths
- Parse GENCODE v19 GTF
- Calculate exonic union length per gene (merge overlapping exons, sum widths)
- Handle the 17 genes that failed initial matching (15 rescued by symbol, 2 dropped)
- **Validate:** Compare computed gene lengths against the `gene_length_bp` column in Jason's final 12,544-gene CSV (join on Ensembl ID)

### 04: Calculate TPM
- Implement TPM formula from raw read counts + gene lengths
- Calculate per-replicate TPMs for all 6 replicates
- Calculate PB_TPM, Cytosol_TPM averages and NEW_Log2FC_TPM
- **Critical test:** Calculate TPM both before and after join, quantify difference
- **Validate:** Compare against Jason's TPM columns

### 05: Process FANTOM5 CAGE-seq
- Extract HEK293 untreated data
- Annotate peaks with Ensembl IDs, filter protein-coding, sum per gene
- **Critical test:** Compare annotate-first vs sum-first approaches
- **Validate:** Compare per-gene CAGE TPM sums against the `FANTOM_Total_CAGE_TPM` column in Jason's final CSV

### 06: Inner Join Hub + FANTOM5
- Join on Ensembl ID
- **Validate:** Confirm 12,544 genes, check overlap_status distribution

### 07: Calculate Cap Index
- Compute both ratio and log-delta forms
- **Critical test:** Verify which formula matches Jason's `Cap_index` column
- **Validate:** Gene-by-gene comparison of Cap Index values

### 08: Integrate DepMap + IMPC
- Add DepMap essential/selective classifications
- Add IMPC embryonic lethality orthologs
- **Validate:** Compare against Jason's augmented CSV

### 09: Filter Decapping Candidates
- Apply filtering criteria (need to determine exact thresholds from Jason's data)
- **Validate:** Compare 374 candidates gene-for-gene

### 10: Reproduce Plots (Rmd)
- Recreate the three scatter plots from Jason's PDFs
- Visual comparison against reference PDFs

## Resolved Questions

- **Language:** Python for processing, R for plotting
- **Data sourcing:** Automatic download in scripts
- **Validation rigor:** Quantitative diffs at every step
- **Project structure:** Notebook-based, one per pipeline step, inline validation
- **Approach:** Single pipeline with co-located validation (Approach A)
