---
title: "feat: Independent validation of CAGE x P-body analysis pipeline"
type: feat
date: 2026-03-30
---

## Enhancement Summary

**Deepened on:** 2026-03-30
**Sections enhanced:** 10 (all notebook phases + infrastructure)
**Research agents used:** architecture-strategist, kieran-python-reviewer, performance-oracle, code-simplicity-reviewer, best-practices-researcher (data sources), best-practices-researcher (pyranges/GTF patterns)

### Key Improvements
1. **Concrete download URLs** for all 5 data sources, including FANTOM5 expression matrix file and DepMap gene lists
2. **pyranges code patterns** for exonic union gene lengths and CAGE peak annotation with validated gotchas
3. **Floating point tolerance strategy** — use `np.allclose(rtol=1e-5)` not exact equality when comparing Python vs R outputs
4. **GTF parsing optimization** — cache parsed exons as parquet, saves 30-90s per pipeline run
5. **PAR gene ID edge case** — Ensembl version stripping must preserve `_PAR_Y` suffixes
6. **Deferred utils.py** — don't pre-build; extract only if duplication becomes painful across 3+ notebooks

### New Considerations Discovered
- Reference CSVs have inconsistent quoting: main CSV uses double-quoted headers, Lethality CSV does not. Reference loader must handle both.
- pyranges 0.1.4 is installed (v0 line, stable). Pin this version — the v1.x rewrite has a different API.
- DepMap release version is unknown — may need to test 24Q2/24Q4 to match Jason's exact distributions.
- FANTOM5 expression matrix is already in TPM — no conversion needed, but must grep for HEK293 sample columns by library ID (e.g., `CNhs12328`).

---

# Independent Validation of CAGE x P-body Analysis Pipeline

## Overview

Build a 10-notebook pipeline (9 Python Jupyter + 1 R Markdown) that independently reproduces Jason's analysis integrating FANTOM5 CAGE-seq with Hubstenberger et al. P-body RNA-seq data. Each notebook performs one processing step, exports intermediates, and quantitatively validates against Jason's reference CSVs. Three critical methodology decisions are tested with both-ways comparisons.

## Problem Statement

A colleague (Jason) produced a merged dataset (12,544 genes) and candidate list (374 decapping P-body genes) for paper inclusion. Before publishing, we need independent verification that the processing pipeline produces the same results. The reference outputs exist but no intermediate files are available — only the final CSVs.

## Proposed Solution

A linear notebook pipeline where each step is self-contained (given its prerequisites in `data/processed/`), self-documenting, and validates its output against columns extracted from Jason's final reference CSVs.

## Technical Approach

### Directory Structure

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
│   ├── raw/              # Downloaded source files (gitignored)
│   ├── cache/            # Parsed GTF exons as parquet (gitignored)
│   └── processed/        # Intermediate outputs per step
├── reference/            # Jason's outputs (moved from root)
├── venv/                 # Python virtual environment
├── run_pipeline.sh       # Sequential notebook executor
├── requirements.txt
├── method_summary.md
└── CLAUDE.md
```

### Intermediate File Contracts

Each notebook produces a named output. Downstream notebooks reference these exact filenames:

| Notebook | Output File | Key Columns |
|----------|------------|-------------|
| 01 | `data/raw/*` (multiple files) | N/A |
| 02 | `data/processed/02_protein_coding_genes.csv` | `ensembl_id`, `Associated Gene Name`, raw count columns, gene_biotype |
| 03 | `data/processed/03_gene_lengths.csv` | `ensembl_id`, `gene_length_bp`, `rescue_ensembl_gene_id` |
| 03 | `data/cache/gencode_v19_exons.parquet` | Parsed exons for reuse in NB05 |
| 04 | `data/processed/04_hub_with_tpm.csv` | All from 02 + `PB_rep1-3_TPM`, `Cytosol_rep1-3_TPM`, `PB_TPM`, `Cytosol_TPM`, `NEW_Log2FC_TPM` |
| 05 | `data/processed/05_fantom5_cage_per_gene.csv` | `ensembl_id`, `FANTOM_Total_CAGE_TPM`, `FANTOM_n_peaks` |
| 06 | `data/processed/06_hub_cage_merged.csv` | All from 04 + all from 05 + `overlap_status` |
| 07 | `data/processed/07_hub_cage_with_cap_index.csv` | All from 06 + `Cap_index`, `log_Hub`, `log_FANTOM`, `delta_log10_Hub_minus_FANTOM` |
| 08 | `data/processed/08_hub_cage_depmap_impc.csv` | All from 07 + 8 DepMap/IMPC columns (48 total) |
| 09 | `data/processed/09_decapping_candidates.csv` | Subset of 08 columns for 374 genes |

### Pipeline Runner

A simple shell script to execute all notebooks sequentially:

```bash
#!/bin/bash
# run_pipeline.sh
set -e
for nb in notebooks/0[1-9]_*.ipynb; do
    echo "Running $nb..."
    jupyter nbconvert --to notebook --execute "$nb" --output "$(basename $nb)" \
        --ExecutePreprocessor.kernel_name=cage-pbody
done
echo "Pipeline complete. Run 10_plots.Rmd manually in RStudio."
```

### Validation Approach

**Standard two-metric check** at each step: Pearson r + max absolute difference. Gene-by-gene drill-down only when these metrics flag an issue (r < 0.999 or max diff exceeds tolerance).

**Floating point tolerance strategy** (Python vs R comparison):
- Use `numpy.allclose(rtol=1e-5, atol=1e-8)` — never use `==` for numeric validation
- TPM scaling factor involves summation over ~14,690 values; floating-point summation order differences between Python and R produce ~1e-12 level discrepancies
- If tighter agreement needed, use `math.fsum` for compensated summation in TPM scaling factor

**Reference CSV loading caveat:** Jason's CSVs have inconsistent formats — the main 12,544-gene CSV uses double-quoted headers while the Lethality CSV does not, and column ordering differs (`ensembl_id` first vs `Associated Gene Name` first). The reference loader must handle both.

**Do not pre-build a shared utils.py.** Start with inline validation in each notebook. If identical comparison logic appears in 3+ notebooks, extract at that point.

### Implementation Phases

---

#### Phase 1: Data Acquisition and Gene Filtering (Notebooks 01–02)

##### Notebook 01: Download Source Data

**Purpose:** Fetch all raw inputs programmatically. No validation step (data sourcing only).

**Data sources with concrete URLs:**

| Source | File | URL | Destination |
|--------|------|-----|-------------|
| Hubstenberger et al. 2017 | `mmc3.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S1097276517306512-mmc3.xlsx` | `data/raw/mmc3.xlsx` |
| GENCODE v19 | `gencode.v19.annotation.gtf.gz` | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz` | `data/raw/` |
| FANTOM5 CAGE expression matrix | `hg19_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz` | `https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg19_latest/extra/CAGE_peaks/hg19_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz` | `data/raw/` |
| FANTOM5 peak annotation | `hg19_fair+new_CAGE_peaks_phase1and2_ann.txt.gz` | `https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg19_latest/extra/CAGE_peaks/hg19_fair+new_CAGE_peaks_phase1and2_ann.txt.gz` | `data/raw/` |
| FANTOM5 sample table | `Human.sample_name2library_id.txt` | `https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/Human.sample_name2library_id.txt` | `data/raw/` |
| DepMap Common Essential | `AchillesCommonEssentialControls.csv` | `https://depmap.org/portal/download/all/?release=DepMap+Public+24Q2&file=AchillesCommonEssentialControls.csv` | `data/raw/` |
| DepMap Strongly Selective | `AchillesStronglySelectiveControls.csv` | `https://depmap.org/portal/download/all/?release=DepMap+Public+24Q2&file=AchillesStronglySelectiveControls.csv` | `data/raw/` |
| DepMap Gene Dependency | `CRISPRGeneDependency.csv` | `https://depmap.org/portal/download/all/?release=DepMap+Public+24Q2&file=CRISPRGeneDependency.csv` | `data/raw/` |
| IMPC Lethality | `genotype-phenotype-assertions-ALL.csv.gz` | `https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/` | `data/raw/` |
| MGI Orthologs | `HOM_MouseHumanSequence.rpt` | `https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt` | `data/raw/` |

**Research Insights:**

- **FANTOM5 data is already in TPM** in the `_tpm.osc.txt.gz` file — no conversion needed. The file is a large TSV with CAGE peak IDs as rows and sample library IDs as columns. Download the sample table first, grep for `HEK293` to identify the correct column(s) (library IDs like `CNhs12328`), then extract those columns.
- **FANTOM5 peak IDs** use the format `chr10:100013403..100013414,-` — will need parsing into chr/start/end/strand.
- **DepMap gene format** is `GENE_NAME (ENTREZ_ID)` (e.g., `RPL3 (6122)`). Will need parsing to extract gene symbol for joining.
- **DepMap release version** matters — Jason's exact release is unknown. Start with 24Q2; if distributions don't match reference (8,727/2,557/1,095/165), try 24Q4 and 25Q1.
- **IMPC FTP structure** reorganizes periodically. If the `latest` symlink breaks, check numbered release directories (e.g., `DR21.0/`).
- **Cell Press CDN URL** for mmc3.xlsx is the most reliable programmatic path but may change. Fallback: scrape supplementary links from the paper's fulltext page.

**Tasks:**
- [ ] `notebooks/01_download_source_data.ipynb` — Download all files listed above
- [ ] Download FANTOM5 sample table first, identify HEK293 library IDs, then download expression matrix
- [ ] Verify all downloads (file size, row counts)
- [ ] Save manifest of downloaded files with URLs and dates to `data/raw/MANIFEST.md`

**Output:** All raw files in `data/raw/`

##### Notebook 02: Extract Protein-Coding Genes

**Purpose:** Load Hubstenberger data, identify protein-coding genes, filter.

**Processing:**
1. Read `data/raw/mmc3.xlsx` with openpyxl/pandas
2. Inspect column names (expected: `Ensembl Gene ID`, raw read columns, CPM columns, log2FC, p-values)
3. Strip Ensembl version suffixes: `ENSG00000123456.7` → `ENSG00000123456`
4. Map Ensembl IDs to gene biotype — use GENCODE v19 GTF `gene_type` attribute (more reliable than org.Hs.eg.db for hg19 data, and independent of Jason's R approach)
5. Filter to `gene_type == "protein_coding"`
6. Report: total genes in mmc3, protein-coding count, genes lost

**Research Insights:**

- **PAR gene ID edge case:** GENCODE v19 uses `ENSG00000182378.10_PAR_Y` for pseudoautosomal region genes. Naive `split(".")` stripping loses the `_PAR_Y` disambiguation — two rows with the same base ID will silently collapse during merges. Use this pattern:
  ```python
  def strip_ensembl_version(eid: str) -> str:
      base, *rest = eid.split(".")
      par_suffix = ""
      if rest:
          par_idx = rest[0].find("_PAR_")
          if par_idx != -1:
              par_suffix = rest[0][par_idx:]
      return base + par_suffix
  ```
  Then validate uniqueness: `assert stripped.nunique() == len(stripped)`

- **Nullable integer dtypes:** When pandas reads integer count columns with NAs, it silently promotes to float64 (`42` → `42.0`). Use `pd.Int64Dtype()` (capital I) for raw count columns to preserve integer type with NA support.

- **Biotype source difference from Jason:** Jason used R's `org.Hs.eg.db` for biotype mapping. We use GENCODE v19 GTF's `gene_type` attribute. These may disagree on edge-case genes — explicitly compare our protein-coding gene list against the 12,544 Ensembl IDs in Jason's reference to catch any discrepancies.

**Validation:**
- [ ] Confirm ~14,692 protein-coding genes (Jason's count)
- [ ] Cross-check: all 12,544 Ensembl IDs from Jason's reference CSV should be present in our set (superset check)
- [ ] Identify any reference Ensembl IDs NOT in our protein-coding set
- [ ] Validate no duplicate Ensembl IDs after version stripping

**Output:** `data/processed/02_protein_coding_genes.csv`

---

#### Phase 2: Gene Lengths and TPM (Notebooks 03–04)

##### Notebook 03: Compute Gene Lengths

**Purpose:** Calculate exonic union length per gene from GENCODE v19 GTF.

**Processing:**
1. Parse `gencode.v19.annotation.gtf.gz` with `pyranges.read_gtf()`
2. Filter to `Feature == "exon"` and `gene_type == "protein_coding"`
3. Strip version suffixes from gene_id (preserving PAR suffixes)
4. Reconstruct PyRanges from filtered DataFrame
5. `merged = exons_pr.merge(by="gene_id")` — equivalent to R's `GenomicRanges::reduce()`
6. Calculate width: `End - Start` (pyranges auto-converts GTF 1-based to 0-based half-open)
7. Sum widths per gene → `gene_length_bp`
8. Cache parsed exon data as parquet for reuse in NB05

**Research Insights:**

- **Validated pyranges pattern:**
  ```python
  import pyranges as pr
  gtf = pr.read_gtf("data/raw/gencode.v19.annotation.gtf.gz")
  exons = gtf[gtf.Feature == "exon"]
  df = exons.df.copy()
  df["gene_id"] = df["gene_id"].str.replace(r"\.\d+$", "", regex=True)
  df = df[df["gene_type"] == "protein_coding"]
  exons_pr = pr.PyRanges(df)
  merged = exons_pr.merge(by="gene_id")  # merge overlapping exons PER GENE
  merged_df = merged.df.copy()
  merged_df["width"] = merged_df["End"] - merged_df["Start"]
  gene_lengths = merged_df.groupby("gene_id")["width"].sum().reset_index()
  gene_lengths.columns = ["ensembl_id", "gene_length_bp"]
  ```

- **Coordinate system already handled:** pyranges `read_gtf()` converts from GTF 1-based fully-closed to 0-based half-open (BED convention). `End - Start` gives the correct width. Do NOT add 1.

- **`merge(by="gene_id")` is critical:** Without the `by` parameter, all overlapping exons on the same chr/strand merge regardless of gene, which corrupts overlapping genes.

- **Strand-awareness is default and correct:** `merge()` operates per-strand automatically, which is right since all exons of a gene share the same strand.

- **GTF parsing is slow (30-90s):** Cache the filtered exon-only DataFrame as parquet after first parse. Downstream notebooks (NB05) can reuse it:
  ```python
  df.to_parquet("data/cache/gencode_v19_exons.parquet", index=False)
  ```

- **Memory optimization:** Drop unused GTF attribute columns immediately after parsing (transcript_type, tag, level, etc.). Convert remaining string columns like gene_type to categorical dtype to reduce memory by 40-60%.

- **Pin pyranges version** to 0.1.4 in requirements.txt. The v1.x rewrite has a different API (`merge_overlaps()` instead of `merge()`).

**Handling the 17 unmatched genes:**
- After joining gene lengths to Hubstenberger data on Ensembl ID, check for NA lengths
- For genes with missing lengths, attempt rescue by gene symbol: look up the Ensembl ID for the gene symbol in GENCODE, use that gene's exonic union length
- Jason rescued 15 of 17; 2 (HEATR9, BAHCC1) were dropped
- Document which genes we rescue and whether our rescue matches Jason's

**Validation:**
- [ ] Join our gene lengths with Jason's reference CSV on `ensembl_id`
- [ ] Compare `gene_length_bp` values using `np.allclose(rtol=1e-5)` — expect exact integer matches since both derive from the same GTF
- [ ] Specifically check the 11 rescued genes visible in reference CSV
- [ ] Confirm HEATR9 and BAHCC1 are absent from our final set

**Output:** `data/processed/03_gene_lengths.csv`, `data/cache/gencode_v19_exons.parquet`

##### Notebook 04: Calculate TPM

**Purpose:** Compute TPM from raw read counts using gene lengths. This is where Critical Question #2 (TPM order) is tested.

**Processing:**
1. Load protein-coding genes (NB02 output) with gene lengths (NB03 output)
2. Identify raw read count columns from mmc3:
   - P-body: `sorted P-body replicate 1/2/3 (mapped read number)`
   - Cytosol: `pre-sorted fraction replicate 1/2/3 (mapped read number)`
3. TPM formula — fully vectorized:
   ```python
   gene_lengths_kb = gene_lengths_bp / 1000
   rpk = counts[sample_columns].div(gene_lengths_kb, axis=0)
   scaling_factors = rpk.sum(axis=0) / 1e6
   tpm = rpk.div(scaling_factors, axis=1)
   ```
4. Calculate for all 6 replicates → `PB_rep1_TPM` through `Cytosol_rep3_TPM`
5. Averages: `PB_TPM = mean(PB_rep1-3)`, `Cytosol_TPM = mean(Cytosol_rep1-3)`
6. `NEW_Log2FC_TPM = log2((PB_TPM + 0.01) / (Cytosol_TPM + 0.01))` — pseudocount 0.01

**Research Insights:**

- **Do NOT iterate over rows or samples.** The pandas `div` with `axis` parameter handles broadcasting. This completes in under 10ms for 14,690 genes.
- **Floating-point summation order:** The TPM scaling factor involves `rpk.sum()` over ~14,690 genes. Python and R may sum in different orders, producing ~1e-12 level discrepancies. Use `rtol=1e-5` tolerance when comparing, not exact equality. If tighter agreement is needed, use `math.fsum` for compensated summation.
- **Count column NA handling:** Use `pd.Int64Dtype()` when loading to prevent silent float promotion. Convert to float only at the point of TPM calculation.

**Critical Test: TPM Before vs After Join**
- Calculate TPM on the full ~14,690 gene set (correct method)
- Also calculate TPM on just the 12,544 genes that survive the inner join
- Compare: Pearson r between the two, max/mean absolute difference, whether any gene's Log2FC sign flips
- Report quantitative impact on downstream Cap Index

**Validation:**
- [ ] Compare our 6 replicate TPMs against Jason's reference columns using `np.allclose(rtol=1e-5)`
- [ ] Compare PB_TPM, Cytosol_TPM, NEW_Log2FC_TPM
- [ ] Report correlation and discrepancies
- [ ] Document the TPM-order impact quantitatively

**Output:** `data/processed/04_hub_with_tpm.csv`

---

#### Phase 3: CAGE-seq Processing and Join (Notebooks 05–06)

##### Notebook 05: Process FANTOM5 CAGE-seq

**Purpose:** Extract HEK293 CAGE data, annotate, aggregate. This is where Critical Question #3 (annotation order) is tested.

**Processing:**
1. Load FANTOM5 expression matrix (already in TPM)
2. Parse peak coordinates from FANTOM5 peak IDs (format: `chr10:100013403..100013414,-`)
3. Extract HEK293 column(s) by library ID (identified in NB01)
4. Load protein-coding gene coordinates from cached GENCODE exon data (`data/cache/gencode_v19_exons.parquet`) — extract gene-level boundaries (Feature == "gene")
5. Annotate CAGE peaks with gene IDs using pyranges overlap
6. Filter to protein-coding genes
7. Sum CAGE TPM per gene → `FANTOM_Total_CAGE_TPM`
8. Count peaks per gene → `FANTOM_n_peaks`

**Research Insights:**

- **CAGE peak annotation pattern with pyranges:**
  ```python
  # Load gene boundaries (not exons — we want gene-level coordinates)
  genes_df = gtf_df[gtf_df["Feature"] == "gene"]
  genes_df = genes_df[genes_df["gene_type"] == "protein_coding"]
  genes_pr = pr.PyRanges(genes_df[["Chromosome", "Start", "End", "Strand", "gene_id"]])

  # Overlap CAGE peaks with genes (strand-matched)
  annotated = cage_pr.join(genes_pr, strandedness="same")

  # Sum TPM per gene (annotate-first, then sum)
  cage_per_gene = annotated.df.groupby("gene_id").agg(
      FANTOM_Total_CAGE_TPM=("tpm_column", "sum"),
      FANTOM_n_peaks=("tpm_column", "count")
  ).reset_index()
  ```

- **Peaks overlapping multiple genes:** `join()` produces one row per peak-gene overlap pair. A peak overlapping two genes gets its TPM counted for both. This is standard practice and likely what Jason did.

- **Validate row counts after join:** pyranges joins can produce cartesian products for overlapping intervals. Check `len(result) / len(peaks)` — if >>1.5x, investigate many-to-many joins.

- **MT pseudogene exclusion:** Filtering genes to `gene_type == "protein_coding"` automatically excludes nuclear-encoded mitochondrial pseudogenes (they have biotypes like `processed_pseudogene`).

- **Intergenic peaks** are automatically dropped by the inner join — peaks with no gene overlap are excluded.

**Critical Test: Annotate-First vs Sum-First**
- **Method A (correct):** Annotate peaks with Ensembl IDs → filter protein-coding → sum per gene
- **Method B (old):** Sum peaks by genomic position first → then annotate summed values with gene IDs
- Compare: gene counts, total CAGE TPM per gene, which genes differ, whether MT pseudogenes appear in Method B but not A

**Validation:**
- [ ] Compare our `FANTOM_Total_CAGE_TPM` per gene against Jason's reference column using `np.allclose(rtol=1e-5)`
- [ ] Compare `FANTOM_n_peaks`
- [ ] Report correlation, exact match count, worst discrepancies
- [ ] Document the annotate-first vs sum-first impact

**Output:** `data/processed/05_fantom5_cage_per_gene.csv`

##### Notebook 06: Inner Join Hub + FANTOM5

**Purpose:** Merge the two datasets on Ensembl ID.

**Processing:**
1. Load NB04 output (Hub with TPM, ~14,690 genes) and NB05 output (FANTOM5 CAGE per gene)
2. Inner join on `ensembl_id`
3. Add `overlap_status` column (`"overlap"` for all rows since it's an inner join)

**Note:** NB04 and NB05 have no mutual dependency — they can genuinely run in parallel, converging here.

**Validation:**
- [ ] Confirm exactly 12,544 genes after join
- [ ] Verify all 12,544 Ensembl IDs match Jason's reference
- [ ] Compare all numeric columns carried through from NB04 and NB05
- [ ] Identify any extra or missing genes vs reference

**Output:** `data/processed/06_hub_cage_merged.csv`

---

#### Phase 4: Cap Index and Annotations (Notebooks 07–08)

##### Notebook 07: Calculate Cap Index

**Purpose:** Compute Cap Index and related metrics. This is where Critical Question #1 (formula) is resolved.

**Processing:**
1. `Cap_index = FANTOM_Total_CAGE_TPM / Cytosol_TPM`
   - Handle division by zero: when Cytosol_TPM = 0, Cap_index = 0 (matching Jason's behavior)
2. Compute `log_Hub` and `log_FANTOM` columns
   - **Investigation needed:** These columns do NOT correspond to simple log10 of visible columns. From the reference data, they appear to derive from a separate "Hub CAGE" dataset or from CPM columns rather than TPM.
   - **Resolve this definition first** before computing anything — try candidate formulas against reference values:
     - `log_Hub = log10(sorted P-body replicate average CPM + 1)` ?
     - `log_FANTOM = log10(FANTOM_Total_CAGE_TPM + 1)` ?
     - Other transforms of CPM or count columns
3. `delta_log10_Hub_minus_FANTOM = log_Hub - log_FANTOM`

**Research Insight:** These `log_Hub`/`log_FANTOM` columns are NOT used in the scatter plots (which use `Cap_index` on the Y-axis). If they cannot be reproduced, document as unresolved — focus validation effort on `Cap_index` which is the biologically meaningful metric.

**Validation:**
- [ ] Gene-by-gene comparison of `Cap_index` using `np.allclose(rtol=1e-5)`
- [ ] Attempt to reproduce `log_Hub`, `log_FANTOM`, `delta_log10_Hub_minus_FANTOM`
- [ ] Document the Cap Index formula definitively

**Output:** `data/processed/07_hub_cage_with_cap_index.csv`

##### Notebook 08: Integrate DepMap + IMPC

**Purpose:** Add cancer dependency and embryonic lethality annotations.

**Processing:**
1. **DepMap integration:**
   - Load `AchillesCommonEssentialControls.csv` — parse gene names from `GENE_NAME (ENTREZ_ID)` format
   - Load `AchillesStronglySelectiveControls.csv` — same format
   - Load `CRISPRGeneDependency.csv` — count columns where dependency probability >= 0.5 per gene
   - Map to `DepMap_Common_Essential` (binary), `DepMap_Strongly_Selective` (binary)
   - Derive `DepMap_Flag_Status`: "not flagged", "common essential", "strongly selective", "both"
   - Join on gene symbol (`Associated Gene Name`)

2. **IMPC integration:**
   - Load IMPC genotype-phenotype data, filter for viability phenotype terms (MP:0011100 for preweaning lethality, embryo stage-specific terms)
   - Load MGI `HOM_MouseHumanSequence.rpt` for human-mouse ortholog mapping
   - Map human genes to mouse orthologs, flag embryonic lethal knockouts
   - `IMPC_Embryonic_Lethal_Ortholog` (binary 0/1)
   - `IMPC_Embryonic_Lethal_Class_Ortholog`: "Early lethal", "Late lethal", "Intermediate lethal"
   - `IMPC_Mouse_Source_Genes`, `IMPC_Human_Ortholog_Symbol`

**Research Insight:** DepMap release version matters — the gene lists change slightly between releases. If distributions don't match reference (not flagged=8,727; selective=2,557; essential=1,095; both=165), try DepMap 24Q4 and 25Q1.

**Validation:**
- [ ] Compare DepMap column distributions against reference
- [ ] Compare IMPC distributions: 908 embryonic lethal orthologs (422 Early, 356 Late, 130 Intermediate)
- [ ] Gene-by-gene comparison of all annotation columns (these are categorical — use exact match, not tolerance)

**Output:** `data/processed/08_hub_cage_depmap_impc.csv` — full 12,544-gene dataset with all 48 columns

---

#### Phase 5: Candidate Filtering and Visualization (Notebooks 09–10)

##### Notebook 09: Filter Decapping Candidates

**Purpose:** Apply filtering criteria to identify decapping P-body candidate genes.

**Processing (reverse-engineered thresholds):**
1. Filter: `Cap_index <= 0.1` (confirmed — max in candidates is 0.0998)
2. Filter: `Cytosol_TPM >= 10` (confirmed — filename "10threshold", min = 10.34)
3. Filter: `Pub_Log2FC > 0` (all candidates are positive; min = 2.02)
4. Additional criteria to verify:
   - Is there a `NEW_Log2FC_TPM` threshold? (min in candidates = 1.10)
   - Is there a p-value or FDR threshold?

**Threshold discovery approach:**
- Start with confirmed thresholds (Cap_index <= 0.1, Cytosol_TPM >= 10)
- Apply progressively and check gene counts
- If count doesn't match 374, test additional threshold combinations systematically
- Compare filtered set against Jason's 374 candidates gene-for-gene

**Validation:**
- [ ] Achieve exactly 374 candidates (or document why the count differs)
- [ ] Gene-for-gene comparison: every gene in our list should be in Jason's list and vice versa
- [ ] If any genes differ, report which and investigate why
- [ ] Document the exact filtering criteria used

**Output:** `data/processed/09_decapping_candidates.csv`

##### Notebook 10: Reproduce Plots (R Markdown)

**Purpose:** Recreate the three reference scatter plots in R for visual comparison.

**Plot specifications (from reference PDFs):**

All three plots share:
- X-axis: `Pub_Log2FC` (published P-body enrichment log2 fold-change)
- Y-axis: `Cap_index` on pseudo-log10 scale (accommodates zeros)
- Horizontal dashed magenta line at `Cap_index = 0.1`
- Vertical dashed black line at `Pub_Log2FC = 0`
- GAPDH labeled with black diamond (~(-4.5, 1.63))
- Gray band at bottom for zero/near-zero Cap_index values

**Plot 1:** `CapIndex_v_PubLog2FC_30bg_100thresh_CommonOnlygreen.pdf`
- Filter: `Cytosol_TPM >= 30`
- Gray = all genes; Green = Common Essential genes

**Plot 2:** `CapIndex_v_PubLog2FC_30bg_100thresh_SelectiveOnly.pdf`
- Filter: `Cytosol_TPM >= 30`
- Gray = all genes; Orchid/magenta = Strongly Selective genes

**Plot 3:** `CapIndex_v_PubLog2FC_5TPM_EmbryonicLethalRed.pdf`
- Filter: `Cytosol_TPM >= 5`
- Gray = all genes; Red = IMPC embryonic lethal orthologs

**R packages:** ggplot2, dplyr, scales (pseudo_log_trans), ggrepel, readr

**Tasks:**
- [ ] `notebooks/10_plots.Rmd` — Load OUR processed data from `data/processed/08_hub_cage_depmap_impc.csv` (not Jason's reference)
- [ ] Reproduce each plot with matching aesthetics
- [ ] Export to PDF for side-by-side comparison with Jason's reference PDFs
- [ ] Document any visual differences

**Output:** Three PDF plots in `figures/`

---

## Acceptance Criteria

### Functional Requirements

- [ ] All raw data sources downloaded automatically in notebook 01
- [ ] Protein-coding gene count within 0.5% of Jason's ~14,692
- [ ] Gene lengths match Jason's `gene_length_bp` (expect exact integer match for most genes)
- [ ] TPM values pass `np.allclose(rtol=1e-5)` against Jason's per-replicate TPMs
- [ ] FANTOM CAGE TPM passes `np.allclose(rtol=1e-5)` against Jason's `FANTOM_Total_CAGE_TPM`
- [ ] Inner join produces exactly 12,544 genes
- [ ] Cap_index passes `np.allclose(rtol=1e-5)` against Jason's values
- [ ] DepMap and IMPC annotation distributions match reference exactly
- [ ] Candidate filter produces exactly 374 genes matching Jason's list
- [ ] Three scatter plots visually match Jason's reference PDFs

### Methodology Investigation Requirements

- [ ] TPM before-vs-after join: quantitative impact documented (correlation, max diff, sign flips)
- [ ] CAGE annotate-first vs sum-first: gene count difference and affected genes documented
- [ ] Cap Index formula: `Cap_index` column formula confirmed, `log_Hub`/`log_FANTOM` columns investigated

### Quality Gates

- [ ] Each notebook runs given its prerequisites in `data/processed/`
- [ ] Each notebook includes a validation section with Pearson r + max diff
- [ ] All discrepancies documented with specific gene names and values
- [ ] `requirements.txt` captures all Python dependencies (done)
- [ ] Reference files moved to `reference/` directory
- [ ] `run_pipeline.sh` executes all notebooks sequentially

## Dependencies & Prerequisites

- [x] Python 3.12 virtual environment with all packages (`venv/`)
- [x] `requirements.txt` frozen
- [x] Jupyter kernel registered ("CAGE P-body Analysis")
- [ ] Pin pyranges==0.1.4 in requirements.txt
- [ ] R installation with ggplot2, dplyr, scales, ggrepel, readr (for notebook 10)
- [ ] Internet access for data downloads (notebook 01)
- [ ] Jason's reference CSVs moved to `reference/` directory

## Risk Analysis & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| FANTOM5 HEK293 sample ID unknown | Blocks NB05 | Download sample table first; grep for HEK293; fallback library IDs: CNhs12328/29/30 |
| DepMap release mismatch | Wrong distributions | Start with 24Q2; compare distributions; try 24Q4/25Q1 if needed |
| Gene biotype mapping differs (GENCODE vs org.Hs.eg.db) | Different gene count | Cross-check against reference 12,544 gene list; document which genes differ |
| `log_Hub`/`log_FANTOM` cannot be reproduced | Incomplete validation | Not used in scatter plots. Document as unresolved; focus on `Cap_index` |
| Python vs R floating-point differences | False validation failures | Use `np.allclose(rtol=1e-5, atol=1e-8)`; document systematic offsets |
| PAR gene ID collision after version stripping | Silent data corruption | Preserve `_PAR_Y` suffixes; validate uniqueness after stripping |
| Large GTF file (1.5GB) slow to parse | Developer friction | Cache parsed exons as parquet; drop unused attribute columns |
| pyranges join produces many-to-many explosion | Wrong CAGE TPM | Validate `len(result) / len(peaks)` ratio after every join |
| Cell Press CDN URL changes | NB01 download fails | Fallback: scrape supplementary links from paper fulltext page |

## References & Research

### Internal References
- Brainstorm: `docs/brainstorms/2026-03-30-cage-pbody-validation-brainstorm.md`
- Method summary: `method_summary.md`
- Correspondence: `email.md`
- Workflow doc: `1- Workflow and HUB dataset processing.docx`

### External Data Sources
- Hubstenberger et al. 2017: DOI 10.1016/j.molcel.2017.09.003 (`ars.els-cdn.com/.../mmc3.xlsx`)
- FANTOM5: `fantom.gsc.riken.jp/5/datafiles/reprocessed/hg19_latest/extra/CAGE_peaks/`
- GENCODE v19: `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/`
- DepMap: `depmap.org/portal/download/all/?release=DepMap+Public+24Q2`
- IMPC: `ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/`
- MGI Orthologs: `informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt`

### Key Verified Facts (from repo research)
- `Cap_index = FANTOM_Total_CAGE_TPM / Cytosol_TPM` — algebraically confirmed (GAPDH: 3238.908/1987.596 = 1.6296)
- Candidate thresholds: `Cap_index <= 0.1` AND `Cytosol_TPM >= 10`
- Reference CSV: 40 columns (main), 48 columns (with DepMap+IMPC)
- 11 rescued genes with non-NA `rescue_ensembl_gene_id` in reference
- `overlap_status` = "overlap" for all 12,544 rows
- DepMap: not flagged=8,727; selective=2,557; essential=1,095; both=165
- IMPC: 908 lethal orthologs (422 Early, 356 Late, 130 Intermediate)
