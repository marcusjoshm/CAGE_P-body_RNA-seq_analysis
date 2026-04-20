# CAGE × P-body RNA-seq Analysis

Independent reproduction of a pipeline that integrates FANTOM5 CAGE-seq with Hubstenberger et al. (2017) P-body RNA-seq data to compute a **Cap Index** for inferring decapping status of P-body-enriched mRNAs. Results are cross-referenced with DepMap cancer-dependency and IMPC embryonic-lethality annotations.

See [`method_summary.md`](method_summary.md) for the full methodology.

## Repository Layout

```
notebooks/        # 01–09 Python pipeline, 10 R Markdown plots
data/
  raw/            # Downloaded source data (gitignored; re-created by notebook 01)
  processed/      # Intermediate + final outputs from notebooks 02–09
  cache/          # Scratch (gitignored)
reference/        # Jason's reference outputs to validate against
figures/          # Reproduced PDF plots
docs/             # Brainstorm + plan documents
run_pipeline.sh   # Runs notebooks 01–09 end-to-end
requirements.txt  # Python pinned dependencies
```

## Prerequisites

- **Python 3.11+** with Jupyter
- **R 4.x** with `ggplot2`, `dplyr`, `readr`, `scales` (for notebook 10)
- ~5 GB free disk space for raw downloads (FANTOM5 matrix alone is ~4 GB compressed)

## Setup

```bash
# Clone
git clone https://github.com/marcusjoshm/CAGE_P-body_RNA-seq_analysis.git
cd CAGE_P-body_RNA-seq_analysis

# Python environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Register the Jupyter kernel used by run_pipeline.sh
python -m ipykernel install --user --name cage-pbody --display-name "cage-pbody"
```

## Running the Pipeline

### Full pipeline (notebooks 01–09)

```bash
./run_pipeline.sh
```

This executes each notebook in order via `jupyter nbconvert --execute`, writing outputs back into the notebook files and producing CSVs in `data/processed/`.

### Plots (notebook 10)

Open `notebooks/10_plots.Rmd` in RStudio and knit, or from the command line:

```bash
Rscript -e "rmarkdown::render('notebooks/10_plots.Rmd')"
```

Plots are written to `figures/`.

### Running notebooks individually

```bash
jupyter lab notebooks/
```

Run notebooks **in order** — each step consumes the output of the previous one.

## Pipeline Steps

| Notebook | Purpose | Key output |
|---|---|---|
| `01_download_source_data.ipynb` | Fetch Hubstenberger, GENCODE v19, FANTOM5, DepMap, IMPC/MGI source files | `data/raw/` |
| `02_extract_protein_coding_genes.ipynb` | Filter Hub dataset to protein-coding genes | `02_protein_coding_genes.csv` |
| `03_compute_gene_lengths.ipynb` | Exonic union gene lengths from GENCODE v19 | `03_hub_with_gene_lengths.csv` |
| `04_calculate_tpm.ipynb` | TPM per replicate on full protein-coding set | `04_hub_with_tpm.csv` |
| `05_process_fantom5_cage.ipynb` | Annotate then sum CAGE peaks per gene (HEK293 untreated) | `05_fantom5_cage_per_gene.csv` |
| `06_inner_join_hub_fantom.ipynb` | Join on Ensembl ID → 12,544 genes | `06_hub_cage_merged.csv` |
| `07_calculate_cap_index.ipynb` | `Cap_Index = log10(Hub_CAGE) − log10(FANTOM_CAGE)` | `07_hub_cage_with_cap_index.csv` |
| `08_integrate_depmap_impc.ipynb` | Add DepMap + IMPC embryonic-lethality annotations | `08_hub_cage_depmap_impc.csv` |
| `09_filter_decapping_candidates.ipynb` | Filter to decapping candidates (374 genes) | `09_decapping_candidates.csv` |
| `10_plots.Rmd` | Reproduce Cap Index vs Pub_Log2FC scatter plots | `figures/*.pdf` |

## Validating Against Reference

Reference outputs from the original analysis live in `reference/`. After running the pipeline, compare `data/processed/08_hub_cage_depmap_impc.csv` against `reference/Lethality_DepMap_12544_Hub_CAGE_MERGE_with_CapIndex.csv`, and `data/processed/09_decapping_candidates.csv` against `reference/Decapping PB candidates 10threshold+DepMap_with_IMPC_EmbryonicLethal_Ortholog.csv`.

## Reproducibility Notes

- **TPM is computed on the full ~14,690 protein-coding set before the inner join** — computing after the join (12,544 genes) distorts the scaling denominator.
- **CAGE peaks are annotated with Ensembl IDs before summing**, not by genomic position. This avoids misattribution and correctly excludes mitochondrial-origin nuclear pseudogenes.
- Genome build: **hg19 / GENCODE v19** throughout.
- Log2FC uses a pseudocount of 0.01.
- 2 genes (HEATR9, BAHCC1) are dropped for lack of a gene-length match.
