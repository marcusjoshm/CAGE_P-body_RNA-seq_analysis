# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Goal

Independent reproduction and review of a colleague's (Jason's) analysis integrating FANTOM5 CAGE-seq with Hubstenberger et al. P-body RNA-seq data. The provided files (CSVs, PDFs, workflow doc) are **reference outputs to validate against**, not source code to extend. We are building analysis scripts from scratch to verify methodology and results for inclusion in a paper.

See `method_summary.md` for a detailed description of the processing pipeline and critical methodological decisions.

## Key Data Files

- **12544_Hub_CAGE_MERGE_with_CapIndex.csv / .rds** — Master dataset (12,544 genes). Merges P-body RNA-seq (3 replicates each for P-body and cytosol fractions, TPM-normalized), published P-body enrichment scores (Pub_Log2FC), FANTOM CAGE expression data, and a computed "Cap Index" (`delta_log10_Hub_minus_FANTOM`). Available as both CSV and R serialized (RDS) format.
- **Decapping PB candidates 10threshold+DepMap.csv** — Filtered candidate list (374 genes) with DepMap cancer dependency annotations (Common Essential, Strongly Selective, dependent cell lines).
- **Lethality_DepMap_12544_Hub_CAGE_MERGE_with_CapIndex.csv** — Full 12,544-gene dataset augmented with DepMap dependency and IMPC (International Mouse Phenotyping Consortium) embryonic lethality ortholog data.
- **PDF plots** — Scatter plots of Cap Index vs Pub_Log2FC under various filtering thresholds, with annotations for common essential genes, selective genes, and embryonic lethal orthologs.

## Key Columns and Metrics

- **PB_TPM / Cytosol_TPM**: Average TPM across 3 replicates for P-body and cytosol fractions
- **NEW_Log2FC_TPM**: Log2 fold-change of P-body vs cytosol (recalculated from TPM)
- **Pub_Log2FC**: Published P-body enrichment log2 fold-change
- **Cap_index**: `log10(Hub_CAGE) - log10(FANTOM_CAGE)` — measures cap status; negative values suggest decapping
- **overlap_status**: Whether gene overlaps between Hub and FANTOM CAGE datasets
- **DepMap annotations**: Cancer dependency classifications from DepMap (Common Essential, Strongly Selective)
- **IMPC embryonic lethality**: Mouse ortholog lethality phenotype data

## Analysis Context

The workflow (described in `1- Workflow and HUB dataset processing.docx`) covers:
1. Hub CAGE dataset processing and merging with FANTOM CAGE data
2. Cap Index calculation to infer decapping status of P-body-enriched mRNAs
3. Integration with DepMap and IMPC lethality data to assess functional importance
4. Visualization of Cap Index vs P-body enrichment with various gene set highlights

## Working with This Data

- The `.rds` file can be loaded in R with `readRDS("12544_Hub_CAGE_MERGE_with_CapIndex.rds")`
- CSV files use standard comma separation with quoted headers
- Gene identifiers use Ensembl IDs (`ensembl_id`) and HGNC symbols (`Associated Gene Name`)
