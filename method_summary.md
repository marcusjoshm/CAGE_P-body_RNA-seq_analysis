# Method and Data Processing Summary

## Objective

Assess the cap status of mRNAs enriched/depleted from P-bodies by integrating FANTOM5 CAGE-seq data with P-body/cytosol RNA-seq data from Hubstenberger et al. (2017). A "Cap Index" is computed to infer whether P-body-enriched transcripts are preferentially decapped.

## Data Sources

1. **Hubstenberger et al. 2017 P-body RNA-seq** — Supplementary file `mmc3.xlsx` from [Molecular Cell paper](https://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30651-2). Contains raw mapped read counts for 3 P-body replicates and 3 pre-sorted (cytosol) replicates, aligned to hg19.
2. **FANTOM5 CAGE-seq** — HEK293 untreated CAGE-seq data from the FANTOM5 mega-dataset. Readout is in TPM. Aligned to hg19.
3. **GENCODE v19 GTF** — `gencode.v19.annotation.gtf` used to calculate exonic union gene lengths (hg19-compatible).
4. **DepMap** — Cancer dependency data classifying genes as Common Essential or Strongly Selective, with dependent cell line counts.
5. **IMPC** — International Mouse Phenotyping Consortium embryonic lethality data linked via human-mouse ortholog mapping.

## Processing Pipeline

### Step 1: Extract protein-coding genes from Hubstenberger dataset
- Load `mmc3.xlsx`, strip Ensembl ID version suffixes (e.g., `ENSG00000123456.7` → `ENSG00000123456`)
- Map Ensembl IDs to gene biotype using `org.Hs.eg.db` (R AnnotationDbi)
- Filter to protein-coding genes only (~14,692 genes)

### Step 2: Compute exonic union gene lengths from GENCODE v19
- Parse `gencode.v19.annotation.gtf` with `GenomicFeatures::makeTxDbFromGFF()`
- Group exons by gene (`exonsBy(txdb, by = "gene")`)
- For each gene: merge overlapping exons (`reduce()`), then sum widths → **exonic union length**
- Strip version suffixes from GTF gene IDs before joining
- 17 genes initially failed to match; 15 were rescued by gene symbol lookup; 2 dropped (HEATR9, BAHCC1)

### Step 3: Calculate TPM for each replicate
- TPM formula per replicate:
  1. RPK = raw_count / (gene_length_bp / 1000)
  2. Scaling factor = sum(RPK) / 1,000,000
  3. TPM = RPK / scaling_factor
- **Critical detail**: TPM is calculated on the full protein-coding population (~14,690 genes) *before* joining with CAGE data. This preserves the correct denominator; calculating after the inner join (12,544 genes) would distort TPM values.
- Six replicate TPMs computed: PB_rep1-3_TPM, Cytosol_rep1-3_TPM
- Averages: PB_TPM = mean(PB reps), Cytosol_TPM = mean(Cytosol reps)
- NEW_Log2FC_TPM = log2((PB_TPM + 0.01) / (Cytosol_TPM + 0.01)) — pseudocount of 0.01

### Step 4: Process FANTOM5 CAGE-seq data
- Extract HEK293 untreated samples from FANTOM5 mega-dataset
- **Critical change**: Annotate CAGE peaks with Ensembl IDs *before* summing peaks per gene (not after). This prevents misattribution of peaks to wrong genes when summing by genomic position alone.
- This also excludes nuclear pseudogenes of mitochondrial origin — these pseudo-genes near promoters get capped and picked up by CAGE but should not be counted. MT-encoded genes should only appear in the RNA-seq data.
- Filter to protein-coding genes
- Sum all CAGE peaks per gene → `FANTOM_Total_CAGE_TPM`

### Step 5: Inner join Hubstenberger and FANTOM5 datasets
- Join on Ensembl gene ID
- Result: 12,544 genes present in both datasets
- The `overlap_status` column tracks whether a gene was found in Hub, FANTOM, or both

### Step 6: Calculate Cap Index
- The Cap Index quantifies cap status relative to transcript abundance:
  - `Cap_Index = FANTOM_Total_CAGE_TPM / Cytosol_TPM` (ratio form in final spreadsheet)
  - Log-transformed version: `Cap_Index = log10(Hub_CAGE) - log10(FANTOM_CAGE)` = `delta_log10_Hub_minus_FANTOM` (used in scatter plots)
- **Interpretation**: A low Cap Index (negative log-delta) suggests decapping — the transcript is abundant in the cytosol but has low CAGE signal, meaning fewer capped 5' ends.

### Step 7: Integrate DepMap annotations
- Map genes to DepMap classifications: Common Essential, Strongly Selective, or both
- Add dependent cell line counts and flag status
- **Finding**: Common essential genes are slightly de-enriched from P-bodies; strongly selective genes are positively enriched in P-bodies.

### Step 8: Integrate IMPC embryonic lethality
- Link human genes to mouse orthologs with IMPC embryonic lethality phenotypes
- Classify as Early Lethal or other lethality class
- **Finding**: Subtle enrichment of embryonic lethal gene orthologs in P-bodies.

### Step 9: Filter decapping P-body candidates
- Selection criteria for candidates (374 genes):
  - P-body enriched (positive Log2FC)
  - Low Cap Index (suggesting decapping)
  - Thresholds applied: appears to use a minimum abundance/TPM cutoff and a Cap Index threshold
- Final candidate list annotated with DepMap and IMPC data

## Key Outputs

| File | Genes | Description |
|------|-------|-------------|
| `12544_Hub_CAGE_MERGE_with_CapIndex.csv/.rds` | 12,544 | Full merged dataset with TPM, Log2FC, CAGE, and Cap Index |
| `Decapping PB candidates 10threshold+DepMap.csv` | 374 | Filtered decapping candidates with DepMap annotations |
| `Lethality_DepMap_...with_CapIndex.csv` | 12,544 | Full dataset + DepMap + IMPC lethality |
| `Decapping PB candidates...IMPC_EmbryonicLethal_Ortholog.csv` | 374 | Candidates + DepMap + IMPC lethality |

## Scatter Plot Visualizations
- **CapIndex_v_PubLog2FC_30bg_100thresh_CommonOnlygreen.pdf** — Cap Index vs published Log2FC, common essential genes highlighted in green, with 30 background / 100 threshold filtering
- **CapIndex_v_PubLog2FC_30bg_100thresh_SelectiveOnly.pdf** — Same axes, strongly selective genes highlighted
- **CapIndex_v_PubLog2FC_5TPM_EmbryonicLethalRed.pdf** — Same axes with 5 TPM threshold, embryonic lethal orthologs in red

## Notes for Reproducibility

1. **TPM calculation order matters**: TPM must be computed on the full ~14,690 protein-coding gene set per replicate, *before* the inner join reduces the set to 12,544. This affects the scaling factor denominator.
2. **Ensembl ID annotation before peak summing**: CAGE peaks must be annotated with Ensembl IDs before aggregation to avoid misattribution. This is a procedural correction from an earlier version.
3. **Pseudocount**: A pseudocount of 0.01 is used in the Log2FC calculation to handle zeros.
4. **Gene length**: Exonic union length (merged non-overlapping exon spans per gene), not transcript length or genomic span.
5. **2 genes dropped**: HEATR9 and BAHCC1 could not be matched to gene lengths and were omitted.
6. **Genome build**: All data aligned to hg19 / GENCODE v19.
