Hi Guys,

See attached for the detailed workflow including scripts (filenames might not match). Also, add the final spreadsheet in .csv and .rds - ill add the buildup files to Box later the weekend.

Today, I re-processed the data from scratch for documentation. Two procedural changes: 
Ensembl IDs were immediately applied to CAGE-seq data. This did two things:
Summing peaks based on genomic position can lead misappropriating peaks to the wrong gene. Annotating peaks beforehand ensures that peaks are attributed to the correct gene.
Excludes the nuclear pseudo-genes that were “stolen” from mitochondria as we evolved. If these pseudo-genes are near a promoter then they will be capped and picked up by CAGE. 
Now these MT-encoded genes are only picked up in Hub’s RNA-seq data and not CAGE-seq (rightly so).
Gene lengths and TPM were added and calculated, respectively, in the Hubstenberger dataset before joining:
TPM should be calculated with the as much of the protein-coding transcript population (14,690) for each replicate. Since CAGE-seq and Hub-seq have non-overlapping populations, calculating TPM after joining (12,544 transcripts) would result in a different TPM. 
Consequences: Accuracy. The negative correlation is still blatantly obvious! The sensitivity curve with ratios is beautiful!

Cheers,
Jason


Hi Guys,

Attached is the full spreadsheet from the previous email but fused with columns from DepMap identifying essential and selective genes. Also, attached are Decapped PB candidate genes to choose from. Two different scatterplots highlighting the distribution of common essential genes and strongly selective genes. Also, stacked box plots for Low-Medium and High abundance genes across four PB partitioning bins.

For a biological theme (like GO), I thought evaluating the distribution of essential and critical genes might be useful. 

Common essential (for most cell lines) is statistically slightly de-enriched from P-bodies. Highly selective (which means essential for a significant number of cell lines) are positively enriched in P-bodies. 

We can select 5-6 transcripts that low CAGE/high abundance transcripts and selective/common essential. 

Jason


Last addition to the spreadsheets. Embryonic lethal genes from mice have been linked onto our main list and decapping candidate list (both attached). Also, a scatterplot with these lethal KO genes in red. 

It seems strongly selective genes are enriched in PBs with a subtle enrichment of embryonic lethal genes. Common essential genes are slightly depleted from PBs.