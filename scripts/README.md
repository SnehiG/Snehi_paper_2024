# Summary of scripts

## General information

- All scripts use [tidyverse](10.21105/joss.01686).
- When the folder is not specified, the file was save to the corresponding species folder. 
- For the plots, when one figure is saved per species, the name is shown as SPECIES...
- Some files are saved with the date of analysis in the file name, here the name is shown as ...DATE.OF.ANALYSIS...
- Figures exported as both pdf and tiff are shown with the extension (pdf|tiff)
- Two script files, not present in the list below, were used. These are useful for more than the current analysis and thus were made public as github gists.
  1. [getSemData.R](https://gist.github.com/KarenGoncalves/e1b37865d5f6bdeff1c8c490ba8df963): collects the GO term information for Arabidopsis for calculation of similarity between GOterms
  2. [FUNCTIONS_GOplots.R](https://gist.github.com/KarenGoncalves/95b6b81ef391d85b4fb5f8866cb448c3): creates functions that are used in the scripts 13 and 14 below.

## Scripts

1. [09-getReferenceGenes.R](./scripts/2-getReferenceGenes.R)
  - Script uses the package CustomSelection<sup>[1](10.1186/s12864-019-6426-2),[2](10.1186/s12864-021-07743-7)</sup> to:
      1. Get the expression levels in TPM for each gene and sample (output saved to ExpressionLevels_TPM.tsv)
      2. Calculate the minimum TPM value to accept a gene as expressed in each sample (shown as $`log_2(TPM)`$, output written to `MinimumEpresion_log2TPM.tsv`)
      3. Get the genes $`mean(TPM) > 2^{mean(cutoff)}`$ and selects the top .5% genes from lowest to highest coefficient of variation. Output written to `resultsCustomSelection.tsv`.
3. [10-DESeq_plantLine.R](./scripts/3-DESeq_plantLine.R)
    - Script filters the counts data frame to:
      1. Remove genes not passing the minimum TPM level threshold from the counts table. Saves output to `counts_filtered.tsv`
      2. Performs differential expression analysis using [DESeq2](10.1186/s13059-014-0550-8) using the selected reference genes (these are used to estimate teh amount of input material for each sample).
        - Alpha: 0.01; Null hypothesis: $`log_2fold.change = 0 `$; Fold change threshold: $`abs(log_2Fold.Change) > 2`$; Adjusted p-value threshold: $`padj < 0.01`$
        - Result saved in long data frame format in: `DESeq_DATE.OF.ANALYSIS.tsv`
      3. Exports PCA result calculated using the regular log transformation on the normalized read counts.
        - Three files saved: plot saved to `plots/SPECIES_PCA_DATE.OF.ANALYSIS.(pdf|tiff)`
4. [11-loadAndPrepareDEGs_Files.R](./scripts/4-loadAndPrepareDEGs_Files.R)
    - Imports the results from the previous script and adds the columns: plantLine, SampleName (the effector expresed) and Deregulation ("Up-regulated", "Down-Regulated" and "Non-regulated").
5. [12-heatmaps.R](./scripts/5-heatmaps.R)
    - Creates a wide table to calculate a gene dendrogram based on the fold change values.
    - Uses the dendrogram to order the genes in the heatmap and plots the heatmaps of Arabidopsis (on top) and Poplar (at the bottom) in the file `plots/Heatmaps_DATE.OF.ANALYSIS.pdf`
6. [13-annotateDeregulatedGenes.R](./scripts/6-annotateDeregulatedGenes.R)
    - Uses biomaRt<sup>[1](10.1093/bioinformatics/bti525),[2](10.1038/nprot.2009.97)</sup> to browse [Ensembl Plants biomarts](doi:10.1093/database/bar030) for annotation of the deregulated genes and for search of the poplar homologs of arabidopsis deregulated genes, and vice-versa.
    - Annotation includes the gene name and description as well as GO ID, GO name and GO ontology (biological process, molecular function and cellular compartment)
    - Exports the files `DEGs_annotation.txt` (genes deregulated by either effector),   `Annotation_commonDEGs.txt` (genes deregulated by both effectors) and `DEGs_homologs.txt`
7. [14-Commonly_deregulated_genes.R](./scripts/7-Commonly_deregulated_genes.R)
    - Creates the tables `(Mlp72983|Mlp52166)_deregulatedGenes.csv` and the plots 
`plots/SPECIES_scatterPlot_DEGs_DATE.OF.ANALYSIS.pdf` (each gene is a dot, its fold change in the two transgenic lines of each effector is used as its XY coordinate in the plot) and `plots/SPECIES_numberDEGs_DATE.OF.ANALYSIS.pdf`.
8. [15-Heatmap_commonDEGs.R](./scritps/8-Heatmap_commonDEGs.R)
    - Finds genes deregulated in both species and uses cluster 2.1.4 [^1] functions `daisy(metric = "gower")` and `diana(diss = T)` to get the dendrogram and order the genes in the heatmap.
    - Exports the figure `plots/Heatmap_homologs_deregulated.pdf`, which only includes genes deregulated in both species (meaning that a deregulated gene is included if one or more of its homologs in the other species is also deregulated).
    - Creates the VennDiagram `plots/SimpleVenn_DEGs.pdf`, which uses the combination of arabidopsis gene id and poplar gene id to find the number of genes deregulated in each condition and in the different combinations, regardless of the species.
9. [16-VennDiagram.R](./scripts/9-VennDiagram.R)
    - Finds the number of genes deregulated in each group (effector-species) and the number of homologs in the other species. In areas that combine the two species, shows the number of homologs also deregulated. Plot saved as `plots/Venn_diagram_deregulatedGenes_homologs.pdf`
10. [17-prepareData_GOenrichment.R](./scripts/10-prepareData_GOenrichment.R)
    - Gets the GOterm annotation for all the genes expressed in each species and exports as RData in `RData_files/GO_annotation.RData`
    - Using the gist `getSemData.R`, obtains and saves (in the file `RData_files/semantic_GO_At.RData`) the go information for arabidopsis for calculation of similarity between terms.
11. [18-get_geneSets_goEnrichment.R](./scripts/11-get_geneSets_goEnrichment.R)
    
    
  - Uses clusterProfiler<sup>[1](10.1016/j.xinn.2021.100141),[2](10.1089/omi.2011.0118)</sup> to perform GO enrichment analysis



[^1]: Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2022).  cluster: Cluster Analysis Basics and Extensions. R package version 2.1.4.
