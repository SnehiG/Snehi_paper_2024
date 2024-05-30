#' loadAndPrepareDEGs_files.R
#' 
#' This script is a helper for others used after the differential expression analysis.
#' It prepares the R environment with the necessary variables, file names and packages
#' It also creates a long version of the result table from DESeq.
#' 

#### Install and load packages ####
pkgs = c("tidyverse", "DESeq2", "devtools", "GenomicFeatures", "ggpubr", "biomaRt", "xlsx")
#BiocManager::install("KEGGREST")

pkgs.To.Install = ! pkgs %in% installed.packages()

if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])

for (curPkg in c(pkgs, "KEGGREST") ) library(curPkg, character.only = T) 

theme_bw() %>% theme_set()


#### VARIABLES ####
SPECIES = c("Arabidopsis", "Poplar")
filename = sapply(SPECIES, \(curSpecies) 
		  paste0(curSpecies, '/DESeq_2023-09-28.tsv'))

designName = sapply(SPECIES, \(curSpecies) 
		    paste0(curSpecies, "/metadata_selectedReplicates.csv"))

patternContrast = c("_vs_Col_0_GFP", "_vs_POPLAR_GFP")
names(patternContrast) = SPECIES

discarded_line = "P1_LINE1_Mlp729"


#### Import datasets ####

differentially_expressed_genes = 
	sapply(simplify = F,
	       SPECIES,
	       \(species) read.delim(filename[species])
	)
design = sapply(simplify = F,
		SPECIES,
		\(species) read.csv(designName[species])
)

# Let's add columns to the differentially_expressed_genes tables
# One to show the correct plantLine name,
# Another to indicate the effector expressed
# And one to indicate the direction of the deregulation

degs_table = sapply(
	simplify = F, SPECIES, \(species) {
		constrastSuffix = patternContrast[species]
		differentially_expressed_genes[[species]] %>%
			mutate(plantLine = gsub(constrastSuffix, "", ContrastName),
			       SampleName = ifelse(grepl("521", ContrastName),
			       		    "Mlp52166", "Mlp72983"),
			       Deregulation = case_when(
			       	(log2FoldChange > 2 & padj < 0.01) ~ "Up-regulated",
			       	(log2FoldChange < -2 & padj < 0.01) ~"Down-regulated",
			       	.default = "Non-regulated")
			)

	})
