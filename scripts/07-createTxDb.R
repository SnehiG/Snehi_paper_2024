bioconductor_pkgs = c("GenomicFeatures", "GenomicAlignments", "biomaRt")
cran_pkgs = c("tidyverse")
for (i in cran_pkgs) {
	if (!i %in% rownames(installed.packages())) {
		install.packages(i)
	}
	library(i, character.only = T)
}

if (!require("BiocManager", quietly = TRUE)){
		install.packages("BiocManager")}
for (i in bioconductor_pkgs) {
	BiocManager::install(i)	
	library(i, character.only = T)
}

# Edit the following based on what you need
biomart = "plants_mart"
prefix = "ensembl_"
host = "https://plants.ensembl.org"
dir = "."
species = c("Arabidopsis", "Poplar")
dataset = c("athaliana_eg_gene", "ptrichocarpa_eg_gene")
outputPath = paste0(dir, species , "/", species, "_TxDb.RData")

for (i in 1:2) {
	dir.create(species[i])
	TxDb <- makeTxDbFromBiomart(biomart = biomart,
				    dataset = dataset[i],
				    id_prefix = prefix,
				    host = host)
	tx <- transcriptsBy(TxDb,"gene")
	save(TxDb, tx, file = outputPath[i])
}
