# Get reference genes #
SPECIES = c("Arabidopsis", "Poplar")

#### Install and load packages ####
cran_pkgs = c("tidyverse", "devtools",  "ggpubr")
bioconductor_pkgs =  c("DESeq2","GenomicFeatures", "biomaRt")
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

# CustomSelection has to be installed from github, so it cannot be grouped with the others yet
if (! "CustomSelection" %in% installed.packages() ) {
	install_github("KarenGoncalves/CustomSelection")
}
library(CustomSelection) 

#### CustomSelection ####

customRefs = sapply(SPECIES, simplify = F, \(curSpecies) {
	countsData = 
		read.csv(file = paste0("counts/", curSpecies, ".csv"), 
			 header = T, row.names = 1)
	metadata = read.csv(paste0("metadata/DESeq_", curSpecies, ".tsv"))

	counts = countsData[, metadata$replicate]
	# get TxDb object
	load(paste0("./", curSpecies, "_TxDb.RData"))
	
	transcript.lengths = width(tx)
	
	length = sapply(names(transcript.lengths), \(geneName) {
		sort(transcript.lengths[[geneName]], decreasing = T)[1]
	})
	tpm = Counts_to_tpm(counts = counts, featureLength = length)[[1]]
	write.table(tpm, 
		    file = paste0("counts/", curSpecies, 
		    	      "ExpressionLevels_TPM.tsv"),
		    row.names = T, 
		    quote = F, sep = "\t", append = F)
	cutv = DAFS(tpm=tpm)
	
	write.table(data.frame(DAFS_result = unname(cutv),
			       Sample = names(cutv)),
		    file = paste0("counts/", curSpecies, 
		    	      "/MinimumExpression_log2TPM.tsv"),
		    row.names = T, 
		    quote = F, sep = "\t", append = F)
	
	resultCustomRefs = gene_selection(tpm, cutv)
	write.table(resultCustomRefs, 
		    file = paste0("counts/", curSpecies,
		    	      "/resultsCustomSelection.tsv"),
		    row.names = T, col.names = c("Mean_TPM", "Coefficient_Variance"),
		    quote = F, sep = "\t", append = F)
	resultCustomRefs
})
dir.create("plots")
#### Plots ####
theme_bw() %>% theme_set()
joinedInfoRefs = lapply(SPECIES, \(curSpecies) {
	data.frame(Mean = customRefs[[curSpecies]]$Mean, 
		   CoefVar = customRefs[[curSpecies]]$Covariance,
		   Species = rep(curSpecies, nrow(customRefs[[curSpecies]]))
	)	
}) %>% do.call(what = "rbind")

(Mean_plot = joinedInfoRefs %>%
		ggplot(aes(y = log2(Mean), x = Species, fill = Species)) +
		geom_boxplot(show.legend = F) +
		scale_y_continuous(breaks = seq(0, 12, 2),
				   limits = c(0, 13)) +
		ylab("Mean TPM (log2)") +
		xlab("")
)

(CoefVar_plot = joinedInfoRefs %>%
		ggplot(aes(y = CoefVar, x = Species, 
			   fill = Species)) +
		geom_boxplot(show.legend = F)	+
		scale_y_continuous(breaks = seq(0, 0.1, 0.02),
				   limits = c(0, 0.11)) +
		ylab("Coefficient of variation") +
		xlab("")
)

# join the two plots into a single figure with ggarrange from the package ggpubr
ggarrange(plotlist = list(Mean_plot, CoefVar_plot),
	  ncol = 2)
ggsave("plots/Mean_coefVar_references.pdf")

### Reference annotation

ensembl = sapply(SPECIES, simplify = F, \(species) {
	dataset = ifelse(species == "Arabidopsis",
			 "athaliana_eg_gene",
			 "ptrichocarpa_eg_gene")
	useMart(dataset = dataset,
		host = "https://plants.ensembl.org",
		biomart = "plants_mart")
})


reference_annotation <- sapply(SPECIES, simplify = F, \(species) {
	attributes = c("ensembl_gene_id",
		       "description")
	if (species == "Arabidopsis") {
		
		attributes = c(attributes,
			       "external_gene_name",
			       paste0("ptrichocarpa_eg_homolog_", 
			              c("ensembl_gene",
			                "subtype",
			                "orthology_confidence")
			       ))
	} else {
		attributes = c(attributes,
			       paste0("athaliana_eg_homolog_", 
			              c("ensembl_gene",
			                "subtype",
			                "orthology_confidence",
			                "associated_gene_name")
			       ))
	}
	
	getBM(attributes = attributes, 
	      filters = c("ensembl_gene_id"),
	      values = rownames(customRefs[[species]]), 
	      mart = ensembl[[species]])
})
