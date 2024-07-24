# DEA #
SPECIES = c("Arabidopsis", "Poplar")

#### Install and load packages ####
pkgs = c("tidyverse", "DESeq2", "GenomicFeatures", "ggpubr")
pkgs.To.Install = ! pkgs %in% installed.packages()
if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
for (curPkg in pkgs) library(curPkg, character.only = T) 

theme_set(theme_bw())

#### Filter weakly expressed genes ####
### Input files ###
for (i in 1:2) {
	countsData = read.csv(paste0(SPECIES[i], "/counts.csv"), row.names = 1)
	metadata = read.csv(paste0(SPECIES[i], "/metadata_selectedReplicates.csv"))
	selectedReps = metadata$replicate[ metadata$Keep_replicate ]
	countsData = countsData[, selectedReps]
	metadata = metadata[metadata$Keep_replicate,]
	
	tpm = read.delim(paste0(SPECIES[i], "/ExpressionLevels_TPM.tsv"))
	cutv = read.delim(paste0(SPECIES[i], "/MinimumExpression_log2TPM.tsv"))
	resultsCustomSelection = read.delim(paste0(SPECIES[i], "/resultsCustomSelection.tsv")) %>% rownames
	
	
	# The line below calculates the mean tpm for each gene
	# then selects the genes that have mean tpm > the cutoff value 
	counts_filtered = countsData[apply(tpm, 1, mean) > (mean(cutv$DAFS_result)^2),]
	write.table(counts_filtered,
		    file = paste0(SPECIES[i], "/counts_filtered.tsv"), 
		    quote = F, row.names = T, col.names = T,
		    sep = "\t", append = F)

	references = rownames(counts_filtered) %in% resultsCustomSelection
	
	#### Differential expression analysis ####
	# Prepare the design table forcing the control to be the first level
	sampleGroups = c("Control", unique(metadata$sampleName)) %>% unique
	design = sapply(names(metadata), \(column) {
		currentColData = sapply(sampleGroups, \(groupName) {
			metadata[metadata$sampleName == groupName, column]
		}) %>% unlist %>% unname
		factor(currentColData, levels = unique(currentColData))
	}, simplify = F) %>% as.data.frame
	
	
	# Order the columns from the counts object 
	counts_filtered = counts_filtered[, levels(design$replicate)]
	
	# Prepare the dds object
	dds = DESeqDataSetFromMatrix(countData = counts_filtered, 
				     colData = design, design = ~ plantLine)
	
	dds = estimateSizeFactors(dds, controlGenes = which(references))
	dds = DESeq(dds)
	result_dds =
		lapply(resultsNames(dds)[-1], \(contrast) {
			columnName = gsub("plantLine_", "", contrast)
			result_dseq = results(dds, name = contrast)@listData
			data.frame(ContrastName = columnName,
				   Gene = rownames(dds),
				   log2FoldChange = result_dseq$log2FoldChange,
				   padj = result_dseq$padj)
		}) %>% do.call(what = "rbind")
	
	result_dds_filtered = 
		result_dds %>%
		filter(abs(log2FoldChange) > 2 & padj < 0.01)
	write.table(result_dds_filtered,
		    file = paste0(SPECIES[i], '/DESeq_', Sys.Date(), ".tsv"),
		    row.names = F, col.names=T, sep = '\t', quote=F)
	
	#### Sample similarity PCA ####
	rlog = rlogTransformation(dds)
	pca_data = plotPCA(rlog, 
			   intgroup = c('sampleName', "plantLine")
	)
	pca_data$data %>%
		ggplot(aes(x = PC1, y = PC2, 
			   color = sampleName, 
			   shape = plantLine)) +
		geom_point(size = 2) +
		labs(x = pca_data$labels$x,
		     y = pca_data$labels$y)
	
	write.table(pca_data$data, 
		    file = paste0(SPECIES[i], '/pca_data_table', Sys.Date(), '.tsv'),
		    row.names = F, col.names = T, sep = '\t', quote = F)
	save(pca_data,
	     file = paste0(SPECIES[i], '/pca_data.RData'))
	
	ggsave(paste0("plots/", SPECIES[i], "_PCA_", Sys.Date(), ".tiff"))
	ggsave(paste0("plots/", SPECIES[i], "_PCA_", Sys.Date(), ".pdf"))
}
