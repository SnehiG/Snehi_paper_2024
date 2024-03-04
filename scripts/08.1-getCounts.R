args <- commandArgs(trailingOnly = T)

if (length(args) == 1) {
	cat("Usage:
                Rscript script_counts.R bam_baiPath TxDbPath outPutDir
                Where:
                script_counts.R must be writen with the full path
                bamFile.ba[mi] is the full path where the bam and bai files are
                TxDbPath is the full path to the RData containing the TxDb created from biomart
                outPutDir where to save the read count file
        ")
} else {
	bamFiles <- paste0(args[1], list.files(path = args[1], pattern = ".bam$"))
	baiFiles <- paste0(bamFiles, ".bai")
	baseNames <- gsub("Aligned.sortedByCoord.out.bam", "", basename(bamFiles))
	TxDbPath <- args[2]
	outPutDir <- args[3]
	
	# IF there are packages to download, R can use any of the sites in the vector "repos" to download from
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
	
	# This RData file was created using another RScript
	load(TxDbPath)
	
	# This line tells R where the alignment files for the sample are and how much of them to read at a time
	bfl <- BamFileList(bamFiles, baiFiles, yieldSize=200000)
	
	# The code below uses summarizedOverlaps to count the number of alignments (bfl object) falling in each transcript (tx object loaded in line 31)
	# In this case we use the union mode (https://htseq.readthedocs.io/en/release_0.11.1/count.html), paired-end reads (singleEnd = F, fragments = T), and we ignore the strand information
	overlaps <- summarizeOverlaps(
		tx, bfl, mode = "Union", singleEnd = F, fragments = T)
	save(overlaps, file = paste0(outPutDir, "overlaps.RData"))
	countAssays = assays(overlaps)$counts  %>% as.data.frame()
	
	# Then we extract the counts slot, transform into a data.frame and export it as a txt file.
	names(countAssays) <- baseNames
	write.csv(countAssays, file = paste0(outPutDir, "counts.csv"), row.names = T, quote = F)
}
