# use R and bioconductor to get read counts for each species
mkdir counts

for SPECIES in Arabidopsis Poplar;
	outputDIR=./counts
	
	TxDbPath=./${SPECIES}_TxDb.RData
	bam_baiPath=./alignments/${SPECIES}/alignments
	
	Rscript $DIR/scripts/getcounts.R\
	 ${bam_baiPath}\
	 $TxDbPath\
	 $outputDIR
done
