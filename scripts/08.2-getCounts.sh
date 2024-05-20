# use R and bioconductor to get read counts for each species


for SPECIES in Arabidopsis Poplar;
	outputDIR=./$SPECIES/counts
	mkdir $outputDIR
	TxDbPath=./$SPECIES/${SPECIES}_TxDb.RData
	bam_baiPath=./$SPECIES/alignments
	
	Rscript $DIR/scripts/getcounts.R\
	 ${bam_baiPath}\
	 $TxDbPath\
	 $outputDIR
done
