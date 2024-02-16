# Aligning clean reads to genomes using STAR version 2.7.9a
NCPUS=32

# make sure to be in the base directory of the project
for SPECIES in Arabidopsis Poplar; do
	# Folders where files are and will be stored
	inDIR=cleanReads/$SPECIES
	mkdir -p alignments/$SPECIES/
	# List of clean reads
	r1s=($(ls $inDIR/*R1.fastq.gz))
	r2s=($(ls $inDIR/*R2.fastq.gz))
	
	# Set the folder containing the info about the genomes
	if [ $SPECIES == "Arabidopsis" ]; then
	        indexDIR=metadata/TAIR10_idx
	else
	        indexDIR=metadata/Poptriv4_idx
	fi
	
	# This is a safeguard in case we ran already the script,
	# but there was not enough time to finish all alignments
	# We get the name of the file that is produced only when the alignment is done (Log.final.out)
	# and make a list of all files with that name in the folder
	finishedAlignments=($(ls $DIR/alignments/$SPECIES/*Log.final.out))
	
	# Then we remove the folder name and the Log.final.out from the name,
	# which gives just the name of the sample
	for i in ${!finishedAlignments[@]}; do 
	 basename ${finishedAlignments[i]/Log.final.out}
	done > samplesFinished
	
	# Repeat the commands for each file in the r1s list
	for i in ${!r1s[@]}; do
	        r1=${r1s[i]}
	        # remove the folder name from the name of the file
	        r1_basename=$(basename ${r1})
	        # If the sample is in the list of the ones that are finished, we skip it
	        if  grep -q ${r1_basename/_1.fastq} samplesFinished ; 
	         then echo $r1 is done; 
	         continue; 
	        fi
	
	        r2=${r2s[i]}
	        r1_basename=$(basename ${r2/.gz})
	        outPrefix=${r1_basename/_1.fastq}
	
	        # Run the alignment
	        STAR --genomeDir $indexDIR/\
	         --runThreadN $NCPUS \
	         --readFilesIn $r1 $r2 \
	         --outFileNamePrefix alignments/$SPECIES/${outPrefix}\
	         --outSAMtype BAM SortedByCoordinate \
	         --outSAMunmapped Within \
	         --outSAMattributes Standard
	done
done