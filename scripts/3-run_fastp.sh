# Softwares needed scipy-stack fastp
N_threads=8 # change as needed

for species in Arabidopsis Poplar; do
	mkdir -p cleanReads/$species fastpReports
	
	# Paths to files
	## Read1 files
	### We create 2 vectors with the paths of the files for R1 and the output R1, by replacing the folder name links to cleanReads
	R1=($(ls $PWD/links/$species/*_1.fastq))
	CleanR1=($(echo ${R1[@]} | sed 's:links:cleanReads:g'))
	## Read2 files
	### R2 files have the same name as the R1s, just replace R1 by R2
	R2=($(echo ${R1[@]} | sed 's:_R1:_R2:g'))
	CleanR2=($(echo ${CleanR1[@]} | sed 's:_R1:_R2:g'))
	
	# We get the indices of the vector R1 and repeat the code for each one
	for rf in ${!R1[@]}; do
	        reportName=$(basename ${R1[$rf]} | sed 's:_R1.fastq:.html:')
	
	         fastp -i ${R1[$rf]}\
	         -I ${R2[$rf]}\
	         -o ${CleanR1[$rf]}\
	         -O ${CleanR2[$rf]}\
	         --qualified_quality_phred 20\
	         --unqualified_percent_limit 30\
	         --adapter_sequence=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG\
	         --adapter_sequence_r2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG\
	         --cut_front --cut_front_window_size 3\
	         --cut_right --cut_right_window_size 4\
	         --cut_right_mean_quality 15\
	         --length_required 25\
	         --html $PWD/fastpReports/$reportName\
	         --thread ${N_threads}
	done
done
	
# In short, fastp will remove reads with more than 30% of bases with quality score < 20
# Will check the mean quality of bases in a window of 3bp at the beginning of the reads and trim the read until the quality of the window is > 20
# Will do the same as above, but from the end of the read, but the window is 4bp and the quality has to be > 15
# Will drop any read of length < 25bp
# Will write the report to the path given in front of --html