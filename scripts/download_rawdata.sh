# Download file from SRA using sra-toolkit
## Arabidopsis ##
for species in Arabidopsis Poplar; do
 mkdir -p RAW_DATA/$species links/$species
 cut -f 1,7 metadata/PRJNA1072118_$species.tsv |
  while read accession library; do
   fasterq-dump --seq-defline '@$sn[_$rn]/$ri'\
    --split-files $accession\
    --outdir RAW_DATA/$species;
   # Create links to have files with correct names
   ln -s RAW_DATA/$species/${accession}_1.fastq links/$species/${library}_1.fastq
   ln -s RAW_DATA/$species/${accession}_2.fastq links/$species/${library}_2.fastq
 done
 
 
done

