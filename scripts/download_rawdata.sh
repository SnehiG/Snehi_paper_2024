# Download file from SRA using sra-toolkit
## Arabidopsis ##
for species in Arabidopsis Poplar; do
 mkdir -p RAW_DATA/$species
 cut -f 1 metadata/PRJNA1072118_$species.tsv |
  while read accession; do
   fasterq-dump --seq-defline '@$sn[_$rn]/$ri'\
    --split-files $accession\
    --outdir RAW_DATA/$species;
 done
done

