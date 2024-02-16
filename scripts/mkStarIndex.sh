# Make STAR indexes for poplar and arabidosis genomes
# Using STAR version 2.7.9a
NCPUS=8

# Start by downloading genomes
cd metadata
# Fasta
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/fasta/populus_trichocarpa/dna/Populus_trichocarpa.Pop_tri_v4.dna.toplevel.fa.gz
# GTF
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz
wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/plants/gtf/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v4.58.gtf.gz

# the three variable below are vectors, so ${genomeIdxDIR[0]} prints TAIR10_idx
genomeIdxDIR=(TAIR10_idx Poptriv4_idx)
genomeFastaFiles=(Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Populus_trichocarpa.Pop_tri_v4.dna.toplevel.fa)
genomeGTFFile=(Arabidopsis_thaliana.TAIR10.57.gtf Populus_trichocarpa.Pop_tri_v4.57.gtf)
readLength=99

for i in 0 1; do
 mkdir ${genomeIdxDIR[i]}

 STAR\
  --runThreadN ${NCPUS}\
  --runMode genomeGenerate\
  --genomeDir ${genomeIdxDIR[i]}\
  --genomeFastaFiles ${genomeFastaFiles[i]}\
  --sjdbGTFfile ${genomeGTFFile[i]}\
  --sjdbOverhang ${readLength}

 done
