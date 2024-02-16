# Separate metadata file into arabidopsis and poplar
awk -F '\t' -v OFS='\t' 'NR == 1 {print} NR > 1\
 && $6 ~ /P.+/ {print}'\
 metadata/PRJNA1072118.tsv >\
 metadata/PRJNA1072118_poplar.tsv

awk -F '\t' -v OFS='\t' 'NR == 1 {print} NR > 1\
 && $6 ~ /^[AC].+/ {print}'\
 metadata/PRJNA1072118.tsv >\
 metadata/PRJNA1072118_arabidopsis.tsv