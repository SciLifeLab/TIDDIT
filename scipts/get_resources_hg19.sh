## get repeat masker
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
## get know genes
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
## get gaps
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz

## gunzip
gunzip rmsk.txt.gz knownGene.txt.gz gap.txt.gz

cat gap.txt | grep -v "#" | awk '{print $2 "\t" $3  "\t" $4  "\t" "gap_" $8  }' > gaps.bed
cat knownGene.txt | awk '{print $2 "\t" $4  "\t" $5  "\t" "gene_" $1}' > genes.bed
cat rmsk.txt | awk '{print $6 "\t" $7  "\t" $8  "\t" "rmsk_" $12}' > rmsk.bed

rm gap.txt knownGene.txt rmsk.txt
