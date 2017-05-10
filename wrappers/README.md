#TIDDIT_multi_sample.nf

a script used to run multiple samples through TIDDIT. The script uses nextflow to run each sample in parallel.

#Installation
download and install nextflow

https://www.nextflow.io/

set the TIDDIT(line 2) executable path, and TIDDIT_multi_sample.py(line 4) script path in the TIDDIT_multi_sample.nf script.
Install the latest version of svdb

#usage

./nextflow TIDDIT_multi_sample.nf --bam sample1.bam,sample2.bam,sample3.bam --output merged.vcf

this command will perform variant calling on the bam files sample1.bam, sample2.bam, and sample3.bam, and generate an output vcf called merged.vcf

#settings/algorithm
The variants are called at a set split read and discordant pair threshold(according to the variables set at line 10 and 11 of 
TIDDIT_multi_sample.nf).
Once variant calling is complete, the vcf files of each sample are merged using svdb. Thereafter, the variant calls are filtered.
For each variant, if there is any sample having a signal strenght exceeding the split read and discordant pair, limit set on lines 14, and 15, the variant is kept. Otherwise, the variant is discarded.
The purpose of running TIDDIT in multi sample mode  is to be able to detect variants whose signal is low due to noise and random errors.

#nextflow settings
See the nextflow website for info on how to configure nextflow.
https://www.nextflow.io/

