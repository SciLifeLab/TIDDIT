assemblator is an assembly pipeline used to filter the TIDDIT output files. 

Run

nextflow assemblatron.nf --vcf tiddit_vcf_file.vcf --bam file.bam --genome reference.fa

The variants of tiddit_vcf_file.vcf are extracted and assembled using abyss. The contigs are then aligned to the reference.
Variants where less than 50% of the contigs map within the regions denoted window A and window B will be flagged as low quality variants.

Dependencies
		nextflow
		abyss
		bwa
		samtools

