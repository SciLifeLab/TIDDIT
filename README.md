DESCRIPTION
==============
TIDDIT: Is a tool to used to identify  chromosomal rearrangements using Mate Pair or Paired End sequencing data. TIDDIT identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions, using supplementary alignments as well as discordant pairs.
TIDDIT searches for discordant reads and split reads (supplementary alignments). TIDDIT also performs local assembly using a custom local de novo assembler.
Next all signals (contigs, split-reads, and discordant pairs) are clustered using DBSCAN. The resulting clusters are filtered and annotated, and reported as SV depending on the statistics.
TIDDIT has two analysis modules. The sv mode, which is used to search for structural variants. And the cov mode that analyse the read depth of a bam file and generates a coverage report.
On a 30X human genome, the TIDDIT SV module typically completetes within 5 hours, and requires less than 10Gb ram.

INSTALLATION
==============
TIDDIT requires python3 (=> 3.8), cython, pysam, and Numpy.


By default, tiddit will require, bwa, fermi2 and ropebwt2 for local assembly; local assembly may be disabled through the "--skip_assembly" parameter.

Installation

Cloning from Git Hub:

```
git clone https://github.com/SciLifeLab/TIDDIT.git
```

To install TIDDIT:
```
cd tiddit
pip install -e .
```

Next install bwa, I recommend using conda:

conda install bwa

You may also compile bwa yourself. Remember to add executables to path, or provide path through the command line parameters.
```
tiddit --help
tiddit --sv --help
tiddit --cov --help
```

Optionally, the assembly calling may be turned off using the "--skip_assembly" option.

TIDDIT may be installed using bioconda:
```
conda install tiddit
```

or using the docker image on biocontainers
```
docker pull quay.io/biocontainers/tiddit:<tag>
```
visit https://quay.io/repository/biocontainers/tiddit?tab=tags for tags.


The SV module
=============
The main TIDDIT module, detects structural variant using discordant pairs, split reads and coverage information

    python tiddit --sv [Options] --bam in.bam --ref reference.fa

Where bam is the input bam or cram file. And reference.fasta is the reference fasta used to align the sequencing data: TIDDIT will crash if the reference fasta is different from the one used to align the reads. The reads of the input bam file must be sorted on genome position.

TIDDIT may be fine-tuned by altering these optional parameters:


	-o	output prefix(default=output)
	-i	paired reads maximum allowed insert size. Pairs aligning on the same chr at a distance higher than this are considered candidates for SV (default= 99.9th percentile of insert size)
	-d	expected reads orientations, possible values "innie" (-> <-) or "outtie" (<- ->). Default: major orientation within the dataset
	-p	Minimum number of supporting pairs in order to call a variant (default 3)
	-r	Minimum number of supporting split reads to call a variant (default 3)
	--threads	Number of threads (default 1)
	-q	Minimum mapping quality to consider an alignment (default 5)
	-n	the ploidy of the organism,(default = 2)
	-e	clustering distance parameter, discordant pairs closer than this distance are considered to belong to the same variant(default = sqrt(insert-size*2)*12)
	-c	average coverage, overwrites the estimated average coverage (useful for exome or panel data)
	-l	min-pts parameter (default=3),must be set >= 2
	-s	Number of reads to sample when computing library statistics(default=25000000)
	-z 	minimum variant size (default=50), variants smaller than this will not be printed ( z < 10 is not recomended)
	--force_ploidy	force the ploidy to be set to -n across the entire genome (i.e skip coverage normalisation of chromosomes)
	 --force_overwrite     force the analysis and overwrite any data in the output folder
	--n_mask	exclude regions from coverage calculation if they contain more than this fraction of N (default = 0.5)
	--skip_assembly	Skip running local assembly, tiddit will perform worse, but wont require fermi2, bwa, ropebwt and bwa indexed ref
	--bwa	path to bwa executable file(default=bwa)
	--fermi2	path to fermi2 executable file (default=fermi2)
	--ropebwt2	path to ropebwt2 executable file (default=ropebwt2)
	--p_ratio	minimum discordant pair/normal pair ratio at the breakpoint junction(default=0.1)
	--r_ratio	minimum split read/coverage ratio at the breakpoint junction(default=0.1)
	--max_coverage	filter call if X times higher than chromosome average coverage (default=4)
	--min_contig	 Skip calling on small contigs (default < 10000 bp)



output:

TIDDIT SV module produces two output files, a vcf file containing SV calls, and a tab file dscribing the estimated ploidy and coverage across each contig.
         
                                
The cov module
==============
Computes the coverge of different regions of the bam file

    python TIDDIT.py --cov [Options] --bam bam
    
optional parameters:

    -o - the prefix of the output files
    -z - compute the coverage within bins of a specified size across the entire genome, default bin size is 500
    -w - generate a wig file instead of bed
 --ref - reference sequence (fasta), required for reading cram file.

Filters
=============
TIDDIT uses four different filters to detect low quality calls. The filter field of variants passing these tests are set to "PASS". If a variant fail any of these tests, the filter field is set to the filter failing that variant. These are the four filters empoyed by TIDDIT:

    Expectedlinks
	Less than <p_ratio> fraction of the spanning pairs or <r_ratio> fraction reads support the variant
    FewLinks
        The number of discordant pairs supporting the variant is too low compared to the number of discordant pairs within that genomic region.
    Unexpectedcoverage
        High coverage

Failed Variants may be removed using tools such as VCFtools or grep. Removing these variants greatly improves the precision of TIDDIT, but may reduce the sensitivity. It is adviced to remove filtered variants or prioritize the variants that have passed the quality checks.
This command may be usedto filter the TIDDIT vcf:

	grep -E "#|PASS" input.vcf > output.filtered.vcf

Quality column
=============
The scores in the quality column are calculated using non parametric sampling: 1000 points/genomic positions are sampled across each chromosome. And the number of read-pairs and reads spanning these points are counted.
The variant support of each call is compared to these values, and the quality column is set to he lowest percentile higher than the (variant support*ploidy).

Note: SVs usually occur in repetetive regions, hence these scores are expected to be relatively low. A true variant may have a low score, and the score itself depends on the input data (mate-pair vs pe for instance).

Merging the vcf files
=====================
I usually merge vcf files using SVDB (https://github.com/J35P312)

svdb --merge --vcf file1.vcf file2.vcf --bnd_distance 500 --overlap 0.6 > merged.vcf

Merging of vcf files could be useful for tumor-normal analysis or for analysing a pedigree. But also to combine the output of multiple callers.

Tumor normal example
===================

run the tumor sample using a lower ratio treshold (to allow for subclonal events, and to account for low purity)

python TIDDIT.py --sv --p_ratio 0.10 --bam tumor.bam -o tumor --ref reference.fasta
grep -E "#|PASS" tumor.vcf > tumor.pass.vcf

run the normal sample

python TIDDIT.py --sv --bam normal.bam -o normal --ref reference.fasta
grep -E "#|PASS" normal.vcf > normal.pass.vcf

merge files:

svdb --merge --vcf tumor.pass.vcf normal.pass.vcf --bnd_distance 500 --overlap 0.6 > Tumor_normal.vcf

The output vcf should be filtered further and annotated (using a local-frequency database for instance)

Annotation
==========
genes may be annotated using vep or snpeff. NIRVANA may be used for annotating CNVs, and SVDB may be used as a frequency database

Algorithm
=========

Discordant pairs, split reads (supplementary alignments), and contigs are extracted. A discordant pair is any pair having a larger insert size than the  -i paramater, or a pair where the reads map to different chromosomes.
supplementary alignments and discordant pairs are only extracted if their mapping quality exceed the -q parameter. Contigs are generated by assembling all reads with supplementary alignment using fermi2

The most recent version of TIDDIT uses an algorithm similar to DBSCAN: A cluster is formed if -l or more signals are located within the -e distance. Once a cluster is formed, more signals may be added if these signals are within the
-e distance of -l signals within a cluster.

A cluster is rejected if it contains less than -r plus -p signals. If the cluster is rejected, it will not be printed to the vcf file.

If the cluster is not rejected, it will be printed to file, even if it fails any quality filter. 

The sensitivity and precision may be controlled using the -q,r,p, and -l parameters. 

LICENSE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



