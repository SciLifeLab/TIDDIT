DESCRIPTION
==============
TIDDIT: Is a tool to used to identify  chromosomal rearrangements using Mate Pair or Paired End sequencing data. TIDDIT identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions, using supplementary alignments as well as discordant pairs.

TIDDIT has two modes of analysing bam files. The sv mode, which is used to search for structural variants. And the cov mode that analyse the read depth of a bam file and generates a coverage report.


INSTALLATION
==============
TIDDIT requires standard c++/c libraries, python 2.7, Numpy and scipy. To compile TIDDIT, cmake must be installed. 


```
git clone https://github.com/SciLifeLab/TIDDIT.git
```

To install TIDDIT:
```
cd TIDDIT
./INSTALL.sh
```

TIDDIT is run via the TIDDIT.py script:
```

python TIDDIT.py --help
python TIDDIT.py  --sv --help
python TIDDIT.py  --cov --help
```

TIDDIT is also distributed with an experimental Singularity environment (http://singularity.lbl.gov/index.html). This environment could be used to solve issues with the c++ libraries. This environment has been tested with Singularity 2.4.1-dist. You may either use the singularity image (TIDDIT.simg), or install it from the Singularity file.

Type the following to enter a session using the TIDDIT.simg container:

    singularity shell TIDDIT.simg

Now you can install and run TIDDIT. type exit to leave the environment:

    exit

 The Singlularity environment is build typing the following command:

    singularity build TIDDIT_env.simg Singularity

you may need sudo permissions

    sudo singularity build TIDDIT_env.simg Singularity



The SV module
=============
The main TIDDIT module, detects structural variant using discordant pairs, split reads and coverage information

    python TIDDIT.py --sv [Options] --bam bam --ref reference.fasta

Where bam is the input bam file. And reference.fasta is the reference fasta used to align the sequencing data: TIDDIT will crash if the reference fasta is different from the one used to align the reads. The reads of the input bam file must be sorted on genome position.
TIDDIT may be fine tuned by altering these optional parameters:

    -o - The prefix of the output files(default = output)
        
    -i - The maximum allowed insert size of a normal pair. Pairs having larger insert 
         than this is treated as discordant pairs. Default is 3*std+mean insert size
                        
    -d - The pair orientation, use this setting to override the automatic orientation selection

    -l - The density parameter, to create a cluster, more than l signals (split reads+ discordant pairs) must be present, signals are added to a cluster if they are neighbouring atleast this  number of signals (defualt 4)
            
    -p - The minimum number of discordant pairs and supplementary alignments used to call large SV. Default is 5
    
    -r - The minimum number of supplementary alignments used to call small SV. Default is 5
            
    -q - The minimum mapping quality of the discordant pairs/supplementary alignments 
         forming a variant. Default value is 10.

    -n - The ploidy of the organism ,(default = 2)

    --force_ploidy - set the ploidy of all chromosomes to -n (including the sex chromosomes), this option will disable the ploidy estimation.
                     This option is meant to be used for low quality data or for species having equal ploidy across all chromosomes

output:

TIDDIT SV module produces three output files, a vcf file containing SV calls, a tab file describing the coverage across the genome in bins of size 100 bp, and a tab file dscribing the estimated ploidy and coverage across each contig.
                                        
The cov module
==============
Computes the coverge of different regions of the bam file

    python TIDDIT.py --cov [Options] --bam bam
    
optional parameters:

    -o - the prefix of the output files
    -z - compute the coverage within bins of a specified size across the entire genome, default bin size is 500

Filters
=============
TIDDIT uses four different filters to detect low quality calls. The filter field of variants passing these tests are set to "PASS". If a variant fail any of these tests, the filter field is set to the filter failing that variant. These are the four filters empoyed by TIDDIT:

    Expectedlinks
        The number of discordant pairs/supplementary alignments supporting
        the variant is less than 60% of the expected number of supporting reads
    FewLinks
        The number of discordant pairs supporting the variant is too low compared to the number of discordant pairs within that genomic region.
    Unexpectedcoverage
        High coverage
    Smear
        The two windows that define the regions next to the breakpoints overlap.

Failed Variants may be removed using tools such as VCFtools or grep. Removing these variants greatly improves the precision of TIDDIT, but may reduce the sensitivity. It is adviced to remove failed variants or prioritize the variants that have passed the quality checks.

Contents of the VCF INFO field
=============
TIDDIT returns the detected variants into two vcf files, one vcf for intrachromosomal variants, and one for interchromosomal variants. The INFO field of the VCF contains the following entries:

    SVTYPE
        Type of structural variant(DEL,DUP,BND,INV,TDUP)
    END
        End position of an intra-chromosomal variant
    LFA
        The number of discordant pairs at the the first breakpoint of the variant
    LFB
	The number of discordant pairs at the the second breakpoint of the variant
    LTE
        The number of discordnat pairs that form the structural variant.
    COVA
        Coverage on window A
    COVM
        The coverage between A and B
    COVB
        Coverage on window B
    OA
        Orientation of the reads in window A
    OB
        Orientation of the mates in window B
    CIPOS
        start and stop positon of window A
    CIEND
        start and stop position of window B
    EL
        Expected links to window B
    ER
        Expected number of split reads
    QUALA
        The average mapping quality of the reads in window A
    QUALB
        The average mapping quality of the reads in window B

The content of the INFO field can be used to filter out false positives and to gain more understanding of the structure of the variant. More info is found in the vcf file

LICENSE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



