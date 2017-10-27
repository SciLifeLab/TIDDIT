DESCRIPTION
==============
TIDDIT: Is a tool to used to identify  chromosomal rearrangements using Mate Pair or Paired End sequencing data. TIDDIT identifies intra and inter-chromosomal translocations, deletions, tandem-duplications, intersperesed duplications and inversions, using supplementary alignments as well as discordant pairs. 

TIDDIT is distributed together with a database software called SVDB. SVDB is used to create structural variant databases, merge structural variants and to use the structural variant databases as a frequency filter.

TIDDIT has two modes of analysing bam files. The sv mode, which is used to search for structural variants. And the cov mode that analyse the read depth of a bam file and generates a coverage report.

TIDDIT is mainly designed to run on whole genome sequencing data. However, TIDDIT is also able to perform variant calling on exome data if the coverage is supplied through the -c parameter.

A nextflow wrapper for running multiple samples at once is available in the wrappers folder

[![DOI](https://zenodo.org/badge/81584907.svg)](https://zenodo.org/badge/latestdoi/81584907)

INSTALLATION
==============
TIDDIT requires only standard c++/c libraries. To compile TIDDIT, cmake must be installed.

Do a recursive clone of TIDDIT in order to get the database software:
```
git clone --recursive https://github.com/SciLifeLab/TIDDIT.git
```

To install TIDDIT:
```
cd TIDDIT
mkdir build
cd build
cmake ..
make
```
The executable is located in the bin folder:
```
cd ..
cd bin
```
run the executable file to view the help message or run TIDDIT:
```
./TIDDIT
./TIDDIT --help
./TIDDIT  --sv --help
./TIDDIT  --cov --help
```

The SV module
=============
The main TIDDIT module, detects structural variant using discordant pairs, split reads and coverage information

    TIDDIT --sv [Options] -b bam 

Where bam is the input bam file. The reads of the input bam file must be sorted on genome position.
TIDDIT may be fine tuned by altering these optional parameters:

    -n - the ploidy of the organism, 2 is default
    -o - the prefix of the output files(default = output)
        
    -i - the maximum allowed insert size of a normal pair. Pairs having larger insert 
                    than this is treated as discordant pairs. Default is 3*std+mean insert size
                        
    -d - the pair orientation, use this setting to override the automatic orientation selection
            
    -p - the minimum number of discordant pairs and supplementary alignments used to call large SV. Default is 4
    
    -r - the minimum number of supplementary alignments used to call small SV. Default is 6
            
    -q - the minimum mapping quality of the discordant pairs/supplementary alignments 
            forming a variant. Default value is 10.
                                        
    -c - the library coverage. Default is calculated from average genomic coverage.
        

The cov module
==============
Computes the coverge of different regions of the bam file

    TIDDIT --cov [Options] -b bam
    
optional parameters:

    -o - the prefix of the output files
    -z - compute the coverage within bins of a specified size across the entire genome, default bin size is 500

Filters
=============
TIDDIT uses four different filters to detect low quality calls. The filter field of variants passing these tests are set to "PASS". If a variant fail any of these tests, the filter field is set to the filter failing that variant. These are the four filters empoyed by TIDDIT:

    Expectedlinks
        The number of discordant pairs/supplementary alignments supporting
        the variant is less than 40% of the expected number of supporting reads
    FewLinks
        The number of discordant pairs supporting the variant is less than 20% of the 
        discordant pairs within that genomic region.
    Unexpectedcoverage
        The coverage across the variant is more than 10* the mean coverage.
    Smear
        The two windows that define the regions next to the breakpoints overlap.

Failed Variants may be removed using tools such as VCFtools or grep. Removing these variants greatly improves the precision of TIDDIT, but may reduce the sensitivity. It is adviced to remove failed variants or prioritize the variants that have passed the quality checks.

Contents of the VCF INFO field
=============
TIDDIT returns the detected variants into two vcf files, one vcf for intrachromosomal variants, and one for interchromosomal variants. The INFO field of the VCF contains the following entries:

    SVTYPE
        Type of structural variant(DEL,DUP,BND,INV,TDUP,IDUP)
    END
        End position of an intra-chromosomal variant
    LFW
        The number of discordant pairs close to the breakpoints of the variant
    LCB
        The number of discordant pairs close to the breakpoints of the variant, that map to the same chromosome pair as the pairs defining the variant 
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
    CHRA
        The chromosome of window A
    CHRB
        The chromosome of window B
    WINA
        start and stop positon of window A
    WINB
        start and stop position of window B
    EL
        Expected links to window B
    ER
        Expected number of split reads
    RATIO
        The number of links divided by the expected number of links
    QUALA
        The average mapping quality of the reads in window A
    QUALB
        The average mapping quality of the reads in window B

The content of the INFO field can be used to filter out false positives and to gain more understanding of the structure of the variant.

Algorithm
=============
TIDDIT detects structural variants using supplementary alignments as well as discordant pairs. A discordant pair is defined as any pair of reads having a larger distance than the --insert parameter(which is set to 3*std+average library distance as default). Supplementary aligments are produced by reads where one part of the read aligns to a certain part of the reference, while another part of the same read aligns to a distant region.

TIDDIT performs a linear search for SV signatures within the input Bam file. These signatures are sets, or clusters of Discordant pairs/supplementary alignments, having similar patterns regarding coverage and position within the genome.
Any set of discordant pairs or supplementary alignments larger or equal to the -p parameter will be analysed and later returned as a structural variant by printing it to the vcf file. TIDDIT will only consider reads that fullfill the -q parameter: reads having lower mapping quality will not be added to a set, and thus will not contribute to the detection of SV.

TIDDIT detects a wide spectra of strutural variants, and is able to classify deletions, duplications(tandem and interspersed), inversions and translocations(intrachromosomal and interchromosomal). Variants are classified based on the pair orientation of the reads defining a structural variant, as well as the coverage across the structural variant and the regions where the read pairs are aligned. If TIDDIT is unnable to classify a variant, it will be returned as a break end event.

LICENSE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



