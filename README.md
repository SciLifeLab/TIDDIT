DESCRIPTION
==============
TIDDIT: Is a tool to used to identify  chromosomal rearrangements using Mate Pair or Pair End data. TIDDIT identifies intra and inter-chromosomal translocations, deletions, tandem-duplications, intersperesed duplications and inversions. 
TIDDIT is distributed together with a database software called SVDB. SVDB is used to create structural variant databases, merge strctural variants and to use the structural variant databases as a frequency filter.
TIDDIT has two modes of analysing bam files. The sv mode, which is used to search for structural variants. And the cov mode that analyse the read depth of a bam file and generates a coverage report.

INSTALLATION
==============
Do a recursive clone of TIDDIT in order to get the database software:
```
git clone --recursive https://github.com/SciLifeLab/TIDDIT.git
```

Do the following to install TIDDIT:
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
    TIDDIT --sv [Options] -b inputfile 

    options:
    ploidy - the ploidy of the organism, 2 is default
    -o - the prefix of the output files
        
    --insert - the maximum allowed insert size of a normal pair. Pairs having larger insert 
                    than this is treated as discordant pairs. Default is 2*std+mean insert size for PE 
                    data or 4std+ mean on mp data
                        
    --orientation - the pair orientation, use this setting to override the automatic orientation selection
            
    -p - the minimum number of discordant pairs used to all a variant. Default is 4
            
        -q - the minimum mapping quality of the discordant pairs 
            forming a variant. Default value is 10.
                                        
        -c - the library coverage. Default is calculated from average genomic coverage.

The cov module
==============
Computes the coverge of different regions of the bam file
    TIDDIT --cov [Options] -b inputfile
    
    options:
    bin_size - compute the coverage within bins of a specified size across the entire genome, default bin size is 500


LICENCE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



