INSTALLATION
==============

From the TIDDIT directory run
```
mkdir build
cd build
cmake ..
make
```

DESCRIPTION
==============
TIDDIT: tool to identify  chromosomal rearrangements using Mate Pair or Pair End data. The idea is to identify areas having clusters of discordant pairs or abnormal coverage.

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
    TIDDIT --cov [Mode] -b inputfile
    
    options:
    bin_size - compute the coverage within bins of a specified size across the entire genome, default bin size is 500


LICENCE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



