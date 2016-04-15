INSTALLATION
==============

From the FindTranslocation directory run
```
mkdir build
cd build
cmake ..
make
```

DESCRIPTION
==============
FindTranslocations: tool to identify  chromosomal rearrangements using Mate Pair or Pair End data. The idea is to identify areas having clusters of discordant pairs or abnormal coverage.

The SV module
=============
The main FT modules, detects structural variant using discordant pairs, split reads and coverage information
    FindTransloctions --sv [Options] --bam inputfile 

    options:
    ploidy - the ploidy of the organism, 2 is default
    output - the prefix of the output files
        
    max-insert - the maximum allowed insert size of a normal pair. Pairs having larger insert 
                    than this is treated as discordant pairs. Default is 1.5*std+mean insert size for PE 
                    data or 4std+ mean on mp data
                        
    orientation - the pair orientation, use this setting to override the automatic orientation selection
            
    pairs - the minimum number of discordant pairs used to all a variant. Default is 3
            
        q - the minimum mapping quality of the discordant pairs 
            forming a variant. Default value is 0.
                                        
        coverage - the library coverage. Default is calculated from average genomic coverage.

The cov module
==============
Computes the coverge of different regions of the bam file
    FindTranslocations --cov [Mode] --bam inputfile
    
    options:
    bin - compute the coverage within bins of a specified size across the entire genome
            , outputs a tab file of the format chromosome    start  stop coverage
            
            
            light - compute the coverage within bins of a specified size across the 
                    entire genome, outputs a tab file of the format chromosome  coverage
            

LICENCE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



