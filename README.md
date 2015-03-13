INSTALLATION
==============

From the FindTranslocation directory run:
- mkdir build
- cd build
- cmake ..
- make


for Uppmax change step three into 
- cmake .. -DBoost_NO_BOOST_CMAKE=ON

You will find the binaries in the main directory under bin. In case of problems the majority of the times there is a problem
with the local installation of boost.


DESCRIPTION
==============
FindTranslocations: tool to identify  chromosomal rearrangements using Mate Pair or Pair End data. The idea is to identify areas 

Options:
* ``--bam`` alignment file in bam format, sorted by coordinate. If bwa mem is used option -M MUST be specified in order to map as secondary the splitted reads
* ``--min-insert`` paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goes from beginning of first read to end of second read
* ``--max-insert`` paired reads maximum allowed insert size. pairs aligning on the same chr at a distance higher than this are considered candidates for SV.
* ``--orientation`` expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default outtie
* ``--output`` Header of the output file names
* ``--minimum-supporting-pairs`` Minimum number of supporting pairs in order to call a variation event (default 10) 
* ``--minimum-mapping-quality`` Minimum mapping quality to consider an alignment (default 20) 
* ``--window-size`` Size of the sliding window (default 1000) 
* ``--window-step`` Size of the step in overlapping window (must be lower than window-size) (default 100)


LICENCE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



