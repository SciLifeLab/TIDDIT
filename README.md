INSTALLATION
==============

From the FRCurve directory run:
- mkdir build
- cd build
- cmake ..
- make

You will find the binaries in the main directory under bin. In case of problems the majority of the times there is a problem
with the local installation of boost.


DESCRIPTION
==============
 FRCcurve is a package containing tools to process bam files in order to evaluate and analyze de novo assembly/assemblers and identify Structural Variations 
 suspicious genomics regions. The tools have been already successfully applied in several de novo and resequencing projects.
 
 This package contains two tools:
 
1. FRCbam: tool to compute Feature Response Curves in order to validate and rank assemblies and assemblers
2. FindTranslocations: tool to identify  chromosomal rearrangements using Mate Pairs


FRCbam
--------------
 **USAGE: basic, no CE-stats tuning**

1. Assemble your data (n PE libraries and m MP libraries) with your favorite tools. Let us call the assemblies A_tool1, A_tool2, etc.
2. Align one PE library and one MP library against each of your assemblies (e.g., A_tool1)       
  1. Use the same parameters
  2. PE library is mandatory, MP library is highly recommended
  3. sort and index the generated bam files by coordinate. We will call them A_tool1_PE_lib.bam and A_tool1_MP_lib.bam
  4. use PE library with largest read coverage (i.e., vertical coverage) and MP with largest spanning coverage (i.e., horizontal coverage)
3.Run FRCurve for each assembly: 
```FRC --pe-sam A_tool1_PE_lib.bam --pe-min-insert MIN_PE_INS --pe-max-insert MAX_PE_INS -mp-sam A_tool1_MP_lib.bam  --mp-min-insert MIN_MP_INS --mp-max-insert MAX_MP_INS 
		--genome-size ESTIMATED_GENOME_SIZE --output OUTPUT_HEADER```

where:

* ```--pe-sam``` A_tool1_PE_lib.bam```: sorted bam file obtained aligning PE library against assembly obtained with tool A;
* ```--pe-min-insert MIN_PE_INS```: estimated min insert length
* ```--pe-max-insert MAX_PE_INS``` : estimated max insert length
* ```--mp-sam A_tool1_MP_lib.bam```: sorted bam file obtained aligning MP library against assembly obtained with tool A;
* ```--mp-min-insert MIN_MP_INS``` : estimated min insert length
* ```--mp-max-insert MAX_MP_INS``` : estimated max insert length
* ```--genome-size ESTIMATED_GENOME_SIZE```: estimated genome size;
* ```--output OUTPUT_HEADER```: output header;
	
**IMPORTANT**:
If ```--genome-size``` is not specified the assembly length is used to compute FRCurve. In order to be able to compare FRCurves
obtained with different tools (and hence producing slightly different assembly sizes) the same ```ESTIMATED_GENOME_SIZE```
must be specified.
		
OUTPUT:

* ```OUTPUT_HEADER_Features.txt```: human readable description of features: contig start end feature_type
* ```OUTPUT_HEADER_FRC.txt```: FRCurve computed with all the features (to be plotted)
* ```OUTPUT_HEADER_FEATURE.txt```: FRCurve for the corresponding feature
* ```OUTPUT_HEADER_featureType.txt```: for each featureType the specific FRCurve
* ```Features.gff```: features description in GFF format (for visualization)
* ```OUTPUT_HEADER_CEstats_PE.txt```: CEvalues distribution (for CE_stats tuning)
* ```OUTPUT_HEADER_CEstats_MP.txt```: CEvalues distribution (for CE_stats tuning)
		
**USAGE: advanced, CE-stats tuning**

CE-stats are able to identify the presence of insertion and deletion events. Different insert sizes give the possibility to
identify different events. In order to avoid too many False Positives (or too many False Negatives) a tuning phase is 
highly recommended.

Once step 3 of USAGE is done, the user can already plot the FRCurves (for all the features or for only some of them).
CE_stats based features have been computed with default (i.e., not optimal) parameters. Each run of FRCbam produces two
files: ```OUTPUT_HEADER_CEstats_PE.txt``` and ```OUTPUT_HEADER_CEstats_MP.txt```. These files contain the distribution of the ```CE_values```
on each assembly. These values must be plotted as suggested in page 3 of Supplementary Material 
(see http://www.nada.kth.se/~vezzi/publications/supplementary.pdf) to estimate the optimal ```CE_min``` and ```CE_max``` values for 
 the PE and MP library respectively.
 Once the optimal parameters are estimated FRCurves must be recomputed for all assemblies (only ```COMPR``` and ```STRECH``` features
 will change) specifying the following extra parameters:
 
* ```--CEstats-PE-min CE_PE_MIN```: all position with CE values computed with PE library lower than this are considered compressions
* ```--CEstats-PE-max CE_PE_MAX```: all position with CE values computed with PE library higher than this are considered expansions 
* ```--CEstats-MP-min CE_MP_MIN```: all position with CE values computed with MP library lower than this are considered compressions
* ```--CEstats-MP-max CE_MP_MAX```: all position with CE values computed with MP library higher than this are considered expansions
 
FindTranslocations
--------------
The tool is under constant development. THere will be soon a detailed user guide, for now run the tool with ```--help``` to discover the options.




LICENCE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



