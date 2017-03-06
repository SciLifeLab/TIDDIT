/*
   Francesco Vezzi
   Jesper Eisfeldt
*/

#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <string>
#include <queue>

#include <sstream>
#include <iostream>
#include <fstream>

#include "data_structures/ProgramModules.h"

//converts a string to int
int convert_str(string to_be_converted,string option){
	int converted;
	
	istringstream convert(to_be_converted);
	if ( !(convert >> converted) ){
		cout << "ERROR: non-integer option, " << option << " " << to_be_converted << endl;
		exit(0);
	}

	return(converted);
}

int main(int argc, char **argv) {
	//MAIN VARIABLE

	bool outtie 				    = true;	 // library orientation
	uint32_t minimumSupportingPairs = 3;
	uint32_t minimumSupportingReads = 3;
	int min_insert				    = 100;      // min insert size
	int max_insert				    = 100000;  // max insert size
	int minimum_mapping_quality		=10;
	float coverage;
	float coverageStd;
	float meanInsert;
	float insertStd;
	int min_variant_size= 100;
	string outputFileHeader ="output";
	string version = "1.0.2";
	
	//collect all options as a vector
	vector<string> arguments(argv, argv + argc);
	arguments.erase (arguments.begin());

	map<string,string> vm;
	//general options
	vm["-b"]="";
	vm["-o"]="";
	vm["--help"] ="store";
	
	string general_help="TIDDIT-" + version  +  " - a structural variant caller\nUsage: TIDDIT <module> [options] \n";
	general_help+="modules\n\t--help\tproduce help message\n";
	general_help+="\t--sv\tselect the sv module to find structural variations\n";
	general_help+="\t--cov\tselect the cov module to analyse the coverage of the genome using bins of a specified size\n";
	
	//the sv module
	vm["--sv"]="store";
	vm["-i"]="";
	vm["-d"]="";
	vm["-p"]="";
	vm["-r"]="";
	vm["-q"]="";
	vm["-c"]="";
	vm["-n"]="";
	vm["-m"] = "";
	
	string sv_help ="\nUsage: TIDDIT --sv -b inputfile [-o prefix] \nOther options\n";
	sv_help +="\t-b\tcoordinate sorted bam file(required)\n";
	sv_help +="\t-i\tpaired reads maximum allowed insert size. Pairs aligning on the same chr at a distance higher than this are considered candidates for SV (=3std + mean_insert_size)\n";
	sv_help +="\t-d\texpected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default: major orientation within the dataset\n";
	sv_help +="\t-p\tMinimum number of supporting pairs in order to call a variation event (default 3)\n";
	sv_help +="\t-r\tMinimum number of supporting split reads to call a small variant (default 3)\n";
	sv_help +="\t-q\tMinimum mapping quality to consider an alignment (default 10)\n";
	sv_help +="\t-c\taverage coverage, (default= computed from the bam file)\n";
	sv_help +="\t-n\tthe number of sets of chromosomes,(default = 2)\n";
	sv_help +="\t-m\tminimum variant size,(default = 100)\n";
			
	//the coverage module
	vm["--cov"]="store";
	vm["-z"]="";
	
	string coverage_help="\nUsage: TIDDIT --cov [Mode] -b inputfile [-o prefix]\n";
	coverage_help +="\t-b\tcoordinate sorted bam file(required)\n";
	coverage_help +="\n\t-z\tuse bins of specified size(default = 500bp) to measure the coverage of the entire bam file, set output to stdout to print to stdout\n";
	
	
	//store the options in a map
	for(int i = 0; i < arguments.size(); i += 2){
		if(vm.count( arguments[i] ) ){
			if(vm[arguments[i]] == ""){
				if(i+1 < arguments.size()){
					vm[arguments[i]]=arguments[i+1];
				}else{
					cout << "Missing argument: " << arguments[i] << endl;
					return(1);
				}
			}else{
				vm[arguments[i]]="found";
				i += -1;
			}
		}else{
			cout << "invalid option: " << arguments[i] << endl;
			cout << "type --help for more information" << endl;
			return(1);
		}
	}
	
	//print help message
	if(vm["--help"] == "found" or arguments.size() == 0){
		if( vm["--sv"] == "found"){
			cout << sv_help;
		}else if(vm["--cov"] == "found"){
			cout << coverage_help;
		}else{
			cout << general_help;
		}
		return(0);
	}
	
	//if no module is chosen
	if(vm["--sv"] == "store" and vm["--cov"] == "store"){
		cout << general_help;
		cout << "ERROR: select a module" << endl;
		return(0);
	}else if(vm["--sv"] == "found" and vm["--cov"] == "found"){
		cout << general_help;
		cout << "ERROR: only one module may be selected" << endl;
		return(0);
	
	}
	
	//the bam file is required by all modules
	if(vm["-b"] == ""){
		cout << general_help;
		cout << "ERROR: select a bam file using the -b option" << endl;
		return(1);
	}

	string alignmentFile	= vm["-b"];
	uint64_t genomeLength = 0;
	uint32_t contigsNumber = 0;

	// Now parse BAM header and extract information about genome lenght
	BamReader bamFile;
	bamFile.Open(alignmentFile);
	SamHeader head = bamFile.GetHeader();
	if (head.HasSortOrder()) {
		string sortOrder = head.SortOrder;
		if (sortOrder.compare("coordinate") != 0) {
			cout << "sort order is " << sortOrder << ": please sort the bam file by coordinate\n";
			return 1;
		}
	} else {
		cout << "BAM file has no @HD SO:<SortOrder> attribute: it is not possible to determine the sort order\n";
		return 1;
	}

	SamSequenceDictionary sequences  = head.Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		genomeLength += StringToNumber(sequence->Length);
		contigsNumber++;
	}	
	bamFile.Close();

	//if the find structural variations module is chosen collect the options
	if(vm["--sv"] == "found"){
		if(vm["-o"] != ""){
		    outputFileHeader=vm["-o"];
		}
	
	
		if(vm["-q"] != ""){
			minimum_mapping_quality=convert_str( vm["-q"] ,"-q");
		}
		if(vm["-p"] != ""){
			minimumSupportingPairs=convert_str( vm["-p"] , "-p");
		}
		if(vm["-r"] != ""){
			minimumSupportingReads=convert_str( vm["-r"] , "-r");
		}

		if(vm["-d"] != ""){
			if (vm["-d"] == "outtie"){
				outtie=true;
			}else if (vm["-d"] == "innie"){
				outtie=false;
			}else{
				cout << "ERROR: invalid orientation " << vm["-d"] << endl;
				return(0);
			}
		}

		int ploidy = 2;
		if(vm["-n"] != ""){
			ploidy = convert_str( vm["-n"],"-n" );
		}
		
		//now compute library stats
		LibraryStatistics library;
		size_t start = time(NULL);
		library = computeLibraryStats(alignmentFile, genomeLength, max_insert, 50 , outtie,minimum_mapping_quality,outputFileHeader);
		printf ("library stats time consumption= %lds\n", time(NULL) - start);
		
        
		coverage   = library.C_A;
		if(vm["-c"] != ""){
			coverage    = convert_str( vm["-c"],"-c" );
		}
		
		if(vm["-d"] == ""){
			outtie=library.mp;
			if(outtie == true){
				cout << "auto-config orientation: outtie" << endl;
			}else{
				cout << "auto-config orientation: innie" << endl;
			}
		}
		
		meanInsert = library.insertMean;
		insertStd  = library.insertStd;
		max_insert=meanInsert+3*insertStd;
		if(outtie == true){
			max_insert=meanInsert+4*insertStd;
		}
		if(vm["-i"] != ""){
			max_insert  = convert_str( vm["-i"], "-i");
		}else{
			cout << "insert size threshold:" << max_insert << endl;
		}


        
        if (vm["-m"] != ""){
        	min_variant_size = convert_str( vm["-m"],"-m" );
        }
        
        map<string,int> SV_options;
		SV_options["max_insert"]      = max_insert;
		SV_options["pairs"]           = minimumSupportingPairs;
		SV_options["mapping_quality"] = minimum_mapping_quality;
		SV_options["readLength"]      = library.readLength;
		SV_options["ploidy"]          = ploidy;
		SV_options["contigsNumber"]   = contigsNumber;
        SV_options["meanInsert"]      = meanInsert;
        SV_options["STDInsert"]       = insertStd;
        SV_options["min_variant_size"]	  = min_variant_size;
		SV_options["splits"] = minimumSupportingReads;
        
		StructuralVariations *FindTranslocations;
		FindTranslocations = new StructuralVariations();		
		FindTranslocations -> findTranslocationsOnTheFly(alignmentFile, outtie, coverage,outputFileHeader, version ,SV_options);


	//the coverage module
	}else if(vm["--cov"] == "found"){
		Cov *calculateCoverage;
		int binSize=500;
		if(vm["-z"] != ""){
			binSize =convert_str( vm["-z"], "-z");
		}
		if(vm["-o"] != ""){
		    outputFileHeader=vm["-o"];
		}

		calculateCoverage = new Cov(binSize,alignmentFile,outputFileHeader);
		BamReader bam;
		bam.Open(alignmentFile);
		BamAlignment currentRead;
    	while ( bam.GetNextAlignmentCore(currentRead) ) {
	    	calculateCoverage -> bin(currentRead);
		}
		bam.Close();
		calculateCoverage -> printCoverage();
	}

}
