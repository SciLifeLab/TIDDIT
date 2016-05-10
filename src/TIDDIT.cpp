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
	int max_insert				    = 100000;  // max insert size
	int minimum_mapping_quality		= 20;
	float coverage;
	float coverageStd;
	float meanInsert;
	float insertStd;
	string outputFileHeader ="output";
	
	
	//collect all options as a vector
	vector<string> arguments(argv, argv + argc);
	arguments.erase (arguments.begin());

	map<string,string> vm;
	//general options
	vm["--bam"]="";
	vm["--output"]="";
	vm["--help"] ="store";
	
	string general_help="TIDDIT - a structural variant caller\nUsage: TIDDIT <module> [options] \n";
	general_help+="modules\n\t--help\tproduce help message\n";
	general_help+="\t--sv\tselect the sv module to find structural variations\n";
	general_help+="\t--cov\tselect the cov module to analyse the coverage of the genome using bins of a specified size\n";
	
	//the sv module
	vm["--sv"]="store";
	vm["--insert"]="";
	vm["--orientation"]="";
	vm["--supporting-pairs"]="";
	vm["--mapping-quality"]="";
	vm["--coverage"]="";
	vm["--plody"]="";
	
	string sv_help ="\nUsage: TIDDIT --sv --bam inputfile [--output prefix] \nOther options\n";
	sv_help +="\t--insert\tpaired reads maximum allowed insert size. Pairs aligning on the same chr at a distance higher than this are considered candidates for SV (if not specified default=1.5std + mean_insert_size)\n";
	sv_help +="\t--orientation\texpected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default: major orientation within the dataset\n";
	sv_help +="\t--supporting-pairs\tMinimum number of supporting pairs in order to call a variation event (default 3)\n";
	sv_help +="\t--mapping-quality\tMinimum mapping quality to consider an alignment (default 0)\n";
	sv_help +="\t--coverage\tdo not compute coverage from bam file, use the one specified here\n";
	sv_help +="\t--plody\tthe number of sets of chromosomes,(default = 2)\n";
			
	//the coverage module
	vm["--cov"]="store";
	vm["--bin_size"]="";
	
	string coverage_help="\nUsage: TIDDIT --cov [Mode] --bam inputfile [--output prefix] \nOptions:only one mode may be selected\n";
	coverage_help +="\n\t-bin_size\tuse bins of specified size(default = 500bp) to measure the coverage of the entire bam file, set output to stdout to print to stdout\n";
	
	
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
	if(vm["--bam"] == ""){
		cout << general_help;
		cout << "ERROR: select a bam file using the --bam option" << endl;
		return(1);
	}

	string alignmentFile	= vm["--bam"];
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
		if(vm["--output"] != ""){
		    outputFileHeader=vm["--output"];
		}
	
	
		if(vm["--mapping-quality"] != ""){
			minimum_mapping_quality=convert_str( vm["--mapping-quality"] ,"--mapping-quality");
		}
		if(vm["-p"] != ""){
			minimumSupportingPairs=convert_str( vm["--supporting-pairs"] , "--supporting-pairs");
		}
		if(vm["--insert"] != ""){
			int insert_test  = convert_str( vm["--insert"], "--insert");
		}
		
		if(vm["--orientation"] != ""){
			if (vm["--orientation"] == "outtie"){
				outtie=true;
			}else if (vm["--orientation"] == "innie"){
				outtie=false;
			}else{
				cout << "ERROR: invalid orientation " << vm["--orientation"] << endl;
				return(0);
			}
		}
		if(vm["--plody"] != ""){
            		int ploidy = convert_str( vm["--plody"],"--plody" );
        	}
		
		LibraryStatistics library;
		size_t start = time(NULL);
		library = computeLibraryStats(alignmentFile, genomeLength, max_insert,40, outtie,minimum_mapping_quality,outputFileHeader);
		printf ("library stats time consumption= %lds\n", time(NULL) - start);
		coverage   = library.C_A;
		if(vm["--coverage"] != ""){
			coverage    = convert_str( vm["--coverage"],"--coverage" );
		}
		
		if(vm["--orientation"] == ""){
            outtie=library.mp;
			if(outtie == true){
				cout << "auto-config orientation: outtie" << endl;
			}else{
				cout << "auto-config orientation: innie" << endl;
			}
		}
		
		meanInsert = library.insertMean;
		insertStd  = library.insertStd;
		if(outtie == false){
			max_insert =meanInsert+3*insertStd;
		} else {
			max_insert =meanInsert+4*insertStd;
		}

		if(vm["--insert"] != ""){
			max_insert  = convert_str( vm["--insert"], "--insert");
		}

        	int ploidy = 2;
        	if(vm["--plody"] != ""){
            		ploidy = convert_str( vm["--plody"],"--plody" );
        	}
        
        map<string,int> SV_options;
		SV_options["max_insert"]=max_insert;
		SV_options["pairs"]=minimumSupportingPairs;
		SV_options["mapping_quality"]=minimum_mapping_quality;
		SV_options["readLength"]=library.readLength;
		SV_options["ploidy"]=ploidy;
		SV_options["contigsNumber"]=contigsNumber;
        SV_options["meanInsert"]=meanInsert;
        SV_options["STDInsert"]=insertStd;
        
		StructuralVariations *FindTranslocations;
		FindTranslocations = new StructuralVariations();		
		FindTranslocations -> findTranslocationsOnTheFly(alignmentFile, outtie, coverage,outputFileHeader,SV_options);


	//the coverage module
	}else if(vm["--cov"] == "found"){
		Cov *calculateCoverage;
		int binSize=500;
		if(vm["--bin_size"] != ""){
			binSize =convert_str( vm["--bin_size"], "--bin_size");
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
	}

}
