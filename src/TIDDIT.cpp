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
	//MAIN VARIABLES

	bool outtie = true;	 // library orientation
	uint32_t minimumSupportingPairs = 4;
	uint32_t minimumSupportingReads = 6;
	int min_insert = 100;      // min insert size
	int max_insert = 80000;  // max insert size
	int minimum_mapping_quality = 10;
	float coverage;
	float coverageStd;
	float meanInsert;
	float insertStd;
	int min_variant_size= 100;
	int sample = 100000000;
	string outputFileHeader ="output";
	string version = "2.7.1";
	
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
	general_help+="\t--sv\tcollect SV signals\n";
	general_help+="\t--cov\tselect the cov module to analyse the coverage of the genome using bins of a specified size\n";
        general_help+="\t--gc\tselect the gc module to compute the gc content across the genome using bins of a specified size(accepts a fasta through stdin)\n";
	
	//the sv module
	vm["--sv"]="store";
	vm["-i"]="";
	vm["-d"]="";
	vm["-p"]="";
	vm["-r"]="";
	vm["-q"]="";
	vm["-c"]="";
	vm["-s"]="";
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
	sv_help +="\t-s\tNumber of reads to sample when computing library statistics, (default= 100000000)\n";
	sv_help +="\t-n\tthe number of sets of chromosomes,(default = 2)\n";
	sv_help +="\t-m\tminimum variant size,(default = 100)\n";
			
	//the coverage module
	vm["--cov"]="store";
	vm["-z"]="";
	vm["-w"] = "store";
	vm["-u"] = "store";
	vm["-a"] = "store";

	//the gc module
        vm["--gc"]="store";
	
	string coverage_help="\nUsage: TIDDIT --cov [Mode] -b inputfile [-o prefix]\n";
	coverage_help +="\t-b\tcoordinate sorted bam file(required)\n";
	coverage_help +="\n\t-z\tuse bins of specified size(default = 500bp) to measure the coverage of the entire bam file, set output to stdout to print to stdout\n";
	coverage_help +="\n\t-w\tOutput wig instead of bed\n";
	coverage_help +="\n\t-u\tSkip quality values\n";
	coverage_help +="\n\t-a\tSkip print breadth of coverage (reads and spanning pairs, only possible with wig)\n";
	
	string gc_help="\nUsage: cat in.fa | TIDDIT --gc [Mode] [-o prefix]\n";
	gc_help+="\t-r\treference fasta file(required)\n";
	gc_help+="\n\t-z\tuse bins of specified size(default = 500bp) to measure the coverage of the entire bam file, set output to stdout to print to stdout\n";


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

	string argString = "";

	//store the options in a map
	for(int i = 0; i < arguments.size(); i ++){
		argString += " " + arguments[i];
	}

	//print help message
	if(vm["--help"] == "found" or arguments.size() == 0){
		if( vm["--sv"] == "found"){
			cout << sv_help;
		}else if(vm["--cov"] == "found"){
			cout << coverage_help;
		}else if(vm["--cov"] == "found"){
			cout << gc_help;
		}else{
			cout << general_help;
		}
		return(0);
	}
	
	//if no module is chosen
	if(vm["--sv"] == "store" and vm["--cov"] == "store" and vm["--gc"] == "store" ){
		cout << general_help;
		cout << "ERROR: select a module" << endl;
		return(0);
	}else if( (vm["--sv"] == "found" and vm["--cov"] == "found") or (vm["--sv"] == "found" and vm["--gc"] == "found") or (vm["--cov"] == "found" and vm["--gc"] == "found")  ){
		cout << general_help;
		cout << "ERROR: only one module may be selected" << endl;
		return(0);
	
	}

	if (vm["--gc"] == "found"){
		int binSize=500;
		if(vm["-z"] != ""){
			binSize =convert_str( vm["-z"], "-z");
		}

		string fasta = vm["-r"];
		string outputfile="output";
		if(vm["-o"] != ""){
		    outputfile=vm["-o"];
		}
	
		ifstream in(fasta.c_str());

		string chromosome="";
		bool first = true;
		long pos = 0;
		long element =0;
		vector <string> chromosomes;
		map<string, vector< vector<int> > > bins;

		int A =0;
		int G =0;
		int other=0;
		for (string line; getline(cin, line);) {
			if (line[0] == '>'){

				if (A > 0 or G > 0 or other > 0 ){
					vector<int> bin;
					bin.push_back(A);
					bin.push_back(G);
					bin.push_back(other);
					bins[chromosome].push_back(bin);
				}

				chromosome="";
				for (int i=1;i < line.size();i++){
					if (line[i] == ' '){
						break;
					}
					chromosome += line[i];
				}
				chromosomes.push_back(chromosome);
				vector< vector<int> > row;
				bins[chromosome]=row;
				pos=0;
				element=0;
				A=0;
				G=0;
				other=0;
			}else{
				for(int i=0; i< line.size();i++){
					pos+=1;
					if (line[i] == 'A' or line[i] =='a'){
						A+=1;
					}else if (line[i] == 'C' or line[i] =='c'){
						G+=1;
					}else if (line[i] == 'g' or line[i] == 'G'){
						G+=1;
					}else if (line[i] == 't' or line[i] == 'T'){
						A+=1;
					}else{
						other+=1;
					}
			
					if (pos >= (element+1)*binSize){
						vector<int> bin;
						bin.push_back(A);
						bin.push_back(G);
						bin.push_back(other);
						bins[chromosome].push_back(bin);
						element+=1;
						A=0;G=0;other=0;
					}
				}
			}

		}

		if (A > 0 or G > 0 or other > 0 ){
			vector<int> bin;
			bin.push_back(A);
			bin.push_back(G);
			bin.push_back(other);
			bins[chromosome].push_back(bin);
		}


		if (outputfile == "stdout"){
			cout << "track type=wiggle_0 name=\"GC\" description=\"Per bin GC values\"" << endl;
			for(int i=0;i< chromosomes.size();i++){
				cout << "fixedStep chrom=" << chromosomes[i] << " start=1 step=50" << endl;  
	 			for (int j=0;j<bins[chromosomes[i]].size();j++){
					float gc =0;
					if (bins[chromosomes[i]][j][1]+bins[chromosomes[i]][j][0] > 0){
						gc=(float)bins[chromosomes[i]][j][1]/(bins[chromosomes[i]][j][1]+bins[chromosomes[i]][j][0]);
					}
					cout << gc << "\n";
				}
			}

			cout << "track type=wiggle_0 name=\"N-count\" description=\"Per bin fraction of N\"" << endl;
			for(int i=0;i< chromosomes.size();i++){
				cout << "fixedStep chrom=" << chromosomes[i] << " start=1 step=50" << endl;  
	 			for (int j=0;j<bins[chromosomes[i]].size();j++){
					float n = 0;
					if (bins[chromosomes[i]][j][0]+bins[chromosomes[i]][j][1] > 0){
						n=(float)bins[chromosomes[i]][j][2]/( bins[chromosomes[i]][j][0]+bins[chromosomes[i]][j][1]+bins[chromosomes[i]][j][2] );
					}else{
						n=1;
					}
					cout << n << endl;
					
				}
			}

		}else{
			ofstream gcOutput;
			gcOutput.open((outputfile+".gc.wig").c_str());
			ostream& gcout=gcOutput;
			gcout << "track type=wiggle_0 name=\"GC\" description=\"Per bin GC values\"" << endl;
			for(int i=0;i< chromosomes.size();i++){
				gcout << "fixedStep chrom=" << chromosomes[i] << " start=1 step=50" << endl;  
	 			for (int j=0;j<bins[chromosomes[i]].size();j++){
					float gc =0;
					if (bins[chromosomes[i]][j][1]+bins[chromosomes[i]][j][0] > 0){
						gc=(float)bins[chromosomes[i]][j][1]/(bins[chromosomes[i]][j][1]+bins[chromosomes[i]][j][0]);
					}
					gcout << gc << "\n";
				}
			}

			gcout << "track type=wiggle_0 name=\"N-count\" description=\"Per bin fraction of N\"" << endl;
			for(int i=0;i< chromosomes.size();i++){
				gcout << "fixedStep chrom=" << chromosomes[i] << " start=1 step=50" << endl;  
	 			for (int j=0;j<bins[chromosomes[i]].size();j++){
					float n = 0;
					if (bins[chromosomes[i]][j][0]+bins[chromosomes[i]][j][1] > 0){
						n=(float)bins[chromosomes[i]][j][2]/( bins[chromosomes[i]][j][0]+bins[chromosomes[i]][j][1]+bins[chromosomes[i]][j][2] );
					}else{
						n=1;
					}
					gcout << n << endl;
					
				}
			}

		}
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
		if(vm["-s"] != ""){
			sample=convert_str( vm["-s"] , "-s");
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
		library = computeLibraryStats(alignmentFile, genomeLength, max_insert, 50 , outtie,minimum_mapping_quality,outputFileHeader,sample);
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
		max_insert=meanInsert+2.1*insertStd;
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
		FindTranslocations -> findTranslocationsOnTheFly(alignmentFile, outtie, coverage,outputFileHeader, version, "TIDDIT" + argString,SV_options, genomeLength);

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

		bool wig = false;
		if(vm["-w"] == "found"){
		    wig = true;
		}

		bool skipQual= false;
		if(vm["-u"] == "found"){
		     skipQual = true;
		}

		bool span = false;
		if(vm["-a"] == "found"){
		     span = true;
		}

		calculateCoverage = new Cov(binSize,alignmentFile,outputFileHeader,0,wig,skipQual,span);
		BamReader bam;
		bam.Open(alignmentFile);
		BamAlignment currentRead;
		while ( bam.GetNextAlignmentCore(currentRead) ) {
			readStatus alignmentStatus = computeReadType(currentRead, 100000,100, true);
			calculateCoverage -> bin(currentRead, alignmentStatus);
		}
		bam.Close();
		calculateCoverage -> printCoverage();
	}

}
