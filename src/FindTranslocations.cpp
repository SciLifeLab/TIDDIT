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
#include  <boost/program_options.hpp>

namespace po = boost::program_options;

//Function used to extract the names of the bam file, index file and output
queue<string> checkFiles(po::variables_map vm){
	queue<string> files;

	// PARSE BAM file
	if (!vm.count("bam")) {
		files.push("error");
		return(files);
	}
	files.push(vm["bam"].as<string>());
	//alignmentFile = vm["bam"].as<string>();
	if (vm.count("bai")){
		files.push(vm["bai"].as<string>());
	}else{
		files.push("NOINDEX");
	}
	if (vm.count("output")) {
		//string header = vm["output"].as<string>();
		//outputFileHeader = header ;
		files.push(vm["output"].as<string>());
	}else{
		files.push("output"); // default output name
	}

	return(files);
}

int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string alignmentFile		    = "";       // alignment file name
	bool outtie 				    = true;	 // library orientation
	uint32_t minimumSupportingPairs = 3;
	int min_insert				    = 100;      // min insert size
	int max_insert				    = 100000;  // max insert size
	int minimum_mapping_quality		=20;
	float coverage;
	float coverageStd;
	float meanInsert;
	float insertStd;
	string roi;
	string indexFile="";
	string outputFileHeader ="";

	//The module menu, this module contains every module in the program
	stringstream ssModule;
	ssModule << endl << "Program: FindTranslocations" << endl <<"Version: beta" << endl << "Usage: FindTranslocations <module> [options]" << endl << "modules";
	po::options_description modules(ssModule.str().c_str());
	modules.add_options()
		("help", "produce help message")
		("cov", "select the cov module to analyse the coverage of regions inside a bam file")
		("extract", "select the extract module to extract regions of interest from a bam file")
		("bamStatistics", "select the bamStatistics module to output statistics of the library")
		("sv", "select the sv module to find structural variations");
	
	//The general menu, this menu contains options that all modules share, such as bam file and output
	stringstream ssGeneral;
	ssGeneral << "Generic Options";
	po::options_description general(ssGeneral.str().c_str());
	general.add_options()
		("bam",po::value<string>(),      "alignment file in bam format, sorted by coordinate. If bwa mem is used option -M MUST be specified in order to map as secondary the splitted reads")
		("output",po::value<string>(),      "Set the name of the output file/folders");
				
	//The sv options menu, this menu is used if the find structural variations menu is selected.
	po::options_description desc("\nUsage: FindTranslocations --sv [Options] --bam inputfile --output prefix(optional) \nOptions");
	desc.add_options()
        ("ploidy",po::value<int>(), "the number of sets of chromosomes,(default = 2)")
		("max-insert",   po::value<int>(),         "paired reads maximum allowed insert size. pairs aligning on the same chr at a distance higher than this are considered candidates for SV(default=1.5std+mean insert size).")
		("orientation",  po::value<string>(),      "expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default: major orientation within the dataset")				
		("minimum-supporting-pairs",  po::value<unsigned int>(), "Minimum number of supporting pairs in order to call a variation event (default 3)")
		("minimum-mapping-quality",  po::value<int>(), "Minimum mapping quality to consider an alignment (default 0)")
		("coverage",     po::value<float>(), "do not compute coverage from bam file, use the one specified here (must be used in combination with --insert --insert-std)");

	po::options_description svModule(" ");
	svModule.add_options()
		("help", "produce help message")
		("sv", "select the sv module to find structural variations");

	//extraction parameters
	po::options_description desc2("\nUsage: FindTranslocations --extract [Options] --bam inputfile --bai index --output outputFolder(optional) \nOptions");
	desc2.add_options()	
		("input-file",       po::value<string>(), "VCF or bed file containing the events(required)")
		("no-singlets"	, "the files containing single events/regions will not be generated")
		("no-total" 	, "the file contining all events will not be generated")
		("filter",       po::value<string>(), "filter the reads that are to be extracted(TODO)")
		("extract-events",	"Select this option if the input file is the vcf output of FindTranslocations(default)")
		("extract-bed",		"select this option if the input file is a bed file");
	po::options_description extractModule(" ");
	extractModule.add_options()
		("help", "produce help message")
		("extract", "select the extract module to extract regions of interest from a bam file");

	//Analyse coverage
	po::options_description desc4("\nUsage: FindTranslocations --cov [Mode] --bam inputfile --output outputFile(optional) \nOptions:only one mode may be selected");
	desc4.add_options()
		("bin",		po::value<int>(), "use bins of specified size to measure the coverage of the entire bam file, set output to stdout to print to stdout")
		("light",		po::value<int>(), "use bins of specified size to measure the coverage of the entire bam file, only prints chromosome and coverage for each bin, set output to stdout to print to stdout");
		
	po::options_description coverageModule(" ");
	coverageModule.add_options()
		("help", "produce help message")
		("cov", "select the coverage module to analyse the coverage of a bam file");



	stringstream ss;
	ss << package_description() << endl;
	po::options_description menu_union(ss.str().c_str());
	menu_union.add(modules).add(general).add(desc).add(desc2).add(desc4);



//TODO: add minimum number of reads to support a translocation, window size etc.

	po::variables_map vm;

	//parse the command line arguments
	try {
		po::store(po::parse_command_line(argc, argv, menu_union), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		return 2;
	}
	
	//if no module is selected a help message is shown
	if( !vm.count("extract") and !vm.count("sv") and !vm.count("cov")){
		DEFAULT_CHANNEL << modules << endl;
		return 2;
	}


	//retrieve the name of the bam file, the bai file and the output
	queue<string> fileQueue=checkFiles(vm);
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;
	uint64_t genomeLength = 0;
	uint32_t contigsNumber = 0;

	if(fileQueue.back() != "error"){
		alignmentFile = fileQueue.front();
		fileQueue.pop();
		indexFile= fileQueue.front();
		fileQueue.pop();
		outputFileHeader = fileQueue.front();

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
			contig2position[sequence->Name] = contigsNumber; // keep track of contig name and position in order to avoid problems when processing two libraries
	
			position2contig[contigsNumber] = sequence->Name;
			contigsNumber++;
		}	
		bamFile.Close();;
	}
	//if the bam extraction module is chosen;
	if (vm.count("extract")){
		
		po::options_description menu_union(ss.str().c_str());
		po::variables_map vm;
		menu_union.add(general).add(desc2).add(extractModule);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
			return 2;
		}

		if(vm.count("help") or fileQueue.back()== "error"){
			DEFAULT_CHANNEL << desc2 << endl;
			return 2;
		}


	if (vm.count("input-file")){
		roi=vm["input-file"].as<string>();
		Extract *BamExtraction;
		BamExtraction = new Extract();
		int exclusion =0;
		if (vm.count("no-singlets")){
			exclusion+=1;
			cout << "yes" << endl;
		}if (vm.count("no-total")){
			exclusion+=2;
		}

		if (vm.count("extract-events")){
			 BamExtraction -> extract(alignmentFile,outputFileHeader,roi,contig2position,indexFile,1,exclusion);
		}
		if (vm.count("extract-bed")){
			 BamExtraction -> extract(alignmentFile,outputFileHeader,roi,contig2position,indexFile,2,exclusion);
		}else{
			 BamExtraction -> extract(alignmentFile,outputFileHeader,roi,contig2position,indexFile,1,exclusion);
		}
		return 1;

	}else{
		DEFAULT_CHANNEL << desc2 << endl;
		return 0;
	}	

	//if the find structural variations module is chosen
	}else if(vm.count("sv")){

		po::options_description menu_union(" ");
		po::variables_map vm;
		menu_union.add(general).add(desc).add(svModule);



		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
			return 2;
		}

		if(vm.count("help") or fileQueue.back()== "error"){
			DEFAULT_CHANNEL << desc << endl;
			return 2;
		}


	
		if (vm.count("minimum-supporting-pairs")) {
			minimumSupportingPairs = vm["minimum-supporting-pairs"].as<unsigned int>();
		}

		LibraryStatistics library;
		size_t start = time(NULL);
		library = computeLibraryStats(alignmentFile, genomeLength, max_insert,40, outtie,minimum_mapping_quality,outputFileHeader);
		printf ("library stats time consumption= %lds\n", time(NULL) - start);
		coverage   = library.C_A;
		if(vm.count("coverage")){
			coverage    = vm["coverage"].as<float>();
		}
		outtie=library.mp;
		if(outtie == true){
			cout << "auto-config orientation: outtie" << endl;
		}else{
			cout << "auto-config orientation: inne" << endl;
		}
		meanInsert = library.insertMean;
		insertStd  = library.insertStd;
		if(outtie == false){
			max_insert =meanInsert+1.5*insertStd;
		}else{
			max_insert =meanInsert+4*insertStd;
		}

        min_insert = meanInsert/2; 
		if(vm.count("max-insert")){
			max_insert  = vm["insert"].as<float>();
		}

        int ploidy = 2;
        if(vm.count("ploidy")){
            ploidy = vm["ploidy"].as<int>();
        }
		StructuralVariations *FindTranslocations;
		FindTranslocations = new StructuralVariations();
		FindTranslocations -> findTranslocationsOnTheFly(alignmentFile, min_insert, max_insert, outtie, minimum_mapping_quality, minimumSupportingPairs, coverage, meanInsert, insertStd, outputFileHeader, indexFile,contigsNumber,ploidy,library.readLength);


	//the coverage module
	}else if(vm.count("cov")){

		po::options_description menu_union(ss.str().c_str());
		po::variables_map vm;
		menu_union.add(general).add(desc4).add(coverageModule);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
			return 2;
		}

		if(vm.count("help") or fileQueue.back()== "error"){
			DEFAULT_CHANNEL << desc4 << endl;
			return 2;
		}

		
		Cov *calculateCoverage;
		
		//calculate the coverage of the entre library
		if(not vm.count("light") and not vm.count("bin")){
			cout << "no mode selected, shuting down" << endl;
			return(0);
		}
		if(vm.count("bin") or vm.count("light") ){
			int option;
			int binSize;
			if(vm.count("bin")){
				option=1;
				binSize =vm["bin"].as<int>();
			}else{
				option=4;
				binSize =vm["light"].as<int>();
			}
			calculateCoverage = new Cov(binSize,alignmentFile,outputFileHeader);
			calculateCoverage -> coverageMain(alignmentFile,outputFileHeader,contig2position,option,binSize);
		}

	}
}
