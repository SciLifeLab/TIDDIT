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
	if (!vm.count("bai")){
		files.push("error");
		return(files);
	}
	//string indexFile=vm["bai"].as<string>();
	files.push(vm["bai"].as<string>());
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
	uint32_t minimumSupportingPairs = 10;
	int min_insert				    = 100;      // min insert size
	int max_insert				    = 1000000;  // max insert size
	int minimum_mapping_quality     = 20;
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
		("cnv", "Select the cnv module to find copy number variations ")
		("cov", "select the cov module to analyse the coverage of regions inside a bam file")
		("extract", "select the extract module to extract regions of interest from a bam file")
		("region", "select the region module to search for a certain variation given in a region file, read the README file for more information on the format of the region file")
		("sv", "select the sv module to find structural variations");
	
	//The general menu, this menu contains options that all modules share, such as bam file and output
	stringstream ssGeneral;
	ssGeneral << "Generic Options";
	po::options_description general(ssGeneral.str().c_str());
	general.add_options()
		("bai",po::value<string>(),	"the reference file")
		("bam",po::value<string>(),      "alignment file in bam format, sorted by coordinate. If bwa mem is used option -M MUST be specified in order to map as secondary the splitted reads")
		("output",po::value<string>(),      "Set the name of the output file/folders");
				
	//The sv options menu, this menu is used if the find structural variations menu is selected.
	po::options_description desc("\nUsage: FindTranslocations --sv [Options] --bam inputfile --bai indexfile --output outputFile(optional) \nOptions");
	desc.add_options()
		("auto","find min-insert,max-insert and orientation based on the statistics of the bam input file")
        ("ploidy",po::value<int>(), "the number of sets of chromosomes,(default = 2)")
		("min-insert",   po::value<int>(),         "paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goes from beginning of first read to end of second read")
		("max-insert",   po::value<int>(),         "paired reads maximum allowed insert size. pairs aligning on the same chr at a distance higher than this are considered candidates for SV.")
		("orientation",  po::value<string>(),      "expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default outtie")				
		("minimum-supporting-pairs",  po::value<unsigned int>(), "Minimum number of supporting pairs in order to call a variation event (default 10)")
		("minimum-mapping-quality",  po::value<int>(), "Minimum mapping quality to consider an alignment (default 20)")
		("coverage",     po::value<float>(), "do not compute coverage from bam file, use the one specified here (must be used in combination with --insert --insert-std)")
		("insert",       po::value<float>(), "do not compute insert size from bam file, use the one specified here (must be used in combination with --coverage --insert-std)")
		("insert-std",   po::value<float>(), "do not compute insert size standard deviation from bam file, use the one specified here (must be used in combination with --insert --coverage)");

	po::options_description svModule(" ");
	svModule.add_options()
		("help", "produce help message")
		("sv", "select the sv module to find structural variations");

	//extraction parameters
	po::options_description desc2("\nUsage: FindTranslocations --extract [Options] --bam inputfile --bai indexfile --output outputFolder(optional) \nOptions");
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


	//find copy number variation parameters
	po::options_description desc3("\nUsage: FindTranslocations --cnv [Options] --bam inputfile --bai indexfile --output outputFile(optional) \nOptions");
	desc3.add_options()
		("TODO",       po::value<string>(), "analyse copy number variations"" ");
	po::options_description cnvModule(" ");
	cnvModule.add_options()
		("help", "produce help message")
		("cov", "select the cov module to analyse the coverage of regions inside a bam file");

	//Analyse coverage
	po::options_description desc4("\nUsage: FindTranslocations --cov [Mode] --bam inputfile --bai indexfile --output outputFile(optional) \nOptions:only one mode may be selected");
	desc4.add_options()
		("bin",		po::value<int>(), "use bins of specified size to measure the coverage of the entire bam file")
		("light",		po::value<int>(), "use bins of specified size to measure the coverage of the entire bam file, only prints chromosome and coverage for each bin")
		("intra-vcf", 	po::value<string>() ,"Select this option if the input file is the output intra-chromosomal vcf of Findtranslocations. coverage will be calculated between, before and after the event")
		("inter-vcf", 	po::value<string>() ,"Select this option if the input file is the output inter-chromosomal vcf of Findtranslocations. The coverage before and after window A and B will be calculated")
		("bed", 	po::value<string>() ,"Select this option if the input file is a bedfile");
		
	po::options_description coverageModule(" ");
	coverageModule.add_options()
		("help", "produce help message")
		("cov", "select the coverage module to analyse the coverage of a bam file");
	//used to control the input files

	//analyse a certain variation within a given region
	po::options_description desc5("\nUsage: FindTranslocations --region [options] --bam inputfile --bai indexfile --output outputFile(optional)");
	desc5.add_options()
		("region-file", 	po::value<string>() ,"Select this option if the input file is a bedfile")
		("vcf", 	po::value<string>() ,"Select this option if the input file is a vcf file")
		("ploidy",po::value<int>(), "the number of sets of chromosomes,(default = 2)");
	po::options_description regionModule(" ");
	regionModule.add_options()
		("help", "produce help message")
		("region", "select the region module to search for a certain variation given in a region file, read the README file for more information on the format of the region file");
	//used to control the input files




	stringstream ss;
	ss << package_description() << endl;
	po::options_description menu_union(ss.str().c_str());
	menu_union.add(modules).add(general).add(desc).add(desc2).add(desc3).add(desc4).add(desc5);



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
	if( !vm.count("extract") and !vm.count("sv") and !vm.count("cnv") and !vm.count("cov") and !vm.count("region") ){
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
			cout << sortOrder << "\n";
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
		bamFile.Close();
	
		cout << "total number of contigs " 	<< contigsNumber << endl;
		cout << "assembly length " 			<< genomeLength << "\n";
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

		

		if (vm.count("min-insert")) {
			min_insert = vm["min-insert"].as<int>();
			if(min_insert <= 0) {
				DEFAULT_CHANNEL << desc << endl;
				DEFAULT_CHANNEL << "minimum insert should be at least 1\n";
				return 2;
			}
		}

		if (vm.count("max-insert")) {
			max_insert = vm["max-insert"].as<int>();
		}

		if (vm.count("orientation")) {
			string userOrientation = vm["orientation"].as<string>();
			if(userOrientation.compare("innie") == 0) {
				outtie = false;
			} else if (userOrientation.compare("outtie") == 0) {
				outtie = true;
			} else {
				DEFAULT_CHANNEL << "outtie (<- ->) or innie (-> <-) only allowed orientations\n";
				DEFAULT_CHANNEL << desc << endl;
				return 2;
			}
		}

		if (vm.count("minimum-mapping-quality")) {
			minimum_mapping_quality = vm["minimum-mapping-quality"].as<int>();
		}
	
		if (vm.count("minimum-supporting-pairs")) {
			minimumSupportingPairs = vm["minimum-supporting-pairs"].as<unsigned int>();
		}

		if (vm.count("coverage") or vm.count("insert") or vm.count("insert-std") and not vm.count("auto") ) {
			if ( !(vm.count("coverage") and vm.count("insert") and vm.count("insert-std")) ){
				DEFAULT_CHANNEL << "--coverage , --insert, and --insert-std must be all specified at the same time. Please specify all three or let the program compute the numbers (Bam file wil be read twice).\n";
				return 2;
			}
		}

		if(vm.count("auto")){
			vector<int> libraryStats;
			//compute max_insert,mean insert and outtie
			autoSettings *autoStats;
			autoStats = new autoSettings();
			size_t start = time(NULL);
			libraryStats=autoStats->autoConfig(alignmentFile,minimum_mapping_quality);
			
			printf ("auto config time consumption= %lds\n", time(NULL) - start);
			max_insert=libraryStats[0];
			min_insert=libraryStats[1];
			outtie=libraryStats[2] != 0;

		}


		if (vm.count("coverage") or vm.count("insert") or vm.count("insert-std") and not vm.count("auto")) {
			coverage    = vm["coverage"].as<float>();
			meanInsert  = vm["insert"].as<float>();
			insertStd   = vm["insert-std"].as<float>();

		 }else {
			LibraryStatistics library;
			size_t start = time(NULL);
			library = computeLibraryStats(alignmentFile, genomeLength, max_insert,min_insert, outtie);
			printf ("library stats time consumption= %lds\n", time(NULL) - start);
			coverage   = library.C_A;
			meanInsert = library.insertMean;
			insertStd  = library.insertStd;
			//update the max_insert
			if(vm.count("auto")){
				max_insert =meanInsert+4*insertStd;
			}

		}
		
        int ploidity = 2;
        if(vm.count("ploidiy")){
            ploidity = vm["ploidiy"].as<int>();
        }
		size_t start = time(NULL);
		Cov *calculateCoverage;
		calculateCoverage = new Cov();
		calculateCoverage -> coverageMain(alignmentFile,indexFile,"",outputFileHeader,coverage,contig2position,4,200);
		printf ("time consumption of the coverage computation= %lds\n", time(NULL) - start);

		StructuralVariations *FindTranslocations;
		FindTranslocations = new StructuralVariations();
		FindTranslocations -> findTranslocationsOnTheFly(alignmentFile, min_insert, max_insert, outtie, minimum_mapping_quality, minimumSupportingPairs, coverage, meanInsert, insertStd, outputFileHeader, indexFile,contigsNumber,ploidity);


	//if the find copy number variation module is chosen
	}else if(vm.count("cnv")){

		po::options_description menu_union(ss.str().c_str());
		po::variables_map vm;
		menu_union.add(general).add(desc3).add(cnvModule);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			DEFAULT_CHANNEL << desc3 << endl;
			return 2;
		}
		DEFAULT_CHANNEL << desc3 << endl;

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
		calculateCoverage = new Cov();
		LibraryStatistics library;
		//calculate the coverage of the entre library
		if(vm.count("bin") or vm.count("bed") or vm.count("intra-vcf") or vm.count("inter-vcf")){
			library = computeLibraryStats(alignmentFile, genomeLength, max_insert,min_insert ,outtie);
			coverage   = library.C_A;
		}else if(not vm.count("light")){
			cout << "no mode selected, shuting down" << endl;
			return(0);
		}
				
		int option=0;
		int binSize=0;
		if(vm.count("bin") or vm.count("light") ){
			if(vm.count("bin")){
				option=1;
				binSize =vm["bin"].as<int>();
			}else{
				option=4;
				binSize =vm["light"].as<int>();
			}
			calculateCoverage -> coverageMain(alignmentFile,indexFile,"",outputFileHeader,coverage,contig2position,option,binSize);
		}else if(vm.count("bed")){
			option=2;
			string inFile=vm["bed"].as<string>();
			calculateCoverage -> coverageMain(alignmentFile,indexFile,inFile,outputFileHeader,coverage,contig2position,option,binSize);
		}else if(vm.count("intra-vcf")){
			option=0;
			string inFile=vm["intra-vcf"].as<string>();
			calculateCoverage -> coverageMain(alignmentFile,indexFile,inFile,outputFileHeader,coverage,contig2position,option,binSize);
		}else if(vm.count("inter-vcf")){
			string inFile=vm["inter-vcf"].as<string>();
			option=3;
			calculateCoverage -> coverageMain(alignmentFile,indexFile,inFile,outputFileHeader,coverage,contig2position,option,binSize);
		}
	}else if( vm.count("region") ){

		po::options_description menu_union(" ");
		po::variables_map vm;
		menu_union.add(general).add(desc5).add(regionModule);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
			return 2;
		}

		if(vm.count("help") or fileQueue.back()== "error"){
			DEFAULT_CHANNEL << desc5 << endl;
			return 2;
		}

        int ploidy = 2;
        if(vm.count("ploidiy")){
            ploidy = vm["ploidiy"].as<int>();
        }

		Region *AnalyseRegion;
		AnalyseRegion = new Region();

		string regionFile;
		//use either the bed like input or a vcf file
		if(vm.count("region-file") or vm.count("vcf") and not ( vm.count("region-file") and vm.count("vcf") ) ){
			if(vm.count("region-file")){
				AnalyseRegion -> input = "bed";
				regionFile= vm["region-file"].as<string>();
			}else{
				AnalyseRegion -> input = "vcf";
				regionFile= vm["vcf"].as<string>();
			}
		}else{
			ERROR_CHANNEL <<  "please enter a region file" << endl;
			return 2;
		}

		cout << "searching the regions...." << endl;



		AnalyseRegion -> region(alignmentFile,indexFile,regionFile,outputFileHeader,contig2position,position2contig,ploidy,minimum_mapping_quality,genomeLength,contigsNumber);

	}


}



