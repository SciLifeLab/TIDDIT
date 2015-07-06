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
		("extract", "select the extract module to extract regions of interest from a bam file")
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
		("input-file",       po::value<string>(), "tab file containing the events(required)")
		("filter",       po::value<string>(), "filter the reads that are to be extracted(TODO)");
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
		("extract", "select the extract module to extract regions of interest from a bam file");




	stringstream ss;
	ss << package_description() << endl;
	po::options_description menu_union(ss.str().c_str());
	menu_union.add(modules).add(general).add(desc).add(desc2).add(desc3);



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
	if( !vm.count("extract") and !vm.count("sv") and !vm.count("cnv") ){
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
	
			position2contig[contigsNumber] = contig2position[sequence->Name];
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
		 BamExtraction -> extract(alignmentFile,outputFileHeader,roi,contig2position,indexFile);
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

		if (vm.count("coverage") or vm.count("insert") or vm.count("insert-std")) {
			if ( !(vm.count("coverage") and vm.count("insert") and vm.count("insert-std")) ){
				DEFAULT_CHANNEL << "--coverage , --insert, and --insert-std must be all specified at the same time. Please specify all three or let the program compute the numbers (Bam file wil be read twice).\n";
				return 2;
			}
		}

		if (vm.count("coverage") or vm.count("insert") or vm.count("insert-std")) {
			coverage    = vm["coverage"].as<float>();
			meanInsert  = vm["insert"].as<float>();
			insertStd   = vm["insert-std"].as<float>();

		} else {
			LibraryStatistics library;
			library = computeLibraryStats(alignmentFile, genomeLength, max_insert, outtie);

			coverage   = library.C_A;
			meanInsert = library.insertMean;
			insertStd  = library.insertStd;
		}

		StructuralVariations *FindTranslocations;
		FindTranslocations = new StructuralVariations();
		FindTranslocations -> findTranslocationsOnTheFly(alignmentFile, min_insert, max_insert, outtie, minimum_mapping_quality, minimumSupportingPairs, coverage, meanInsert, insertStd, outputFileHeader, indexFile,contigsNumber);


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
	}

}



