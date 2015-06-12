/*
   Francesco Vezzi
 */

#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <string>

#include <sstream>
#include <iostream>
#include <fstream>

#include "data_structures/Translocation.h"
#include  <boost/program_options.hpp>
namespace po = boost::program_options;



void findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality,
		uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs, float coverage, float meanInsertSize, float StdInsertSize, string outputFileHeader);

int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string alignmentFile		    = "";       // alignment file name
	bool outtie 				    = true;	 // library orientation
	uint32_t windowSize 		    = 1000;     // Window size
	uint32_t windowStep 		    = 100;     // Window step
	uint32_t minimumSupportingPairs = 10;
	int min_insert				    = 100;      // min insert size
	int max_insert				    = 1000000;  // max insert size
	string outputFileHeader         = "output"; // default output name
	int minimum_mapping_quality     = 20;
	float coverage;
	float meanInsert;
	float insertStd;

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
						("bam",          po::value<string>(),      "alignment file in bam format, sorted by coordinate. If bwa mem is used option -M MUST be specified in order to map as secondary the splitted reads")
						("min-insert",   po::value<int>(),         "paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goes from beginning of first read to end of second read")
						("max-insert",   po::value<int>(),         "paired reads maximum allowed insert size. pairs aligning on the same chr at a distance higher than this are considered candidates for SV.")
						("orientation",  po::value<string>(),      "expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default outtie")
						("output",       po::value<string>(),      "Header output file names")
						("minimum-supporting-pairs",  po::value<unsigned int>(), "Minimum number of supporting pairs in order to call a variation event (default 10)")
						("minimum-mapping-quality",  po::value<int>(), "Minimum mapping quality to consider an alignment (default 20)")
						("window-size",  po::value<unsigned int>(),    "Size of the sliding window (default 1000)")
						("window-step",  po::value<unsigned int>(),    "size of the step in overlapping window (must be lower than window-size) (default 100)")
						("coverage",     po::value<float>(), "do not compute coverage from bam file, use the one specified here (must be used in combination with --insert --insert-std)")
						("insert",       po::value<float>(), "do not compute insert size from bam file, use the one specified here (must be used in combination with --coverage --insert-std)")
						("insert-std",   po::value<float>(), "do not compute insert size standard deviation from bam file, use the one specified here (must be used in combination with --insert --coverage)")
						;
//TODO: add minimum number of reads to support a translocation, window size etc.

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		return 2;
	}

	if (vm.count("help")) {
		DEFAULT_CHANNEL << desc << endl;
		return 2;
	}

	// PARSE BAM file
	if (!vm.count("bam")) {
		DEFAULT_CHANNEL << "Please specify --bam " << endl;
		return 2;
	}
	alignmentFile = vm["bam"].as<string>();

	if (!vm.count("max-insert")) {
		DEFAULT_CHANNEL << "Please specify max-insert " << endl;
		return 2;
	}

	if (vm.count("min-insert")) {
		min_insert = vm["min-insert"].as<int>();
		if(min_insert <= 0) {
			DEFAULT_CHANNEL << "minimum insert should be at least 1\n";
			DEFAULT_CHANNEL << desc << endl;
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

	if (vm.count("output")) {
		string header = vm["output"].as<string>();
		outputFileHeader = header ;
	}

	if (vm.count("minimum-mapping-quality")) {
		minimum_mapping_quality = vm["minimum-mapping-quality"].as<int>();
	}

	if (vm.count("window-size")) {
		windowSize = vm["window-size"].as<unsigned int>();
	}

	if (vm.count("window-step")) {
		windowStep = vm["window-step"].as<unsigned int>();
		if (windowStep > windowSize) {
			DEFAULT_CHANNEL << "window-step cannot be larger than window-size\n";
			return 2;
		}
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


	// Now parse BAM header and extract information about genome lenght
	uint64_t genomeLength = 0;
	uint32_t contigsNumber = 0;
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
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;

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

	findTranslocationsOnTheFly(alignmentFile, min_insert, max_insert, outtie, minimum_mapping_quality,
			windowSize, windowStep, minimumSupportingPairs, coverage, meanInsert, insertStd, outputFileHeader);

}


void findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality,
		uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs, float meanCoverage, float meanInsertSize, float StdInsertSize, string outputFileHeader) {
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	// now create Translocation on the fly
	Window *window;

	window = new Window(windowSize, windowStep, max_insert, minimum_mapping_quality,
		outtie,  meanInsertSize,  StdInsertSize,  minimumSupportingPairs,
		 meanCoverage,  outputFileHeader);
	window->initTrans(head);

	//Initialize bam entity
	BamAlignment currentRead;
	//now start to iterate over the bam file
	int counter = 0;
	while ( bamFile.GetNextAlignment(currentRead) ) {
		if(currentRead.IsMapped()) {
			window->insertRead(currentRead);
		}
	}
	if(window->windowOpen) {
		window->computeVariations();
	}
	window->interChrVariations.close();
	window->intraChrVariations.close();

}


