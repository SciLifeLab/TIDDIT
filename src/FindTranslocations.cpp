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



#include "data_structures/Translocation.h"
#include  <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "api/BamWriter.h"

namespace po = boost::program_options;

//returns a bam file contaning the reads located in a given region of the genome
void TranslocationBAM(string BamFileName,string outputFileHeader,string inputFileName,map<string,unsigned int> contig2position, string indexFile){
	//the number of regions involved in each event
	int n=2;
	int j=0;
	string tmpstr;
	string line;
	string eventbam;
	string eventfolder;

	//creates a folder for the outputdata
	boost::filesystem::path dir(outputFileHeader.c_str());
	if(!boost::filesystem::create_directory(dir)) {
		std::cout << "error creating output folder" << "\n";
		return;
	}
	eventfolder=outputFileHeader+"/events";
	boost::filesystem::path eventdir( eventfolder.c_str());
	if(!boost::filesystem::create_directory(eventdir)) {
		std::cout << "error creating output folder" << "\n";
		return;
	}

	//opens one writer for writing one file containing all the output sequences, and one writer to create files containing each separate region
	BamWriter regionwriter;
	BamWriter writer;
	BamReader bamFile;
	eventbam=outputFileHeader+"/all_events_roi.bam";
	//open the bam file
	bamFile.Open(BamFileName);
	BamAlignment currentRead;

	// retrieve 'metadata' from BAM files, these are required by BamWriter
	const SamHeader header = bamFile.GetHeader();
	const RefVector references = bamFile.GetReferenceData();

	// attempt to open our BamWriter{
	cout << eventbam << endl;
	if ( !writer.Open(eventbam.c_str(), header, references) ) {
		cout << "Could not open output BAM file" << endl;
		return;
	}
	//test if an index file is available
	if(bamFile.OpenIndex(indexFile) == 0){
		cout << "Failed to load the index file" << endl;
	}





	//opens the input file and reads each line
	cout << "initiates extraction, please wait" << endl;
	ifstream inputFile( inputFileName.c_str() );
	if (inputFile){
		while (getline( inputFile, line )){
			//skips the input header
			if (j > 0){
				//splits on tab
				vector<std::string> splitline;
				boost::split(splitline, line, boost::is_any_of("\t"));

				queue<int> chr;
				queue<int> startPos;queue<int> endPos;
				//the filename of the file containing single regions
				string regionfile = "";

				try{
					//extracts the startpos,endpos and chromosome of each given region
					for(int i=0;i<=1;i++){
						//extracts the chromosome REFID
						tmpstr=splitline.front();
						splitline.erase( splitline.begin() );
						chr.push(contig2position[tmpstr]);
						regionfile+=tmpstr;
						regionfile+="_";
			
						//extracts the start position on the given chromosome
						tmpstr=splitline.front();
						splitline.erase( splitline.begin() );
						startPos.push( atoi(tmpstr.c_str()));
						regionfile+=tmpstr;
						regionfile+="_";

						//extracts the end position on the given chromosome
						tmpstr=splitline.front();
						splitline.erase( splitline.begin() );
						endPos.push(atoi(tmpstr.c_str()));
						regionfile+=tmpstr;
						if(i == 0){regionfile+="_";}


						if(startPos.back() > endPos.back()){
							cout << "ERROR: chromosome startpos must be located before endpos" << endl;
							return;
						}
			
					}
				}catch(...){
					cout << "error invalid --extract input" << endl;
					return;
				}

				regionfile+=".bam";
				regionfile=eventfolder+"/"+regionfile;
				cout << regionfile << endl;
					// attempt to open our BamWriter{
					if ( !regionwriter.Open(regionfile.c_str(), header, references) ) {
						cout << "Could not open output event BAM file" << endl;
						return;
					}
				





				//collects all reads found in the given intervals and writes them to the file outputFileHeader+"_roi.bam"
				for(int i=0;i<n;i++){
					bamFile.SetRegion(chr.front(),startPos.front(),chr.front(),endPos.front()); 
					chr.pop();startPos.pop();endPos.pop();
					while ( bamFile.GetNextAlignment(currentRead) ) {
						if(currentRead.IsMapped()) {
						//prints the read to bam file
						writer.SaveAlignment(currentRead);
						regionwriter.SaveAlignment(currentRead);
						}
					}
				}
				regionwriter.Close();
				

			}
			j++;
		}
		
		inputFile.close();
	}

	cout << "writing complete" << endl;
	cout << outputFileHeader << endl;
	writer.Close();
	bamFile.Close();
	

}


void findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality,
		uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs, float coverage, float meanInsertSize, float StdInsertSize, string outputFileHeader, string Indexfile);

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
	string roi;

	stringstream ssGeneral;
	ssGeneral << "General" << endl << endl << "Allowed options";
	po::options_description general(ssGeneral.str().c_str());
	general.add_options() ("help", "produce help message")
	("bai",		po::value<string>(),	"the reference file")
	("bam",          po::value<string>(),      "alignment file in bam format, sorted by coordinate. If bwa mem is used option -M MUST be specified in order to map as secondary the splitted reads")
	("module-help", po::value<string>(), "produce the help message of a module, \n --module-help X \n X =sv,extract,cnv or general")
	("output",       po::value<string>(),      "Header output file names");
	
	
	 

	// PROCESS PARAMETERS


	//find structural variations parameters
	stringstream ssSV;
	ssSV << "sv module" << endl << endl << "Allowed options";
	po::options_description desc(ssSV.str().c_str());
	desc.add_options()			("sv", "select the sv module to find structural variations")
						("min-insert",   po::value<int>(),         "paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goes from beginning of first read to end of second read")
						("max-insert",   po::value<int>(),         "paired reads maximum allowed insert size. pairs aligning on the same chr at a distance higher than this are considered candidates for SV.")
						("orientation",  po::value<string>(),      "expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default outtie")
						
						("minimum-supporting-pairs",  po::value<unsigned int>(), "Minimum number of supporting pairs in order to call a variation event (default 10)")
						("minimum-mapping-quality",  po::value<int>(), "Minimum mapping quality to consider an alignment (default 20)")
						("window-size",  po::value<unsigned int>(),    "Size of the sliding window (default 1000)")
						("window-step",  po::value<unsigned int>(),    "size of the step in overlapping window (must be lower than window-size) (default 100)")
						("coverage",     po::value<float>(), "do not compute coverage from bam file, use the one specified here (must be used in combination with --insert --insert-std)")
						("insert",       po::value<float>(), "do not compute insert size from bam file, use the one specified here (must be used in combination with --coverage --insert-std)")
						("insert-std",   po::value<float>(), "do not compute insert size standard deviation from bam file, use the one specified here (must be used in combination with --insert --coverage)");


	//extraction parameters
	stringstream ssExtract;
	ssExtract << "extract module" << endl << endl << "Allowed options";
	po::options_description desc2(ssExtract.str().c_str());
	desc2.add_options()	("extract", "select the extract module to extract regions of interest as bam file")
				("input-file",       po::value<string>(), "tab file containing the events")
				("filter",       po::value<string>(), "filter the reads that are to be extracted");

	//find copy number variation parameters
	stringstream ssCNV;
	ssCNV << "cnv module" << endl << endl << "Allowed options";
	po::options_description desc3(ssCNV.str().c_str());
	desc3.add_options()	("cnv", "Select the cnv module to find copy number variations ")
				("TODO",       po::value<string>(), "analyse copy number variations"" ");




	stringstream ss;
	ss << package_description() << endl;
	po::options_description menu_union(ss.str().c_str());
	menu_union.add(general).add(desc).add(desc2).add(desc3);



//TODO: add minimum number of reads to support a translocation, window size etc.

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, menu_union), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		return 2;
	}

	if(vm.count("help")) {
		DEFAULT_CHANNEL << menu_union << endl;
		return 2;
	}else if(vm.count("module-help")){
		string moduleHelp = vm["module-help"].as<string>();
		if(moduleHelp=="general"){
			DEFAULT_CHANNEL << general << endl;
		}else if(moduleHelp=="sv"){
			DEFAULT_CHANNEL << desc << endl;
		}else if(moduleHelp=="extract"){
			DEFAULT_CHANNEL << desc2 << endl;
		}else if(moduleHelp=="cnv"){
			DEFAULT_CHANNEL << desc3 << endl;
		}else{
			cout << "Error: the module was not found" << endl;
		}
		return 2;
	}

	// PARSE BAM file
	if (!vm.count("bam")) {
		DEFAULT_CHANNEL << "Please specify --bam " << endl;
		return 2;
	}
	alignmentFile = vm["bam"].as<string>();
	if (!vm.count("bai")){
		DEFAULT_CHANNEL << "a reference file must be selected" << endl;
		return 2;
	}
	string indexFile=vm["bai"].as<string>();

	if (vm.count("output")) {
		string header = vm["output"].as<string>();
		outputFileHeader = header ;
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


	//if the bam extraction module is chosen;
	if (vm.count("extract")){
		
		po::options_description menu_union(ss.str().c_str());
		po::variables_map vm;
		menu_union.add(general).add(desc2);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
			return 2;
		}

	if (vm.count("input-file")){
		roi=vm["input-file"].as<string>();
		TranslocationBAM(alignmentFile,outputFileHeader,roi,contig2position,indexFile);
		return 1;
	}else{
		cout << "please select an input file" << endl;
		return 0;
	}	

	//if the find structural variations module is chosen
	}else if(vm.count("sv")){

		po::options_description menu_union(ss.str().c_str());
		po::variables_map vm;
		menu_union.add(general).add(desc);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
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
				windowSize, windowStep, minimumSupportingPairs, coverage, meanInsert, insertStd, outputFileHeader, indexFile);


	//if the find copy number variation module is chosen
	}else if(vm.count("cnv")){

		po::options_description menu_union(ss.str().c_str());
		po::variables_map vm;
		menu_union.add(general).add(desc3);

		try {
			po::store(po::parse_command_line(argc, argv, menu_union), vm);
			po::notify(vm);
		} catch (boost::program_options::error & error) {
			ERROR_CHANNEL <<  "Only one module may be selected" << endl;
			return 2;
		}
		cout << "cnv" << endl;
	}

}


void findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality,
		uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs
, float meanCoverage, float meanInsertSize, float StdInsertSize, string outputFileHeader, string indexFile) {
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//open a reade used to jump across the bamfile



	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	// now create Translocation on the fly
	Window *window;

	window = new Window(windowSize, windowStep, max_insert, minimum_mapping_quality,
		outtie,  meanInsertSize,  StdInsertSize,  minimumSupportingPairs,
		 meanCoverage,  outputFileHeader,bamFileName,indexFile);
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




