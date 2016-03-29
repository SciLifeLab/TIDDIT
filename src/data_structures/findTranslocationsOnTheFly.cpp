/*
   Francesco Vezzi
   Jesper Eisfeldt
 */
#include "ProgramModules.h"
#include "data_structures/Translocation.h"
#include <boost/thread.hpp>
#include <boost/algorithm/string.hpp>

//function used to find translocations
StructuralVariations::StructuralVariations() { }

void StructuralVariations::findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality, uint32_t minimumSupportingPairs
, float meanCoverage, float meanInsertSize, float StdInsertSize, string outputFileHeader, string indexFile,int contigsNumber,int ploidy, int readLength) {
	size_t start = time(NULL);
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	// now create Translocation on the fly
	Window *window;

	window = new Window(max_insert, min_insert,minimum_mapping_quality,
		outtie,  meanInsertSize,  StdInsertSize,  minimumSupportingPairs,
		 meanCoverage,  outputFileHeader,bamFileName,indexFile,ploidy,readLength);
	window->initTrans(head);
	//expands a vector so that it is large enough to hold reads from each contig in separate elements
	window->eventReads.resize(contigsNumber);
	window->eventSplitReads.resize(contigsNumber);

	window-> binnedCoverage.resize(contigsNumber);
	window-> linksFromWin.resize(contigsNumber);
	
	window -> numberOfEvents = 0;

	string line;
	string coverageFile=outputFileHeader+".tab";
	ifstream inputFile( coverageFile.c_str() );
	while (getline( inputFile, line )){
		vector<std::string> splitline;
		boost::split(splitline, line, boost::is_any_of("\t"));
		window -> binnedCoverage[window -> contig2position[splitline[0]]].push_back(atof(splitline[1].c_str()));
	}
	inputFile.close();


	//Initialize bam entity
	BamAlignment currentRead;
	//now start to iterate over the bam file
	int counter = 0;
	while ( bamFile.GetNextAlignmentCore(currentRead) ) {
	  if(currentRead.IsMapped()) {
	    window->insertRead(currentRead);
	  }
	}
	for(int i=0;i< window-> eventReads.size();i++){
	  if(window -> eventReads[i].size() >= window -> minimumPairs){
	    window->computeVariations(i);
	  }
	  window->eventReads[i]=queue<BamAlignment>();
	  window->eventSplitReads[i] = vector<BamAlignment>();
	}
	  
	window->interChrVariationsVCF.close();
	window->intraChrVariationsVCF.close();
	printf ("variant calling time consumption= %lds\n", time(NULL) - start);
}

