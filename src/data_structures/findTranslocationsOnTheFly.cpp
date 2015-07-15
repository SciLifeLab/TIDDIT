/*
   Francesco Vezzi
   Jesper Eisfeldt
 */
#include "ProgramModules.h"
#include "data_structures/Translocation.h"


//function used to find translocations
StructuralVariations::StructuralVariations() { }

void StructuralVariations::findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality, uint32_t minimumSupportingPairs
, float meanCoverage, float meanInsertSize, float StdInsertSize, string outputFileHeader, string indexFile,int contigsNumber) {
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	// now create Translocation on the fly
	Window *window;

	window = new Window(max_insert, minimum_mapping_quality,
		outtie,  meanInsertSize,  StdInsertSize,  minimumSupportingPairs,
		 meanCoverage,  outputFileHeader,bamFileName,indexFile);
	window->initTrans(head);
	//expands a vector so that it is large enough to hold reads from each contig in separate elements
	window->eventReads.resize(contigsNumber);
	window->covOnChrA.resize(contigsNumber);
	window->tmpCovOnChrA.resize(contigsNumber);
	window->linksFromWin.resize(contigsNumber);
	window->tmpLinksFromWin.resize(contigsNumber);


	//Initialize bam entity
	BamAlignment currentRead;
	//now start to iterate over the bam file
	int counter = 0;
	while ( bamFile.GetNextAlignment(currentRead) ) {
		if(currentRead.IsMapped()) {
			window->insertRead(currentRead);
		}
	}
		for(int i=0;i< window-> eventReads.size();i++){
			if(window -> eventReads[i].size() >= window -> minimumPairs){
			window->computeVariations(i);
			}
			window->eventReads[i]=queue<BamAlignment>();
		}
	window->interChrVariations.close();
	window->intraChrVariations.close();
	window->interChrVariationsVCF.close();
	window->intraChrVariationsVCF.close();

}

