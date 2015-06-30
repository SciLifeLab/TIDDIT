/*
   Francesco Vezzi
   Jesper Eisfeldt
 */




//function used to find translocations
void findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality,
		uint32_t windowSize , uint32_t windowStep, uint32_t minimumSupportingPairs
, float meanCoverage, float meanInsertSize, float StdInsertSize, string outputFileHeader, string indexFile) {
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
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

