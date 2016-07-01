/*
   Francesco Vezzi
   Jesper Eisfeldt
 */

#include "ProgramModules.h"
#include "data_structures/Translocation.h"

//function used to find translocations
StructuralVariations::StructuralVariations() { }

void StructuralVariations::findTranslocationsOnTheFly(string bamFileName, bool outtie, float meanCoverage, string outputFileHeader, map<string,int> SV_options) {
	size_t start = time(NULL);
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	// now create Translocation on the fly
	Window *window;

	window = new Window(bamFileName,outtie,meanCoverage,outputFileHeader,SV_options);
	window->initTrans(head);
	//expands a vector so that it is large enough to hold reads from each contig in separate elements
	window->eventReads.resize(4);
	window->eventReads[0].resize(SV_options["contigsNumber"]);
	window->eventReads[1].resize(SV_options["contigsNumber"]);
	window->eventReads[2].resize(SV_options["contigsNumber"]);
	window->eventReads[3].resize(SV_options["contigsNumber"]);
	
	window->eventSplitReads.resize(SV_options["contigsNumber"]);

	window-> binnedCoverage.resize(SV_options["contigsNumber"]);
	
	window-> linksFromWin.resize(4);
	window-> linksFromWin[0].resize(SV_options["contigsNumber"]);
	window-> linksFromWin[1].resize(SV_options["contigsNumber"]);
	window-> linksFromWin[2].resize(SV_options["contigsNumber"]);
	window-> linksFromWin[3].resize(SV_options["contigsNumber"]);
	
	window -> numberOfEvents = 0;

	string line;
	string coverageFile=outputFileHeader+".tab";
	ifstream inputFile( coverageFile.c_str() );
	int line_number=0;
	while (std::getline( inputFile, line )){
		if(line_number > 0){
			vector<string> splitline;
    		std::stringstream ss(line);
    		std::string item;
    		while (std::getline(ss, item, '\t')) {
        		splitline.push_back(item);
    		}
			window -> binnedCoverage[window -> contig2position[splitline[0]]].push_back(atof(splitline[3].c_str()));
		}
		line_number += 1;
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
	
	
	for(int i=0;i< SV_options["contigsNumber"];i++){
		for (int j=0;j<4;j++){
			if(window -> eventReads[j][i].size() >= window -> minimumPairs){
				window -> O = j;
				window->computeVariations(i);
			}
			window->eventReads[j][i]=queue<BamAlignment>();
		}
		window->eventSplitReads[i] = vector<BamAlignment>();
		
	}
	
	window->interChrVariationsVCF.close();
	window->intraChrVariationsVCF.close();
	printf ("variant calling time consumption= %lds\n", time(NULL) - start);
}

