/*
   Francesco Vezzi
   Jesper Eisfeldt
 */

#include "ProgramModules.h"
#include "data_structures/Translocation.h"

//function used to find translocations
StructuralVariations::StructuralVariations() { }

bool zeroth_columnsort(const vector<int>& lhs,const vector<int>& rhs){
	return(lhs[1] < rhs[1]);
}

void StructuralVariations::findTranslocationsOnTheFly(string bamFileName, bool outtie, float meanCoverage, string outputFileHeader, string version, map<string,int> SV_options) {
	size_t start = time(NULL);
	//open the bam file
	BamReader bamFile;
	bamFile.Open(bamFileName);
	//Information from the header is needed to initialize the data structure
	SamHeader head = bamFile.GetHeader();
	// now create Translocation on the fly
	Window *window;

	window = new Window(bamFileName,outtie,meanCoverage,outputFileHeader,SV_options);
	window-> version = version;
	window->initTrans(head);
	//expands a vector so that it is large enough to hold reads from each contig in separate elements
	window->eventReads.resize(4);
	window->eventReads[0].resize(SV_options["contigsNumber"]);
	window->eventReads[1].resize(SV_options["contigsNumber"]);
	window->eventReads[2].resize(SV_options["contigsNumber"]);
	window->eventReads[3].resize(SV_options["contigsNumber"]);
	
	window->eventSplitReads.resize(4);
	window->eventSplitReads[0].resize(SV_options["contigsNumber"]);
	window->eventSplitReads[1].resize(SV_options["contigsNumber"]);
	window->eventSplitReads[2].resize(SV_options["contigsNumber"]);
	window->eventSplitReads[3].resize(SV_options["contigsNumber"]);

	window-> binnedCoverage.resize(SV_options["contigsNumber"]);
	window-> binnedQuality.resize(SV_options["contigsNumber"]);
	
	window-> linksFromWin.resize(4);
	window-> linksFromWin[0].resize(SV_options["contigsNumber"]);
	window-> linksFromWin[1].resize(SV_options["contigsNumber"]);
	window-> linksFromWin[2].resize(SV_options["contigsNumber"]);
	window-> linksFromWin[3].resize(SV_options["contigsNumber"]);

	for (int i=0;i< SV_options["contigsNumber"];i++){
		window -> SV_calls[i] = vector<string>();
		window -> SV_positions[i]= vector< vector< int> >();
	}
	
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
			window -> binnedQuality[window -> contig2position[splitline[0]]].push_back(atof(splitline[4].c_str()));
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
	
	//call any remaining variant
	for(int i=0;i< SV_options["contigsNumber"];i++){
		for (int j=0;j<4;j++){
			if(window -> eventReads[j][i].size() >= window -> minimumPairs){
				window -> pairOrientation = j;
				window->computeVariations(i);
			}
			window->eventReads[j][i]=queue<BamAlignment>();
			window->eventSplitReads[j][i] = vector<BamAlignment>();
			window->linksFromWin[j][i]=queue<int>();
		}
		
		
	}
	//print calls
	for(int i=0;i< SV_options["contigsNumber"];i++){
		if(window->SV_calls[i].size() == 0){
			continue;
		}
		vector<  vector <int> > sorted_SV_positions;
		sorted_SV_positions = window -> SV_positions[i];
		sort(sorted_SV_positions.begin(), sorted_SV_positions.end(),zeroth_columnsort );
		for (int j=0;j < sorted_SV_positions.size();j++){
			window-> TIDDITVCF <<  window -> SV_calls[i][ sorted_SV_positions[j][0] ];
		}
	}
	window->TIDDITVCF.close();
	printf ("variant calling time consumption= %lds\n", time(NULL) - start);
}

