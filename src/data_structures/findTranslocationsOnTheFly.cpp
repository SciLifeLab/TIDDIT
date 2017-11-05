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

void StructuralVariations::findTranslocationsOnTheFly(string bamFileName, bool outtie, float meanCoverage, string outputFileHeader, string version,string command, map<string,int> SV_options) {
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
	stringstream ss;
	string orientation_string="innie";
	if(outtie == true){
		orientation_string="outtie";
	}
	ss << "##LibraryStats=TIDDIT-" << version <<   " Coverage=" << meanCoverage << " ReadLength=" << SV_options["readLength"] << " MeanInsertSize=" << SV_options["meanInsert"] << " STDInsertSize=" << SV_options["STDInsert"] << " Orientation=" << orientation_string << "\n" << "##TIDDITcmd=\"" << command << "\""; 
	string libraryData=ss.str();
	window->initTrans(head,libraryData);

	for (int i=0;i< SV_options["contigsNumber"];i++){
		window -> SV_calls[i] = vector<string>();
		window -> SV_positions[i]= vector< vector< int> >();
	}
	
	window -> numberOfEvents = 0;
	//Initialize bam entity
	BamAlignment currentRead;
	//now start to iterate over the bam file
	int counter = 0;
	while ( bamFile.GetNextAlignmentCore(currentRead) ) {
	  if(currentRead.IsMapped()) {
	    window->insertRead(currentRead);
	  }
	}
	
	//print calls
	for(int i=0;i< SV_options["contigsNumber"];i++){
		for (int j=0;j < window -> SV_calls[i].size();j++){
			window-> TIDDITVCF <<  window -> SV_calls[i][j];
		}
	}
	window->TIDDITVCF.close();
	printf ("variant calling time consumption= %lds\n", time(NULL) - start);
}

