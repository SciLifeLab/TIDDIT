/*
   Francesco Vezzi
   Jesper Eisfeldt
 */

#include "ProgramModules.h"
#include "data_structures/Translocation.h"

//function used to find translocations
StructuralVariations::StructuralVariations() { }

void StructuralVariations::findTranslocationsOnTheFly(string bamFileName, BamReader bamFile, bool outtie, float meanCoverage, string outputFileHeader, string version,string command, map<string,int> SV_options,uint64_t genomeLength) {
	size_t start = time(NULL);
	//Information from the header is needed to initialize the data structure

	SamHeader head = bamFile.GetHeader();

	Window *window;
	window = new Window(bamFileName,outtie,meanCoverage,outputFileHeader,SV_options);

	window-> version = version;
	window->initTrans(head);

	string orientation_string="innie";
	if(outtie == true){
		orientation_string="outtie";
	}

	for (int i=0;i< SV_options["contigsNumber"];i++){window -> SV_calls[i] = vector<string>();}

        Cov *calculateCoverage;
        calculateCoverage = new Cov(50,head,outputFileHeader,SV_options["mapping_quality"],true,false,true);

	uint64_t mappedReadsLength= 0;
	//analyse the already sampled reads

	BamAlignment currentRead;
        BamReader sampleBamFile;
        sampleBamFile.Open(outputFileHeader+".sample.bam");
	sampleBamFile.Rewind();

	while ( sampleBamFile.GetNextAlignmentCore(currentRead) ){
		if(currentRead.IsMapped()){
			readStatus alignmentStatus = computeReadType(currentRead, window->max_insert,window-> min_insert, window-> outtie);

        		if( alignmentStatus == lowQualty) {
				continue;
		        }
		
			currentRead.BuildCharData();
			window->insertRead(currentRead, alignmentStatus);
			calculateCoverage -> bin(currentRead, alignmentStatus);

			mappedReadsLength+= currentRead.Length;
		}

	}
	sampleBamFile.Close();

	//now start to iterate over the bam file
	while ( bamFile.GetNextAlignmentCore(currentRead) ) {
		if(currentRead.IsMapped()){

			readStatus alignmentStatus = computeReadType(currentRead, window->max_insert,window-> min_insert, window-> outtie);

        		if( alignmentStatus == lowQualty) {
				continue;
		        }

			currentRead.BuildCharData();
			window->insertRead(currentRead, alignmentStatus);
			calculateCoverage -> bin(currentRead, alignmentStatus);
			mappedReadsLength+= currentRead.Length;
		}

	}
	//print header
	stringstream ss;
	ss << "##LibraryStats=TIDDIT-" << version <<   " Coverage=" << float(mappedReadsLength)/genomeLength << " ReadLength=" << SV_options["readLength"] << " MeanInsertSize=" << SV_options["meanInsert"] << " STDInsertSize=" << SV_options["STDInsert"] << " Orientation=" << orientation_string << "\n" << "##TIDDITcmd=\"" << command << "\""; 
	string libraryData=ss.str();
	window->printHeader(head,libraryData);

	//print calls
	for(int i=0;i< SV_options["contigsNumber"];i++){
		for (int j=0;j < window -> SV_calls[i].size();j++){
			window-> TIDDITVCF <<  window -> SV_calls[i][j];
		}
	}

	for ( map<string,string>::iterator it =  window -> SV_calls_discordant.begin(); it != window -> SV_calls_discordant.end(); ++it  ){
		window-> TIDDITVCF <<  it-> second;
	} 

	window->TIDDITVCF.close();

	calculateCoverage -> printCoverage();
	printf ("signal extraction time consumption= %lds\n", time(NULL) - start);
	exit(0);

}
