/*
 * ExtractModule.cpp
 *
 *  Created on: Jul 1, 2015
 *      Author: vezzi, Eisfeldt
 */

#include "ProgramModules.h"


#include <sstream>
#include <iostream>
#include <fstream>

#include  <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "api/BamWriter.h"
#include <queue>


Extract::Extract() { }

//This function is used to extract the event regions from a bam file using the output vcf of FindTranslocations
void Extract::extract(string BamFileName,string outputFileHeader,string inputFileName,map<string,unsigned int> contig2position, string indexFile, int BedOrVCF,int exclusion){
	int VCFheader=0;
	string line;
	string eventbam;
	string eventfolder;

	if(exclusion == 3){
		cout << "All output excluded!" << endl;
		return;
	}

	//creates a folder for the outputdata
	boost::filesystem::path dir(outputFileHeader.c_str());
	if(!boost::filesystem::create_directory(dir)) {
		std::cout << "error creating output folder" << "\n";
		return;
	}
	
	//only generate this folder if the user wants files of each separate event
	if(exclusion != 1){
		eventfolder=outputFileHeader+"/events";
		boost::filesystem::path eventdir( eventfolder.c_str());
		if(!boost::filesystem::create_directory(eventdir)) {
			std::cout << "error creating output folder(is a folder with the output name already existing?)" << "\n";
			return;
		}
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

	//open the bam file writer if the user wants to output a file containing all the events
	if( exclusion != 2){
		// attempt to open our BamWriter{
		cout << eventbam << endl;
		if ( !writer.Open(eventbam.c_str(), header, references) ) {
			cout << "Could not open output BAM file" << endl;
			return;
		}
	}

	//test if an index file is available
	if(bamFile.OpenIndex(indexFile) == 0){
		cout << "Failed to load the index file" << endl;
	}
	

	//opens the input file and reads each line
	cout << "initiates extraction, please wait" << endl;

	VcfBedInput *openInputFile;
	openInputFile = new VcfBedInput();
	queue<string> regions= openInputFile -> findRegions(inputFileName, BedOrVCF);
	int numberOfRegions;
	int regionsPerLine;
	//if vcf is chosen
	if(BedOrVCF ==1){
		regionsPerLine=2;
		numberOfRegions=regions.size()/6;
	}else{
		regionsPerLine=1;
		numberOfRegions=regions.size()/3;
	}
	for(int j =0;j < numberOfRegions;j++){
		queue<int> chr;
		queue<int> startPos;queue<int> endPos;
		//the filename of the file containing single regions
		string regionfile = "";
		for(int i =0;i<regionsPerLine;i++){	
			//extract the chromosome
			chr.push(contig2position[regions.front()]);
			regionfile+=regions.front();
			regionfile+="_";
			regions.pop();

			//extracts the start position on the given chromosome
			startPos.push( atoi(regions.front().c_str()));
			regionfile+=regions.front();
			regionfile+="_";
			regions.pop();

			//extracts the end position on the given chromosome
			endPos.push(atoi(regions.front().c_str()));
			regionfile+=regions.front();
			regions.pop();

			if(i ==0){
				regionfile+="_";
			}
		}

		//only created if the user wants to output every event as separate files
		if(exclusion != 1){
			regionfile+=".bam";
			regionfile=eventfolder+"/"+regionfile;
			// attempt to open our BamWriter{
			if ( !regionwriter.Open(regionfile.c_str(), header, references) ) {
				cout << "Could not open output event BAM file" << endl;
				return;
			}
		}		





		//collects all reads found in the given intervals and writes them to the file outputFileHeader+"_roi.bam"
		for(int i=0;i<regionsPerLine;i++){
			bamFile.SetRegion(chr.front(),startPos.front(),chr.front(),endPos.front()); 
			while ( bamFile.GetNextAlignment(currentRead) ) {
				if(currentRead.IsMapped()) {
         	           		if(currentRead.Position >= startPos.front() and currentRead.Position <= endPos.front()){
						//prints the read to bam file

						//only created if the user wants to output every event as separate files
						if(exclusion != 2){	
							writer.SaveAlignment(currentRead);
						}
						//only saved if the user wants to output one file containing all the selected reads
						if(exclusion !=1){
							regionwriter.SaveAlignment(currentRead);
						}
					}			         
				}
			}
			chr.pop();startPos.pop();endPos.pop();                    
		}
		regionwriter.Close();	
	}

	cout << "writing complete" << endl;
	cout << outputFileHeader << endl;
	writer.Close();
	bamFile.Close();
}


