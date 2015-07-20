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
void Extract::extract(string BamFileName,string outputFileHeader,string inputFileName,map<string,unsigned int> contig2position, string indexFile){
	int VCFheader=0;
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
		std::cout << "error creating output folder(is a folder with the output name already existing?)" << "\n";
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
			//splits on tab
			vector<std::string> splitline;
			boost::split(splitline, line, boost::is_any_of("\t"));


			//make sure that the reader is not collecting values from the metadata of the vcf file
			if(VCFheader != 0){
				queue<int> chr;
				queue<int> startPos;queue<int> endPos;
				//the filename of the file containing single regions
				string regionfile = "";


				string INFOLine;
				//extract the info line
				INFOLine=splitline[7];
			

				//get the region from the info line
				vector<std::string> splitINFOLine;
				boost::split(splitINFOLine, INFOLine, boost::is_any_of(";"));
				//
				for(int i =0;i<2;i++){
					vector<std::string> chromosomeVector;
					boost::split(chromosomeVector,splitINFOLine[1+2*i], boost::is_any_of("="));
					//extract the chromosome
					chr.push(contig2position[chromosomeVector[1]]);
					regionfile+=chromosomeVector[1];
					regionfile+="_";
			
					//extracts the start position on the given chromosome
					vector<std::string> positionVector;
					boost::split(positionVector,splitINFOLine[2+2*i], boost::is_any_of("="));
					boost::split(positionVector,positionVector[1], boost::is_any_of(","));
					startPos.push( atoi(positionVector[0].c_str()));
					regionfile+=positionVector[0];
					regionfile+="_";

						//extracts the end position on the given chromosome
					endPos.push(atoi(positionVector[1].c_str()));
					regionfile+=positionVector[1];
					if(i ==0){
						regionfile+="_";
					}
				}


				regionfile+=".bam";
				regionfile=eventfolder+"/"+regionfile;
					// attempt to open our BamWriter{
				if ( !regionwriter.Open(regionfile.c_str(), header, references) ) {
					cout << "Could not open output event BAM file" << endl;
					return;
				}
			





				//collects all reads found in the given intervals and writes them to the file outputFileHeader+"_roi.bam"
				for(int i=0;i<2;i++){
					bamFile.SetRegion(chr.front(),startPos.front(),chr.front(),endPos.front()); 
					chr.pop();startPos.pop();endPos.pop();
					while ( bamFile.GetNextAlignment(currentRead) ) {
						if(currentRead.IsMapped()) {
							if(currentRead.Position >= startPos.front() and currentRead.Position <= endPos.front()){
							//prints the read to bam file
							writer.SaveAlignment(currentRead);
							regionwriter.SaveAlignment(currentRead);
							}
						}
					}
				}
				regionwriter.Close();
				

			}else{
				if(splitline[0] == "#CHROM"){
					VCFheader =1;
				}
			
			}
		}
		
		inputFile.close();
	}

	cout << "writing complete" << endl;
	cout << outputFileHeader << endl;
	writer.Close();
	bamFile.Close();
}


