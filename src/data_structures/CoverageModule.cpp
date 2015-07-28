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

//constructor
Cov::Cov() { }
//The main function
void Cov::coverageMain(string bamFile,string baiFile,string inputFileName,string output, double averageCoverage,map<string,unsigned int> contig2position,int selection,int binSize){
	//the vcf mode outputs two files, one tab/bed file like the two other modes, and one vcf file, similar to the inout vcf but with added coverage data.
	ofstream coverageOutput;
	coverageOutput.open((output+".tab").c_str());

	ofstream coverageOutputVcf;
	coverageOutputVcf.open((output+".vcf").c_str());
	//if vcf was chosen
	if(selection == 0){
		cout << "analysing the coverage between the findtranslocation events" << endl;
		coverageOutput << "CHR" << "\t" << "endA" << "\t" << "startB" << "\t" << "coverage" <<"\t" << "libraryCoverage"<< "\t" << "coverage/librarycoverage" << endl;
		//read the input file containing regions of interest and return them as a vector
		VcfBedInput *openInputFile;
		openInputFile = new VcfBedInput();
		queue<string> regions= openInputFile -> findRegions(inputFileName, 1);
		ifstream inputFile( inputFileName.c_str() );

		int VCFheader =0;
		if (inputFile){
			string line;
			while (getline( inputFile, line )){
				//splits on tab
				vector<std::string> splitline;
				boost::split(splitline, line, boost::is_any_of("\t"));
				if(VCFheader != 0){

					//remove the first position
					//read through each region and calculate its coverage

					regions.pop();regions.pop();
					int chr;int startPos;int endPos;
					string chromosome;string startPosition;string endPosition;
					//extracts the end position of the curent region
					endPosition=regions.front();
					endPos=atoi(regions.front().c_str());
					regions.pop();
					//extract the chromosome
					chr=contig2position[regions.front()];
					chromosome=regions.front();
					regions.pop();
					//remove the start position of the curent region, we are only interested to calculate the coverage between region an and n+1
					startPos=atoi(regions.front().c_str());
					startPosition = regions.front();
					regions.pop();
					//remove the end position of the last region
					regions.pop();
					double coverage=findCoverage(bamFile,baiFile,chr,endPos,startPos);
					//print the results
					//the tab file
					coverageOutput << chromosome << "\t" << endPosition << "\t" << startPosition << "\t" << coverage
					 << "\t" << averageCoverage << "\t" << (float)coverage/(float)averageCoverage << endl;
					
					//the vcf file
					string INFOLine=splitline[7];
					vector<std::string> splitINFOLine;
					boost::split(splitINFOLine, INFOLine, boost::is_any_of(";"));

					coverageOutputVcf << splitline[0] << "\t" << splitline[0] << "\t" << splitline[1] << "\t" << splitline[2] << "\t"
							  << splitline[3] << "\t" << splitline[4] << "\t" << splitline[5] << "\t" << splitline[6] << "\t";

					//The coverage data is always added behind COVB, no matter what statistics have been added to the file
					for(int i =0;i< splitINFOLine.size();i++){
						string INFO=splitINFOLine[i];
						vector<std::string> splitINFO;
						boost::split(splitINFO,INFO, boost::is_any_of("="));
						if(splitINFO[0] == "COVB"){
							coverageOutputVcf << INFO << ";COVAB=" << coverage << ";LIBCOV=" << averageCoverage << ";COVRATIO="<< (float)coverage/(float)averageCoverage;

						}else{
							coverageOutputVcf << INFO;
						}

						//at the last INFO column, an end line is added as separator, for the other a ";" is added.
						if(i == splitINFOLine.size() -1){
							coverageOutputVcf << "\n";
						}else{
							coverageOutputVcf << ";";
						}
					}

				}else{
					if(splitline[0] == "#CHROM"){
						VCFheader =1;
						coverageOutputVcf << line << endl;

					}else if(splitline[0] == "##INFO=<ID=COVB,Number=1,Type=Integer,Description=\"Coverage on window B\">"){
						coverageOutputVcf << line << endl;
						coverageOutputVcf << "##INFO=<ID=COVAB,Number=1,Type=Integer,Description=\"The coverage between window A and B\">\n";
						coverageOutputVcf << "##INFO=<ID=LIBCOV,Number=1,Type=Integer,Description=\"The average library coverage\">\n";
						coverageOutputVcf << "##INFO=<ID=COVRATIO,Number=1,Type=Integer,Description=\"COVAB/LIBCOV\">\n";

					}else{
						coverageOutputVcf << line << endl;
					}
				}
			}
		}
	//if bin was chosen
	}else if(selection == 1){
		cout << "analysing the coverage in bins of size " << binSize << endl;
		coverageOutput << "CHR" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" <<"\t" << "libraryCoverage"<< "\t" << "coverage/librarycoverage" << endl;
		int binStart =0;
		int binEnd=binSize+binStart;
		int currentChr=-1;
		vector<int> sequencedBases;
		sequencedBases.push_back(0);

		BamReader alignmentFile;
		//open the bam file
		alignmentFile.Open(bamFile);
		BamAlignment currentRead;

		//get which refID belongs to which chromosome
		map<unsigned int,string> position2contig;
		uint32_t contigsNumber = 0;
		SamSequenceDictionary sequences  = alignmentFile.GetHeader().Sequences;
		for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
			position2contig[contigsNumber]  = sequence->Name;
			contigsNumber++;
		}


		
		while ( alignmentFile.GetNextAlignment(currentRead) ) {
			if(currentRead.IsMapped()) {
				//initialise the chromosome ID
				if(currentChr == -1){
					currentChr=currentRead.RefID;
				}

				//check if we have left the chromosome
				if(currentChr == currentRead.RefID){
					//if the entire read is inside the region, add all the bases to sequenced bases
         	           		if(currentRead.Position >= binStart and currentRead.Position+currentRead.Length-1 <= binEnd){
						sequencedBases[0]+=currentRead.Length;

					}		
					//if the read starts within the region but reaches outside it, add only those bases that fit inside the region.
         	           		if(currentRead.Position >= binStart and currentRead.Position+currentRead.Length-1 > binEnd and currentRead.Position <= binEnd){
	
						sequencedBases[0]+=binEnd-currentRead.Position+1;
						//the part of the read hanging out of the bin is added to the bins following the currentbin
						int remainingRead=currentRead.Length-(binEnd-currentRead.Position+1);
						int binNumber=1;

						while (remainingRead > 0){
							//if there are not enough bins to store the bases, add another bin
							if(binNumber > sequencedBases.size()){
								sequencedBases.push_back(0);
							}
							//check if the entire read fits in the next bin or not
							if(binSize >= remainingRead){	
								sequencedBases[binNumber]+=remainingRead;
								remainingRead=0;		
							//if the amount of bases exceeds an entire bin, take a bin of bases and add it to the next bin, 
							//itterate once more and see how many bases that will fit into the next				
							}else{
								sequencedBases[binNumber]+=binSize;
								remainingRead=remainingRead-binSize;
							}
							binNumber++;
						}
					}
				}
					
				//when there are no more reads in the current bin, calculate the coverage and start analysing the next bin
				if(binEnd < currentRead.Position or currentChr != currentRead.RefID){
					double coverage=(double)sequencedBases[0]/(double)binSize;
					unsigned int position=currentRead.RefID;
					cout << position2contig[position] << "\t" << binStart << "\t" << binEnd << "\t" << coverage << endl;
					coverageOutput << position2contig[position] << "\t" << binStart << "\t" << binEnd << "\t" << coverage << "\t" << averageCoverage << "\t" << (float)coverage/(float)averageCoverage << endl;
					binStart=binEnd;
					binEnd=binStart+binSize;
					sequencedBases.erase(sequencedBases.begin());
					if(sequencedBases.size() == 0){
						sequencedBases.push_back(0);
					}
					if(currentChr != currentRead.RefID){
						currentChr = currentRead.RefID;
						binStart=0;
						binEnd=binSize;
					}

					
				}	
		         
			}
		}
		







	//if bed was chosen
	}else{
		cout << "analysing the coverage of the bed file regions" << endl;
		coverageOutput << "CHR" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" <<"\t" << "libraryCoverage"<< "\t" << "coverage/librarycoverage" << endl;
		ofstream coverageOutput;
		coverageOutput.open(output.c_str());

		//read the input file containing regions of interest and return them as a vector
		VcfBedInput *openInputFile;
		openInputFile = new VcfBedInput();
		queue<string> regions= openInputFile -> findRegions(inputFileName, 2);
		//read through each region and calculate its coverage
		int numberOfregions=regions.size()/3;
		for(int i =0; i < numberOfregions;i++){
			int chr;int startPos;int endPos;
			string chromosome;string startPosition;string endPosition;

			chr=contig2position[regions.front()];
			chromosome=regions.front();
			regions.pop();

			startPos=atoi(regions.front().c_str());
			startPosition = regions.front();
			regions.pop();

			endPosition=regions.front();
			endPos=atoi(regions.front().c_str());
			regions.pop();
			double coverage=findCoverage(bamFile,baiFile,chr,startPos,endPos);
			coverageOutput << chromosome << "\t" << startPosition << "\t" << endPosition << "\t" << coverage << "\t" << averageCoverage << "\t" << (float)coverage/(float)averageCoverage << endl;
		}
	}

		
}

//This function is used to calculate the coverage of a given region
double Cov::findCoverage(string bamFileName, string baiFile,int chr, int start, int end){
	double coverage = 0;
	long sequencedBases;
	double regionLength = end-start+1;

	BamReader bamFile;
	//open the bam file
	bamFile.Open(bamFileName);
	BamAlignment currentRead;


	//test if an index file is available
	if(bamFile.OpenIndex(baiFile) == 0){
		cout << "Failed to load the index file" << endl;
	}

	bamFile.SetRegion(chr,start,chr,end); 
		while ( bamFile.GetNextAlignment(currentRead) ) {
			if(currentRead.IsMapped()) {

				//if the entire read is inside the region, add all the bases to sequenced bases
         	           	if(currentRead.Position >= start and currentRead.Position+currentRead.Length-1 <= end){
					sequencedBases+=currentRead.Length;

				}	
				//if the read start outside the region, but reaches inside it, add those bases that enter the region to the sequenced bases
         	           	if(currentRead.Position < start and currentRead.Position+currentRead.Length-1 <= end and currentRead.Position+currentRead.Length-1 >= start){
					sequencedBases+=currentRead.Length-(start-currentRead.Position)+1;

				}	
				//if the read starts within the region but reaches outside it, add only those bases that fit inside the region.
         	           	if(currentRead.Position >= start and currentRead.Position+currentRead.Length-1 > end and currentRead.Position <= end){
					sequencedBases+=end-currentRead.Position+1;

				}	
		         
			}
		}


	coverage=sequencedBases/regionLength;
	return(coverage);
}
