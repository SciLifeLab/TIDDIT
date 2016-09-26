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
#include <queue>

ostream& test(ofstream &coverageOutput,string output){
	if(output != "stdout"){
		ostream& covout=coverageOutput;
		return(covout);
	}else{
		ostream& covout=cout;
		return(covout);
	}

}
//constructor
Cov::Cov(int binSize,string bamFile,string output){
	ostream& covout=test(coverageOutput,output);
	if(output != "stdout"){
		coverageOutput.open((output+".tab").c_str());
		static ostream& covout = coverageOutput;
		
	}else{
		static ostream& covout=cout;
	}

	covout << "#CHR" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" <<"\t" << "quality" << endl;

	//initialize the function
	this -> binStart =0;
	this -> binEnd=binSize+binStart;
	this -> binSize = binSize;
	this -> currentChr=-1;
	
    this -> contigsNumber = 0;
    this -> bamFile = bamFile;

    BamReader alignmentFile;
	//open the bam file
	alignmentFile.Open(bamFile);
	
	//get which refID belongs to which chromosome
	SamSequenceDictionary sequences  = alignmentFile.GetHeader().Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		position2contig[contigsNumber]  = sequence->Name;
		contig2position[sequence->Name] = contigsNumber;
		contigLength.push_back(StringToNumber(sequence->Length));
		contigsNumber++;
	}
	this -> coverageStructure.resize(contigsNumber);
	this -> qualityStructure.resize(2);
	qualityStructure[0].resize(contigsNumber);
	qualityStructure[1].resize(contigsNumber);
	for(int i=0;i<contigsNumber;i++){
		//cout << i << endl;
		//cout << ceil(contigLength[i]/double(binSize)) << endl;
		coverageStructure[i].resize(ceil(contigLength[i]/double(binSize)),0);
		qualityStructure[0][i].resize(ceil(contigLength[i]/double(binSize)),0);
		qualityStructure[1][i].resize(ceil(contigLength[i]/double(binSize)),0);
	}
	
	
	
    alignmentFile.Close();
}



//this function calculates the coverage in bins of size binSize across the entire bamfile
void Cov::bin(BamAlignment currentRead){
	//initialise the chromosome ID
	if(currentChr == -1){
		currentChr=currentRead.RefID;
	}
	
	//if the quality of the read is high enough, it will be added to the data structure
	readStatus alignmentStatus = computeReadType(currentRead, 100000,100, true);
	if(alignmentStatus != unmapped and alignmentStatus != lowQualty) {
		int element=floor(double(currentRead.Position)/double(binSize));

		//if the entire read is inside the region, add all the bases to sequenced bases
		if(currentRead.Position >= element*binSize and currentRead.Position+currentRead.Length-1 <= (element+1)*binSize){
			coverageStructure[currentRead.RefID][element]+=currentRead.Length;
			qualityStructure[0][currentRead.RefID][element] += currentRead.MapQuality;
			qualityStructure[1][currentRead.RefID][element] += 1;
		}else{
		//if the read starts within the region but reaches outside it, add only those bases that fit inside the region.
			coverageStructure[currentRead.RefID][element]+=(element+1)*binSize-currentRead.Position+1;
			qualityStructure[0][currentRead.RefID][element] += currentRead.MapQuality;
			qualityStructure[1][currentRead.RefID][element] += 1;
		
			//the part of the read hanging out of the bin is added to the bins following the currentbin
			int remainingRead=currentRead.Length-((element+1)*binSize-currentRead.Position+1);
			while (remainingRead >= binSize){
				coverageStructure[currentRead.RefID][element]+=binSize;
				qualityStructure[0][currentRead.RefID][element] += currentRead.MapQuality;
				qualityStructure[1][currentRead.RefID][element] += 1;
				remainingRead=remainingRead-binSize;
			}
			if (remainingRead > 0){
				element++;
				coverageStructure[currentRead.RefID][element]+=remainingRead;
				qualityStructure[0][currentRead.RefID][element] += currentRead.MapQuality;
				qualityStructure[1][currentRead.RefID][element] += 1;
			}

		}
	}
}

//prints the results
void Cov::printCoverage(){
    ostream& covout = test(coverageOutput,output);
    for(int i=0;i<contigsNumber;i++){
        for(int j=0;j<coverageStructure[i].size();j++){
            int binStart = j*binSize;
            int binEnd=binStart+binSize;
            
            if(binEnd > contigLength[i]){
						binEnd=contigLength[i];
			}
            double coverage=double(coverageStructure[i][j])/double(binEnd-binStart);
            double quality = 0;
            if(double(qualityStructure[1][i][j]) > 0){
            	quality=double(qualityStructure[0][i][j])/double(qualityStructure[1][i][j]);
            }
            covout << position2contig[i] << "\t" << binStart << "\t" << binEnd << "\t" << coverage << "\t" << quality << endl;
        }
	}
	
}
