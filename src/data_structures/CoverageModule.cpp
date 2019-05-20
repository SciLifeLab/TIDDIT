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
Cov::Cov(int binSize,string bamFile,string output,int minQ,bool wig, bool skipQual, bool span){
	ostream& covout=test(coverageOutput,output);
	if(output != "stdout"){
		if (wig == false){
		coverageOutput.open((output+".tab").c_str());
		}else{
		coverageOutput.open((output+".wig").c_str());
		}
		static ostream& covout = coverageOutput;
		
	}else{
		static ostream& covout=cout;
	}



	//initialize the function
	this -> binStart =0;
	this -> binEnd=binSize+binStart;
	this -> binSize = binSize;
	this -> wig = wig;
	this -> skipQual = skipQual;
	this -> currentChr=-1;
	this -> minQ = minQ;
	this -> contigsNumber = 0;
	this -> bamFile = bamFile;
	this -> output = output;
	this -> span = span;
	if (wig == false){
		if (skipQual == true){
			covout << "#CHR" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" << endl;
		}else{
			covout << "#CHR" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" <<"\t" << "quality" << endl;
		}
	}
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
	this -> spanCoverageStructure.resize(2);

	qualityStructure[0].resize(contigsNumber);
	qualityStructure[1].resize(contigsNumber);

	if (this -> span){
		spanCoverageStructure[0].resize(contigsNumber);
		spanCoverageStructure[1].resize(contigsNumber);
	}

	for(int i=0;i<contigsNumber;i++){
		coverageStructure[i].resize(ceil(contigLength[i]/double(binSize)),0);

		qualityStructure[1][i].resize(ceil(contigLength[i]/double(binSize)),0);
		qualityStructure[0][i].resize(ceil(contigLength[i]/double(binSize)),0);
		if (this -> span){
			spanCoverageStructure[1][i].resize(ceil(contigLength[i]/double(binSize)),0);
			spanCoverageStructure[0][i].resize(ceil(contigLength[i]/double(binSize)),0);
		}
	}
	
	
	
    alignmentFile.Close();
}



//this function calculates the coverage in bins of size binSize across the entire bamfile
void Cov::bin(BamAlignment currentRead, readStatus alignmentStatus){
	//initialise the chromosome ID
	if(currentChr == -1){
		currentChr=currentRead.RefID;
	}
	
	//if the quality of the read is high enough, it will be added to the data structure
	if(alignmentStatus != lowQualty and alignmentStatus != unmapped) {
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
			while (remainingRead >= binSize and  coverageStructure[currentRead.RefID].size() > element+1 ){
				element++;
				coverageStructure[currentRead.RefID][element]+=binSize;
				qualityStructure[0][currentRead.RefID][element] += currentRead.MapQuality;
				qualityStructure[1][currentRead.RefID][element] += 1;
				remainingRead=remainingRead-binSize;
			}
			if (remainingRead > 0 and coverageStructure[currentRead.RefID].size() > element+1){
				element++;
				coverageStructure[currentRead.RefID][element]+=remainingRead;
				qualityStructure[0][currentRead.RefID][element] += currentRead.MapQuality;
				qualityStructure[1][currentRead.RefID][element] += 1;
			}

		}
		//span coverage based on discordants
		if ( alignmentStatus != pair_wrongDistance and alignmentStatus != pair_wrongOrientation and this -> span and currentRead.MapQuality > this-> minQ ){
			element=floor(double(currentRead.Position)/double(binSize));
			int elements=ceil(double(currentRead.InsertSize-binSize+ currentRead.Position-element*binSize)/double(binSize));
			if ( currentRead.Position-element*binSize < 20){
				element++;
				elements=elements-1;
			}

			for (int i=0;i < elements;i++){
				spanCoverageStructure[0][currentRead.RefID][element+i] +=1;
			}

			if (not currentRead.HasTag("SA")){
				elements=ceil(double(currentRead.Length-binSize+currentRead.Position-element*binSize)/double(binSize));
				for (int i=0;i < elements;i++){
					spanCoverageStructure[1][currentRead.RefID][element+i] +=1;
				}
			}
		}
	}

}

//prints the results
void Cov::printCoverage(){
	ostream& covout=test(coverageOutput,this -> output);
	if ( this -> output == "stdout") {
		static ostream& covout=cout;
	}
	if (this -> wig == false){
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
				if (this -> skipQual == false){
					covout << position2contig[i] << "\t" << binStart << "\t" << binEnd << "\t" << coverage << "\t" << quality << endl;
				}else{
					covout << position2contig[i] << "\t" << binStart << "\t" << binEnd << "\t" << coverage << endl;
				}

			}
		}
	}else{
		covout << "track type=wiggle_0 name=\"Coverage\" description=\"Per bin average coverage\"" << endl;
		for(int i=0;i<contigsNumber;i++){
			covout << "fixedStep chrom=" << position2contig[i] << " start=1 step=" << binSize << endl;

			for(int j=0;j<coverageStructure[i].size();j++){
				int binStart = j*binSize;
				int binEnd=binStart+binSize;
            
				if(binEnd > contigLength[i]){
					binEnd=contigLength[i];
				}
				double coverage=double(coverageStructure[i][j])/double(binEnd-binStart);
				covout << coverage << endl;

			}
		}
		if (this -> skipQual == false){
			covout << "track type=wiggle_0 name=\"MapQ\" description=\"Per bin average mapping quality\"" << endl;
			for(int i=0;i<contigsNumber;i++){
				covout << "fixedStep chrom=" << position2contig[i] << " start=1 step=" << binSize << endl;

				for(int j=0;j<coverageStructure[i].size();j++){
					double quality = 0;
					if(double(qualityStructure[1][i][j]) > 0){
						quality=double(qualityStructure[0][i][j])/double(qualityStructure[1][i][j]);
					}
					covout << quality << endl;
				}
			}
		}
		if (this -> span){
			covout << "track type=wiggle_0 name=\"SpanPairs\" description=\"Spanning pairs per bin\"" << endl;
			for(int i=0;i<contigsNumber;i++){
				covout << "fixedStep chrom=" << position2contig[i] << " start=1 step=" << binSize << endl;
				for(int j=0;j<coverageStructure[i].size();j++){
					covout << spanCoverageStructure[0][i][j] << endl;
				}
			}
			covout << "track type=wiggle_0 name=\"SpanReads\" description=\"Spanning reads per bin\"" << endl;
			for(int i=0;i<contigsNumber;i++){
				covout << "fixedStep chrom=" << position2contig[i] << " start=1 step=" << binSize << endl;
				for(int j=0;j<coverageStructure[i].size();j++){
					covout << spanCoverageStructure[1][i][j] << endl;
				}
			}
		}

	}
	
}
