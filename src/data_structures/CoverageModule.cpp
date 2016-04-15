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


vector<unsigned int> addBases(vector<unsigned int> sequencedBases,BamAlignment currentRead,int binStart,int binEnd,int currentChr,int binSize){
	//check the quality of the alignment, the hardcoded data is only used to initialise the function, the numbers themselves wont affect the quality measurement.
	readStatus alignmentStatus = computeReadType(currentRead, 100000,100, true);
	if(alignmentStatus != unmapped and alignmentStatus != lowQualty ) {
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
	}
	return(sequencedBases);
}

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

	covout << "#CHR" << "\t" << "start" << "\t" << "end" << "\t" << "coverage" <<"\t" << endl;

	//initialize the function
	this -> binStart =0;
	this -> binEnd=binSize+binStart;
	this -> binSize = binSize;
	this -> currentChr=-1;
	this -> sequencedBases.push_back(0);
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
	
    alignmentFile.Close();
}



//this function calculates the coverage in bins of size binSize across the entire bamfile
void Cov::bin(BamAlignment currentRead){
    ostream& covout = test(coverageOutput,output);
	//initialise the chromosome ID
	if(currentChr == -1){
		currentChr=currentRead.RefID;
	}

	sequencedBases=addBases(sequencedBases,currentRead,binStart,binEnd,currentChr,binSize);
	//when there are no more reads in the current bin, calculate the coverage and start analysing the next bin
	if(binEnd < currentRead.Position or currentChr != currentRead.RefID){
		//check if the read will fit into the next bin, otherwise, keep adding bins
		while(binEnd < currentRead.Position or currentChr != currentRead.RefID){
			
			if(binEnd > contigLength[currentChr]){
				binEnd=contigLength[currentChr];
			}
			double coverage=double(sequencedBases[0])/double(binEnd-binStart);
			covout << position2contig[currentChr] << "\t"<< binStart << "\t" << binEnd << "\t"<< coverage << endl;
			
			binStart=binEnd;
			binEnd=binStart+binSize;
			sequencedBases.erase(sequencedBases.begin());
			
			if(sequencedBases.size() == 0){sequencedBases.push_back(0);}

			if(currentChr != currentRead.RefID){
				//create bins to hold the last bases, and create empty bins to span the entire contig
				while (binEnd < contigLength[currentChr]){
					double coverage=double(sequencedBases[0])/double(binSize);
					covout << position2contig[currentChr] << "\t" << binStart << "\t" << binEnd << "\t" << coverage << endl;

					binStart=binEnd;
					binEnd=binStart+binSize;

					sequencedBases.erase(sequencedBases.begin());
					if(sequencedBases.size() == 0){
						sequencedBases.push_back(0);
					}
				}
				//clear the base container and start analysing the next chromosome		
				sequencedBases.clear();
				sequencedBases.push_back(0);
				currentChr = currentRead.RefID;
				binStart=0;
				binEnd=binSize;
			}
		}
			//now we are at the right position, then add the read
		sequencedBases=addBases(sequencedBases,currentRead,binStart,binEnd,currentChr,binSize);	
	}		
}
