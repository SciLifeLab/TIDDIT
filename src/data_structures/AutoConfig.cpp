#include "ProgramModules.h"
//The auto config module. This module finds the input parameters based on the input file

autoSettings::autoSettings(){}

//online calculation of the average value
double onlineAverage(double formerAVG,double currentVal,long long n){
	double delta = currentVal-formerAVG;
	double average= formerAVG+delta/double(n);
	return(average);
}

//returns 1 if a read pair is outtie, -1 if read pair in innie, 0 otherwise
int findOrientation(int firstPos,bool firstOrientation,int matePos,bool mateOrientation){
	int orientation=0;
	if(!firstOrientation and mateOrientation){
		orientation = -1;
	}else if(firstOrientation and !mateOrientation){
		orientation = 1;
	}	

	return(orientation);
} 



//this function is used to estimate the orientation,max insert size and min insert size.
vector<int> autoSettings::autoConfig(string bamFile,int Quality){
	vector<int> configuration;
	long orientation=0;
	double meanInsert=0;
	unsigned long long MInsert=0;
	int readLen=0;
	double stdReadLen=0;
	long long n=1;

	BamReader alignmentFile;
	//open the bam file
	alignmentFile.Open(bamFile);
	BamAlignment currentRead;

	int currentPos=0;
    int startPos=0;
	while ( alignmentFile.GetNextAlignmentCore(currentRead) ) {
			if(currentRead.IsMapped() and currentRead.MapQuality > Quality) {
				if(currentRead.RefID == currentRead.MateRefID){
					if(currentRead.Position <= currentRead.MatePosition){
						double insertSize=(currentRead.MatePosition-currentRead.Position+1)/1000.0;
						meanInsert=	onlineAverage(meanInsert,insertSize,n);	
						MInsert=meanInsert;
						readLen=onlineAverage(readLen,currentRead.Length,n);				
						orientation+= findOrientation(currentRead.Position,currentRead.IsReverseStrand(),currentRead.MatePosition,currentRead.IsMateReverseStrand());
						n++;
					}
				}
			}		
	}
	cout << "Auto Configuration:" << endl;
	

	int stdInsert=sqrt(MInsert)*1000;	
	meanInsert=meanInsert*1000;
	cout << meanInsert << endl;
	int maxInsert=2*meanInsert;
	//the min insert size is set to 3*avgread len
	int minInsert = 3*readLen;
	//outie if orientation > 0, innie otherwise
	if (orientation > 0){
		orientation=1;
		cout << "outtie" << endl;
	}else{
		orientation=0;
		cout << "innie" << endl;
	}


	configuration.push_back(maxInsert);configuration.push_back(minInsert);configuration.push_back(orientation);
	
	return(configuration);
}
