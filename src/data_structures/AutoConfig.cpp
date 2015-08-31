#include "ProgramModules.h"


autoSettings::autoSettings(){}

//online calculation of the average value
double onlineAverage(double formerAVG,double currentVal,long long n){
	double delta = currentVal-formerAVG;
	double average= formerAVG+delta/double(n);
	return(average);
}

//returns 1 if a read pair is outtie, -1 if read pair in inttie, 0 otherwise
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

	int i =0;
	while ( alignmentFile.GetNextAlignmentCore(currentRead) ) {
		//analyse the mapped reads exceeding the user set minimum value
			if(currentRead.IsMapped() and currentRead.MapQuality > Quality) {
				if(currentRead.RefID == currentRead.MateRefID){
					if(currentRead.Position <= currentRead.MatePosition){
						double insertSize=(currentRead.MatePosition-currentRead.Position+1)/1000.0;
						//no insert should be 100 mb large
						//MInsert = onlinevarriance(MInsert,meanInsert,insertSize,n);
						meanInsert=	onlineAverage(meanInsert,insertSize,n);	
						MInsert=meanInsert;
						readLen=onlineAverage(readLen,currentRead.Length,n);				
						orientation+= findOrientation(currentRead.Position,currentRead.IsReverseStrand(),currentRead.MatePosition,currentRead.IsMateReverseStrand());
						//cout << insertSize<< " " << meanInsert<< " " << MInsert <<  " " << n << endl;
						n++;
					}
				}
			}
		//}
		//i++;		
	}
	cout << "Auto Configuration:" << endl;
	
	//assuming normal distribution, the ~1% most extreme readpairs have a distance of avgInsert + 2.58*stdInsert
	int stdInsert=sqrt(MInsert)*1000;
	
	meanInsert=meanInsert*1000;
	cout << meanInsert << endl;
	int maxInsert=meanInsert;
	//the min insert size is set to 3*avgread len
	int minInsert = 3*readLen;
	cout << "Max insert size: " << maxInsert << endl;
	cout << "Min insert size: " << minInsert << endl;
	//outie if orientation > 0, inttie otherwise
	if (orientation > 0){
		orientation=1;
		cout << "outtie" << endl;
	}else{
		orientation=0;
		cout << "intie" << endl;
	}


	configuration.push_back(maxInsert);configuration.push_back(minInsert);configuration.push_back(orientation);
	
	return(configuration);
}
