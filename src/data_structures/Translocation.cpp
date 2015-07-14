/*
 * Translocations.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi, Eisfeldt
 */

#include "Translocation.h"


bool sortMate(long i, long  j) {
	return (i < j);
}

string VCFHeader(){
	string headerString ="";
	//Define fileformat and source
	headerString+="##fileformat=VCFv4.1\n";
	headerString+="##source=FindTranslocations\n";
	//define the alowed events
	headerString+="##ALT=<ID=DEL,Description=\"Deletion\">\n";
	headerString+="##ALT=<ID=DUP,Description=\"Duplication\">\n";
	headerString+="##ALT=<ID=BND,Description=\"Duplication\">\n";
	//Define the info fields
	headerString+="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
	headerString+="##INFO=<ID=LFW,Number=1,Type=Integer,Description=\"Links from window\">\n";
	headerString+="##INFO=<ID=LCB,Number=1,Type=Integer,Description=\"Links to chromosome B\">\n";
	headerString+="##INFO=<ID=LTE,Number=1,Type=Integer,Description=\"Links to event\"\b>\n";
	headerString+="##INFO=<ID=COVA,Number=1,Type=Integer,Description=\"Coverage on window A\">\n";
	headerString+="##INFO=<ID=COVB,Number=1,Type=Integer,Description=\"Coverage on window B\">\n";
	headerString+="##INFO=<ID=OA,Number=1,Type=Integer,Description=\"Orientation of the reads in window A\">\n";
	headerString+="##INFO=<ID=OB,Number=1,Type=Integer,Description=\"Orientation of the mates in window B\">\n";
	headerString+="##INFO=<ID=CHRA,Number=1,Type=String,Description=\"The chromosome of window A\">\n";
	headerString+="##INFO=<ID=CHRB,Number=1,Type=String,Description=\"The chromosome of window B\">\n";
	headerString+="##INFO=<ID=WINA,Number=2,Type=Integer,Description=\"start and stop positon of window A\">\n";
	headerString+="##INFO=<ID=WINB,Number=2,Type=Integer,Description=\"start and stop position of window B\">\n";
	headerString+="##INFO=<ID=EL,Number=1,Type=Integer,Description=\"Expected links to window B\">\n";
	headerString+="##INFO=<ID=RATIO,Number=1,Type=Integer,Description=\"The number of links divided by the expected number of links\">\n";
	headerString+="##INFO=<ID=ED,Number=1,Type=Integer,Description=\"The average estimated distance between paired ends within the window\">\n";
	//set filters
	headerString+="##FILTER=<ID=BelowExpectedLinks,Description=\"The number of links between A and B is less than 40\% of the expected value\">\n";
	headerString+="=##FILTER=<ID=FewLinks,Description=\"Fewer than 40% of the links in window A link to window B\">\n";
	headerString+="=##FILTER=<ID=UnexpectedDistance,Description=\"The average paired reads distance is deviating\">\n";
	//Header
	headerString+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";		

	return(headerString);		
}

string filterFunction(double RATIO, int linksToB, int LinksFromWindow,float mean_insert, float std_insert,int estimatedDistance){
	string filter = "PASS";
	double linkratio= (double)linksToB/(double)LinksFromWindow;
	if(RATIO < 0.4){
		filter="BelowExpectedLinks";
	}else if(linkratio < 0.4){
		filter="FewLinks";
	}else if(estimatedDistance > 2*std_insert + mean_insert or estimatedDistance < -2*std_insert+mean_insert){
		filter="UnexpectedDistance";
	}


	return(filter);
}

Window::Window(int max_insert, uint16_t minimum_mapping_quality,
		bool outtie, float mean_insert, float std_insert, int minimumPairs,
		float meanCoverage, string outputFileHeader,string bamFileName, string indexFile) {
	this->max_insert		 = max_insert;
	this->minimum_mapping_quality = minimum_mapping_quality;
	this->outtie			 = outtie;
	this->mean_insert		 = mean_insert;
	this ->std_insert		 = std_insert;
	this->minimumPairs		 = minimumPairs;
	this->meanCoverage		 = meanCoverage;
	this->bamFileName		=bamFileName;
	this -> indexFile		=indexFile;

	this->outputFileHeader   = outputFileHeader;
	string inter_chr_events = outputFileHeader + "_inter_chr_events.tab";
	this->interChrVariations.open(inter_chr_events.c_str());
	this->interChrVariations <<"chrA\tstartOnA\tendOnA\tchrB\tstartOnB\tendOnB\tLinksFromWindow\tLinksToChrB\tLinksToEvent\tCoverageOnChrA\tCoverageOnChrB\tOrientationA\tOrientationB\tExpectedLinksToEvent\tLinksToEvent/ExpectedLinksToEvent\tEstimatedDistance\n";
	string intra_chr_events = outputFileHeader + "_intra_chr_events.tab";
	this->intraChrVariations.open(intra_chr_events.c_str());
	this->intraChrVariations << "chrA\tstartOnA\tendOnA\tchrB\tstartOnB\tendOnB\tLinksFromWindow\tLinksToChrB\tLinksToEvent\tCoverageOnChrA\tCoverageOnChrB\tOrientationA\tOrientationB\tExpectedLinksToEvent\tLinksToEvent/ExpectedLinksToEvent\tEstimatedDistance\n";

	string inter_chr_eventsVCF = outputFileHeader + "_inter_chr_events.vcf";
	this->interChrVariationsVCF.open(inter_chr_eventsVCF.c_str());
	this->interChrVariationsVCF << VCFHeader();

	string intra_chr_eventsVCF = outputFileHeader + "_intra_chr_events.vcf";
	this->intraChrVariationsVCF.open(intra_chr_eventsVCF.c_str());
	this->intraChrVariationsVCF << VCFHeader();

	this->coverage          	= 0;
	this->chr			=-1;

}

void Window::insertRead(BamAlignment alignment) {
	
	readStatus alignmentStatus = computeReadType(alignment, this->max_insert, this->outtie);
	if(alignmentStatus == unmapped or alignmentStatus == lowQualty ) {
		return; // in case the alignment is of no use discard it
	}

	if(this->chr == -1) { //sets chr to the first chromosome at startup
		cout << "working on sequence " << position2contig[alignment.RefID] << "\n";
		chr=alignment.RefID;
	}

	if(alignment.RefID != chr) { // I am moving to a new chromosomes, need to check if the current window can be used or not		
		for(int i=0;i<eventReads.size();i++){
			//check for events
			if( eventReads[i].size() >= minimumPairs){
				computeVariations(i);
			}
			//empty the alignmentqueues
			eventReads[i]=queue<BamAlignment>();

		}
		cout << "working on sequence " << position2contig[alignment.RefID] << "\n";
	}

	if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
		if(alignment.RefID < alignment.MateRefID or (alignment.RefID == alignment.MateRefID and alignment.Position < alignment.MatePosition)) {  // insert only "forward" variations

			int alignmentNumber= eventReads[alignment.MateRefID].size();


			//if there currently is no sign of a translocation event between the current chromosome and the chromosome of the current alignment
			if( alignmentNumber == 0){
				eventReads[alignment.MateRefID].push(alignment);
				
				
			}else{
				//if there already are alignments of the current chromosome and chromosome of the current alignment:

				//Check if the distance between the current and the previous alignment is larger than max distance
				int currrentAlignmentPos=alignment.Position;
				int pastAlignmentPos= eventReads[alignment.MateRefID].back().Position;
				int distance= currrentAlignmentPos - pastAlignmentPos;

				//If the distance between the two reads is less than the maximum allowed distace, add it to the other reads of the event
				int MaxDistance=1000;
				if(distance <= 1000){
					//add the read to the current window
					eventReads[alignment.MateRefID].push(alignment);
				}else{
					if(eventReads[alignment.MateRefID].size() >= minimumPairs){
						computeVariations(alignment.MateRefID);
					}
					//Thereafter the alignments are removed
					eventReads[alignment.MateRefID]=queue<BamAlignment>();	

					//the current alignment is inserted
					eventReads[alignment.MateRefID].push(alignment);
							
				}
			}
		}

	}

	chr=alignment.RefID;



}

float Window::computeCoverageB(int chrB, int start, int end, int32_t secondWindowLength){
	int bases=0;
	float coverage;
	BamReader bamFile;
	BamAlignment currentRead;


	//test if an index file is available

	if(!bamFile.Open(bamFileName)){
		return -1;
	}else{
		if(bamFile.OpenIndex(indexFile) == 0){
			cout << "warning no index file found, extraction will proceed in slow mode" << endl;
		}
		//moves to a region and itterates through every read inside that region
		bamFile.SetRegion(chrB,start,chrB,end+1); 
		while ( bamFile.GetNextAlignment(currentRead) ) {
			if(start <= currentRead.Position and end >= currentRead.Position){
				readStatus alignmentStatus = computeReadType(currentRead, this->max_insert, this->outtie);
				if(alignmentStatus != unmapped and alignmentStatus != lowQualty ) {
					//calculates the length of a read
					bases+=currentRead.Length;
				}
			}
		}
	bamFile.Close();
		//calculates the coverage and returns the coverage within the window
		coverage=bases/float(secondWindowLength+1);
		return(coverage);
	}	
	
}

//This function accepts a queue with alignments that links to chrB from chrA, the function returns the number of reads that are positioned inside the region on A and have mates that are linking anywhere on B
vector<int> findLinksToChr2(queue<BamAlignment> ReadQueue,long startChrA,long stopChrA,long startChrB,long endChrB, int pairsFormingLink){
	int linkNumber=0;
	int QueueSize=ReadQueue.size();
	int estimatedDistance=0;
	int lengthWindowA=stopChrA-startChrA+1;

	for(int i=0;i<QueueSize;i++){
		if(startChrA <= ReadQueue.front().Position and stopChrA >= ReadQueue.front().Position ){
			linkNumber+=1;
			//If the reads are bridging the two windows, calculate the mate and read distance
			if(startChrB <= ReadQueue.front().MatePosition and endChrB >= ReadQueue.front().MatePosition){
				estimatedDistance+=lengthWindowA-(ReadQueue.front().Position-startChrA)+(ReadQueue.front().MatePosition-startChrB);
			}
		}
		ReadQueue.pop();
	}
	vector<int> output;
	output.push_back(linkNumber);
	if(pairsFormingLink > 0){
		output.push_back(estimatedDistance/pairsFormingLink);	
	}else{
		//if no links exists, the region will not count as an event, thus any number may be returned as the distance(the statistics wont be printed anyway)
		output.push_back(-1);
	}
	return(output);
}

vector<double> Window::computeStatisticsA(string bamFileName, int chrB, int start, int end, int32_t WindowLength, string indexFile){
	vector<double> statisticsVector;
	int bases=0;
	double coverage;
	double AverageReadLength=0;
	BamReader bamFile;
	BamAlignment currentRead;
	int linksFromWindow=0;
	int nreads=0;

	if(!bamFile.Open(bamFileName)){
		statisticsVector.push_back(-1);statisticsVector.push_back(-1);
		return(statisticsVector);
	}else{
		if(bamFile.OpenIndex(indexFile) == 0){
			cout << "warning no index file found, extraction will proceed in slow mode" << endl;
		}
		//moves to a region and itterates through every read inside that region
		bamFile.SetRegion(chrB,start,chrB,end+1); 
		while ( bamFile.GetNextAlignment(currentRead) ) {
			//makes sure that we are inside the ragion
			if(start <= currentRead.Position and end >= currentRead.Position){
				readStatus alignmentStatus = computeReadType(currentRead, this->max_insert, this->outtie);
				//only mapped and high quality reads are used
				if(alignmentStatus != unmapped and alignmentStatus != lowQualty ) {
					//calculates the length of a read
					bases+=currentRead.Length;
					AverageReadLength+=currentRead.Length;
					nreads++;
				
					//if the distance between a pair is too long, the amount of links from the winddoe is increased
					if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
						if(currentRead.RefID < currentRead.MateRefID or (currentRead.RefID == currentRead.MateRefID and currentRead.Position < currentRead.MatePosition)) {
							linksFromWindow+=1;
						}
					}
				}
			}
		}
	AverageReadLength=AverageReadLength/nreads;
	bamFile.Close();
		//calculates the coverage and returns the coverage within the window
		coverage=bases/float(WindowLength+1);
		statisticsVector.push_back(coverage);statisticsVector.push_back(linksFromWindow);statisticsVector.push_back(AverageReadLength);
		return(statisticsVector);
	}	
	
}



bool Window::computeVariations(int chr2) {
	bool found = false;
		
	//by construction I have only forward links, i.e., from chr_i to chr_j with i<j
	int realFirstWindowStart=this->eventReads[chr2].front().Position;
	int realFirstWindowEnd=this->eventReads[chr2].back().Position;
	
	vector<long> chr2regions= findRegionOnB( this-> eventReads[chr2] ,minimumPairs,1000);
	for(int i=0;i < chr2regions.size()/3;i++){
		long startSecondWindow=chr2regions[i*3];
		long stopSecondWindow=chr2regions[i*3+1];	
		long pairsFormingLink=chr2regions[i*3+2];

		

		//resize region so that it is just large enough to cover the reads that have mates in the present cluster
		vector<long> chrALimit=newChrALimit(this-> eventReads[chr2],startSecondWindow,stopSecondWindow);
		

		long startchrA=chrALimit[0];
		long stopchrA=chrALimit[1];

		vector<int> statsOnB=findLinksToChr2(eventReads[chr2],startchrA, stopchrA,startSecondWindow,stopSecondWindow,pairsFormingLink);
		int numLinksToChr2=statsOnB[0];
		int estimatedDistance=statsOnB[1];

		vector<string> orientation=computeOrientation(eventReads[chr2],startchrA ,stopchrA,startSecondWindow,stopSecondWindow);
		string read1_orientation= orientation[0];
		string read2_orientation= orientation[1];
	

		vector<double> statisticsFirstWindow =computeStatisticsA(bamFileName, eventReads[chr2].front().RefID, startchrA, stopchrA, (stopchrA-startchrA), indexFile);
		double coverageRealFirstWindow	= statisticsFirstWindow[0];
		int linksFromWindow=int(statisticsFirstWindow[1]);
		int averageReadLength=int(statisticsFirstWindow[2]);
		double coverageRealSecondWindow=computeCoverageB(chr2, startSecondWindow, stopSecondWindow, (stopSecondWindow-startSecondWindow+1) );

		if(pairsFormingLink >= minimumPairs  ) { //ration between coverage  // and coverageRealFirstWindow < 5*this->meanCoverage
		//I need to find the window that maximise this
			int secondWindowLength=(stopSecondWindow-startSecondWindow+1);
			int firstWindowLength=stopchrA-startchrA+1;
			float expectedLinksInWindow = ExpectedLinks(firstWindowLength, secondWindowLength, estimatedDistance, mean_insert, std_insert, coverageRealFirstWindow, averageReadLength);
			if(this->chr == chr2) {
				intraChrVariations << position2contig[this->chr]  << "\t" <<     startchrA   << "\t" <<       stopchrA               << "\t"  ;
				intraChrVariations << position2contig[chr2]       << "\t" <<         startSecondWindow  << "\t" <<        stopSecondWindow                << "\t"  ;
				intraChrVariations <<      linksFromWindow        << "\t" <<        numLinksToChr2      << "\t" <<          pairsFormingLink              << "\t";
				intraChrVariations <<    coverageRealFirstWindow   << "\t" << coverageRealSecondWindow << "\t";
    				intraChrVariations <<	  read1_orientation    << "\t" <<         read2_orientation << "\t"
				<< expectedLinksInWindow << "\t" << pairsFormingLink/expectedLinksInWindow << "\t" << estimatedDistance << "\n";

				string filter=filterFunction(pairsFormingLink/expectedLinksInWindow, pairsFormingLink, linksFromWindow,mean_insert, std_insert,estimatedDistance);
				string svType="BND";

				intraChrVariationsVCF << position2contig[this -> chr]  << "\t" <<     stopchrA   << "\t.\t"  ;
				intraChrVariationsVCF << "N"       << "\t"	<< "N[" << position2contig[chr2] << ":" << startSecondWindow << "[";
     				intraChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
				intraChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
				intraChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
				intraChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
				intraChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow <<";ED=" << estimatedDistance << "\n";

	
			} else {
				interChrVariations << position2contig[this -> chr]  << "\t" <<     startchrA   << "\t" <<       stopchrA               << "\t"  ;
				interChrVariations << position2contig[chr2]       << "\t" <<         startSecondWindow  << "\t" <<         stopSecondWindow               << "\t"  ;
				interChrVariations <<     linksFromWindow        << "\t" <<        numLinksToChr2      << "\t" <<          pairsFormingLink              << "\t";
				interChrVariations <<     coverageRealFirstWindow	<< "\t" << coverageRealSecondWindow << "\t";
				interChrVariations <<       read1_orientation    << "\t" <<         read2_orientation              << "\t"
				 << expectedLinksInWindow << "\t" << pairsFormingLink/expectedLinksInWindow << "\t" << estimatedDistance << "\n";

				string filter=filterFunction(pairsFormingLink/expectedLinksInWindow, pairsFormingLink, linksFromWindow,mean_insert, std_insert,estimatedDistance);
				string svType="BND";

				interChrVariationsVCF << position2contig[this -> chr]  << "\t" <<     stopchrA   << "\t.\t"  ;
				interChrVariationsVCF << "N"       << "\t"	<< "N[" << position2contig[chr2] << ":" << startSecondWindow << "[";
     				interChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
				interChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
				interChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
				interChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
				interChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow <<";ED=" << estimatedDistance << "\n";

			}
			found = true;
				
			
		
		}
	}
		return(found);

}

void Window::initTrans(SamHeader head) {
	uint32_t contigsNumber = 0;
	SamSequenceDictionary sequences  = head.Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		this->contig2position[sequence->Name] = contigsNumber; // keep track of contig name and position in order to avoid problems when processing two libraries
		this->position2contig[contigsNumber]  = sequence->Name;
		contigsNumber++;
	}

}

//accepts two queues as argument, the second queue is appended onto the first, the total queue is returned.
queue<BamAlignment> Window::queueAppend(queue<BamAlignment> queueOne,queue<BamAlignment> queueTwo){
	for(int i=0; i<queueTwo.size();i++){
	queueOne.push(queueTwo.front());
	queueTwo.pop();
	}
	
	
	return(queueOne);
}

//this is a function that takes a queue containing reads that are involved in an event and calculates the position of the event on chromosome B(the chromosome of the of the mates). The total number of reads
//inside the region is reported, aswell as the start and stop position on B
vector<long> Window::findRegionOnB( queue<BamAlignment> alignmentQueue, int minimumPairs,int maxDistance){
	vector<long> mate_positions;
	queue<long> linksToRegionQueue;
	//The positions of the mates are transfered to a vector
	mate_positions.resize(alignmentQueue.size());
	int QueueSize=alignmentQueue.size();
	for(int i=0; i< QueueSize;i++){
		mate_positions[i]=alignmentQueue.front().MatePosition;
		alignmentQueue.pop();
	}
	
	vector<long> eventRegionVector;
	//sort the mate positions
	sort(mate_positions.begin(), mate_positions.end(), sortMate);
	//finds the region of the event
	linksToRegionQueue.push(mate_positions[0]);
	mate_positions.erase(mate_positions.begin());
	for(int i=0; i< QueueSize-1;i++){

		//if the current mate is close enough to the previous mate
		if( maxDistance >= (mate_positions[0]-linksToRegionQueue.back()) ){
			linksToRegionQueue.push(mate_positions[0]);
			mate_positions.erase(mate_positions.begin());

		}else{
			//if there are enough links beetween two regions, the region and the number of pairs are saved, and later on returned
			if(linksToRegionQueue.size() >= minimumPairs  ) {
				//the front of the queue is the start position
				eventRegionVector.push_back(linksToRegionQueue.front());
				//the back is the end position
				eventRegionVector.push_back(linksToRegionQueue.back());
				//the length is the number of links to this event
				eventRegionVector.push_back(linksToRegionQueue.size());
			}
			//the queue is reset so that the next event may be found(if any)
			linksToRegionQueue=queue<long>();
			linksToRegionQueue.push(mate_positions[0]);
			mate_positions.erase(mate_positions.begin());
		}
		
	}
	//if the queue contain any event 
	if(linksToRegionQueue.size() >= minimumPairs  ) {
		eventRegionVector.push_back(linksToRegionQueue.front());
		eventRegionVector.push_back(linksToRegionQueue.back());
		eventRegionVector.push_back(linksToRegionQueue.size());
			}
	

	return(eventRegionVector);

}

//this function resizes the region of CHRA so that it fits only around reads that are linking to an event on chrB
vector<long> Window::newChrALimit(queue<BamAlignment> alignmentQueue,long Bstart,long Bend){
	vector<long> startStopPos;
	long startA=-1;
	long endA=-1;
	int len=alignmentQueue.size();
	for(int i=0;i< len;i++){
		//Checks if the alignment if the mate of the current read is located inside the chrB region	
		if(alignmentQueue.front().MatePosition >= Bstart and alignmentQueue.front().MatePosition <= Bend){
		//resize the chr A region so that it fits around all the reads that have a mate inside region 
			if(alignmentQueue.front().Position <= startA or startA ==-1){
				startA=alignmentQueue.front().Position;
			}
			if(alignmentQueue.front().Position >= endA or endA ==-1){
				endA=alignmentQueue.front().Position;
			}

		}
		alignmentQueue.pop();
	

	}
	
	startStopPos.push_back(startA);
	startStopPos.push_back(endA);
	return(startStopPos);
}

//This function computes the orientation of the pairs
vector<string> Window::computeOrientation(queue<BamAlignment> alignmentQueue,long Astart,long Aend,long Bstart,long Bend){
	vector<string> orientationVector;
	int numberOfReads=alignmentQueue.size();
	int numberOfReverseReads=0;
	int numberOfReverseMates=0;
	int readsInRegion=0;

	for(int i=0;i<numberOfReads;i++){
		//Checks if the current read is located inside the region
		if(alignmentQueue.front().Position >=Astart and alignmentQueue.front().Position <= Aend){
			if(alignmentQueue.front().MatePosition >=Bstart and alignmentQueue.front().MatePosition <= Bend){
				//check the direction of the read on chrA
				if(alignmentQueue.front().IsReverseStrand()){
					numberOfReverseReads+=1;
				}

				//check the direction of the mate
				if(alignmentQueue.front().IsMateReverseStrand()){
					numberOfReverseMates+=1;
				}
				readsInRegion+=1;

			}
		}
		alignmentQueue.pop();
	}

	double FractionReverseReads=(double)numberOfReverseReads/(double)readsInRegion;
	string readOrientation="";
	ostringstream convertRead;
	

	if(FractionReverseReads > 0.5 ){
		convertRead << round(FractionReverseReads*100);
		readOrientation="- (" + convertRead.str() + "%)";
	}else{
		convertRead << 100-round(FractionReverseReads*100);
		readOrientation="+ (" + convertRead.str() + "%)";
	}

	string mateOrientation="";
	ostringstream convertMate;
	double FractionReverseMates=(double)numberOfReverseMates/(double)readsInRegion;

	if(FractionReverseMates > 0.5 ){
		convertMate << round(FractionReverseMates*100);
		mateOrientation="- (" + convertMate.str() + "%)";
	}else{
		convertMate << 100-round(FractionReverseMates*100);
		mateOrientation="+ (" + convertMate.str() + "%)";
	}

	orientationVector.push_back(readOrientation);orientationVector.push_back(mateOrientation);
	return(orientationVector);
}




