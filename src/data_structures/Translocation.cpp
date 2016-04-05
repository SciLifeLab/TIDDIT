/*
 * Translocations.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi, Eisfeldt
 */

#include "Translocation.h"
#include <string>

string int2str(int to_be_converted){
	string converted= static_cast<ostringstream*>( &(ostringstream() << to_be_converted) )->str();
	return(converted);
}

//converts a string to int
int str2int(string to_be_converted){
	int converted;
	istringstream convert(to_be_converted);
	return(converted);
}


//this function tries to classify intrachromosomal events
vector<string> Window::classification(int chr, int startA,int endA,double covA,int startB,int endB,double covB,int meanInsert,int STDInsert,bool outtie,vector<double> isReverse){
	string svType="BND";
	string start=int2str(startA);
	string end= int2str(endB);
	
	string GT="./1";
	int CN=1;

	double coverage= this -> meanCoverage;
	//the tolerance is dependent on read depth
	float coverageToleranceFraction = 0.6*pow((coverage+1)/coverage,2);
	float coverageTolerance=this -> meanCoverage*coverageToleranceFraction*double(this -> ploidy);
	float covAB;
	if(endA <= startB){
		covAB=computeCoverageB(chr,endA,startB, startB-endA);
	}else{
		covAB=computeCoverageB(chr,startB,endA, endA-startB);
	}
	double averageInsert=isReverse[2];
	

	// for heterozygous events
	CN = round(covAB / (coverage/double(this->ploidy)) ) - 1; // CN sounds like ratio to coverage per chr, minus 1? at least for het vars with the other allele at count 1.
	if (CN < 0) { CN = 0; }
        else if ( (CN+1 < 1.5 + 0.5*coverageToleranceFraction) && (CN+1 > 1 + 0.5*coverageToleranceFraction)) { 
	  // dup genotype is 0/1
          GT="0/1";
	}

	//if the orientation of one of the reads is normal, but the other read is inverted, the variant is an inversion
	//both reads are pointing to the left
	if(isReverse[0] > 0.5 and isReverse[1] > 0.5){
			svType= "INV";
	}else if(isReverse[0] < 0.5 and isReverse[1] < 0.5){
			svType= "INV";
	}

	if(svType == "BND"){
		//Find large tandemduplications
		if(outtie == true){
			if(isReverse[0] < 0.5 and isReverse[1] > 0.5){
    		    if( covA > coverage+coverageTolerance and covB > coverage+coverageTolerance ){
    				svType = "TDUP";
    			}
			}

		}else{
			if(isReverse[0] > 0.5 and isReverse[1] < 0.5){
        	    //if the coverage of the two windows are too high aswell, the event is a tandem dulplication
    		    if( covA > coverage+coverageTolerance and covB > coverage+coverageTolerance){
    				svType = "TDUP";

    			}
			}
		}
	}

	//if the insert size is longer than expected 
	if(svType == "BND"){
		if(averageInsert > this->max_insert){
			if(covAB < coverage-coverageTolerance){
				svType= "DEL";
				if (covAB < coverageTolerance) {
				  GT = "1/1";
				  //				  CN = 0; 
				} else {
				  GT = "0/1";
				  // CN = 1;
				}
			}//if the coverage of all regions is normal and region A and B are disjunct, the event is an insertion(mobile element)
			else if( ( covAB < coverage+coverageTolerance and covAB > coverage-coverageTolerance ) and 
				( covB < coverage+coverageTolerance and covB > coverage-coverageTolerance )	and
				( covA < coverage+coverageTolerance and covA > coverage-coverageTolerance )){
				svType = "INS";
			}
			//if A and B are disjunct and only one of them have coverage above average, the event is an interspersed duplication
			else if( (covA > coverage+coverageTolerance and covB < coverage+coverageTolerance ) or (covB > coverage+coverageTolerance and covA < coverage+coverageTolerance ) ){
				svType = "IDUP";
				// it would be a good idea to split this into one DUP and one BND, indicating the insertion point
				//label the region marked as high coverage as the duplications
				if(covA > coverage+coverageTolerance){
					start=int2str(startA);
					end = int2str(endA);
				}else if(covB > coverage+coverageTolerance){
					start=int2str(startB);
					end=int2str(endB);			
				}
			}
		}
	}
		
	//if the event is a duplication, but the exact type cannot be specified
	if(svType == "BND"){
		if(covA > coverage+coverageTolerance or covB > coverage+coverageTolerance or covAB > coverage+coverageTolerance){
			svType = "DUP";
		}
		//if a deletion is small enough, the coverage will be low on the sides of the variation
		if(covA < coverage-coverageTolerance and covB < coverage-coverageTolerance){
			svType = "DEL";
			if (covA < coverageTolerance or covB < coverageTolerance) {
			  GT = "1/1";
			  CN = 0;
			}
		}
	}
	vector<string> svVector;
	svVector.push_back(svType);svVector.push_back(start);svVector.push_back(end);svVector.push_back(GT);svVector.push_back(int2str(CN));
	return(svVector);
}

string filterFunction(double RATIO, int linksToB, int LinksFromWindow,float mean_insert, float std_insert,int estimatedDistance,double coverageA,double coverageB,double coverageAVG){
	string filter = "PASS";
	double linkratio= (double)linksToB/(double)LinksFromWindow;
	if(RATIO < 0.4){
		filter="BelowExpectedLinks";
	} else if(linkratio < 0.4){
		filter="FewLinks";
	} else if(coverageA > 10*coverageAVG or coverageB > coverageAVG*10){
		filter="UnexpectedCoverage";
	}


	return(filter);
}

bool sortMate(long i, long  j) {
	return (i < j);
}

string Window::VCFHeader(){
	string headerString ="";
	//Define fileformat and source
	headerString+="##fileformat=VCFv4.1\n";
	headerString+="##source=FindTranslocations\n";
	//define the alowed events
	headerString+="##ALT=<ID=DEL,Description=\"Deletion\">\n";
	headerString+="##ALT=<ID=DUP,Description=\"Duplication\">\n";
	headerString+="##ALT=<ID=TDUP,Description=\"Tandem duplication\">\n";
	headerString+="##ALT=<ID=IDUP,Description=\"Interspersed duplication\">\n";
	headerString+="##ALT=<ID=INV,Description=\"Inversion\">\n";
	headerString+="##ALT=<ID=INS,Description=\"Insertion\">\n";
	headerString+="##ALT=<ID=BND,Description=\"Break end\">\n";
	
	//Define the info fields
	headerString+="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
	headerString+="##INFO=<ID=END,Number=1,Type=String,Description=\"End of an intra-chromosomal variant\">\n";
	headerString+="##INFO=<ID=LFW,Number=1,Type=Integer,Description=\"Links from window\">\n";
	headerString+="##INFO=<ID=LCB,Number=1,Type=Integer,Description=\"Links to chromosome B\">\n";
	headerString+="##INFO=<ID=LTE,Number=1,Type=Integer,Description=\"Links to event\">\n";
	headerString+="##INFO=<ID=COVA,Number=1,Type=Float,Description=\"Coverage on window A\">\n";
	headerString+="##INFO=<ID=COVB,Number=1,Type=Float,Description=\"Coverage on window B\">\n";
	headerString+="##INFO=<ID=OA,Number=1,Type=String,Description=\"Orientation of the reads in window A\">\n";
	headerString+="##INFO=<ID=OB,Number=1,Type=String,Description=\"Orientation of the mates in window B\">\n";
	headerString+="##INFO=<ID=CHRA,Number=1,Type=String,Description=\"The chromosome of window A\">\n";
	headerString+="##INFO=<ID=CHRB,Number=1,Type=String,Description=\"The chromosome of window B\">\n";
	headerString+="##INFO=<ID=WINA,Number=2,Type=Integer,Description=\"start and stop positon of window A\">\n";
	headerString+="##INFO=<ID=WINB,Number=2,Type=Integer,Description=\"start and stop position of window B\">\n";
	headerString+="##INFO=<ID=EL,Number=1,Type=Float,Description=\"Expected links to window B\">\n";
	headerString+="##INFO=<ID=RATIO,Number=1,Type=Float,Description=\"The number of links divided by the expected number of links\">\n";
	//set filters
	headerString+="##FILTER=<ID=BelowExpectedLinks,Description=\"The number of links between A and B is less than 40\% of the expected value\">\n";
	headerString+="##FILTER=<ID=FewLinks,Description=\"Fewer than 40% of the links in window A link to chromosome B\">\n";
	headerString+="##FILTER=<ID=UnexpectedCoverage,Description=\"The coverage of the window on chromosome B or A is higher than 10*average coverage\">\n";
	//set format 
	headerString+="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	headerString+="##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n";
	headerString+="##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-ends that support the event\">\n";
	headerString+="##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads that support the event\">\n";
	//Header
	headerString+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"; // now just add the individuals!
	return(headerString);
}



Window::Window(string bamFileName, bool outtie, float meanCoverage,string outputFileHeader, map<string,int> SV_options) {
	this->max_insert		 = SV_options["max_insert"];
	this->minimum_mapping_quality = SV_options["mapping_quality"];
	this->outtie			 = outtie;
	this->mean_insert		 = SV_options["meanInsert"];
	this ->std_insert		 = SV_options["STDInsert"];
	this->minimumPairs		 = SV_options["pairs"];
	this->meanCoverage		 = meanCoverage;
	this->bamFileName		=bamFileName;
	this -> ploidy        =SV_options["ploidy"];
	this -> readLength = SV_options["readLength"];

	this->outputFileHeader   = outputFileHeader;
	string inter_chr_eventsVCF = outputFileHeader + "_inter_chr_events.vcf";
	this->interChrVariationsVCF.open(inter_chr_eventsVCF.c_str());
	this->interChrVariationsVCF << VCFHeader() << outputFileHeader << "\n";

	string intra_chr_eventsVCF = outputFileHeader + "_intra_chr_events.vcf";
	this->intraChrVariationsVCF.open(intra_chr_eventsVCF.c_str());
	this->intraChrVariationsVCF << VCFHeader() << outputFileHeader << "\n";

}

void Window::insertRead(BamAlignment alignment) {
	readStatus alignmentStatus = computeReadType(alignment, this->max_insert,this->min_insert, this->outtie);

	if( not alignment.IsMateMapped()  or alignment.MapQuality < minimum_mapping_quality ) {
		return; // in case the alignment is of no use discard it
	}

	vector <int > clipSizes;
	vector< int > readPositions;
	vector<int> genomePos;
	bool test=alignment.GetSoftClips(clipSizes,readPositions,genomePos);
	bool alignment_split = false;	
	if (test){
		alignment.BuildCharData(); // Somewhat costly, but needed to check for any SA tag.			
		alignment_split = alignment.HasTag("SA");
	}

	if( not alignment_split and (not alignment.IsMateMapped()) or alignment.MapQuality < minimum_mapping_quality) { 
	  return; // discard alignment - no use for the following. But keep splits for now.
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
			eventSplitReads[i]=vector<BamAlignment>();
			linksFromWin[alignment.MateRefID]=queue<int>();
		}
		cout << "working on seqence " << position2contig[alignment.RefID] << "\n";
	}

	if (alignment_split) {
	  
	  // parse split read to get the other segment position, akin to a mate.

	  string SA;
	  alignment.GetTag("SA",SA);
	  /* From the VCF documentation:
	  // Other canonical alignments in a chimeric alignment, formatted as a semicolon-delimited list: 
	  // (rname,pos,strand,CIGAR,mapQ,NM;)+. 
	  // Each element in the list represents a part of the chimeric alignment. 
	  // Conventionally, at a supplementary line, the first element points to the primary line.
	  */

	  vector <string> SA_elements;
      stringstream ss(SA);
	  std::string item;
      while (std::getline(ss, item, ';')) {
        stringstream SS(item);
        string SA_data;
        while (std::getline(SS, SA_data, ',')) {
      	    SA_elements.push_back(SA_data);
      	}
      }

	  //vector <string> SA_elements;
	  //Now split on ,
	  //
 	  //
      //
      //	SA_elements.push_back(item);
      //}
	  // there will be an "extra" empty string element at the end of the list, corresponding to the split of the final ;
	  
	  int contigNr = contig2position[SA_elements[0]];
	  //for(vector<string>::iterator it = SA_elements.begin(), it_end = --SA_elements.end(); it != it_end; it += 6) {
	    //string contig = *it;
	   // cout << contig << endl;
	   // contigNr = this -> ;
	  //}
	  
	  //Check if the distance between the split read and the previous alignment is larger than max distance


	  bool add = false;
	  // we may or may not have a past event
	  if (eventReads[contigNr].size() > 0) {	    
	    int currrentAlignmentPos=alignment.Position;
	    int pastAlignmentPos= eventReads[contigNr].back().Position;
	    int distance= currrentAlignmentPos - pastAlignmentPos;
	    // If the distance between the two reads is less than the maximum allowed distace, add it to the other reads of the event
	    if(distance <= 1000){
	      // parse tag and add it to the destination chr
	      add = true;
	    }
	  }
	  if (add) {;
	    eventSplitReads[contigNr].push_back(alignment);    
	  } else {
	    // for now, we don't allow split reads to start a new window.
	    //	    cout << "Warning: split read outside current window found. Ignored.\n";
	  }
	}
	if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
		if(alignment.RefID < alignment.MateRefID or (alignment.RefID == alignment.MateRefID and alignment.Position < alignment.MatePosition)) {  // insert only "forward" variations

			int alignmentNumber= eventReads[alignment.MateRefID].size();

			//if there currently is no sign of a translocation event between the current chromosome and the chromosome of the current alignment
			if( alignmentNumber == 0) {
				eventReads[alignment.MateRefID].push(alignment);
			} else {
				//if there already are alignments of the current chromosome and chromosome of the current alignment:

				//Check if the distance between the current and the previous alignment is larger than max distance
				int currrentAlignmentPos=alignment.Position;
				int pastAlignmentPos= eventReads[alignment.MateRefID].back().Position;
				int distance= currrentAlignmentPos - pastAlignmentPos;

				//If the distance between the two reads is less than the maximum allowed distace, add it to the other reads of the event
				if(distance <= 1000){
					//add the read to the current window
					eventReads[alignment.MateRefID].push(alignment);
				}else{
					if(eventReads[alignment.MateRefID].size() >= minimumPairs){
						computeVariations(alignment.MateRefID);
					}
					//Thereafter the alignments are removed
					eventReads[alignment.MateRefID]=queue<BamAlignment>();
					linksFromWin[alignment.MateRefID]=queue<int>();										
					eventSplitReads[alignment.MateRefID]=vector<BamAlignment>();

					//the current alignment is inserted
					eventReads[alignment.MateRefID].push(alignment);
				}
			}
			for(int i=0;i< eventReads.size();i++){
				linksFromWin[i].push(alignment.Position);
			}

		}

	}

	chr=alignment.RefID;
}



float Window::computeCoverageB(int chrB, int start, int end, int32_t secondWindowLength){
	int bins=0;
	float coverageB=0;
	int element=0;
	unsigned int pos =floor(double(start)/500.0)*500;
	bool on = false;

	unsigned int nextpos =pos+500;
	string line;
	string coverageFile=this ->outputFileHeader+".tab";
	ifstream inputFile( coverageFile.c_str() );
	while(nextpos >= start and pos <= end and pos/500 < binnedCoverage[chrB].size() ){
			bins++;
			element=pos/500;
			coverageB+=double(binnedCoverage[chrB][element]-coverageB)/double(bins);
			pos+=500;
			nextpos=pos+500;
	}
	return(coverageB);
	
}

//This function accepts a queue with alignments that links to chrB from chrA, the function returns the number of reads that are positioned inside the region on A and have mates that are linking anywhere on B
vector<int> Window::findLinksToChr2(queue<BamAlignment> ReadQueue,long startChrA,long stopChrA,long startChrB,long endChrB, int pairsFormingLink){
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

//Computes statstics such as coverage on the window of chromosome A
vector<double> Window::computeStatisticsA(string bamFileName, int chrB, int start, int end, int32_t WindowLength, string indexFile){
	queue<int> linksFromA=linksFromWin[chrB];
	vector<double> statisticsVector;
	int currentReadPosition=0;
	int linksFromWindow=0;
	while ( linksFromA.size() > 0 ) {
		currentReadPosition=linksFromA.front();
		linksFromA.pop();
		if(start <= currentReadPosition and end >= currentReadPosition){
				linksFromWindow+=1;
			}
	}
	//calculates the coverage and returns the coverage within the window
	double coverageA=computeCoverageB(this -> chr,start,end,WindowLength);
	//todo change the 100 to the read length
	statisticsVector.push_back(coverageA);statisticsVector.push_back(linksFromWindow);statisticsVector.push_back(100);
	return(statisticsVector);
	
}



//Compute statistics of the variations and print them to file
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
		if(pairsFormingLink >= minimumPairs) {

		      //resize region so that it is just large enough to cover the reads that have mates in the present cluster
		      vector<long> chrALimit=newChrALimit(this-> eventReads[chr2],startSecondWindow,stopSecondWindow);					

			long startchrA=chrALimit[0];
			long stopchrA=chrALimit[1];




			vector<int> statsOnB=findLinksToChr2(eventReads[chr2],startchrA, stopchrA,startSecondWindow,stopSecondWindow,pairsFormingLink);
			int numLinksToChr2=statsOnB[0];
			int estimatedDistance=statsOnB[1];
            if(pairsFormingLink/double(numLinksToChr2) > 0.15){
	      // event accepted on discordant pairs. any split reads?
	      int splitsFormingLink = 0;
	      if (this->eventSplitReads[chr2].size() > 0) {
		  	for (vector<BamAlignment>::iterator alnit = eventSplitReads[chr2].begin(); alnit != eventSplitReads[chr2].end(); ++alnit) {
		  	// parse split read to get the other segment position, akin to a mate.
		  	string SA;
		  	(alnit)->GetTag("SA",SA);
		  	/* From the VCF documentation: Other canonical
		  	// alignments in a chimeric alignment, formatted as
		  	// a semicolon-delimited list:
		  	// (rname,pos,strand,CIGAR,mapQ,NM;)+.  Each element
		  	// in the list represents a part of the chimeric
		  	// alignment.  Conventionally, at a supplementary
		  	// line, the first element points to the primary
		  	// line.
		  	*/

		  	//Old boost code
		  	//vector <string> SA_elements;
		  	//boost::split(SA_elements,SA,boost::is_any_of(",;"));
		  	// 6 elements per SA entry, but there will be an
		  	// "extra" empty string element at the end of the
		  	// list, corresponding to the split of the final ;
		  	

		  
			//First split at ;, keep only the 0th elemnt
			vector <string> SA_elements;
			stringstream ss(SA);
	    	std::string item;
			while (std::getline(ss, item, ';')) {
				stringstream SS(item);
				string SA_data;
				while (std::getline(SS, SA_data, ',')) {
					SA_elements.push_back(SA_data);
				}
			}	
		
			for(int j=0; j< SA_elements.size(); j+=6) {
				string contig = SA_elements[j];
				string chr2name = this -> position2contig[chr2];
			    long pos = atol(SA_elements[j+1].c_str());
				if (contig == chr2name && pos >= startSecondWindow && pos <= stopSecondWindow) {
				      ++splitsFormingLink;
			    }
			  }
			}		
		}

	      
	      VCFLine(chr2,startSecondWindow,stopSecondWindow,startchrA,stopchrA,pairsFormingLink,splitsFormingLink,numLinksToChr2,estimatedDistance);
			    	found = true;
			}
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
	for(int i=1; i< QueueSize;i++){

		//if the current mate is close enough to the previous mate
		if( maxDistance >= (mate_positions[i]-linksToRegionQueue.back()) ){
			linksToRegionQueue.push(mate_positions[i]);

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
			linksToRegionQueue.push(mate_positions[i]);
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
	int average =0;
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
vector<double> Window::computeOrientation(queue<BamAlignment> alignmentQueue,long Astart,long Aend,long Bstart,long Bend){
	vector<double> orientationVector;
	int numberOfReads=alignmentQueue.size();
	int numberOfReverseReads=0;
	int numberOfReverseMates=0;
	int readsInRegion=0;
	double average=0;
	int n=0;

	for(int i=0;i<numberOfReads;i++){
		//Checks if the current read is located inside the region
		if(alignmentQueue.front().Position >=Astart and alignmentQueue.front().Position <= Aend){
			if(alignmentQueue.front().MatePosition >=Bstart and alignmentQueue.front().MatePosition <= Bend){
				//check the direction of the read on chrA
				average+=alignmentQueue.front().MatePosition-alignmentQueue.front().Position+1;
				n++;
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
	double FractionReverseMates=(double)numberOfReverseMates/(double)readsInRegion;

	orientationVector.push_back(FractionReverseReads);orientationVector.push_back(FractionReverseMates);orientationVector.push_back(average/n);
	return(orientationVector);
}

//function that prints the statistics to a vcf file
void Window::VCFLine(int chr2,int startSecondWindow, int stopSecondWindow,int startchrA,int stopchrA,int pairsFormingLink, int splitsFormingLink, int numLinksToChr2,int estimatedDistance){
	double coverageRealSecondWindow=computeCoverageB(chr2, startSecondWindow, stopSecondWindow, (stopSecondWindow-startSecondWindow+1) );
	vector<double> statisticsFirstWindow =computeStatisticsA(bamFileName, chr2, startchrA, stopchrA, (stopchrA-startchrA), indexFile);
	double coverageRealFirstWindow	= statisticsFirstWindow[0];
	int linksFromWindow=int(statisticsFirstWindow[1]);
	int averageReadLength=int(statisticsFirstWindow[2]);

	//calculate the expected number of reads
	int secondWindowLength=(stopSecondWindow-startSecondWindow+1);
	int firstWindowLength=stopchrA-startchrA+1;

	if(estimatedDistance > firstWindowLength) {
		estimatedDistance = estimatedDistance - firstWindowLength;
	}else{
		estimatedDistance = 1;
	}
	float expectedLinksInWindow = ExpectedLinks(firstWindowLength, secondWindowLength, estimatedDistance, mean_insert, std_insert, coverageRealFirstWindow, this -> readLength);
	vector<double> orientation=computeOrientation(eventReads[chr2],startchrA ,stopchrA,startSecondWindow,stopSecondWindow);
	string read1_orientation;
	string read2_orientation;
	string readOrientation="";
	ostringstream convertRead;
	//calculate the orientation of the reads
	if(orientation[0] > 0.5 ){
		convertRead << round(orientation[0]*100);
		read1_orientation="-(" + convertRead.str() + "%)";
	}else{
		convertRead << 100-round(orientation[0]*100);
		read1_orientation="+(" + convertRead.str() + "%)";
	}

		string mateOrientation="";
		ostringstream convertMate;


	if(orientation[1] > 0.5 ){
		convertMate << round(orientation[1]*100);
		read2_orientation="-(" + convertMate.str() + "%)";
	}else{
		convertMate << 100-round(orientation[1]*100);
		read2_orientation="+(" + convertMate.str() + "%)";
	}
	string filter;
	filter=filterFunction(pairsFormingLink/expectedLinksInWindow,numLinksToChr2,linksFromWindow,mean_insert,std_insert,estimatedDistance,coverageRealFirstWindow,coverageRealSecondWindow,meanCoverage);	
	if(this->chr == chr2) {

		vector<string> svVector=classification(chr2,startchrA, stopchrA,coverageRealFirstWindow,startSecondWindow, stopSecondWindow,coverageRealSecondWindow,this -> mean_insert,this -> std_insert,this -> outtie,orientation);
		string svType=svVector[0];
		string start=svVector[1];
		string end = svVector[2];
		string GT = svVector[3];
		string CN = svVector[4];
		this -> numberOfEvents++;

		//TODO generate the info field as a string, instead of printing the separate variable directly to the file
		if(svType != "INS" and svType !="BND"){
			intraChrVariationsVCF << this -> position2contig[this -> chr]  << "\t" <<     start   << "\tSV_" << this -> numberOfEvents << "_1" <<  "\t"  ;
			intraChrVariationsVCF << "N"       << "\t"	<< "<" << svType << ">";
			intraChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
			intraChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";END="<< end <<";LFW=" << linksFromWindow;
			intraChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
			intraChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
			intraChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow ;
			intraChrVariationsVCF << "\tGT:CN:PE:SR\t" << GT << ":" << CN << ":" <<pairsFormingLink << ":" << splitsFormingLink << "\n";
		}else{
				intraChrVariationsVCF << position2contig[this -> chr]  << "\t" <<     stopchrA   << "\tSV_" << this -> numberOfEvents << "_1" << "\t";
				intraChrVariationsVCF << "N"       << "\t"	<< "N[" << position2contig[chr2] << ":" << startSecondWindow << "[";
				intraChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
				intraChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
				intraChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
				intraChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
				intraChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow;
				intraChrVariationsVCF << "\tGT:PE:SR\t" << "./1" << ":" << pairsFormingLink << ":" << splitsFormingLink << "\n";

				//print the second breakend
				intraChrVariationsVCF <<   position2contig[chr2]<< "\t" <<    startSecondWindow    << "\tSV_" << this -> numberOfEvents <<  "_2" << "\t";
				intraChrVariationsVCF << "N"       << "\t"	<< "N]" << position2contig[this -> chr]  << ":" << stopchrA << "]";
				intraChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
				intraChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
				intraChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
				intraChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
				intraChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow;
				intraChrVariationsVCF << "\tGT:PE:SR\t" << "./1" << ":" << pairsFormingLink << ":" << splitsFormingLink << "\n";

		}
		if(svType == "IDUP"){
				intraChrVariationsVCF << position2contig[this -> chr]  << "\t" <<     stopchrA   << "\tSV_" << this -> numberOfEvents << "_2" << "\t";
				intraChrVariationsVCF << "N"       << "\t"	<< "N[" << position2contig[chr2] << ":" << startSecondWindow << "[";
				intraChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
				intraChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
				intraChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
				intraChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
				intraChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow;
				intraChrVariationsVCF << "\tGT:PE:SR\t" << "./1" << ":" << pairsFormingLink << ":" << splitsFormingLink << "\n";

				//print the second breakend
				intraChrVariationsVCF <<   position2contig[chr2]<< "\t" <<    startSecondWindow    << "\tSV_" << this -> numberOfEvents <<  "_3" << "\t";
				intraChrVariationsVCF << "N"       << "\t"	<< "N]" << position2contig[this -> chr]  << ":" << stopchrA << "]";
				intraChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
				intraChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
				intraChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
				intraChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
				intraChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow;
				intraChrVariationsVCF << "\tGT:PE:SR\t" << "./1" << ":" << pairsFormingLink << ":" << splitsFormingLink << "\n";
		}


	} else {
		string svType="BND";
		this -> numberOfEvents++;
		//print the first breakend
		interChrVariationsVCF << position2contig[this -> chr]  << "\t" <<     stopchrA   << "\tSV_" << this -> numberOfEvents << "_1" << "\t";
		interChrVariationsVCF << "N"       << "\t"	<< "N[" << position2contig[chr2] << ":" << startSecondWindow << "[";
		interChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
		interChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
		interChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
		interChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
		interChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow;
		interChrVariationsVCF << "\tGT:PE:SR\t" << "./1" << ":" << pairsFormingLink << ":" << splitsFormingLink << "\n";

		//print the second breakend
		interChrVariationsVCF <<   position2contig[chr2]<< "\t" <<    startSecondWindow    << "\tSV_" << this -> numberOfEvents <<  "_2" << "\t";
		interChrVariationsVCF << "N"       << "\t"	<< "N]" << position2contig[this -> chr]  << ":" << stopchrA << "]";
		interChrVariationsVCF << "\t.\t"  << filter << "\tSVTYPE="+svType <<";CHRA="<<position2contig[this->chr]<<";WINA=" << startchrA << "," <<  stopchrA;
		interChrVariationsVCF <<";CHRB="<< position2contig[chr2] <<";WINB=" <<  startSecondWindow << "," << stopSecondWindow << ";LFW=" << linksFromWindow;
		interChrVariationsVCF << ";LCB=" << numLinksToChr2 << ";LTE=" << pairsFormingLink << ";COVA=" << coverageRealFirstWindow;
		interChrVariationsVCF << ";COVB=" << coverageRealSecondWindow << ";OA=" << read1_orientation << ";OB=" << read2_orientation;
		interChrVariationsVCF << ";EL=" << expectedLinksInWindow << ";RATIO="<< pairsFormingLink/expectedLinksInWindow;
		interChrVariationsVCF << "\tGT:PE:SR\t" << "./1" << ":" << pairsFormingLink << ":" << splitsFormingLink << "\n";
	}
	return;

}
