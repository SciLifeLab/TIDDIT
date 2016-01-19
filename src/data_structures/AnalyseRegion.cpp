#include "ProgramModules.h"
#include "data_structures/Translocation.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
//The region module, searches for variations given by a region file, and outputs the most similar variation/variations as a VCF.

Region::Region(){}
//main function
void Region::region(string bamFile,string baiFile,string regionFile,string output,map<string,unsigned int> contig2position, map<unsigned int,string> position2contig,int ploidy,int minimum_mapping_quality,uint64_t genomeLength,uint32_t contigsNumber){
	this -> contig2position = contig2position;
	this -> position2contig= position2contig;
	this -> chrPrefix="";

	string testPrefix=position2contig[1];
	if(testPrefix.find("chr")  != string::npos){
		this -> chrPrefix="chr";
	}else if(testPrefix.find("Chr")  != string::npos){
		this -> chrPrefix="Chr";
	}else if(testPrefix.find("CHR")  != string::npos){
		this -> chrPrefix="CHR";
	}
	//open the region file if it is available, otherwise, crash
	ifstream inputFile( regionFile.c_str() );
	cout << regionFile << endl;
	if(inputFile){
		cout << "region file found, calculating bam file statistics..." << endl;
	}else{
		cout << "Error, unnable to load the region file" << endl;
		return;
	}

	bool outtie = true;
	int min_insert = 100;
	int max_insert = 200000;
	
	LibraryStatistics library;
	library = computeLibraryStats(bamFile, genomeLength, max_insert,min_insert, outtie,minimum_mapping_quality);
	double coverage   = library.C_A;
	int meanInsert = library.insertMean;
	int insertStd  = library.insertStd;
	max_insert =meanInsert+2*insertStd;
	this -> max_insert = max_insert;
	outtie=library.mp;

	Cov *calculateCoverage;
	calculateCoverage = new Cov();
	calculateCoverage -> coverageMain(bamFile,baiFile,"",output,coverage,contig2position,4,500);

	//open the bam file
	BamReader alignmentFile;
	alignmentFile.Open(bamFile);
	//Information from the header is needed to initialize the data structure
	SamHeader head = alignmentFile.GetHeader();

	Window *window;

	window = new Window(max_insert,min_insert, minimum_mapping_quality,outtie,  meanInsert,  insertStd,  0,coverage,  output,bamFile,baiFile,ploidy);
	
	window->initTrans(head);
	//expands a vector so that it is large enough to hold reads from each contig in separate elements
	window->eventReads.resize(contigsNumber);
	window-> binnedCoverage.resize(contigsNumber);
	window-> linksFromWin.resize(contigsNumber);

	window -> numberOfEvents=1;

	string line;
	string coverageFile=output+".tab";
	ifstream coverageInputFile( coverageFile.c_str() );
	while (getline( coverageInputFile, line )){
		vector<std::string> splitline;
		boost::split(splitline, line, boost::is_any_of("\t"));
		window -> binnedCoverage[window -> contig2position[splitline[0]]].push_back(atof(splitline[1].c_str()));
	}
	coverageInputFile.close();
	//read the region file
	if (inputFile){
		while (getline( inputFile, line )){
			int pairsFormingLink=0;
			vector<std::string> splitline;
			boost::split(splitline, line, boost::is_any_of("\t"));
			queue<int> chr;queue<int> startPos;queue<int> endPos;
			vector< queue<int> > regionVector;
			if(this -> input == "bed"){
				regionVector=readRegionFile(splitline);

			}else{
				if(splitline[0].find("#")  != string::npos)
					continue;
				else{
					regionVector=readVCF(splitline);
				}
			}
			chr=regionVector[0];
			startPos=regionVector[1];
			endPos=regionVector[2];

			BamAlignment currentRead;
			

			//test if an index file is available
			if(alignmentFile.OpenIndex(baiFile) == 0){
				cout << "Failed to load the index file" << endl;
				return;
			}

			//Define the chromosomes
			
			int chrA=chr.front();int startA=startPos.front(); int endA=endPos.front();
			chr.pop();startPos.pop();endPos.pop();
			int chrB=chr.front();int startB=startPos.front(); int endB=endPos.front();

			//sort the events
			if(chrB < chrA){
				int tmpChr;int tmpStart;int tmpEnd;
				tmpChr=chrB;tmpStart=startB;tmpEnd=endB;
				chrB=chrA;startB=startA;endB=endA;
				chrA=tmpChr;startA=tmpStart;endA=tmpEnd;
			}else if(chrA == chrB and startB < startA){
				int tmpChr;int tmpStart;int tmpEnd;
				tmpChr=chrB;tmpStart=startB;tmpEnd=endB;
				chrB=chrA;startB=startA;endB=endA;
				chrA=tmpChr;startA=tmpStart;endA=tmpEnd;

			}
			//read through the file
			window ->chr =chrA;
			bool first=true;
			int newChrBstart=endB;
			int newChrBend=startB;
			int newChrAEnd=endA;
			alignmentFile.SetRegion(chrA,startA,chrA,endA+1); 
			while ( alignmentFile.GetNextAlignmentCore(currentRead) ) {
				//makes sure that we are inside the region
				if(startA <= currentRead.Position and endA >= currentRead.Position){
					readStatus alignmentStatus = computeReadType(currentRead, max_insert,min_insert, outtie);
					if(alignmentStatus != unmapped and alignmentStatus != lowQualty ) {
						if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
							//store the read if its mate is reaching region B
							//store the read in the discordant read queue if it is discordant, but do not reach region B
							if(currentRead.MatePosition >= startB and currentRead.MatePosition <= endB and currentRead.MateRefID == chrB) {
								pairsFormingLink++;
								//set the first position of A to the first read in the A window
								if (first == true){
									first = false;
									startA=currentRead.Position;
								}
								//the end position is continously updated, the last read will become the end
								newChrAEnd=currentRead.Position;
								if(currentRead.MatePosition < newChrBstart){
									newChrBstart=currentRead.MatePosition;
								}
	
								if(currentRead.MatePosition >newChrBend){
									newChrBend=currentRead.MatePosition;
								}
	
							}
							//only add reads to the vector if we are inside the event window
							if(first == false){
								window->eventReads[currentRead.MateRefID].push(currentRead);
								window -> linksFromWin[chrB].push(currentRead.Position);	
							}
		         		}
					}
				}
			}
			if(newChrBend >= newChrBstart){
				endB=newChrBend;
				startB=newChrBstart;
			}
			endA=newChrAEnd;
			vector<int> statsOnB=window->findLinksToChr2(window->eventReads[chrB],startA, endA,startB,endB,pairsFormingLink);
			int numLinksToChr2=statsOnB[0];
			int estimatedDistance=statsOnB[1];
			if(estimatedDistance < 0){
				estimatedDistance=1;
			}
			window->VCFLine(chrB,startB,endB,startA,endA,pairsFormingLink,numLinksToChr2,estimatedDistance);
		}
		window->interChrVariationsVCF.close();
		window->intraChrVariationsVCF.close();

	}
}


vector< queue<int> > Region::readRegionFile(vector<string> splitline){
	queue<int> chr;queue<int> startPos;queue<int> endPos;
	if(splitline.size() == 6){
		for(int i =0;i<2;i++){	
			//extract the chromosome
			
			chr.push(contig2position[matchContig(splitline[i*3])]);

			//extracts the start position on the given chromosome
			startPos.push( atoi(splitline[i*3+1].c_str()));

			//extracts the end position on the given chromosome
			endPos.push(atoi(splitline[i*3+2].c_str()));
		}

	//the region is given by a chr,start,stop. This is common in CNV events
	}else if(splitline.size() == 3){
		for(int i =0;i<2;i++){	
					
			//extract the chromosome
			chr.push(contig2position[matchContig(splitline[0])]);

			//extracts the start position on the given chromosome
			int pos =atoi(splitline[(i+1)].c_str())-max_insert;
			if (pos > 0){
				startPos.push(pos);
			}else{
				startPos.push(1);
			}

			//extracts the end position on the given chromosome
			endPos.push(atoi(splitline[(i+1)].c_str())+max_insert);
		}		
	//the region is given by two position, each position is also accompanied by its chromosome, exact translocations are given this way
	}else if(splitline.size() == 4){
		for(int i =0;i<2;i++){	
			//extract the chromosome
			chr.push(contig2position[matchContig(splitline[i*2])]);

			//extracts the start position on the given chromosome
			int pos =atoi(splitline[(2*i+1)].c_str())-max_insert;
			if (pos > 0){
				startPos.push(pos);
			}else{
				startPos.push(1);
			}
			//extracts the end position on the given chromosome
			endPos.push(atoi(splitline[(2*i+1)].c_str())+max_insert);
		}	
	}
	vector< queue<int> > regionVector;	
	regionVector.push_back(chr);regionVector.push_back(startPos);regionVector.push_back(endPos);
	return(regionVector);
}

//VCF input
vector< queue<int> > Region::readVCF(vector<string> splitline){
	queue<int> chr;queue<int> startPos;queue<int> endPos;

	//for all kinds of variations, the first column is chrA and the second column is the start position of a variation
	chr.push( this -> contig2position[matchContig(splitline[0])] );
	int pos =atoi(splitline[1].c_str())-max_insert;
	if (pos > 0){
		startPos.push(pos);
	}else{
		startPos.push(1);
	}
	endPos.push(atoi(splitline[1].c_str())+max_insert);

	//if the variation is a simple intra chromosomal variation
	if(splitline[4].find("<")  != string::npos){
		chr.push( this ->contig2position[matchContig(splitline[0])] );
		vector<std::string> infoField;
		boost::split(infoField, splitline[7], boost::is_any_of(";"));
		//itterate through the info field, find entry with subfield END, split on = use element 1 as end
		for(int i =0; i < infoField.size();i++){
			if(infoField[i].find("END")  != string::npos){
				vector<std::string> endPosition;
				boost::split(endPosition, infoField[i], boost::is_any_of("="));

				pos =atoi(endPosition[1].c_str())-max_insert;
				if (pos > 0){
					startPos.push(pos);
				}else{
					startPos.push(1);
				}
				endPos.push(atoi(endPosition[1].c_str())+max_insert);

				break;
			}
		}
	//if the variation is given as a complex break end
	}else{
		vector<std::string> positionColumn;
		//split the alt field
		boost::split(positionColumn,splitline[4],boost::is_any_of("]["));
		for(int i =0; i < positionColumn.size();i++){
			if(positionColumn[i].find(":")  != string::npos){
				vector<std::string> altField;
				boost::split(altField, positionColumn[i], boost::is_any_of(":"));
				chr.push( this ->contig2position[matchContig(altField[0])] );
				pos =atoi(altField[1].c_str())-max_insert;
				if (pos > 0){
					startPos.push(pos);
				}else{
					startPos.push(1);
				}
				endPos.push(atoi(altField[1].c_str())+max_insert);

				break;
			}
		}
		
	}
	vector< queue<int> > regionVector;	
	regionVector.push_back(chr);regionVector.push_back(startPos);regionVector.push_back(endPos);
	return(regionVector);
}

//solves chromosome prefix differences in caller output and reference data 
string Region::matchContig(string inputChr){
	string chr;
	string inputPrefix="";
	this -> position2contig= position2contig;
	if(inputChr.find("chr")  != string::npos){
		inputPrefix="chr";
	}else if(inputChr.find("Chr")  != string::npos){
		inputPrefix="Chr";
	}else if(inputChr.find("CHR")  != string::npos){
		inputPrefix="CHR";
	}
	
	if(inputPrefix == this-> chrPrefix){
		chr=inputChr;
	}else if(inputPrefix==""){
		chr=this->chrPrefix+inputChr;
	}else{
		boost::replace_all(inputChr, inputPrefix, this -> chrPrefix);
		chr= inputChr;

	}

	return(chr);
}
