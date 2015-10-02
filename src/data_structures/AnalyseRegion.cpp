#include "ProgramModules.h"
#include "data_structures/Translocation.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
//The region module, searches for variations given by a region file, and outputs the most similar variation/variations as a VCF.

Region::Region(){}
//main function
void Region::region(string bamFile,string baiFile,string regionFile,string output,map<string,unsigned int> contig2position, int ploidy,int minimum_mapping_quality,uint64_t genomeLength,uint32_t contigsNumber){
	//open the region file if it is available, otherwise, crash
	ifstream inputFile( regionFile.c_str() );
	cout << regionFile << endl;
	if(inputFile){
		cout << "region file found, calculating bam file statistics..." << endl;
	}else{
		cout << "Error, unnable to load the region file" << endl;
		return;
	}

	vector<int> libraryStats;
	//compute max_insert,mean insert and outtie
	autoSettings *autoStats;
	autoStats = new autoSettings();
	size_t start = time(NULL);
	libraryStats=autoStats->autoConfig(bamFile,minimum_mapping_quality);
			
	printf ("auto config time consumption= %lds\n", time(NULL) - start);
	int max_insert=libraryStats[0];
	int min_insert=libraryStats[1];
	bool outtie=libraryStats[2] != 0;

	LibraryStatistics library;
	start = time(NULL);
	library = computeLibraryStats(bamFile, genomeLength, max_insert,min_insert, outtie);
	printf ("library stats time consumption= %lds\n", time(NULL) - start);
	double coverage   = library.C_A;
	int meanInsert = library.insertMean;
	int insertStd  = library.insertStd;
	max_insert =meanInsert+4*insertStd;


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
	window -> numberOfEvents=1;

	//read the region file
	string line;
	if (inputFile){
		cout << "input file" << endl;
		while (getline( inputFile, line )){
			int pairsFormingLink=0;
			queue<int> chr;
			queue<int> startPos;queue<int> endPos;

			vector<std::string> splitline;
			boost::split(splitline, line, boost::is_any_of("\t"));
			//if the variation is given by two regions, such as a translocation or FIndTranslocations intrachromosomal event
			if(splitline.size() == 6){
				for(int i =0;i<2;i++){	
					//extract the chromosome
					chr.push(contig2position[splitline[i*3]]);

					//extracts the start position on the given chromosome
					startPos.push( atoi(splitline[i*3+1].c_str()));

					//extracts the end position on the given chromosome
					endPos.push(atoi(splitline[i*3+2].c_str()));
				}

			//the region is given by a chr,start,stop. This is common in CNV events
			}else if(splitline.size() == 3){
				for(int i =0;i<2;i++){	
					
					//extract the chromosome
					chr.push(contig2position[splitline[0]]);

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
					chr.push(contig2position[splitline[i*2]]);

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
			


			


			BamAlignment currentRead;
			

			//test if an index file is available
			if(alignmentFile.OpenIndex(baiFile) == 0){
				cout << "Failed to load the index file" << endl;
				return;
			}

			//Define the chromosomes
			alignmentFile.SetRegion(chr.front(),startPos.front(),chr.front(),endPos.front()+1);
			int chrA=chr.front();int startA=startPos.front(); int endA=endPos.front();
			chr.pop();startPos.pop();endPos.pop();
			int chrB=chr.front();int startB=startPos.front(); int endB=endPos.front();

			//read through the file
			while ( alignmentFile.GetNextAlignmentCore(currentRead) ) {
				
				readStatus alignmentStatus = computeReadType(currentRead, max_insert,min_insert, outtie);
				if(alignmentStatus != unmapped or alignmentStatus != lowQualty ) {
					if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance or alignmentStatus == pair_wrongOrientation) {
						//store the read if its mate is reaching region B
						//store the read in the discordant read queue if it is discordant, but do not reach region B
						if(currentRead.RefID < currentRead.MateRefID or (currentRead.RefID == currentRead.MateRefID and currentRead.Position < currentRead.MatePosition)) {
							if(currentRead.MatePosition >= startB and currentRead.MatePosition <= endB and currentRead.MateRefID == chrB) {
								pairsFormingLink++;
							}
							window->eventReads[currentRead.MateRefID].push(currentRead);
						}
		         	}
				}
			}
		
			//calculate statistics and print the vcf file
			vector<long> chrALimit=window->newChrALimit(window->eventReads[chrB],startB,endB);
			vector<int> statsOnB=window->findLinksToChr2(window->eventReads[chrB],startA, endA,startB,endB,pairsFormingLink);
			int numLinksToChr2=statsOnB[0];
			int estimatedDistance=statsOnB[1];
			window -> chr =chrA;
			if(estimatedDistance < 0){
				estimatedDistance=1;
			}
			window->VCFLine(chrB,startB,endB,startA,endA,pairsFormingLink,numLinksToChr2,estimatedDistance);

		}
		window->interChrVariationsVCF.close();
		window->intraChrVariationsVCF.close();

	}
}
