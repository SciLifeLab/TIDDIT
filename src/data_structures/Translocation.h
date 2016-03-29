/*
 * Translocations.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi, Eisfeldt
 */

#ifndef TRANSLOCATION_H_
#define TRANSLOCATION_H_

#include "common.h"


class Window {
public:
	//object fields that change as the window moves
	float coverage;
	int chr;

	vector< queue<BamAlignment> >	eventReads;
	
	vector< vector<BamAlignment> > eventSplitReads;

	vector<long> covOnChrA;
	vector<long> tmpCovOnChrA;

	vector< queue<int> > linksFromWin;


	//More static parts initialized by constructor
	int max_insert;
	int min_insert;
	uint16_t minimum_mapping_quality;
	bool outtie;
	float mean_insert;
	float std_insert;
	int minimumPairs;
	float meanCoverage;
	int ploidy;
	int readLength;
	//the file name of the bamfile
	string bamFileName;
	string indexFile;


	//More static parts initialized by init method
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;
	vector< vector<float> > binnedCoverage;
	//Output part
	string outputFileHeader;
	ofstream interChrVariationsVCF;
	ofstream intraChrVariationsVCF;

        //A counter used to keep track of the number of events
        long numberOfEvents;


	Window(int max_insert, int min_insert, uint16_t minimum_mapping_quality,
			bool outtie, float mean_insert, float std_insert, int minimumPairs,
			float meanCoverage, string outputFileHeader, string bamFileName, string indexFile,int ploidity, int readLength); // constructor
	void initTrans(SamHeader head);				   // initialise the contig to position array
	void insertRead(BamAlignment alignment);	   // inserts a new read
	queue<BamAlignment> queueAppend(queue<BamAlignment> queueOne,queue<BamAlignment> queueTwo); //append queues;
	vector<long> findRegionOnB( queue<BamAlignment> alignmentQueue, int minimumPairs,int maxDistance); //Finds the region of the event on chromosome B
	vector<long> newChrALimit(queue<BamAlignment> alignmentQueue,long Bstart,long Bend); //resizes the window on CHRA
	vector<double> computeStatisticsA(string bamFileName, int chrB, int start, int end, int32_t WindowLength, string indexFile); //compute coverage and number of links from window on the chrA
	vector<double> computeOrientation(queue<BamAlignment> alignmentQueue,long Astart,long Aend,long Bstart,long Bend);//compute the orientation of the read and the mate
	vector<string> classification(int chr, int startA,int endA,double covA,int startB,int endB,double covB,int meanInsert,int STDInsert,bool outtie,vector<double> isReverse);
	string VCFHeader();
	void VCFLine(int chr2,int startSecondWindow, int stopSecondWindow,int startchrA,int stopchrA,int pairsFormingLink,int splitsFormingLink,int numLinksToChr2,int estimatedDistance);
	vector<int> findLinksToChr2(queue<BamAlignment> ReadQueue,long startChrA,long stopChrA,long startChrB,long endChrB, int pairsFormingLink);
	float computeCoverageB(int chrB, int start, int end, int32_t secondWindowLength); //computes the coverage of the window of chromosome B
	bool computeVariations(int chr2);

};

#endif /* TRANSLOCATIONS_H_ */
