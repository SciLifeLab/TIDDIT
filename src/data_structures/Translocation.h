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
    vector< queue<BamAlignment> > readsRegionA;
	//queue<int> readsRegionA;

	vector<long> covOnChrA;
	vector<long> tmpCovOnChrA;

	vector<int> linksFromWin;
	vector<int> tmpLinksFromWin;



	//More static parts initialized by constructor
	int max_insert;
	uint16_t minimum_mapping_quality;
	bool outtie;
	float mean_insert;
	float std_insert;
	int minimumPairs;
	float meanCoverage;
    int ploidity;

	//the file name of the bamfile
	string bamFileName;
	string indexFile;


	//More static parts initialized by init method
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;

	//Output part
	string outputFileHeader;
	ofstream interChrVariationsVCF;
	ofstream intraChrVariationsVCF;

        //A counter used to keep track of the number of events
        long numberOfEvents;


	Window(int max_insert, uint16_t minimum_mapping_quality,
			bool outtie, float mean_insert, float std_insert, int minimumPairs,
			float meanCoverage, string outputFileHeader, string bamFileName, string indexFile,int ploidity); // constructor
	void initTrans(SamHeader head);				   // initialise the contig to position array
	void insertRead(BamAlignment alignment);	   // inserts a new read
	queue<BamAlignment> queueAppend(queue<BamAlignment> queueOne,queue<BamAlignment> queueTwo); //append queues;
	vector<long> findRegionOnB( queue<BamAlignment> alignmentQueue, int minimumPairs,int maxDistance); //Finds the region of the event on chromosome B
	vector<long> newChrALimit(queue<BamAlignment> alignmentQueue,long Bstart,long Bend); //resizes the window on CHRA
	vector<double> computeStatisticsA(queue<BamAlignment> readsRegionA, int chrB, int start, int end, int32_t WindowLength); //compute coverage and number of links from window on the chrA
	vector<double> computeOrientation(queue<BamAlignment> alignmentQueue,long Astart,long Aend,long Bstart,long Bend);//compute the orientation of the read and the mate
	vector<string> classification(int chr, int startA,int endA,int covA,int startB,int endB,int covB,int meanInsert,int STDInsert,bool outtie,vector<double> isReverse);

	float computeCoverageB(int chrB, int start, int end, int32_t secondWindowLength); //computes the coverage of the window of chromosome B
	bool computeVariations(int chr2);

};

#endif /* TRANSLOCATIONS_H_ */
