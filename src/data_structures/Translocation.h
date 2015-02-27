/*
 * Translocations.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi
 */

#ifndef TRANSLOCATION_H_
#define TRANSLOCATION_H_

#include "common.h"
#include <list>





class Window {
public:
	//object fields that change as the window moves
	float coverage;
	int currentWindowStart;
	int currentWindowEnd;
	bool windowOpen;
	int chr;
	list<BamAlignment> TranslocationEvents;
	list<BamAlignment> alignmentsOnWindow;

	//More static parts initialized by constructor
	int windowSize;
	int windowStep;
	int max_insert;
	uint16_t minimum_mapping_quality;
	bool outtie;
	float mean_insert;
	float std_insert;
	int minimumPairs;
	float meanCoverage;


	//More static parts initialized by init method
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;

	//Output part
	string outputFileHeader;
	ofstream interChrVariations;
	ofstream intraChrVariations;


	Window(int windowSize, int windowStep, int max_insert, uint16_t minimum_mapping_quality,
			bool outtie, float mean_insert, float std_insert, int minimumPairs,
			float meanCoverage, string outputFileHeader); // constructor
	void initTrans(SamHeader head);				   // initialise the contig to position array
	void insertRead(BamAlignment alignment);	   // inserts a new read
	void goToNextWindow(int position);			   // moves to next window
	void resetWindow(int position, uint32_t chr);  // resets window when moving to next chr

	float computeCoverage(); // computes coverage of the area memorised in the area
	float computeCoverage(uint32_t start, uint32_t end); // computes coverage of the area memorised in the area

	bool computeIntraChr(ofstream & OutputFileDescriptor, uint16_t minimum_mapping_quality, float mean_insert, float std_insert, int minimumPairs, float meanCoverage);
	bool computeInterChr(ofstream & OutputFileDescriptor, uint16_t minimum_mapping_quality, float mean_insert, float std_insert, int minimumPairs, float meanCoverage);
	bool computeVariations();

// the two must return a bool and I need to adjust the window accordingly
// it is likely that a double passage must be performed or I need to doble the data structure

};

class Link {
public:
	uint32_t chr1_start;
	uint32_t chr2_start;
	uint32_t chr2_end;
	uint32_t supportingPairs;
	bool	ischr1_rev;
	bool	ischr2_rev;
};

class Translocations {
public:
	map<uint32_t, vector<Link> > Connections ;
	Translocations();
	void insertConnection(uint32_t chr2, uint32_t pos2  );
	void insertConnection(uint32_t chr1_start, uint32_t chr2, uint32_t pos2 ,  bool reversed1, bool reversed2 );


//	void findEvents(ofstream & OutputFileDescriptor, uint32_t chr1, uint32_t chr2, uint32_t minimumPairs, float minCov, float maxCov, uint32_t windowSize, uint32_t windowStep);

};



#endif /* TRANSLOCATIONS_H_ */
