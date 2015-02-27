/*
 * Translocations.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi
 */

#ifndef TRANSLOCATIONS_H_
#define TRANSLOCATIONS_H_

#include "common.h"



class Link {
public:
	uint32_t chr1_start;
	uint32_t chr1_end;
	uint32_t chr2_start;
	uint32_t chr2_end;
	uint32_t supportingPairs;
};



class Translocations {
public:
	uint32_t chromosomesNum;
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;
	map<uint32_t, map<uint32_t, vector<Link> > > Connections ;

	Translocations(uint32_t numChromosomes);

	void initTrans(SamHeader head);
	void insertConnection(uint32_t chr1, uint32_t pos1, uint32_t chr2, uint32_t pos2  );

	void printConnections();

	void compressConnections(uint32_t chr1, uint32_t chr2);
	void findEvents(ofstream & OutputFileDescriptor, uint32_t chr1, uint32_t chr2, uint32_t minimumPairs, float minCov, float maxCov, uint32_t windowSize, uint32_t windowStep);

};


#endif /* TRANSLOCATIONS_H_ */
