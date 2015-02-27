/*
 * Contig.h
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#ifndef CONTIG_H_
#define CONTIG_H_



//#include "common.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "common.h"
#include "Features.h"
//using namespace BamTools;


class ContigsFeat{
public:
	unsigned int feature;
	vector<pair<unsigned int, unsigned int> > startStopPositions;
	ContigsFeat();

	~ContigsFeat();

};


enum data {readCov, insertCov, cmCov, woCov, singCov, mdcCov};

class Position {
public:
	unsigned int ReadCoverage;
	unsigned int StratingInserts;
	unsigned int InsertCoverage;
	unsigned int CorrectlyMated;
	unsigned int WronglyOriented;
	unsigned int Singleton;
	unsigned int MatedDifferentContig;
	unsigned long int insertsLength;

	Position();

	~Position();

};


#define MIN(x,y) \
  ((x) < (y)) ? (x) : (y)




class Contig{
	unsigned int contigLength;

	float lowCoverageFeat;
	float highCoverageFeat;
	float lowNormalFeat;
	float highNormalFeat;
	float highSingleFeat;
	float highSpanningFeat;
	float highOutieFeat;



	float MINUM_COV;

	void updateCov(unsigned int strat, unsigned int end, data type);

public:
	Position *CONTIG;

	Contig();
	Contig(unsigned int contigLength);
	~Contig();

	void updateContig(BamAlignment b, int max_insert,  bool is_mp); // given an alignment it updates the contig situation

	float getCoverage();

	unsigned int getContigLength();

	unsigned int getLowCoverageAreas(float C_A, unsigned int windowSize, unsigned int windowStep);
	unsigned int getHighCoverageAreas(float C_A, unsigned int windowSize, unsigned int windowStep);
	unsigned int getLowNormalAreas(float C_M, unsigned int windowSize, unsigned int windowStep);
	unsigned int getHighNormalAreas(float C_M, unsigned int windowSize, unsigned int windowStep);
	unsigned int getHighSingleAreas(unsigned int windowSize, unsigned int windowStep, float C_A);
	unsigned int getHighSpanningAreas(unsigned int windowSize, unsigned int windowStep, float C_A);
	unsigned int getHighOutieAreas(unsigned int windowSize, unsigned int windowStep, float C_A);
	unsigned int getCompressionAreas(float insertionMean, float insertionStd, float Zscore, unsigned int windowSize, unsigned int windowStep);
	unsigned int getExpansionAreas(float insertionMean, float insertionStd, float Zscore, unsigned int windowSize, unsigned int windowStep);

	void print();



	vector<pair<unsigned int, unsigned int> > lowCoverageAreas;
	vector<pair<unsigned int, unsigned int> > highCoverageAreas;
	vector<pair<unsigned int, unsigned int> > lowNormalAreas;
	vector<pair<unsigned int, unsigned int> > highNormalAreas;
	vector<pair<unsigned int, unsigned int> > highSingleAreas;
	vector<pair<unsigned int, unsigned int> > highSpanningAreas;
	vector<pair<unsigned int, unsigned int> > highOutieAreas;
	vector<pair<unsigned int, unsigned int> > compressionAreas;
	vector<pair<unsigned int, unsigned int> > expansionAreas;

};





#endif /* CONTIG_H_ */
