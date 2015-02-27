/*
 * FRC.h
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#ifndef FRC_H_
#define FRC_H_



#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

//#include "common.h"
//#include "Features.h"
#include "Contig.h"


class contigFeatures {
	string contigID;
	unsigned long int contigLength;

	unsigned int TOTAL;

public:


	Features PE;
	Features MP;

	contigFeatures();
	~contigFeatures();

	void setID(string ID);
	string getID();

	void setContigLength(unsigned int contigLength);
	unsigned long int getContigLength();


	unsigned int getTotal();


	unsigned int getLOW_COV_PE();
	unsigned int getHIGH_COV_PE();
	unsigned int getLOW_NORM_COV_PE();
	unsigned int getHIGH_NORM_COV_PE();
	unsigned int getHIGH_SINGLE_PE();
	unsigned int getHIGH_OUTIE_PE();
	unsigned int getHIGH_SPAN_PE();
	unsigned int getCOMPR_PE();
	unsigned int getSTRECH_PE();

	unsigned int getHIGH_SINGLE_MP();
	unsigned int getHIGH_OUTIE_MP();
	unsigned int getHIGH_SPAN_MP();
	unsigned int getCOMPR_MP();
	unsigned int getSTRECH_MP();

	vector<ternary> SUSPICIOUS_AREAS;

	void printFeatures(ofstream &file);
	void printFeaturesGFF3(ofstream &file);

};



class FRC {

    vector<contigFeatures> CONTIG;
    unsigned int contigs;

    float C_A; // total read coverage
    float S_A; // total span coverage

    float C_M; // coverage induced by correctly aligned pairs
    float C_W; // coverage induced by wrongly mated pairs
    float C_S; // coverage induced by singletons
    float C_D; // coverage induced by reads with mate on a diferent contif
    float Expansion;
    float Compression;

    float insertMean;
    float insertStd;

public:

	FRC();
	FRC(unsigned int contigs);
	~FRC();

	map<float, unsigned int> CEstatistics;

	unsigned int returnContigs();
	void setContigLength(unsigned int ctg, unsigned int contigLength);
	unsigned int getContigLength(unsigned int ctg);

	void setC_A(float C_A);
	void setS_A(float S_A);
	void setC_M(float C_M);
	void setC_W(float C_W);
	void setC_S(float C_S);
	void setC_D(float C_D);
	void setInsertMean(float insertMean);
	void setInsertStd(float insertStd);
	void setID(unsigned int i, string ID);
	string  getID(unsigned int i);

	void sortFRC();
	float obtainCoverage(unsigned int ctg, Contig *contig);

	void computeCEstats(Contig *contig, unsigned int WindowSize, unsigned int WindowStep, float mean, float std);

	void computeLowCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int WindowSize, unsigned int WindowStep);
	void computeHighCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeLowNormalArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighNormalArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighSingleArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighSpanningArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeHighOutieArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep);
	void computeCompressionArea(string type, unsigned int ctg, Contig *contig, float Zscore, unsigned int windowSize, unsigned int windowStep);
	void computeStrechArea(string type, unsigned int ctg, Contig *contig, float Zscore, unsigned int windowSize, unsigned int windowStep);

	unsigned int getTotal(unsigned int ctg);


	unsigned int getFeatures(FeatureTypes type, int contig);

	unsigned int getLOW_COV_PE(unsigned int ctg);
	unsigned int getHIGH_COV_PE(unsigned int ctg);
	unsigned int getLOW_NORM_COV_PE(unsigned int ctg);
	unsigned int getHIGH_NORM_COV_PE(unsigned int ctg);
	unsigned int getHIGH_SINGLE_PE(unsigned int ctg);
	unsigned int getHIGH_OUTIE_PE(unsigned int ctg);
	unsigned int getHIGH_SPAN_PE(unsigned int ctg);
	unsigned int getCOMPR_PE(unsigned int ctg);
	unsigned int getSTRECH_PE(unsigned int ctg);
	unsigned int getHIGH_SINGLE_MP(unsigned int ctg);
	unsigned int getHIGH_OUTIE_MP(unsigned int ctg);
	unsigned int getHIGH_SPAN_MP(unsigned int ctg);
	unsigned int getCOMPR_MP(unsigned int ctg);
	unsigned int getSTRECH_MP(unsigned int ctg);



	void printFeatures(unsigned int ctg, ofstream &f);

	void printFeaturesGFF3(unsigned int ctg, ofstream &f);


};




#endif /* FRC_H_ */
