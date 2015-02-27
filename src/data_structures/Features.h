/*
 * Features.h
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#ifndef FEATURES_H_
#define FEATURES_H_

//#include "common.h"
#include <string>
using namespace std;


struct ternary{
	string feature;
	unsigned int start;
	unsigned int end;
};




class Features {
	unsigned int LOW_COVERAGE_AREA;
	unsigned int HIGH_COVERAGE_AREA;
	unsigned int LOW_NORMAL_AREA;
	unsigned int HIGH_NORMAL_AREA;
	unsigned int HIGH_SINGLE_AREA;
	unsigned int HIGH_SPANNING_AREA;
	unsigned int HIGH_OUTIE_AREA;
	unsigned int COMPRESSION_AREA;
	unsigned int STRECH_AREA;

public:
	Features();
	~Features();

	void setLOW_COVERAGE_AREA(unsigned int numFeat);
	void setHIGH_COVERAGE_AREA(unsigned int numFeat);
	void setLOW_NORMAL_AREA(unsigned int numFeat);
	void setHIGH_NORMAL_AREA(unsigned int numFeat);
	void setHIGH_SINGLE_AREA(unsigned int numFeat);
	void setHIGH_SPANNING_AREA(unsigned int numFeat);
	void setHIGH_OUTIE_AREA(unsigned int numFeat);
	void setCOMPRESSION_AREA(unsigned int numFeat);
	void setSTRECH_AREA(unsigned int numFeat);

	void updateLOW_COVERAGE_AREA(unsigned int numFeat);
	void updateHIGH_COVERAGE_AREA(unsigned int numFeat);
	void updateLOW_NORMAL_AREA(unsigned int numFeat);
	void updateHIGH_NORMAL_AREA(unsigned int numFeat);
	void updateHIGH_SINGLE_AREA(unsigned int numFeat);
	void updateHIGH_SPANNING_AREA(unsigned int numFeat);
	void updateHIGH_OUTIE_AREA(unsigned int numFeat);
	void updateCOMPRESSION_AREA(unsigned int numFeat);
	void updateSTRECH_AREA(unsigned int numFeat);

	unsigned int getLOW_COVERAGE_AREA();
	unsigned int getHIGH_COVERAGE_AREA();
	unsigned int getLOW_NORMAL_AREA();
	unsigned int getHIGH_NORMAL_AREA();
	unsigned int getHIGH_SINGLE_AREA();
	unsigned int getHIGH_SPANNING_AREA();
	unsigned int getHIGH_OUTIE_AREA();
	unsigned int getCOMPRESSION_AREA();
	unsigned int getSTRECH_AREA();

	unsigned int returnTotal();

	unsigned int returnLOW_COV();
	unsigned int returnHIGH_COV();
	unsigned int returnLOW_NORM_COV();
	unsigned int returnHIGH_NORM_COV();
	unsigned int returnHIGH_SINGLE();
	unsigned int returnHIGH_OUTIE();
	unsigned int returnHIGH_SPAN();
	unsigned int returnCOMPR();
	unsigned int returnSTRECH();



};



#endif /* FEATURES_H_ */
