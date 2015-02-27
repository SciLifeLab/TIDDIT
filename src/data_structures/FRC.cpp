/*
 * FRC.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#include "FRC.h"
#include <fstream>

bool sortContigs(contigFeatures i, contigFeatures j) {
	return (i.getContigLength() > j.getContigLength());
}

FRC::FRC() {

}

FRC::~FRC() {
	if (CONTIG.empty() != 0) {
		CONTIG.~vector();
	}
}

FRC::FRC(unsigned int contigs) {
	this->contigs = contigs;
	this->CONTIG.resize(contigs);
}


void FRC::sortFRC() {
	sort(CONTIG.begin(), CONTIG.end(), sortContigs);
}


unsigned int FRC::returnContigs() {
	return this->contigs;
}

void FRC::setContigLength(unsigned int ctg, unsigned int contigLength) {
	CONTIG[ctg].setContigLength(contigLength);
}

unsigned int FRC::getContigLength(unsigned int ctg) {
	return this->CONTIG[ctg].getContigLength();
}


float FRC::obtainCoverage(unsigned int ctg, Contig *contig) {


	float coverage = contig->getCoverage();
	return coverage;

}


void FRC::computeCEstats(Contig *contig, unsigned int windowSize, unsigned int windowStep, float insertionMean, float insertionStd ) {

	unsigned int contigLength = contig->getContigLength();
	unsigned long int spanningCoverage = 0; // total insert length
	unsigned int inserts = 0; // number of inserts
	float localMean;
	float Z_stats = 0;

	unsigned int minInsertNum = 5;
	if(contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < contigLength ; i++ ) {
			if(contig->CONTIG[i].StratingInserts > 0) {
				inserts += contig->CONTIG[i].StratingInserts;
				spanningCoverage += contig->CONTIG[i].insertsLength;
			}
			//spanningCoverage += contig->CONTIG[i].InsertCoverage;
		}
		if(inserts > minInsertNum) {
			localMean = spanningCoverage/(float)inserts;
			Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
			Z_stats = floorf(Z_stats * 10) / 10;
			if(this->CEstatistics.count(Z_stats) == 1) {
				this->CEstatistics[Z_stats]++;
			} else {
				this->CEstatistics[Z_stats] = 1;
			}
			//cout << Z_stats << "\n";

		}
	} else { //otherwise compute features on sliding window
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		for(unsigned int i=startWindow; i < endWindow; i++ ) {
			if(contig->CONTIG[i].StratingInserts > 0) {
				inserts += contig->CONTIG[i].StratingInserts;
				spanningCoverage += contig->CONTIG[i].insertsLength;
			}
			//spanningCoverage += contig->CONTIG[i].InsertCoverage;
		}
		if(inserts > minInsertNum) {
			localMean = spanningCoverage/(float)inserts;
			Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
			Z_stats = floorf(Z_stats * 10) / 10;
			if(this->CEstatistics.count(Z_stats) == 1) {
				this->CEstatistics[Z_stats]++;
			} else {
				this->CEstatistics[Z_stats] = 1;
			}
			//cout << Z_stats << "\n";

		}
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > contigLength) {
			endWindow = contigLength;
		}
		while(endWindow < contigLength) {
			inserts = 0; //reset window coverage
			spanningCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow; i++ ) {
				if(contig->CONTIG[i].StratingInserts > 0) {
					inserts += contig->CONTIG[i].StratingInserts;
					spanningCoverage += contig->CONTIG[i].insertsLength;
				}
				//spanningCoverage += contig->CONTIG[i].InsertCoverage;
			}
			if(inserts > minInsertNum) {
				localMean = spanningCoverage/(float)inserts;
				Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
				Z_stats = floorf(Z_stats * 10) / 10;
				if(this->CEstatistics.count(Z_stats) == 1) {
					this->CEstatistics[Z_stats]++;
				} else {
					this->CEstatistics[Z_stats] = 1;
				}
				//cout << Z_stats << "\n";
			}
			startWindow += windowStep;
			endWindow += windowStep;
			if(endWindow > contigLength) {
				endWindow = contigLength;
			}
		}
	}
}



void FRC::computeLowCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getLowCoverageAreas(C_A,windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateLOW_COVERAGE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateLOW_COVERAGE_AREA(feat);
	}
	for(unsigned int i=0; i< contig->lowCoverageAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "LOW_COV_"+type;
		tmp.start = contig->lowCoverageAreas.at(i).first;
		tmp.end = contig->lowCoverageAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}

}

void FRC::computeHighCoverageArea(string type, unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighCoverageAreas(this->C_A, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_COVERAGE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_COVERAGE_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highCoverageAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_COV_"+type;
		tmp.start = contig->highCoverageAreas.at(i).first;
		tmp.end = contig->highCoverageAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
		//cout << "\t" << tmp.feature << " " << tmp.start << " " << tmp.end << "\n";
	}
}

void FRC::computeLowNormalArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getLowNormalAreas(this->C_M, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateLOW_NORMAL_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateLOW_NORMAL_AREA(feat);
	}
	for(unsigned int i=0; i < contig->lowNormalAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "LOW_NORM_COV_"+type;
		tmp.start = contig->lowNormalAreas.at(i).first;
		tmp.end = contig->lowNormalAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
}

void FRC::computeHighNormalArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighNormalAreas(this->C_M, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_NORMAL_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_NORMAL_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highNormalAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_NORM_COV_"+type;
		tmp.start = contig->highNormalAreas.at(i).first;
		tmp.end = contig->highNormalAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
}

void FRC::computeHighSingleArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighSingleAreas( windowSize, windowStep, this->C_A);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_SINGLE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_SINGLE_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highSingleAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_SINGLE_"+type;
		tmp.start = contig->highSingleAreas.at(i).first;
		tmp.end = contig->highSingleAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);

	}

}

void FRC::computeHighSpanningArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighSpanningAreas( windowSize, windowStep, this->C_A);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_SPANNING_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_SPANNING_AREA(feat);
	}

	for(unsigned int i=0; i< contig->highSpanningAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_SPAN_"+type;
		tmp.start = contig->highSpanningAreas.at(i).first;
		tmp.end = contig->highSpanningAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);

	}
}

void FRC::computeHighOutieArea(string type,unsigned int ctg, Contig *contig, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getHighOutieAreas( windowSize, windowStep, this->C_A);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateHIGH_OUTIE_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateHIGH_OUTIE_AREA(feat);
	}

	for(unsigned int i=0; i < contig->highOutieAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "HIGH_OUTIE_"+type;
		tmp.start = contig->highOutieAreas.at(i).first;
		tmp.end = contig->highOutieAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}

}

void FRC::computeCompressionArea(string type,unsigned int ctg, Contig *contig, float Zscore , unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getCompressionAreas(this->insertMean, this->insertStd, Zscore, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateCOMPRESSION_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateCOMPRESSION_AREA(feat);
	}

	for(unsigned int i=0; i< contig->compressionAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "COMPR_"+type;
		tmp.start = contig->compressionAreas.at(i).first;
		tmp.end = contig->compressionAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}

}

void FRC::computeStrechArea(string type,unsigned int ctg, Contig *contig, float Zscore, unsigned int windowSize, unsigned int windowStep) {
	unsigned int feat = contig->getExpansionAreas(this->insertMean, this->insertStd, Zscore, windowSize, windowStep);
	if(type.compare("PE") == 0) {
		this->CONTIG[ctg].PE.updateSTRECH_AREA(feat);
	} else {
		this->CONTIG[ctg].MP.updateSTRECH_AREA(feat);
	}

	for(unsigned int i=0; i < contig->expansionAreas.size(); i++) {
		ternary tmp;
		tmp.feature = "STRECH_"+type;
		tmp.start = contig->expansionAreas.at(i).first;
		tmp.end = contig->expansionAreas.at(i).second;
		this->CONTIG[ctg].SUSPICIOUS_AREAS.push_back(tmp);
	}
}



void FRC::setC_A(float C_A) {
	this->C_A = C_A;
}
void FRC::setS_A(float S_A) {
	this->S_A = S_A;
}
void FRC::setC_M(float C_M) {
	this->C_M = C_M;
}
void FRC::setC_W(float C_W) {
	this->C_W = C_W;
}
void FRC::setC_S(float C_S) {
	this->C_S = C_S;
}
void FRC::setC_D(float C_D) {
	this->C_D = C_D;
}
void FRC::setInsertMean(float insertMean) {
	this->insertMean = insertMean;
}
void FRC::setInsertStd(float insertStd) {
	this->insertStd = insertStd;
}

void FRC::setID(unsigned int i, string ID) {
	CONTIG[i].setID(ID);
}

string  FRC::getID(unsigned int i) {
	return CONTIG[i].getID();
}


void FRC::printFeaturesGFF3(unsigned int ctg, ofstream &file) {
	this->CONTIG[ctg].printFeaturesGFF3(file);

}


void FRC::printFeatures(unsigned int ctg, ofstream &file) {
	this->CONTIG[ctg].printFeatures(file);

}


unsigned int FRC::getTotal(unsigned int ctg) {
	return this->CONTIG[ctg].getTotal();
}


unsigned int FRC::getFeatures(FeatureTypes type, int contig) {
	switch (type) {
		 case FRC_TOTAL:
			 return getTotal(contig);
			 break;
		 case LOW_COV_PE:
			 return getLOW_COV_PE(contig);
			 break;
		 case HIGH_COV_PE:
			 return getHIGH_COV_PE(contig);
			 break;
		 case LOW_NORM_COV_PE:
			 return getLOW_NORM_COV_PE(contig);
			 break;
		 case HIGH_NORM_COV_PE:
			 return getHIGH_NORM_COV_PE(contig);
			 break;
		 case HIGH_SINGLE_PE:
			 return getHIGH_SINGLE_PE(contig);
			 break;
		 case HIGH_SPAN_PE:
			 return getHIGH_SPAN_PE(contig);
			 break;
		 case HIGH_OUTIE_PE:
			 return getHIGH_OUTIE_PE(contig);
			 break;
		 case COMPR_PE:
			 return getCOMPR_PE(contig);
			 break;
		 case STRECH_PE:
			 return getSTRECH_PE(contig);
			 break;
		 case HIGH_SINGLE_MP:
			 return getHIGH_SINGLE_MP(contig);
			 break;
		 case HIGH_OUTIE_MP:
			 return getHIGH_OUTIE_MP(contig);
			 break;
		 case HIGH_SPAN_MP:
			 return getHIGH_SPAN_MP(contig);
			 break;
		 case COMPR_MP:
			 return getCOMPR_MP(contig);
			 break;
		 case STRECH_MP:
			 return getSTRECH_MP(contig);
			 break;
		 default:
			 cout << "THis whould never happen\n";
	}
}


unsigned int FRC::getLOW_COV_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getLOW_COV_PE();
}
unsigned int FRC::getHIGH_COV_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_COV_PE();
}
unsigned int FRC::getLOW_NORM_COV_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getLOW_NORM_COV_PE();
}
unsigned int FRC::getHIGH_NORM_COV_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_NORM_COV_PE();
}
unsigned int FRC::getHIGH_SINGLE_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_SINGLE_PE();
}
unsigned int FRC::getHIGH_OUTIE_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_OUTIE_PE();
}
unsigned int FRC::getHIGH_SPAN_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_SPAN_PE();
}
unsigned int FRC::getCOMPR_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getCOMPR_PE();
}
unsigned int FRC::getSTRECH_PE(unsigned int ctg) {
	return this->CONTIG[ctg].getSTRECH_PE();
}

unsigned int FRC::getHIGH_SINGLE_MP(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_SINGLE_MP();
}
unsigned int FRC::getHIGH_OUTIE_MP(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_OUTIE_MP();
}
unsigned int FRC::getHIGH_SPAN_MP(unsigned int ctg) {
	return this->CONTIG[ctg].getHIGH_SPAN_MP();
}
unsigned int FRC::getCOMPR_MP(unsigned int ctg) {
	return this->CONTIG[ctg].getCOMPR_MP();
}
unsigned int FRC::getSTRECH_MP(unsigned int ctg) {
	return this->CONTIG[ctg].getSTRECH_MP();
}






/////////////////////

contigFeatures::contigFeatures() {
	contigLength = 0;
	TOTAL = 0;
	SUSPICIOUS_AREAS.clear();
}

contigFeatures::~contigFeatures() {

}

void contigFeatures::setContigLength(unsigned int contigLength) {
	this->contigLength = contigLength;
}

unsigned long int contigFeatures::getContigLength() {
	return this->contigLength;
}




void contigFeatures::setID(string ID) {
	this->contigID = ID;
}

string contigFeatures::getID() {
	return this->contigID;
}

unsigned int contigFeatures::getTotal() {
	TOTAL = PE.returnTotal() + MP.returnTotal();
	return TOTAL;
}

unsigned int contigFeatures::getLOW_COV_PE() {
	unsigned int features = PE.returnLOW_COV();
	return features;
}
unsigned int contigFeatures::getHIGH_COV_PE() {
	unsigned int features = PE.returnHIGH_COV();
	return features;
}
unsigned int contigFeatures::getLOW_NORM_COV_PE() {
	unsigned int features = PE.returnLOW_NORM_COV();
	return features;
}
unsigned int contigFeatures::getHIGH_NORM_COV_PE() {
	unsigned int features = PE.returnHIGH_NORM_COV();
	return features;
}
unsigned int contigFeatures::getHIGH_SINGLE_PE() {
	unsigned int features = PE.returnHIGH_SINGLE();
	return features;
}
unsigned int contigFeatures::getHIGH_OUTIE_PE() {
	unsigned int features = PE.returnHIGH_OUTIE();
	return features;
}
unsigned int contigFeatures::getHIGH_SPAN_PE() {
	unsigned int features = PE.returnHIGH_SPAN();
	return features;
}
unsigned int contigFeatures::getCOMPR_PE() {
	unsigned int features = PE.returnCOMPR();
	return features;
}
unsigned int contigFeatures::getSTRECH_PE() {
	unsigned int features = PE.returnSTRECH();
	return features;
}

unsigned int contigFeatures::getHIGH_SINGLE_MP() {
	unsigned int features = MP.returnHIGH_SINGLE();
	return features;
}
unsigned int contigFeatures::getHIGH_OUTIE_MP() {
	unsigned int features =  MP.returnHIGH_OUTIE();
	return features;
}
unsigned int contigFeatures::getHIGH_SPAN_MP() {
	unsigned int features = MP.returnHIGH_SPAN();
	return features;
}
unsigned int contigFeatures::getCOMPR_MP() {
	unsigned int features = MP.returnCOMPR();
	return features;
}
unsigned int contigFeatures::getSTRECH_MP() {
	unsigned int features = MP.returnSTRECH();
	return features;
}




bool sortTernary(ternary t1, ternary t2) {return (t1.start < t2.start);}

void contigFeatures::printFeatures(ofstream &file) {
	unsigned int total = 0;
	sort(SUSPICIOUS_AREAS.begin(), SUSPICIOUS_AREAS.end(), sortTernary);
	for(unsigned int i=0; i < SUSPICIOUS_AREAS.size(); i++) {
		file << this->contigID << " " << SUSPICIOUS_AREAS[i].feature << " " << SUSPICIOUS_AREAS[i].start << " " << SUSPICIOUS_AREAS[i].end << "\n";
		total += floor((SUSPICIOUS_AREAS[i].end - SUSPICIOUS_AREAS[i].start + 1)/(float)1000 + 0.5);
		//total++;
	}
	this->TOTAL = total;


}

void contigFeatures::printFeaturesGFF3(ofstream &file) {

	file << "##sequence-region\t" << this->contigID <<  "\t" << 1 << "\t" << this->contigLength << "\n";

	sort(SUSPICIOUS_AREAS.begin(), SUSPICIOUS_AREAS.end(), sortTernary);
	for(unsigned int i=0; i < SUSPICIOUS_AREAS.size(); i++) {
		file << this->contigID << "\t" << "." << "\t" << SUSPICIOUS_AREAS[i].feature << "\t";
		file << SUSPICIOUS_AREAS[i].start + 1 << "\t" <<  SUSPICIOUS_AREAS[i].end -1 << "\t";
		file << "." << "\t" << "+" << "\t" << "." << "\t" << "Name=" << SUSPICIOUS_AREAS[i].feature << "\n";
	}


}



