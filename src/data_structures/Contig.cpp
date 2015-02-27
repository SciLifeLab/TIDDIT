/*
 * Contig.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */


#include "Contig.h"



ContigsFeat::ContigsFeat() {
	this->feature=0;
	this->startStopPositions.clear();
}

ContigsFeat::~ContigsFeat() {

}



Position::Position() {
	ReadCoverage = 0;
	StratingInserts = 0;
	InsertCoverage = 0;
	CorrectlyMated = 0;
	WronglyOriented = 0;
	Singleton = 0;
	MatedDifferentContig = 0;
	insertsLength = 0;
}

Position::~Position() {

}


Contig::Contig() {
	contigLength = 0;
	MINUM_COV = 0;

	lowCoverageFeat = 1/(float)3.5;
	highCoverageFeat = 2.5;
	lowNormalFeat = 1/(float)3.5;
	highNormalFeat = 2.5;
	highSingleFeat = 0.4;
	highSpanningFeat = 0.51;
	highOutieFeat = 0.51;
}

Contig::Contig(unsigned int contigLength) {
	this->contigLength = contigLength;
	this->CONTIG =  new Position[contigLength];
	MINUM_COV = 2;

	lowCoverageFeat = 1/(float)2;
	highCoverageFeat = 2;
	lowNormalFeat = 1/(float)2;
	highNormalFeat = 2;
	highSingleFeat = 0.41;
	highSpanningFeat = 0.41;
	highOutieFeat = 0.41;


}

Contig::~Contig() {
	if(CONTIG != NULL) {
		delete [] CONTIG;
	}
}


void Contig::updateCov(unsigned int start, unsigned int end, data type) {
	if(end > this->contigLength) {
//		cout << "hoops, end longer than contig length when updating CONTIG " << type << "\n";
//		cout << "\tcontig length " << this->contigLength << " starting point " << start << " ending point " << end << "\n";
		end = this->contigLength;
	}
	// now update
	if(type == insertCov) {
		CONTIG[start].StratingInserts++; // a new inserts starts in position start
		CONTIG[start].insertsLength += (end - start + 1); // save total length of inserts starting at start
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].InsertCoverage++;
	} else if(type == readCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].ReadCoverage++;
	} else if(type ==  cmCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].CorrectlyMated++;
	} else if(type == woCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].WronglyOriented++;
	} else if(type == singCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].Singleton++;
	} else if(type == mdcCov) {
		for(unsigned int i = start; i< end; i++)
			CONTIG[i].MatedDifferentContig++;
	} else {
		cout << "hoops, unknown type " << type << " there must be something wrong!!!\n";
	}
}

unsigned int Contig::getContigLength() {
	return this->contigLength;
}




void Contig::updateContig(BamAlignment b, int max_nsert, bool is_mp) {
	readStatus read_status 	= computeReadType(b, max_nsert, is_mp);
	uint32_t readLength     = b.Length;
	uint32_t iSize 			= abs(b.InsertSize);
	uint32_t startRead 		= b.Position;
	uint32_t endRead 		= startRead + readLength ; // position where reads ends
	uint32_t startMateRead  = b.MatePosition;
	if (read_status == unmapped or read_status == lowQualty) {
		return;
	}
	if (read_status != unmapped and read_status != lowQualty) { //if the read is aligned and is not duplicated or low quality use it in cov computation
		updateCov(startRead, endRead, readCov); // update coverage
	}
	if (b.IsFirstMate() && read_status == pair_proper) {
		int iSize = abs(b.InsertSize);
		if(startRead < startMateRead) {
			updateCov(startRead, startRead + iSize, insertCov);
		} else {
			updateCov(startMateRead, startMateRead + iSize, insertCov);
		}
	}
	switch (read_status) {
	case singleton:
		updateCov(startRead, endRead, singCov);
		break;
	case pair_wrongChrs:
		updateCov(startRead, endRead, mdcCov);
		break;
	case pair_wrongDistance:
		updateCov(startRead, endRead, woCov); //
		break;
	case pair_wrongOrientation:
		updateCov(startRead, endRead, woCov);
		break;
	case pair_proper:
		updateCov(startRead, endRead, cmCov);
		break;
	default:
		cout << read_status << " --> This should never be printed\n";
		break;

	}
}



void Contig::print() {
	cout << "Contig size " << this->contigLength << "\n";
	for(unsigned int i= 0; i < this->contigLength; i++) {
		if(i % 6 == 0 && i > 0) {
			cout << "\n";
		}
		cout << "(" << i << ":" << this->CONTIG[i].ReadCoverage  << "," << this->CONTIG[i].InsertCoverage << "," << this->CONTIG[i].CorrectlyMated << "," <<
				this->CONTIG[i].MatedDifferentContig  << "," << "," <<   this->CONTIG[i].Singleton << "," <<   this->CONTIG[i].WronglyOriented
				<< "," << this->CONTIG[i].StratingInserts << ") " ;
	}
	cout << "\n\n";
}


float Contig::getCoverage() {
	//compute length of low coverage areas in the contig
	//use a 1K sliding window

	unsigned int totalCoverage = 0;
	float meanCov;
	for(unsigned int i=0; i < this->contigLength ; i++ ) {
		totalCoverage += CONTIG[i].ReadCoverage;
	}
	meanCov = totalCoverage/(float)this->contigLength; // this is the "window" coverage
	return meanCov;

}


unsigned int Contig::getLowCoverageAreas(float C_A, unsigned int windowSize, unsigned int windowStep) {
	//compute length of low coverage areas in the contig
	//use a 1K sliding window

	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" coverage
		if(meanCov < lowCoverageFeat*C_A and meanCov > MINUM_COV ) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->lowCoverageAreas.push_back(SS);
			//Areas->feature=1;
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov < lowCoverageFeat*C_A and meanCov > MINUM_COV ) { // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}

		//cout << feat << " " << startWindow << " " << meanCov << " " << lowCoverageFeat*C_A << "\n";

		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			//cout << feat << " " << startWindow << " " << meanCov << " " << lowCoverageFeat*C_A << "\n";
			if(meanCov < lowCoverageFeat*C_A and meanCov > MINUM_COV) { // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->lowCoverageAreas.push_back(SS);
					//Areas->feature+=floor((endFeat - startFeat + 1)/(float)windowSize + 0.5);


					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		//need to check if last window was a feature
		if(feat) { //reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->lowCoverageAreas.push_back(SS);
//			Areas->feature+=floor((endFeat - startFeat + 1)/(float)windowSize + 0.5);
		}
	}


	return features;
}




unsigned int  Contig::getHighCoverageAreas(float C_A, unsigned int windowSize, unsigned int windowStep) {
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		if(meanCov > highCoverageFeat*C_A  ) { // this is a feature
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->highCoverageAreas.push_back(SS);
			features = 1; // one feature found (in one window)
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov > highCoverageFeat*C_A ) { // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}

		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			if(meanCov > highCoverageFeat*C_A ) { // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->highCoverageAreas.push_back(SS);
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		if(feat) { // a feature  reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->highCoverageAreas.push_back(SS);
		}
	}


	return features;
}




unsigned int Contig::getLowNormalAreas(float C_M, unsigned int windowSize, unsigned int windowStep) {

	//compute length of low coverage areas in the contig
	//use a 1K sliding window
	unsigned int totalCoverageRead = 0;
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov, thr;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].CorrectlyMated ;
			totalCoverageRead += CONTIG[i].ReadCoverage ;

		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		thr = totalCoverageRead/(float)this->contigLength;
		if(meanCov < lowNormalFeat*C_M and thr > MINUM_COV ) { // this is a feature
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->lowNormalAreas.push_back(SS);
			features = 1; // one feature found (in one window)

		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
				totalCoverageRead += CONTIG[i].ReadCoverage ;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		thr = totalCoverageRead/(float)winSize;
		if(meanCov < lowNormalFeat*C_M and thr > MINUM_COV) { // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			totalCoverageRead = 0;
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
				totalCoverageRead += CONTIG[i].ReadCoverage ;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			thr = totalCoverageRead/(float)(endWindow - startWindow);

			if(meanCov < lowNormalFeat*C_M and thr > MINUM_COV) { // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->lowNormalAreas.push_back(SS);
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		if(feat) { // a feature reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->lowNormalAreas.push_back(SS);
		}
	}

	return features;
}



unsigned int Contig::getHighNormalAreas(float C_M, unsigned int windowSize, unsigned int windowStep) {

	//compute length of low coverage areas in the contig
	//use a 1K sliding window
	unsigned int totalCoverage = 0;
	unsigned int features = 0;
	float meanCov;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].CorrectlyMated ;
		}
		meanCov = totalCoverage/(float)this->contigLength; // this is the "window" covrage
		if(meanCov > highNormalFeat*C_M ) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->highNormalAreas.push_back(SS);
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
		}
		meanCov = totalCoverage/(float)winSize; // first window's covrage
		if(meanCov > highNormalFeat*C_M ) {  // in the first window already present a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].CorrectlyMated;
			}
			meanCov = totalCoverage/(float)(endWindow - startWindow); // compute window coverage
			if(meanCov > highNormalFeat*C_M ) {  // in the first window already present a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->highNormalAreas.push_back(SS);
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		if(feat) { // a feature  reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->highNormalAreas.push_back(SS);
		}
	}
	return features;
}

unsigned int Contig::getHighSingleAreas( unsigned int windowSize, unsigned int windowStep, float C_A) {

	unsigned int totalCoverage = 0;
	unsigned int singleReadCoverage = 0;
	unsigned int features = 0;
	float meanTotalCov;
	float meanSingleCov;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			singleReadCoverage += CONTIG[i].Singleton;
		}
		meanTotalCov = totalCoverage/(float)this->contigLength; // this is the "window" total coverage
		meanSingleCov = singleReadCoverage/(float)this->contigLength; // this is the "window" single read coverage
		if( meanSingleCov > highSingleFeat*meanTotalCov and meanTotalCov > MINUM_COV ) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->highSingleAreas.push_back(SS);
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			singleReadCoverage += CONTIG[i].Singleton;
		}
		meanTotalCov = totalCoverage/(float)winSize; // this is the "window" total coverage
		meanSingleCov = singleReadCoverage/(float)winSize; // this is the "window" single read coverage

		if(meanTotalCov > 0.5 * C_A) {
			if( meanSingleCov > highSingleFeat*meanTotalCov  and meanTotalCov > MINUM_COV ) { //first window's covrage
				startFeat = 0;
				endFeat = windowSize;
				feat = true; // there is an open feature
			}
		}
//		cout << feat << " " << meanTotalCov << " " << meanSingleCov << "\n";

		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			singleReadCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage ;
				singleReadCoverage += CONTIG[i].Singleton;
			}
			meanTotalCov = totalCoverage/(float)(endWindow - startWindow); // compute window total coverage
			meanSingleCov = singleReadCoverage/(float)(endWindow - startWindow); // compute window single read coverage
			if( meanSingleCov > highSingleFeat*meanTotalCov and meanTotalCov > MINUM_COV) {
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->highSingleAreas.push_back(SS);
					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
			//cout << feat << " " << startWindow << " " << meanTotalCov << " " << meanSingleCov << "\n";
		}
		if(feat) { // a feature  reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->highSingleAreas.push_back(SS);
		}
	}
	return features;

}



unsigned int Contig::getHighSpanningAreas( unsigned int windowSize, unsigned int windowStep, float C_A ) {

	unsigned int totalCoverage = 0;
	unsigned int matedDifferentContigCoverage = 0;
	unsigned int features = 0;
	float meanTotalCov;
	float meanMatedDifferentContigCoverage;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			matedDifferentContigCoverage += CONTIG[i].MatedDifferentContig;
		}
		meanTotalCov = totalCoverage/(float)this->contigLength; // this is the "window" total coverage
		meanMatedDifferentContigCoverage = matedDifferentContigCoverage/(float)this->contigLength; // this is the "window" single read coverage
		if( meanMatedDifferentContigCoverage > highSpanningFeat*meanTotalCov  and meanTotalCov > MINUM_COV ) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->highSpanningAreas.push_back(SS);
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			matedDifferentContigCoverage += CONTIG[i].MatedDifferentContig;
		}
		meanTotalCov = totalCoverage/(float)winSize; //
		meanMatedDifferentContigCoverage = matedDifferentContigCoverage/(float)winSize; //
		if( meanMatedDifferentContigCoverage > highSpanningFeat*meanTotalCov  and meanTotalCov > MINUM_COV ) { // this is a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			totalCoverage = 0; //reset window coverage
			matedDifferentContigCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage ;
				matedDifferentContigCoverage += CONTIG[i].MatedDifferentContig;
			}
			meanTotalCov = totalCoverage/(float)(endWindow - startWindow); // compute window total coverage
			meanMatedDifferentContigCoverage = matedDifferentContigCoverage/(float)(endWindow - startWindow); // compute window single read coverage
			if( meanMatedDifferentContigCoverage > highSpanningFeat*meanTotalCov  and meanTotalCov > MINUM_COV) { // this is a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->highSpanningAreas.push_back(SS);

					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		if(feat) { // a feature  reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->highSpanningAreas.push_back(SS);
		}
	}

	return features;


}



unsigned int Contig::getHighOutieAreas( unsigned int windowSize, unsigned int windowStep, float C_A ) {
	unsigned int totalCoverage = 0;
	unsigned int outieCoverage = 0;
	unsigned int features = 0;
	float meanTotalCov;
	float meanOutieCoverage;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			outieCoverage += CONTIG[i].WronglyOriented;
		}
		meanTotalCov = totalCoverage/(float)this->contigLength; // this is the "window" total coverage
		meanOutieCoverage = outieCoverage/(float)this->contigLength; // this is the "window" single read coverage
		if( meanOutieCoverage > highOutieFeat*meanTotalCov  and meanTotalCov > MINUM_COV) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->highOutieAreas.push_back(SS);
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		unsigned int winSize     = windowSize;
		for(unsigned int i=startWindow; i < endWindow ; i++ ) {
			totalCoverage += CONTIG[i].ReadCoverage ;
			outieCoverage +=  (CONTIG[i].WronglyOriented);
		}
		meanTotalCov = totalCoverage/(float)winSize; //
		meanOutieCoverage = outieCoverage/(float)winSize; //
		if(  meanOutieCoverage > highOutieFeat*meanTotalCov   and meanTotalCov > MINUM_COV) { // this is a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}


		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			meanTotalCov = 0;
			meanOutieCoverage = 0;
			totalCoverage = 0;
			outieCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow ; i++ ) {
				totalCoverage += CONTIG[i].ReadCoverage ;
				outieCoverage += CONTIG[i].WronglyOriented;
			}
			meanTotalCov = totalCoverage/(float)(endWindow - startWindow); // compute window total coverage
			meanOutieCoverage = outieCoverage/(float)(endWindow - startWindow); // compute window single read coverage
			if(  meanOutieCoverage > highOutieFeat*meanTotalCov   and meanTotalCov > MINUM_COV) { // this is a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->highOutieAreas.push_back(SS);

					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		if(feat) { // a feature  reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->highOutieAreas.push_back(SS);
		}
	}
	return features;

}




unsigned int Contig::getCompressionAreas(float insertionMean, float insertionStd, float Zscore, unsigned int windowSize, unsigned int windowStep) {

	unsigned long int spanningCoverage = 0; // total insert length
	unsigned int inserts = 0; // number of inserts
	unsigned int features = 0;
	float Z_stats;
	float localMean;
	unsigned int minInsertNum = 5;

	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			if(CONTIG[i].StratingInserts > 0) {
				inserts += CONTIG[i].StratingInserts;
				spanningCoverage += CONTIG[i].insertsLength;
			}
		}
		if(inserts > minInsertNum) {
			localMean = spanningCoverage/(float)inserts;
			Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
		} else {
			Z_stats = 0;
		}
		if( Z_stats <  Zscore ) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->compressionAreas.push_back(SS);
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		for(unsigned int i=startWindow; i < endWindow; i++ ) {
			if(CONTIG[i].StratingInserts > 0) {
				inserts += CONTIG[i].StratingInserts;
				spanningCoverage += CONTIG[i].insertsLength;
			}
		}
		if(inserts > minInsertNum) {
			localMean = spanningCoverage/(float)inserts;
			Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
		} else {
			Z_stats = 0;
		}

		//cout << feat << " " << startWindow << " " << Z_stats << " " << inserts << " " << insertionMean << " " << insertionStd << " " << localMean << " " << localMean - insertionMean << "\n";

		if( Z_stats <  Zscore ) { // this is a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			inserts = 0; //reset window coverage
			spanningCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow; i++ ) {
				if(CONTIG[i].StratingInserts > 0) {
					inserts += CONTIG[i].StratingInserts;
					spanningCoverage += CONTIG[i].insertsLength;
				}
			}
			if(inserts > minInsertNum) {
				localMean = spanningCoverage/(float)inserts;
				Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
			} else {
				Z_stats = 0;
			}

			//cout << feat << " " << startWindow << " " << Z_stats << " " << inserts << " " << insertionMean << " " << insertionStd << " " << localMean << " " << localMean - insertionMean << "\n";
			if( Z_stats <  Zscore ) { // this is a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->compressionAreas.push_back(SS);

					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
		}
		if(feat) { // a feature  reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->compressionAreas.push_back(SS);
		}
	}

	return features;

}


unsigned int Contig::getExpansionAreas(float insertionMean, float insertionStd, float Zscore, unsigned int windowSize, unsigned int windowStep) {
	unsigned long int spanningCoverage = 0; // total insert length
	unsigned int inserts = 0; // number of inserts
	unsigned int features = 0;
	float localMean;
	float Z_stats = 0;
	unsigned int minInsertNum = 5;
	if(this->contigLength < windowSize) { // if contig less than window size, only one window
		for(unsigned int i=0; i < this->contigLength ; i++ ) {
			if(CONTIG[i].StratingInserts > 0) {
				inserts += CONTIG[i].StratingInserts;
				spanningCoverage += CONTIG[i].insertsLength;
			}
		}
		if(inserts > minInsertNum) {
			localMean = spanningCoverage/(float)inserts;
			Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
		} else {
			Z_stats = 0;
		}
		if( Z_stats >  Zscore ) { // this is a feature
			features = 1; // one feature found (in one window)
			pair<unsigned int , unsigned int > SS (0, this->contigLength);
			this->expansionAreas.push_back(SS);
		}
	} else { //otherwise compute features on sliding window of 200 bp
		unsigned int startFeat, endFeat;
		bool feat = false;
		unsigned int startWindow = 0;
		unsigned int endWindow   = windowSize;
		for(unsigned int i=startWindow; i < endWindow; i++ ) {
			if(CONTIG[i].StratingInserts > 0) {
				inserts += CONTIG[i].StratingInserts;
				spanningCoverage += CONTIG[i].insertsLength;
			}
		}
		if(inserts > minInsertNum) {
			localMean = spanningCoverage/(float)inserts;
			Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
		} else {
			Z_stats = 0;
		}
		if( Z_stats > Zscore ) { // this is a feature
			startFeat = 0;
			endFeat = windowSize;
			feat = true; // there is an open feature
		}
	//	cout << feat << " " << startWindow << " " << Z_stats << " " << inserts << " " << insertionMean << " " << insertionStd << " " << localMean << " " << localMean - insertionMean << "\n";
		//now update
		startWindow += windowStep;
		endWindow += windowStep;
		if(endWindow > this->contigLength) {
			endWindow = this->contigLength;
		}

		while(endWindow < this->contigLength) {
			inserts = 0; //reset window coverage
			spanningCoverage = 0;
			for(unsigned int i=startWindow; i < endWindow; i++ ) {
				if(CONTIG[i].StratingInserts > 0) {
					inserts += CONTIG[i].StratingInserts;
					spanningCoverage += CONTIG[i].insertsLength;
				}
			}
			if(inserts > minInsertNum) {
				localMean = spanningCoverage/(float)inserts;
				Z_stats   = (localMean - insertionMean)/(float)(insertionStd/sqrt(inserts)); // CE statistics
			} else {
				Z_stats = 0;
			}
			if( Z_stats > Zscore ) { // this is a feature
				if(feat) { // if we are already inside a feature area
					endFeat = endWindow; // simply extend the feature area
				} else {
					startFeat = startWindow;
					endFeat = endWindow;
					feat = true; // open feature area
				}
				startWindow += windowStep;
				endWindow += windowStep;
				if(endWindow > this->contigLength) {
					endWindow = this->contigLength;
				}
			} else { // this window is not affected by feature
				if(feat) { // if before a
					features += floor((endFeat - startFeat + 1)/(float)windowSize + 0.5) ; // compute number of features
					pair<unsigned int , unsigned int > SS (startFeat, endFeat);
					this->expansionAreas.push_back(SS);

					startWindow = endWindow;
					endWindow = startWindow + windowSize;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
					feat = false; //close feature area
				} else { // no feature was present in the window before
					startWindow += windowStep;
					endWindow += windowStep;
					if(endWindow > this->contigLength) {
						endWindow = this->contigLength;
					}
				}
			}
	//		cout << feat << " " << startWindow << " " << Z_stats << " " << inserts << " " << insertionMean << " " << insertionStd << " " << localMean << " " << localMean - insertionMean << "\n";
		}
		//compute last window statistics

		if(feat) { // a feature reached contig end
			features +=   floor((endFeat - startFeat)/(float)windowSize + 0.5); // compute number of features
			pair<unsigned int , unsigned int > SS (startFeat, endFeat);
			this->expansionAreas.push_back(SS);
		}
	}
	return features;

}



