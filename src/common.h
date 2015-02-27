#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <sstream>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"


#ifdef INLINE_DISABLED
#define INLINE
#else
#define INLINE inline
#endif


using namespace BamTools;
using namespace std;

#define DEFAULT_CHANNEL std::cout
#define ERROR_CHANNEL std::cerr

#define VERBOSE_CHANNEL std::cerr
#define DEBUG_CHANNEL std::cerr
#define DEFAULT_CHANNEL std::cout

static inline std::string package_description() {
	std::string line("FRC");
	line.append(" version ");
	line.append("1.2.0");
	return line;
}

enum readStatus {unmapped, lowQualty, singleton, pair_wrongChrs,
	pair_proper, pair_wrongDistance, pair_wrongOrientation};


enum Feature {LOW_COVERAGE_AREA, HIGH_COVERAGE_AREA, LOW_NORMAL_AREA, HIGH_NORMAL_AREA, HIGH_SINGLE_AREA, HIGH_SPANNING_AREA, HIGH_OUTIE_AREA, COMPRESSION_AREA, STRECH_AREA, TOTAL};


enum FeatureTypes {
	FRC_TOTAL,
	LOW_COV_PE,
	HIGH_COV_PE,
	LOW_NORM_COV_PE,
	HIGH_NORM_COV_PE,
	HIGH_SINGLE_PE,
	HIGH_SPAN_PE,
	HIGH_OUTIE_PE,
	COMPR_PE,
	STRECH_PE,
	HIGH_SINGLE_MP,
	HIGH_OUTIE_MP,
	HIGH_SPAN_MP,
	COMPR_MP,
	STRECH_MP,
};

static string returnFeatureName(FeatureTypes type) {
	switch (type) {
	case FRC_TOTAL:
		return "TOTAL";
		break;
	case LOW_COV_PE:
		return "LOW_COV_PE";
		break;
	case HIGH_COV_PE:
		return "HIGH_COV_PE";
		break;
	case LOW_NORM_COV_PE:
		return "LOW_NORM_COV_PE";
		break;
	case HIGH_NORM_COV_PE:
		return "HIGH_NORM_COV_PE";
		break;
	case HIGH_SINGLE_PE:
		return "HIGH_SINGLE_PE";
		break;
	case HIGH_SPAN_PE:
		return "HIGH_SPAN_PE";
		break;
	case HIGH_OUTIE_PE:
		return "HIGH_OUTIE_PE";
		break;
	case COMPR_PE:
		return "COMPR_PE";
		break;
	case STRECH_PE:
		return "STRECH_PE";
		break;
	case HIGH_SINGLE_MP:
		return "HIGH_SINGLE_MP";
		break;
	case HIGH_OUTIE_MP:
		return "HIGH_OUTIE_MP";
		break;
	case HIGH_SPAN_MP:
		return "HIGH_SPAN_MP";
		break;
	case COMPR_MP:
		return "COMPR_MP";
		break;
	case STRECH_MP:
		return "STRECH_MP";
		break;
	default:
		cout << "THis whould never happen\n";
	}

}

static int StringToNumber ( string Text ) {
	stringstream ss(Text);
	int result;
	return ss >> result ? result : 0;
}


struct LibraryStatistics{
	float C_A;
	float S_A;
	float C_D;
	float C_M;
	float C_S;
	float C_W;
	float insertMean;
	float insertStd;
};






static readStatus computeReadType(BamAlignment al, uint32_t max_insert, bool is_mp) {
	if (!al.IsMapped()) {
		return unmapped;
	}
	if((al.IsDuplicate()) || (al.IsFailedQC()) || (!al.IsPrimaryAlignment())) {
		return lowQualty;
	}
	if(!(al.IsMateMapped())) {
		return singleton;
	}
	if (al.IsMateMapped() && al.RefID != al.MateRefID) {
		return pair_wrongChrs;
	}
	//If I am here the read must be aligned, with the pair/mate aligned on the same contig/scaffold
	uint32_t startRead   = al.Position; // start position on the contig
	uint32_t startPaired = al.MatePosition;
	int iSize = al.InsertSize;
	if (iSize < 0) { iSize = -1 * iSize;}
	//Now check if reads belong to a proper pair: both reads aligned on the same contig at the expected distance and orientation
	if (iSize > max_insert) {
		return pair_wrongDistance;
	}
	if (! is_mp) { // I have a paired end
		if(startRead < startPaired) { //
			if(!(al.IsReverseStrand()) && (al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		} else {
			if((al.IsReverseStrand()) && !(al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		}
	} else {
		if(startRead < startPaired) { //
			if((al.IsReverseStrand()) && !(al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		} else {
			if(!(al.IsReverseStrand()) && (al.IsMateReverseStrand())) {
				return pair_proper;
			} else {
				return pair_wrongOrientation;
			}
		}
	}
}




static float normcdf(float x,float mu,float sigma) {
	float t = x - mu;
	float y = 0.5 * erfc(-t / (sigma * sqrt(2.0)));
	if (y > 1.0) {
		y = 1.0;
	}
	return y;
}

static float normpdf(float x,float mu,float sigma) {
    float u = (x - mu) / abs(sigma);
    float y = (1 /(sqrt(2 * M_PI) * abs(sigma))) * exp((-u * u)/2) ;
    return y;
}


static float Part(float a, float b,  uint32_t sizeA, uint32_t sizeB, uint32_t gap, float insert_mean, float insert_stddev, float coverage, uint32_t readLength ) {
	float readfrequency = 2 * readLength / coverage;
	float expr1 = (min(sizeA, sizeB) - (readLength - 0)) / readfrequency * normcdf(a, 0, 1);
	float expr2 = -(- 0 ) / readfrequency * normcdf(b, 0, 1);
	float expr3 = (b * insert_stddev) / readfrequency * (normcdf(b, 0, 1) - normcdf(a, 0, 1));
	float expr4 = (insert_stddev / readfrequency) * (normpdf(b, 0, 1) - normpdf(a, 0, 1));
	float value = expr1 + expr2 + expr3 + expr4;
	return value;
}


static float ExpectedLinks(uint32_t sizeA, uint32_t sizeB, uint32_t gap, float insert_mean, float insert_stddev, float coverage, uint32_t readLength) {
	float b1 = (sizeA + sizeB + gap - insert_mean) / insert_stddev;
	float a1 = (max(sizeA, sizeB) + gap + readLength  - insert_mean) / insert_stddev;
	float b2 = (min(sizeA, sizeB) + gap + readLength  - insert_mean) / insert_stddev;
	float a2 = (gap + 2 * readLength - insert_mean) / insert_stddev;

	float E_links = Part(a1, b1, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength ) - Part(a2, b2, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength);
	return E_links;
}


static LibraryStatistics computeLibraryStats(string bamFileName, uint64_t genomeLength, uint32_t max_insert, bool is_mp) {
	BamReader bamFile;
	bamFile.Open(bamFileName);
	LibraryStatistics library;

	//All var declarations
	uint32_t reads 		 	  = 0;
	uint32_t unmappedReads 	  = 0;
	uint32_t lowQualityReads  = 0;
	uint32_t mappedReads 	  = 0;
	uint64_t mappedReadsLength= 0;

	uint64_t insertsLength = 0; // total inserts length
	float insertMean;
	float insertStd;
	// mated reads (not necessary correctly mated)
	uint32_t matedReads 	  = 0;        // reads that align on a contig with the mate
	uint64_t matedReadsLength = 0;  // total length of mated reads
	// wrongly distance
	uint32_t wrongDistanceReads 		= 0;  // number of paired reads too far away
	uint64_t wrongDistanceReadsLength   = 0; // length  of paired reads too far away
	// wrongly oriented reads
	uint32_t wronglyOrientedReads 		= 0;       // number of wrongly oriented reads
	uint64_t wronglyOrientedReadsLength = 0; // length of wrongly oriented reads
	// singletons
	uint32_t singletonReads 	  = 0; // number of singleton reads
	uint64_t singletonReadsLength = 0;     // total length of singleton reads
	// mates on different contigs
	uint32_t matedDifferentContig 		= 0; // number of contig placed in a different contig
	uint64_t matedDifferentContigLength = 0; // total number of reads placed in different contigs

	float C_A = 0; // total read coverage
	float S_A = 0; // total span coverage
	float C_M = 0; // coverage induced by proper pairs (same contig and correct orientation)
	float C_W = 0; // coverage induced by wrongly mated pairs
	float C_S = 0; // coverage induced by singletons
	float C_D = 0; // coverage induced by reads with mate on a different contigs

	// compute mean and std on the fly
	float Mk = 0;
	float Qk = 0;
	uint32_t counterK = 1;
	//Keep header for further reference
	int32_t currentTid = -1;
	int32_t iSize;

	BamAlignment al;
	while ( bamFile.GetNextAlignmentCore(al) ) {
		reads ++;
		readStatus read_status = computeReadType(al, max_insert, is_mp);
		if (read_status != unmapped and read_status != lowQualty) {
			mappedReads ++;
			mappedReadsLength += al.Length;
		}

		if (al.IsFirstMate() && read_status == pair_proper) {
			iSize = abs(al.InsertSize);
			if(counterK == 1) {
				Mk = iSize;
				Qk = 0;
				counterK++;
			} else {
				float oldMk = Mk;
				float oldQk = Qk;
				Mk = oldMk + (iSize - oldMk)/counterK;
				Qk = oldQk + (counterK-1)*(iSize - oldMk)*(iSize - oldMk)/(float)counterK;
				counterK++;
			}
			insertsLength += iSize;
		}

		switch (read_status) {
		  case unmapped:
			  unmappedReads ++;
		     break;
		  case lowQualty:
			  lowQualityReads ++;
		     break;
		  case singleton:
			  singletonReads ++;
			  singletonReadsLength += al.Length ;
			  break;;
		  case pair_wrongChrs:
			  matedDifferentContig ++;
			  matedDifferentContigLength += al.Length ;
			  break;
		  case pair_proper:
			  matedReads ++;
			  matedReadsLength += al.Length ;
			  break;
		  case pair_wrongDistance:
			  wrongDistanceReads ++;
			  wrongDistanceReadsLength += al.Length ;
			  break;
		  case pair_wrongOrientation:
			  wronglyOrientedReads ++;
			  wronglyOrientedReadsLength += al.Length ;
			  break;
		  default:
		     cout << "This should never be printed\n";
		     break;
		}

	}



	cout << "LIBRARY STATISTICS\n";
	cout << "\t total reads number "	<< reads << "\n";
	cout << "\t total mapped reads " 	<< mappedReads << "\n";
	cout << "\t total unmapped reads " 	<< unmappedReads << "\n";
	cout << "\t proper pairs " 			<< matedReads << "\n";
	cout << "\t wrong distance "		<< wrongDistanceReads << "\n";
	cout << "\t zero quality reads " 	<< lowQualityReads << "\n";
	cout << "\t wrongly oriented "		<< wronglyOrientedReads << "\n";
	cout << "\t wrongly contig "		<< matedDifferentContig << "\n";
	cout << "\t singletons " 			<< singletonReads << "\n";

	uint32_t total = matedReads + wrongDistanceReads +  wronglyOrientedReads +  matedDifferentContig + singletonReads  ;
	cout << "\ttotal " << total << "\n";
	cout << "\tCoverage statistics\n";

	library.C_A = C_A = mappedReadsLength/(float)genomeLength;
	library.S_A = S_A = insertsLength/(float)genomeLength;
	library.C_M = C_M = matedReadsLength/(float)genomeLength;
	library.C_W = C_W = wronglyOrientedReadsLength/(float)genomeLength;
	library.C_S = C_S = singletonReadsLength/(float)genomeLength;
	library.C_D = C_D = matedDifferentContigLength/(float)genomeLength;
	library.insertMean = insertMean = Mk;
	Qk = sqrt(Qk/counterK);
	library.insertStd = insertStd = Qk;

	cout << "\tC_A = " << C_A << endl;
	cout << "\tS_A = " << S_A << endl;
	cout << "\tC_M = " << C_M << endl;
	cout << "\tC_W = " << C_W << endl;
	cout << "\tC_S = " << C_S << endl;
	cout << "\tC_D = " << C_D << endl;
	cout << "\tMean Insert length = " << Mk << endl;
	cout << "\tStd Insert length = " << Qk << endl;
	cout << "----------\n";

	bamFile.Close();
	return library;
}








#endif /*TYPES_H_*/
