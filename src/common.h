

#ifndef TYPES_H_
#define TYPES_H_


#include <string>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <cstring>
#include <queue>


#include <climits>
#include <cstdlib>
#include <sstream>

#include "api/BamAux.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"



#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.141592654
#endif

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
	std::string line("TIDDIT");
	line.append(" version ");
	line.append("1.0.0");
	return line;
}

enum readStatus {unmapped, lowQualty, singleton, pair_wrongChrs,
	pair_proper, pair_wrongDistance, pair_wrongOrientation};



class Cov{
public:
    //the vraiables used in the coverage calculation function
    ofstream coverageOutput;
	int binSize;
    int binStart;
	int binEnd;
	int currentChr;
	vector< vector<unsigned int> > coverageStructure;
	map<unsigned int,string> position2contig;
	map<string,unsigned int> contig2position;
	vector<int> contigLength;
	uint32_t contigsNumber;
	string bamFile;
	string output;

	//constructor
	Cov(int binSize,string bamFile,string output);
	//module used to calculate the coverage of the genome
	void bin(BamAlignment currentRead);
	void printCoverage();
};


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
	int readLength;
	bool mp;
};






static readStatus computeReadType(BamAlignment al, uint32_t max_insert, uint32_t min_insert,bool is_mp) {
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
	float y = 0.5 * erfc(float(-t / (sigma * sqrt(2.0))));
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

static LibraryStatistics computeLibraryStats(string bamFileName, uint64_t genomeLength, uint32_t max_insert, uint32_t min_insert,bool is_mp,int quality,string outputFileHeader) {
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
	// split reads
	uint32_t splitReads = 0;
	uint32_t readLength=0;

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
	
    //initiate the coverage computation class
    Cov *calculateCoverage;
	calculateCoverage = new Cov(300,bamFileName,outputFileHeader);
	BamAlignment al;
	

	while ( bamFile.GetNextAlignmentCore(al) ) {
		reads ++;
		readStatus read_status = computeReadType(al, max_insert, min_insert,is_mp);
		if (read_status != unmapped and read_status != lowQualty) {
			mappedReads ++;
			mappedReadsLength += al.Length;
		}
		
		vector <int > clipSizes;
		vector< int > readPositions;
		vector<int> genomePos;
		if ( al.GetSoftClips(clipSizes,readPositions,genomePos) ){
			splitReads += 1;
		}

        //calculate the coverage in bins of size 400 bases
        calculateCoverage -> bin(al);
        
		if (al.IsFirstMate() && al.IsMateMapped()) {
			if( al.IsReverseStrand() != al.IsMateReverseStrand() ){
				if(al.IsMapped() and al.MapQuality > quality and al.RefID == al.MateRefID and al.MatePosition-al.Position+1 < max_insert ){
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
			}
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
	cout << "\t split reads "               << splitReads << "\n";
	cout << "\t low quality reads " 	<< lowQualityReads << "\n";
	cout << "\t wrongly contig "		<< matedDifferentContig << "\n";
	cout << "\t singletons " 		<< singletonReads << "\n";

	uint32_t total = matedReads + wrongDistanceReads +  wronglyOrientedReads +  matedDifferentContig + singletonReads  ;
	cout << "\ttotal " << total << "\n";
	cout << "\tCoverage statistics\n";

	library.C_A = C_A = mappedReadsLength/(float)genomeLength;
	library.S_A = S_A = insertsLength/(float)genomeLength;
	library.C_M = C_M = matedReadsLength/(float)genomeLength;
	library.C_W = C_W = wronglyOrientedReadsLength/(float)genomeLength;
	library.C_S = C_S = singletonReadsLength/(float)genomeLength;
	library.C_D = C_D = matedDifferentContigLength/(float)genomeLength;
	library.readLength= mappedReadsLength/mappedReads;
	library.insertMean = insertMean = Mk;
	if(reads-wronglyOrientedReads > wronglyOrientedReads){
		library.mp=true;
	}else{
		library.mp=false;
	}
	Qk = sqrt(Qk/counterK);
	library.insertStd = insertStd = Qk;

	cout << "\tC_A = " << C_A << endl;
	cout << "\tS_A = " << S_A << endl;
	cout << "\tC_M = " << C_M << endl;
	cout << "\tC_W = " << C_W << endl;
	cout << "\tC_S = " << C_S << endl;
	cout << "\tC_D = " << C_D << endl;
	cout << "\tRead length = " << library.readLength << endl;
	cout << "\tMean Insert length = " << Mk << endl;
	cout << "\tStd Insert length = " << Qk << endl;
	cout << "----------\n";

	bamFile.Close();
	calculateCoverage -> printCoverage();
	return library;
}

#endif /*TYPES_H_*/
