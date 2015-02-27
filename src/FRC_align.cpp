/*
    FRC: computes the FRC curve starting from alignments
    Copyright (C) 2011  F. Vezzi(vezi84@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h> 
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <string>

#include <sstream>
#include <iostream>
#include <fstream>

#include  <boost/program_options.hpp>
namespace po = boost::program_options;


#include "data_structures/Features.h"
#include "data_structures/FRC.h"

#include "common.h"

//LibraryStatistics computeLibraryStats(string bamFileName, uint64_t estimatedGenomeSize, uint32_t max_insert, bool is_mp);
void computeFRC(FRC &  frc, string bamFileName, LibraryStatistics library,int max_insert, bool is_mp, float CE_min, float CE_max);
void printFRCurve(string outputFile, int totalFeatNum, FeatureTypes type, uint64_t estimatedGenomeSize, FRC frc);


int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string PEalignmentFile = "";
	string MPalignmentFile = "";
	uint32_t max_pe_insert = 5000;
	uint32_t max_mp_insert = 20000;
	uint64_t estimatedGenomeSize;
	float CEstats_PE_min = -5;
	float CEstats_PE_max = +5;
	float CEstats_MP_min = -7;
	float CEstats_MP_max = +7;
	string outputFile = "FRC.txt";
	string featureFile = "Features.txt";

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());
	desc.add_options() ("help", "produce help message")
	("pe-sam"       , po::value<string>(), "paired end alignment file (in sam or bam format). Orientation must be -> <-")
	("pe-max-insert", po::value<int>()   , "maximum allowed insert size for PE (to filter out outleyers)")
	("mp-sam"       , po::value<string>(), "mate pairs alignment file. (in sam or bam format). Orientation must be -> <-")
	("mp-max-insert", po::value<int>()   , "maximum allowed insert size for MP (to filter out outleyers)")
	("genome-size"  , po::value<unsigned long int>(), "estimated genome size (if not supplied genome size is believed to be assembly length")
	("output"       ,  po::value<string>(), "Header output file names (default FRC.txt and Features.txt)")
	("CEstats-PE-min", po::value<float>() , "minimum allowed CE_stats in PE library")
	("CEstats-PE-max", po::value<float>() , "maximum allowed CE_stats in PE library")
	("CEstats-MP-min", po::value<float>() , "minimum allowed CE_stats in MP library")
	("CEstats-MP-max", po::value<float>() , "maximum allowed CE_stats in MP library")
	;

	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (boost::program_options::error & error) {
		ERROR_CHANNEL <<  error.what() << endl;
		ERROR_CHANNEL << "Try \"--help\" for help" << endl;
		exit(2);
	}
	if (vm.count("help")) {
		DEFAULT_CHANNEL << desc << endl;
		exit(0);
	}

	//PARSE CE STATS (if present)
	if (vm.count("CEstats-PE-min")) {
		CEstats_PE_min = vm["CEstats-PE-min"].as<float>();
	}
	if (vm.count("CEstats-PE-max")) {
		CEstats_PE_max = vm["CEstats-PE-max"].as<float>();
	}
	if (vm.count("CEstats-MP-min")) {
		CEstats_MP_min = vm["CEstats-MP-min"].as<float>();
	}
	if (vm.count("CEstats-MP-max")) {
		CEstats_MP_max = vm["CEstats-MP-max"].as<float>();
	}

	cout << "CEstats-PE-min " << CEstats_PE_min << "\n";
	cout << "CEstats-PE-max " << CEstats_PE_max << "\n";
	cout << "CEstats-MP-min " << CEstats_MP_min << "\n";
	cout << "CEstats-MP-max " << CEstats_MP_max << "\n";

	// PARSE PE
	if (!vm.count("pe-sam") && !vm.count("mp-sam")) {
		DEFAULT_CHANNEL << "At least one library must be present. Please specify at least one between pe-sam and mp-sam" << endl;
		exit(0);
	}

	if (vm.count("pe-sam")) {
		PEalignmentFile = vm["pe-sam"].as<string>();
	}
	if (vm.count("max-pe-insert")) {
		max_pe_insert = vm["max-pe-insert"].as<int>();
	}

	// NOW PARSE MP
	if (vm.count("mp-sam")) {
		MPalignmentFile = vm["mp-sam"].as<string>();
	}

	if (vm.count("mp-max-insert")) {
		max_mp_insert = vm["mp-max-insert"].as<int>();
	}
	if(vm.count("pe-sam")){
		cout << "pe-sam file name is " << PEalignmentFile << endl;
	}
	if(vm.count("mp-sam")){
		cout << "mp-sam file name is " << MPalignmentFile << endl;
	}


	string header = "";
	if (vm.count("output")) {
		header = vm["output"].as<string>();
		outputFile = header + "_FRC.txt";
		featureFile = header + "_Features.txt";
	}

	if (vm.count("genome-size")) {
		estimatedGenomeSize = vm["genome-size"].as<unsigned long int>();
	} else {
		estimatedGenomeSize = 0;
	}

	//TODO: PARSING ENDED, CREATE A FUNCTION FOR IT
	uint64_t genomeLength = 0;
	uint32_t contigsNumber = 0;
	BamReader bamFile;
	if(vm.count("pe-sam")) { // paired read library is preset, use it to compute basic contig statistics
		bamFile.Open(PEalignmentFile);
	} else {  // Otherwise use MP library that must be provided
		bamFile.Open(PEalignmentFile);
	}



	SamHeader head = bamFile.GetHeader();
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;

	SamSequenceDictionary sequences  = head.Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		genomeLength += StringToNumber(sequence->Length);
		contig2position[sequence->Name] = contigsNumber; // keep track of contig name and position in order to avoid problems when processing two libraries
		position2contig[contigsNumber] = contig2position[sequence->Name];
		contigsNumber++;
	}
	bamFile.Close();

	if (estimatedGenomeSize == 0) {
		estimatedGenomeSize =  genomeLength;
	}
	cout << "total number of contigs " 	<< contigsNumber << endl;
	cout << "assembly length " 			<< genomeLength << "\n";
	cout << "estimated length " 		<< estimatedGenomeSize << "\n";

	LibraryStatistics libraryPE;
	LibraryStatistics libraryMP;
	uint32_t		  peInsertSize;
	uint32_t 		  peStdDeviation;
	uint32_t 		  mpInsertSize;
	uint32_t 		  mpStdDeviation;
	unsigned int      timesStdDev = 3;

	if(vm.count("pe-sam")) { // in this case file is already OPEN
		cout << "computing statistics for PE library\n";
		libraryPE = computeLibraryStats(PEalignmentFile, estimatedGenomeSize, max_pe_insert, false);
	}

	if(vm.count("mp-sam")) {
		cout << "computing statistics for MP library\n";
		if(!vm.count("pe-sam")) { // in this case file is already OPEN
			libraryMP = computeLibraryStats(MPalignmentFile, estimatedGenomeSize, max_mp_insert, true);
		} else { // in this case I need to open file first
			libraryMP = computeLibraryStats(MPalignmentFile, estimatedGenomeSize, max_mp_insert, true);
		}
	}
// Libraries have been parsed, now I need to compute new maximum and minimum for pe and mp and compute features
//TODO: I need to recompute max and min...

	//parse BAM file again to compute FRC curve
	FRC frc = FRC(contigsNumber); // FRC object, will memorize all information on features and contigs
	uint32_t contigCounter = 0;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		uint32_t contigLength = StringToNumber(sequence->Length);
		frc.setContigLength(contigCounter, contigLength);
		frc.setID(contigCounter, sequence->Name);
		contigCounter++;
	}

	int featuresTotal   = 0;
	int featuresTotalPE = 0;
	int featuresTotalMP = 0;

	if(vm.count("pe-sam")) { // in this case file is already OPEN
		cout << "computing FRC for PE library\n";
		computeFRC(frc, PEalignmentFile, libraryPE, max_pe_insert, false, CEstats_PE_min , CEstats_PE_max);
		string PE_CEstats = header + "_CEstats_PE.txt";
		ofstream CEstats;
		CEstats.open(PE_CEstats.c_str());
		map<float, unsigned int>::iterator it;
		map<float, unsigned int>::iterator secondIterator;
		for ( it = frc.CEstatistics.begin() ; it != frc.CEstatistics.end(); it++ ) {
			unsigned int total = 0;
			if( (*it).first < 0) {
				for(secondIterator = it; secondIterator != frc.CEstatistics.begin(); secondIterator --) {
					total += (*secondIterator).second;
				}
			} else {
				for(secondIterator = it; secondIterator != frc.CEstatistics.end(); secondIterator ++) {
					total += (*secondIterator).second;
				}
			}
			CEstats << (*it).first << " " << total << endl;
		}
		CEstats.close();
		frc.CEstatistics.clear();
		for(unsigned int i=0; i< contigsNumber; i++) {
			featuresTotal   += frc.getTotal(i);
			featuresTotalPE += frc.getTotal(i);
		}
		cout << "TOTAL number of PE features " << featuresTotalPE << "\n";
	}
	//NOW MP
	if(vm.count("mp-sam")) {
		cout << "computing FRC for MP library\n";
		computeFRC(frc, MPalignmentFile, libraryMP, max_mp_insert, true, CEstats_MP_min , CEstats_MP_max);
		string MP_CEstats = header + "_CEstats_MP.txt";
		ofstream CEstats;
		CEstats.open(MP_CEstats.c_str());
		map<float, unsigned int>::iterator it;
		map<float, unsigned int>::iterator secondIterator;
		for ( it= frc.CEstatistics.begin() ; it != frc.CEstatistics.end(); it++ ) {
			unsigned int total = 0;
			if( (*it).first <= 0) {
				for(secondIterator = it; secondIterator != frc.CEstatistics.begin(); secondIterator --) {
					total += (*secondIterator).second;
				}
				CEstats << (*it).first << " " << total << endl;
			} else {
				for(secondIterator = it; secondIterator != frc.CEstatistics.end(); secondIterator ++) {
					total += (*secondIterator).second;
				}
				CEstats << (*it).first << " " << total << endl;
			}
		}
		CEstats.close();
		frc.CEstatistics.clear();
		for(unsigned int i=0; i< contigsNumber; i++) {
			featuresTotal   += frc.getTotal(i);
			featuresTotalMP += frc.getTotal(i);
		}
		cout << "TOTAL number of MP features " << featuresTotalMP << "\n";
	}
	//all features have now been computed
	 cout << "TOTAL number of features " << featuresTotal << "\n\n";

	cout << "\n----------\nNow computing FRC \n------\n";
//print all features
	ofstream featureOutFile;
	featureOutFile.open (featureFile.c_str());
	ofstream GFF3_features;
	string GFF3 = header + "Features.gff";
	GFF3_features.open(GFF3.c_str());
	GFF3_features << "##gff-version   3\n";
    for(unsigned int i=0; i< contigsNumber; i++) {
    	frc.printFeatures(i, featureOutFile);
    	frc.printFeaturesGFF3(i, GFF3_features);
    }

    frc.sortFRC();

    //NOW COMPUTE ALL THE FRCurves
    unsigned int LOW_COV_PE_features = 0;
    unsigned int HIGH_COV_PE_features = 0;
    unsigned int LOW_NORM_COV_PE_features = 0;
    unsigned int HIGH_NORM_COV_PE_features = 0;
    unsigned int HIGH_SINGLE_PE_features = 0;
    unsigned int HIGH_SPAN_PE_features = 0;
    unsigned int HIGH_OUTIE_PE_features = 0;
    unsigned int COMPR_PE_features = 0;
    unsigned int STRECH_PE_features = 0;

    unsigned int HIGH_SINGLE_MP_features = 0;
    unsigned int HIGH_OUTIE_MP_features = 0;
    unsigned int HIGH_SPAN_MP_features = 0;
    unsigned int COMPR_MP_features = 0;
    unsigned int STRECH_MP_features = 0;

    for(unsigned int i=0; i< contigsNumber; i++) {
    	featuresTotal += frc.getTotal(i); // update total number of feature seen so far
    	LOW_COV_PE_features += frc.getLOW_COV_PE(i);
    	HIGH_COV_PE_features += frc.getHIGH_COV_PE(i);
    	LOW_NORM_COV_PE_features += frc.getLOW_NORM_COV_PE(i);
    	HIGH_NORM_COV_PE_features += frc.getHIGH_NORM_COV_PE(i);
    	HIGH_SINGLE_PE_features += frc.getHIGH_SINGLE_PE(i);
    	HIGH_SPAN_PE_features += frc.getHIGH_SPAN_PE(i);
    	HIGH_OUTIE_PE_features += frc.getHIGH_OUTIE_PE(i);
    	COMPR_PE_features += frc.getCOMPR_PE(i);
    	STRECH_PE_features += frc.getSTRECH_PE(i);

    	HIGH_SINGLE_MP_features += frc.getHIGH_SINGLE_MP(i);
    	HIGH_OUTIE_MP_features += frc.getHIGH_OUTIE_MP(i);
    	HIGH_SPAN_MP_features += frc.getHIGH_SPAN_MP(i);
    	COMPR_MP_features += frc.getCOMPR_MP(i);
    	STRECH_MP_features += frc.getSTRECH_MP(i);
    }

    cout << "Total amount of features after merging " << featuresTotal << "\n";
    printFRCurve(outputFile, featuresTotal, FRC_TOTAL, estimatedGenomeSize, frc);
    //now all the others
    outputFile = header + "LOW_COV_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, LOW_COV_PE, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_COV_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_COV_PE, estimatedGenomeSize, frc);

    outputFile = header + "LOW_NORM_COV_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, LOW_NORM_COV_PE, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_NORM_COV_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_NORM_COV_PE, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_SINGLE_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_SINGLE_PE, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_OUTIE_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_OUTIE_PE, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_SPAN_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_SPAN_PE, estimatedGenomeSize, frc);

    outputFile = header + "COMPR_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, COMPR_PE, estimatedGenomeSize, frc);

    outputFile = header + "STRECH_PE_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, STRECH_PE, estimatedGenomeSize, frc);


    outputFile = header + "HIGH_SINGLE_MP_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_SINGLE_MP, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_OUTIE_MP_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_OUTIE_MP, estimatedGenomeSize, frc);

    outputFile = header + "HIGH_SPAN_MP_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, HIGH_SPAN_MP, estimatedGenomeSize, frc);

    outputFile = header + "COMPR_MP_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, COMPR_MP, estimatedGenomeSize, frc);

    outputFile = header + "STRECH_MP_FRC.txt";
    printFRCurve(outputFile, LOW_COV_PE_features, STRECH_MP, estimatedGenomeSize, frc);


    return 0;
}



int getFeatures(FRC frc, FeatureTypes type, int contig) {
	switch (type) {
		 case FRC_TOTAL:
			 return frc.getTotal(contig);
			 break;
		 case LOW_COV_PE:
			 return frc.getLOW_COV_PE(contig);
			 break;
		 case HIGH_COV_PE:
			 return frc.getHIGH_COV_PE(contig);
			 break;
		 case LOW_NORM_COV_PE:
			 return frc.getLOW_NORM_COV_PE(contig);
			 break;
		 case HIGH_NORM_COV_PE:
			 return frc.getHIGH_NORM_COV_PE(contig);
			 break;
		 case HIGH_SINGLE_PE:
			 return frc.getHIGH_SINGLE_PE(contig);
			 break;
		 case HIGH_SPAN_PE:
			 return frc.getHIGH_SPAN_PE(contig);
			 break;
		 case HIGH_OUTIE_PE:
			 return frc.getHIGH_OUTIE_PE(contig);
			 break;
		 case COMPR_PE:
			 return frc.getCOMPR_PE(contig);
			 break;
		 case STRECH_PE:
			 return frc.getSTRECH_PE(contig);
			 break;
		 case HIGH_SINGLE_MP:
			 return frc.getHIGH_SINGLE_MP(contig);
			 break;
		 case HIGH_OUTIE_MP:
			 return frc.getHIGH_OUTIE_MP(contig);
			 break;
		 case HIGH_SPAN_MP:
			 return frc.getHIGH_SPAN_MP(contig);
			 break;
		 case COMPR_MP:
			 return frc.getCOMPR_MP(contig);
			 break;
		 case STRECH_MP:
			 return frc.getSTRECH_MP(contig);
			 break;
		 default:
			 cout << "THis whould never happen\n";
	}
}


void printFRCurve(string outputFile, int totalFeatNum, FeatureTypes type, uint64_t estimatedGenomeSize, FRC frc){
	ofstream myfile;
	myfile.open (outputFile.c_str());

	cout << "now computing " << returnFeatureName(type) << "\t";
	if (totalFeatNum == 0 ) {
		myfile << "0 100\n";
		cout << "No features of this kind: DONE\n";
		return;
	}
	float step = totalFeatNum/(float)100;
	float partial=0;
	uint64_t edgeCoverage = 0;

	while(partial <= totalFeatNum) {
		uint32_t featuresStep = 0;
		uint32_t contigStep    = 0;
		featuresStep += frc.getFeatures(type, contigStep);

		while(featuresStep <= partial) {
			contigStep++;
			if(contigStep < frc.returnContigs()) {
				featuresStep +=  frc.getFeatures(type, contigStep);
			} else {
				featuresStep = partial + 1; // I read all the contigs, time to to stop
			}
		}

		edgeCoverage = 0;
		for(unsigned int i=0; i< contigStep; i++) {
			edgeCoverage += frc.getContigLength(i);
		}
		float coveragePartial =  100*(edgeCoverage/(float)estimatedGenomeSize);
		myfile << partial << " " << coveragePartial << "\n";
		partial += step;
		cout << ".";

		if(partial >= totalFeatNum) {
			partial = totalFeatNum + 1;
		}

		if(contigStep >= frc.returnContigs()) {
			partial = totalFeatNum + 1;
		}
	}
	cout << "\n";
	myfile.close();

}




void computeFRC(FRC & frc, string bamFileName, LibraryStatistics library,int max_insert, bool is_mp, float CE_min, float CE_max) {
	frc.setC_A(library.C_A);
	frc.setS_A(library.S_A);
	frc.setC_D(library.C_D);
	frc.setC_M(library.C_M);
	frc.setC_S(library.C_S);
	frc.setC_W(library.C_W);
	frc.setInsertMean(library.insertMean);
	frc.setInsertStd(library.insertStd);

	unsigned int windowStepCE = library.insertMean;

	BamReader bamFile;
	bamFile.Open(bamFileName);
	SamHeader head = bamFile.GetHeader(); // get the sam header
	BamAlignment al;
	int currentContig 	= -1;
	uint32_t contigSize = 0;
	Contig *contig;
	while ( bamFile.GetNextAlignmentCore(al) ) {
		if (al.IsMapped()) {
			if (al.RefID != currentContig) { // another contig or simply the first one
				//cout << "now porcessing contig " << contig << "\n";
				if(currentContig == -1) { // first read that I`m processing
					contigSize 		= frc.getContigLength(al.RefID) ;
					currentContig 	= al.RefID;
					contig =  new Contig(contigSize);
				} else {
					float coverage = frc.obtainCoverage(currentContig, contig);
					frc.computeCEstats(contig, library.insertMean, windowStepCE, library.insertMean, library.insertStd);
					//frc.computeCEstats(contig, 1000, 200, library.insertMean, library.insertStd);
					if(is_mp) {
						//frc.computeLowCoverageArea("MP", currentContig, contig, 1000, 200);
						//frc.computeHighCoverageArea("MP", currentContig, contig, 1000, 200);
						//frc.computeLowNormalArea("MP", currentContig, contig, 1000, 200);
						//frc.computeHighNormalArea("MP", currentContig, contig, 1000, 200);
						frc.computeHighSingleArea("MP", currentContig, contig, 1000, 200);
						frc.computeHighOutieArea("MP", currentContig, contig, 1000, 200);
						frc.computeHighSpanningArea("MP", currentContig, contig, 1000, 200);
						frc.computeCompressionArea("MP", currentContig, contig, CE_min, library.insertMean, library.insertMean);
						frc.computeStrechArea("MP", currentContig, contig, CE_max, library.insertMean, library.insertMean);
					} else {
						frc.computeLowCoverageArea("PE", currentContig, contig, 1000, 200);
						frc.computeHighCoverageArea("PE", currentContig, contig, 1000, 200);
						frc.computeLowNormalArea("PE", currentContig, contig, 1000, 200);
						frc.computeHighNormalArea("PE", currentContig, contig, 1000, 200);
						frc.computeHighSingleArea("PE", currentContig, contig, 1000, 200);
						frc.computeHighOutieArea("PE", currentContig, contig, 1000, 200);
						frc.computeHighSpanningArea("PE", currentContig, contig, 1000, 200);
						frc.computeCompressionArea("PE", currentContig, contig, CE_min, library.insertMean, library.insertMean);
						frc.computeStrechArea("PE", currentContig, contig, CE_max, library.insertMean, library.insertMean);
					}

					delete contig; // delete hold contig
					contigSize = frc.getContigLength(al.RefID) ;
					if (contigSize < 1) {//We can't have such sizes! this can't be right
						fprintf(stderr,"%d has size %d, which can't be right!\nCheck bam header!",al.RefID,contigSize);
					}
					currentContig 	= al.RefID; // update current identifier
					contig 			= new Contig(contigSize);
				}
				contig->updateContig(al, max_insert, is_mp); // update contig with alignment
			} else {
				//add information to current contig
				contig->updateContig(al, max_insert, is_mp);
			}
		}
	}
	//Last contig needs to be processed (I finished o read the file without parsing it)
	float coverage = frc.obtainCoverage(currentContig, contig);
	frc.computeCEstats(contig, library.insertMean, windowStepCE, library.insertMean, library.insertStd);
	if(is_mp) {
		//frc.computeLowCoverageArea("MP", currentContig, contig, 1000, 200);
		//frc.computeHighCoverageArea("MP", currentContig, contig, 1000, 200);
		//frc.computeLowNormalArea("MP", currentContig, contig, 1000, 200);
		//frc.computeHighNormalArea("MP", currentContig, contig, 1000, 200);
		frc.computeHighSingleArea("MP", currentContig, contig, 1000, 200);
		frc.computeHighOutieArea("MP", currentContig, contig, 1000, 200);
		frc.computeHighSpanningArea("MP", currentContig, contig, 1000, 200);
		frc.computeCompressionArea("MP", currentContig, contig, CE_min, library.insertMean, library.insertMean);
		frc.computeStrechArea("MP", currentContig, contig, CE_max, library.insertMean, library.insertMean);
	} else {
		frc.computeLowCoverageArea("PE", currentContig, contig, 1000, 200);
		frc.computeHighCoverageArea("PE", currentContig, contig, 1000, 200);
		frc.computeLowNormalArea("PE", currentContig, contig, 1000, 200);
		frc.computeHighNormalArea("PE", currentContig, contig, 1000, 200);
		frc.computeHighSingleArea("PE", currentContig, contig, 1000, 200);
		frc.computeHighOutieArea("PE", currentContig, contig, 1000, 200);
		frc.computeHighSpanningArea("PE", currentContig, contig, 1000, 200);
		frc.computeCompressionArea("PE", currentContig, contig, CE_min, library.insertMean, library.insertMean);
		frc.computeStrechArea("PE", currentContig, contig, CE_max, library.insertMean, library.insertMean);
	}


	delete contig; // delete hold contig
	bamFile.Close();

}










