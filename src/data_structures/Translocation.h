/*
 * Translocations.h
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi, Eisfeldt
 */

#ifndef TRANSLOCATION_H_
#define TRANSLOCATION_H_

#include "common.h"


class Window {
public:
	//object fields that change as the window moves
	float coverage;
	int chr;

	//More static parts initialized by constructor
	int pairOrientation;
	int max_insert;
	int min_insert;
	uint16_t minimum_mapping_quality;
	bool outtie;
	float mean_insert;
	float std_insert;
	int ploidy;
	int readLength;
	//the file name of the bamfile
	string bamFileName;
	string version;

	//More static parts initialized by init method
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;
	vector<string> contig_ids;
	vector<string> contig_length;
	vector<string> contig_assembly;	
	map<unsigned int, vector<string> > SV_calls;
	map<string, string > SV_calls_discordant;

	//Output part
	string outputFileHeader;
	ofstream TIDDITVCF;

	Window(string bamFileName, bool outtie, float meanCoverage,string outputFileHeader, map<string,int> SV_options); // constructor

	void initTrans(SamHeader head);	// initialise the contig to position array
	void printHeader(SamHeader head,string libraryData); //print header
	void insertRead(BamAlignment alignment, readStatus alignmentStatus);			// inserts a new read
	string VCFHeader(string libraryData);				//print the vcf header

};

#endif /* TRANSLOCATIONS_H_ */
