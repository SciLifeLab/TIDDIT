/*
 * RunModules.h
 *
 *  Created on: Jul 1, 2015
 *      Author: vezzi, Eisfeldt
 */
//This header files contains the main functions of each program/run module of the find translocations software

#ifndef PROGRAMMODULES_H_
#define PROGRAMMODULES_H_

#include "common.h"
#include <list>

//The extraction module
class Extract{
public:
	//Constructor
	Extract();
	//Main function
	void extract(string BamFileName,string outputFileHeader,string inputFileName,map<string,unsigned int> contig2position, string indexFile);
};

//The module used to find structural variants
class StructuralVariations{
public:
	//constructor
	StructuralVariations();
	//main function
	void findTranslocationsOnTheFly(string bamFileName, int32_t min_insert,  int32_t max_insert, bool outtie, uint16_t minimum_mapping_quality,
		uint32_t minimumSupportingPairs, float meanCoverage, float meanInsertSize, float StdInsertSize, string outputFileHeader, string indexFile, int contigsNumber);
};

#endif /* PROGRAMMODULES_H_ */
