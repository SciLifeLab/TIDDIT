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

//The module used to find structural variants
class StructuralVariations{
public:
	//constructor
	StructuralVariations();
	//main function
	void findTranslocationsOnTheFly(string bamFileName, bool outtie, float meanCoverage,string outputFileHeader,string version,string commandline, map<string,int> SV_options);
};

#endif /* PROGRAMMODULES_H_ */
