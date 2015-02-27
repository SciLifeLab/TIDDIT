/*
 * Features.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */


#include "Features.h"

Features::Features() {
	LOW_COVERAGE_AREA = 0;
	HIGH_COVERAGE_AREA = 0;
	LOW_NORMAL_AREA = 0;
	HIGH_NORMAL_AREA = 0;
	HIGH_SINGLE_AREA = 0;
	HIGH_SPANNING_AREA = 0;
	HIGH_OUTIE_AREA = 0;
	COMPRESSION_AREA  = 0;
	STRECH_AREA = 0;
}

Features::~Features() {

}



void Features::setLOW_COVERAGE_AREA(unsigned int numFeat) {
	this->LOW_COVERAGE_AREA = numFeat;
}
void Features::setHIGH_COVERAGE_AREA(unsigned int numFeat) {
	this->HIGH_COVERAGE_AREA = numFeat;
}
void Features::setLOW_NORMAL_AREA(unsigned int numFeat) {
	this->LOW_NORMAL_AREA = numFeat;
}
void Features::setHIGH_NORMAL_AREA(unsigned int numFeat) {
	this->HIGH_NORMAL_AREA = numFeat;
}
void Features::setHIGH_SINGLE_AREA(unsigned int numFeat) {
	this->HIGH_SINGLE_AREA = numFeat;
}
void Features::setHIGH_SPANNING_AREA(unsigned int numFeat) {
	this->HIGH_SPANNING_AREA = numFeat;
}
void Features::setHIGH_OUTIE_AREA(unsigned int numFeat) {
	this->HIGH_OUTIE_AREA = numFeat;
}
void Features::setCOMPRESSION_AREA(unsigned int numFeat) {
	this->COMPRESSION_AREA = numFeat;
}
void Features::setSTRECH_AREA(unsigned int numFeat) {
	this->STRECH_AREA = numFeat;
}


void Features::updateLOW_COVERAGE_AREA(unsigned int numFeat) {
	this->LOW_COVERAGE_AREA += numFeat;
}
void Features::updateHIGH_COVERAGE_AREA(unsigned int numFeat) {
	this->HIGH_COVERAGE_AREA += numFeat;
}
void Features::updateLOW_NORMAL_AREA(unsigned int numFeat) {
	this->LOW_NORMAL_AREA += numFeat;
}
void Features::updateHIGH_NORMAL_AREA(unsigned int numFeat) {
	this->HIGH_NORMAL_AREA += numFeat;
}
void Features::updateHIGH_SINGLE_AREA(unsigned int numFeat) {
	this->HIGH_SINGLE_AREA += numFeat;
}
void Features::updateHIGH_SPANNING_AREA(unsigned int numFeat) {
	this->HIGH_SPANNING_AREA += numFeat;
}
void Features::updateHIGH_OUTIE_AREA(unsigned int numFeat) {
	this->HIGH_OUTIE_AREA += numFeat;
}
void Features::updateCOMPRESSION_AREA(unsigned int numFeat) {
	this->COMPRESSION_AREA += numFeat;
}
void Features::updateSTRECH_AREA(unsigned int numFeat) {
	this->STRECH_AREA += numFeat;
}


unsigned int Features::getLOW_COVERAGE_AREA() {return LOW_COVERAGE_AREA;}
unsigned int Features::getHIGH_COVERAGE_AREA() {return HIGH_COVERAGE_AREA;}
unsigned int Features::getLOW_NORMAL_AREA() {return LOW_NORMAL_AREA;}
unsigned int Features::getHIGH_NORMAL_AREA() {return HIGH_NORMAL_AREA;}
unsigned int Features::getHIGH_SINGLE_AREA() {return HIGH_SINGLE_AREA;}
unsigned int Features::getHIGH_SPANNING_AREA() {return HIGH_SPANNING_AREA;}
unsigned int Features::getHIGH_OUTIE_AREA() {return HIGH_OUTIE_AREA;}
unsigned int Features::getCOMPRESSION_AREA() {return COMPRESSION_AREA;}
unsigned int Features::getSTRECH_AREA() {return STRECH_AREA;}

unsigned int Features::returnTotal() {
	return COMPRESSION_AREA + HIGH_COVERAGE_AREA + HIGH_NORMAL_AREA + HIGH_OUTIE_AREA +
			HIGH_SINGLE_AREA + HIGH_SPANNING_AREA + LOW_COVERAGE_AREA + LOW_NORMAL_AREA + STRECH_AREA;

}


unsigned int Features::returnLOW_COV() {
	return LOW_COVERAGE_AREA;
}
unsigned int Features::returnHIGH_COV() {
	return HIGH_COVERAGE_AREA;
}
unsigned int Features::returnLOW_NORM_COV(){
	return LOW_NORMAL_AREA;
}
unsigned int Features::returnHIGH_NORM_COV(){
	return HIGH_NORMAL_AREA;
}
unsigned int Features::returnHIGH_SINGLE(){
	return HIGH_SINGLE_AREA;
}
unsigned int Features::returnHIGH_OUTIE(){
	return HIGH_OUTIE_AREA;
}
unsigned int Features::returnHIGH_SPAN(){
	return HIGH_SPANNING_AREA;
}
unsigned int Features::returnCOMPR(){
	return COMPRESSION_AREA;
}
unsigned int Features::returnSTRECH(){
	return STRECH_AREA;
}




