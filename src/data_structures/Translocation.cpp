/*
 * Translocations.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: vezzi, Eisfeldt
 */

#include "Translocation.h"
#include <string>
#include <cmath>  

string int2str(int to_be_converted){
	string converted= static_cast<ostringstream*>( &(ostringstream() << to_be_converted) )->str();
	return(converted);
}

void Window::initTrans(SamHeader head) {
	uint32_t contigsNumber = 0;
	SamSequenceDictionary sequences  = head.Sequences;
	for(SamSequenceIterator sequence = sequences.Begin() ; sequence != sequences.End(); ++sequence) {
		this->contig2position[sequence->Name] = contigsNumber; // keep track of contig name and position in order to avoid problems when processing two libraries
		this->position2contig[contigsNumber]  = sequence->Name;
		this-> contig_ids.push_back(sequence->Name);
		this-> contig_length.push_back(sequence->Length);
		if(sequence->HasAssemblyID() == true){
			this-> contig_assembly.push_back(sequence->AssemblyID);
		}
		contigsNumber++;
	}

}

void Window::printHeader(SamHeader head,string libraryData) {

	string intra_chr_eventsVCF = outputFileHeader + ".signals.tab";
	this->TIDDITVCF.open(intra_chr_eventsVCF.c_str());
	if(head.HasReadGroups() == false){
		this->TIDDITVCF << VCFHeader(libraryData) << outputFileHeader << "\n";
	}else{
		SamReadGroupDictionary readGroups = head.ReadGroups;
		string SM = outputFileHeader;
		for(SamReadGroupIterator readGroup = readGroups.Begin() ; readGroup != readGroups.End(); ++readGroup) {
			if( readGroup -> HasSample() == true){
				SM= readGroup -> Sample;
				break;
			}
		}
		this->TIDDITVCF << VCFHeader(libraryData) << SM << "\n";
	}

}

void Window::insertRead(BamAlignment alignment,readStatus alignmentStatus) {

	if( not alignment.IsMateMapped()  or alignmentStatus == lowQualty) {
		return; // in case the alignment is of no use discard it
	}

	if(this->chr == -1) { //sets chr to the first chromosome at startup
		cout << "working on sequence " << position2contig[alignment.RefID] << "\n";
		chr=alignment.RefID;
	}

	if(alignment.RefID != chr) { // I am moving to a new chromosomes, need to check if the current window can be used or not
		cout << "working on seqence " << position2contig[alignment.RefID] << "\n";
	}
	
	bool alignment_split = false;	
	alignment.BuildCharData();			
	alignment_split = alignment.HasTag("SA");
	if (alignment_split and alignment.IsPrimaryAlignment() and not alignment.IsSupplementaryAlignment() and alignment.MapQuality >= minimum_mapping_quality) {
		// parse split read to get the other segment position, akin to a mate.
		string SA;
		alignment.GetTag("SA",SA);
		
		stringstream ss(SA);
		std::string item;
		while (std::getline(ss, item, ';')) {
			stringstream SS(item);
			string SA_data;
			vector <string> SA_elements;
			while (std::getline(SS, SA_data, ',')) {
				SA_elements.push_back(SA_data);
			}

			int contigNr = contig2position[SA_elements[0]];
			int currrentAlignmentPos=alignment.Position;
			int discordantDistance=0;
			int splitDistance = 0;

			long splitPos = atol(SA_elements[1].c_str());
			string cigar_SA = SA_elements[3];
			string cigar ="";
			// iterate over cigar operations
			vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
 			vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
			std::stringstream sscigar;
 			for ( ; cigarIter != cigarEnd; ++cigarIter) {
				const CigarOp& op = (*cigarIter);
				sscigar << op.Length << op.Type;
			}

			cigar=sscigar.str();

			std::stringstream ss;

			string orientationA="-";
			if( alignment.IsReverseStrand() == false  ){
				orientationA="+";
			}
			//chrA posA orientation cigar q chrB endB orientation cigar qualB resolution
			if(alignment.RefID == contigNr and alignment.Position < splitPos){
				ss << alignment.Name << "\t"  << position2contig[alignment.RefID] << "\t" << alignment.Position +1  << "\t" << orientationA << "\t" << cigar << "\t" << alignment.MapQuality << "\t" << SA_elements[0] << "\t" << splitPos+1 << "\t"<< SA_elements[2] << "\t" << cigar_SA << "\t" << SA_elements[4] <<  "\t" << 1 << "\n";
			}else if(alignment.RefID < contigNr){
				ss << alignment.Name << "\t" << position2contig[alignment.RefID] << "\t" << alignment.Position +1 << "\t" << orientationA << "\t" << cigar << "\t" << alignment.MapQuality << "\t" << SA_elements[0] << "\t" << splitPos+1 << "\t"<< SA_elements[2] << "\t" << cigar_SA << "\t" << SA_elements[4] << "\t" << 1 <<"\n";
			}else{
				ss << alignment.Name << "\t" << SA_elements[0] << "\t" << splitPos+1 << "\t"<< SA_elements[2] << "\t"<< cigar_SA << "\t" << SA_elements[4] << "\t" << position2contig[alignment.RefID] << "\t" << alignment.Position +1 << "\t" << orientationA << "\t" << cigar << "\t" << alignment.MapQuality << "\t" << 1 << "\n";
			}

			SV_calls[alignment.MateRefID].push_back(ss.str());
		}
	}
	
	if( alignment.IsPrimaryAlignment() ){
		if(alignmentStatus == pair_wrongChrs or alignmentStatus ==  pair_wrongDistance) {
			if(alignment.RefID < alignment.MateRefID or (alignment.RefID == alignment.MateRefID and alignment.Position < alignment.MatePosition)) {  // insert only "forward" variations
				if (alignment.MapQuality >= minimum_mapping_quality){
					string orientationA="-";
					if( alignment.IsReverseStrand() == false  ){
						orientationA="+";
					}

					string orientationB="-";
					if(alignment.IsMateReverseStrand() == false  ){
						orientationB="+";
					}

					string cigar ="";
					// iterate over cigar operations
					vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
	 				vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
					std::stringstream sscigar;
	 				for ( ; cigarIter != cigarEnd; ++cigarIter) {
						const CigarOp& op = (*cigarIter);
						sscigar << op.Length << op.Type;
					}

					cigar=sscigar.str();
					std::stringstream ss;

					ss << alignment.Name << "\t" << position2contig[alignment.RefID] << "\t" << alignment.Position +1 << "\t" << orientationA << "\t" << cigar << "\t" << alignment.MapQuality << "\t" << position2contig[alignment.MateRefID] << "\t" << alignment.MatePosition+1 << "\t"<< orientationB << "\t" << "NA" << "\t" << -1 << "\t" << 100 << "\n";
					SV_calls_discordant[alignment.Name]=ss.str();
				}
			}else if (alignment.MapQuality < minimum_mapping_quality){
				if(SV_calls_discordant.count(alignment.Name) == 1){
					SV_calls_discordant.erase(alignment.Name);
				}
			}
		}
	}
	chr=alignment.RefID;
}

//Function for printing the vcf header
string Window::VCFHeader(string libraryData){
	string headerString ="";
	//Define fileformat and source
	headerString+="##fileformat=VCFv4.1\n";
	headerString+="##source=TIDDIT-" + this-> version +  "\n";
	//define the alowed events
	headerString+="##ALT=<ID=DEL,Description=\"Deletion\">\n";
	headerString+="##ALT=<ID=DUP,Description=\"Duplication\">\n";
	headerString+="##ALT=<ID=TDUP,Description=\"Tandem duplication\">\n";
	headerString+="##ALT=<ID=INV,Description=\"Inversion\">\n";
	headerString+="##ALT=<ID=INS,Description=\"Insertion\">\n";
	headerString+="##ALT=<ID=BND,Description=\"Break end\">\n";

	for(int i=0;i < this -> contig_ids.size(); i++){
		headerString+= "##contig=<ID="	+ contig_ids[i] + ",length=" + contig_length[i];
		headerString +=">\n";
	}
	
	//Define the info fields
	headerString+="##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
	headerString+="##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of an intra-chromosomal variant\">\n";
	headerString+="##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
	headerString+="##INFO=<ID=LFA,Number=1,Type=Integer,Description=\"Links from window A\">\n";
	headerString+="##INFO=<ID=LFB,Number=1,Type=Integer,Description=\"Links from window B\">\n";
	headerString+="##INFO=<ID=LTE,Number=1,Type=Integer,Description=\"Links to event\">\n";
	headerString+="##INFO=<ID=COVA,Number=1,Type=Float,Description=\"Coverage on window A\">\n";
	headerString+="##INFO=<ID=COVM,Number=1,Type=Float,Description=\"The coverage between A and B\">\n";
	headerString+="##INFO=<ID=COVB,Number=1,Type=Float,Description=\"Coverage on window B\">\n";
	headerString+="##INFO=<ID=OR,Number=4,Type=Integer,Description=\"Orientation of the pairs (FF,FR,RF,FR)\">\n";
	headerString+="##INFO=<ID=ORSR,Number=2,Type=Integer,Description=\"Orientation of the split reads (inverted,normal)\">\n";
	headerString+="##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n";
	headerString+="##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n";
	headerString+="##INFO=<ID=E1,Number=1,Type=Float,Description=\"Expected links to window B\">\n";
	headerString+="##INFO=<ID=E2,Number=1,Type=Float,Description=\"Expected links to window B(assuming ideal window length)\">\n";
	headerString+="##INFO=<ID=QUALA,Number=1,Type=Float,Description=\"The average mapping quality of the reads in window A\">\n";
	headerString+="##INFO=<ID=QUALB,Number=1,Type=Float,Description=\"The average mapping quality of the reads in window B\">\n";
	//set filters
	headerString+="##FILTER=<ID=BelowExpectedLinks,Description=\"The number of links or reads between A and B is less than 40\% of the expected value\">\n";
	headerString+="##FILTER=<ID=FewLinks,Description=\"Unexpectedly low fraction of discordant reads betwen A and B\">\n";
	headerString+="##FILTER=<ID=UnexpectedCoverage,Description=\"The coverage of the window on chromosome B or A is higher than 4*average coverage\">\n";
	headerString+="##FILTER=<ID=Smear,Description=\"window A and Window B overlap\">\n";
	headerString+="##FILTER=<ID=RegionalQ,Description=\"The mapping quality of the region is lower than the user set limit\">\n";
	headerString+="##FILTER=<ID=MinSize,Description=\"The variant is smaller than the user set limit\">\n";
	headerString+="##FILTER=<ID=Ploidy,Description=\"Intrachromosomal variant on a chromosome having 0 ploidy\">\n";
	headerString+="##FILTER=<ID=SplitsVSDiscs,Description=\"large variant supported mainly by split reads (and not discorant pairs) \">\n";
	headerString+="##FILTER=<ID=Density,Description=\"The discordant reads cluster too tightly\">\n";
	//set format 
	headerString+="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	headerString+="##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n";
	headerString+="##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of paired-ends that support the event\">\n";
	headerString+="##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"Number of split reads that support the event\">\n";
	headerString+="##FORMAT=<ID=DR,Number=2,Type=Integer,Description=\"Number of paired-ends that supporting the reference allele (breakpoint A, and B)\">\n";
	headerString+="##FORMAT=<ID=RR,Number=2,Type=Integer,Description=\"Number of reads supporting the reference allele (breakpoint A, and B)\">\n";
	//this sting contains info such as insert size and mean coverage
	headerString+=libraryData+"\n";
	//Header
	headerString+="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	return(headerString);
}

Window::Window(string bamFileName, bool outtie, float meanCoverage,string outputFileHeader, map<string,int> SV_options) {
	this->max_insert		 = SV_options["max_insert"];
	this->minimum_mapping_quality = SV_options["mapping_quality"];
	this->outtie			 = outtie;
	this->mean_insert		 = SV_options["meanInsert"];
	this ->std_insert		 = SV_options["STDInsert"];
	this->bamFileName		 = bamFileName;
	this -> ploidy           = SV_options["ploidy"];
	this -> readLength       = SV_options["readLength"];
	this -> pairOrientation				 = 0;          
	this->outputFileHeader     = outputFileHeader;
}
