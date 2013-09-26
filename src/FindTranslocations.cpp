/*
   Francesco Vezzi
 */

#include "common.h"
#include "samtools/sam.h"
#include "dataStructures/Translocations.h"


enum pairStatus {pair_Proper, pair_wrongOrientation, pair_wrongDistance, pair_wrongChrs, pair_wrongOptDup, pair_Deletion};


#define MIN(x,y) \
		((x) < (y)) ? (x) : (y)

#define EXIT_IF_NULL(P) \
		if (P == NULL) \
		return 1;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core) {
	if (core->flag&BAM_FUNMAP) {
		return false;
	}
	return true;
}


/**
 * Open a .sam/.bam file.
 * @returns NULL is open failed.
 */
samfile_t * open_alignment_file(std::string path) {
	samfile_t * fp = NULL;
	std::string flag = "r";
	if (path.substr(path.size()-3).compare("bam") == 0) {
		//BAM file!
		flag += "b";
	}
	if ((fp = samopen(path.c_str(), flag.c_str() , 0)) == 0) {
		fprintf(stderr, "Error: Failed to open file %s\n", path.c_str());
	}
	return fp;
}



void computeLibraryStats(samfile_t *fp, uint32_t MinInsert, uint32_t MaxInsert, uint64_t estimatedGenomeSize);
void findTranslocations(ofstream & outputFileDescriptor, samfile_t *fp, int32_t MinInsert,  int32_t MaxInsert, uint64_t estimatedGenomeSize);
pairStatus checkPair(const bam1_core_t *core_read1, const bam1_core_t *core_read2, int32_t MinInsert, int32_t MaxInsert);

// global variables
uint32_t contigsNumber;
float readCoverage;
float spanningCoverage;

int main(int argc, char *argv[]) {
	//MAIN VARIABLE
	string alignmentFile		 = "";

	bool outtie 				 = true;
	unsigned int WINDOW 		 = 1000;
	uint64_t estimatedGenomeSize = 0;
	int MinInsert;
	int MaxInsert;

	// VARIABLES USING DURING COMPUTATION
	string outputFile = "output.bed";
	uint64_t genomeLength = 0;


	ofstream outputFileDescriptor;
	outputFileDescriptor.open (outputFile.c_str());

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());


	desc.add_options() ("help", "produce help message")
							("sam", po::value<string>(), "alignment file in bam format, expected sorted by read name")
							("min-insert",  po::value<int>(), "paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goes from beginning of first read to end of second read")
							("max-insert",  po::value<int>(), "paired reads maximum allowed insert size. Used in order to filter outliers.")
							("orientation", po::value<string>(), "expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default outtie")
							("genome-size", po::value<unsigned long int>(), "estimated genome size (if not supplied genome size is believed to be assembly length")
							("output",  po::value<string>(), "Header output file names")
							("window",  po::value<unsigned int>(), "window size for features computation")
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

	// PARSE SAM/BAM file
	if (!vm.count("sam")) {
		DEFAULT_CHANNEL << "Please specify --sam " << endl;
		exit(0);
	}


	if (vm.count("sam")) {
		alignmentFile = vm["sam"].as<string>();
	}


	if (!vm.count("min-insert") || !vm.count("max-insert")) {
		DEFAULT_CHANNEL << "Please specify min-insert and max-insert " << endl;
		exit(0);
	}

	if (vm.count("min-insert")) {
		MinInsert = vm["min-insert"].as<int>();
		if(MinInsert <= 0) {
			DEFAULT_CHANNEL << "minimum insert should be at least 1\n";
			DEFAULT_CHANNEL << desc << endl;
			exit(2);
		}
	}
	if (vm.count("max-insert")) {
		MaxInsert = vm["max-insert"].as<int>();
	}

	if (vm.count("orientation")) {
		string userOrientation = vm["orientation"].as<string>();
		if(userOrientation.compare("innie") == 0) {
			outtie = false;
		} else if (userOrientation.compare("outtie") == 0) {
			outtie = true;
		} else {
			DEFAULT_CHANNEL << "outtie (<- ->) or innie (-> <-) only allowed orientations\n";
			DEFAULT_CHANNEL << desc << endl;
			exit(2);
		}
	}

	if (vm.count("window")) {
		WINDOW = vm["window"].as<unsigned int>();
	}

	if (vm.count("output")) {
		string header = vm["output"].as<string>();
		outputFile = header + ".txt";
	}

	if (vm.count("genome-size")) {
		estimatedGenomeSize = vm["genome-size"].as<unsigned long int>();
	} else {
		estimatedGenomeSize = 0;
	}

	if(vm.count("sam")){
		cout << "sam file name is " << alignmentFile << endl;
		cout << "library min " << MinInsert << "\n";
		cout << "library max " << MaxInsert << "\n";
		if(outtie) {
			cout << "library orientation <- ->\n";
		} else {
			cout << "library orientation -> <-\n";
		}
	}




	samfile_t *fp;
	fp = open_alignment_file(alignmentFile);
	EXIT_IF_NULL(fp);

	bam_header_t* head = fp->header; // sam header
	map<string,unsigned int> contig2position;
	map<unsigned int,string> position2contig;

	for(int i=0; i< head->n_targets ; i++) {
		genomeLength += head->target_len[i];
		contig2position[head->target_name[i]]=contigsNumber;   // keep track of contig name and position in order to avoid problems when processing two libraries
		position2contig[contigsNumber] = head->target_name[i]; //
		//cout << head->target_name[i] << " " << contigsNumber << " " <<  head->target_len[i]  << " " << genomeLength << "\n";
		contigsNumber++;
	}
	if (estimatedGenomeSize == 0) {
		estimatedGenomeSize = genomeLength; // if user has not provided genome size, approaximated it with assembly size
	}
	cout << "total number of contigs " << contigsNumber << endl;
	//cout << "assembly length " << genomeLength << "\n";
	cout << "estimated Genome length " << estimatedGenomeSize << "\n";

	samclose(fp); // close the file

	//Now compute the BAM file
	fp = open_alignment_file(alignmentFile);
	EXIT_IF_NULL(fp);
	head = fp->header; // sam header
	// obtain library statistics
	//computeLibraryStats(fp,  MinInsert, MaxInsert, estimatedGenomeSize);
	samclose(fp); // close the file


	// now find translocations
	fp = open_alignment_file(alignmentFile);
	EXIT_IF_NULL(fp);
	head = fp->header; // sam header
	findTranslocations(outputFileDescriptor, fp,MinInsert, MaxInsert, estimatedGenomeSize);
	samclose(fp); // close the file


}



void computeLibraryStats(samfile_t *fp, uint32_t MinInsert, uint32_t MaxInsert, uint64_t estimatedGenomeSize) {
	//Initialize bam entity
	bam1_t *b_read1 = bam_init1();
	bam1_t *b_read2 = bam_init1();
	uint32_t reads							 = 0;
	uint64_t properAlignedReadsLength		 = 0; // total length of properly aligned pairs
	uint64_t properAlignedInsertsLength 	 = 0; // total length of proper inserts

	//mated reads (not necessary correctly mated)
	uint32_t matedPairs 					 = 0; // length of reads that align on a contig with the mate
	//correctly aligned mates
	uint32_t correctlyMatedpairs 			 = 0; // total number of correctly mated reads
	//wrongly oriented reads
	uint32_t wronglyOrientedPairs	 		 = 0; // number of wrongly oriented reads
	//wrongly distance reads
	uint32_t wronglyDistancePairs       	 = 0; // number of reads at the wrong distance
	//singletons
	uint32_t singletonPairs 				 = 0; // number of singleton reads
	//mates on different contigs
	uint32_t matesDifferentChromosomesPairs  = 0; // number of contig placed in a different contig
	//optical duplicates
	uint32_t matesOpticalDuplicates			 = 0; // number of contig placed in a different contig


	while (samread(fp, b_read1) >= 0) {
		const bam1_core_t *core_read1 = &b_read1->core;
		if(core_read1->flag&BAM_FMUNMAP) { // mate unmapped
			singletonPairs ++;
			reads ++;
		} else {
			//extract the core from both alignments
			samread(fp, b_read2);
			const bam1_core_t *core_read2 = &b_read2->core;

			if (core_read1 == NULL or core_read2 == NULL) {  //There is something wrong with the read/file
				printf("Input file is corrupt!");
				return ;
			}
			reads += 2; // otherwise one more read is readed
			matedPairs ++;

			pairStatus currentPair = checkPair(core_read1,core_read2, MinInsert, MaxInsert);

			switch (currentPair) {
				case pair_Proper : 				correctlyMatedpairs ++;  break;
				case pair_wrongOrientation : 	wronglyOrientedPairs ++; break;
				case pair_wrongDistance : 		wronglyDistancePairs ++; break;
				case pair_wrongChrs :			matesDifferentChromosomesPairs++; break;
				case pair_wrongOptDup :			matesOpticalDuplicates ++; break;
				default : break;
			}


			if(currentPair ==  pair_Proper) {
				int32_t read1_insertSize  = core_read1->isize;
				int32_t read1_length      = core_read1->l_qseq;
				int32_t read2_length      = core_read2->l_qseq;
				properAlignedReadsLength  += (read1_length + read2_length);
				properAlignedInsertsLength+= read1_insertSize;
			}

		}

	}

	readCoverage 	 = (float)(properAlignedReadsLength)/estimatedGenomeSize;
	spanningCoverage = (float)(properAlignedInsertsLength)/estimatedGenomeSize;;

	cout << "number of reads " 						<< reads 							<< "\n";
	cout << "number of pairs " 						<< matedPairs	 					<< "\n";
	cout << "number of correctly oriented pairs "	<< correctlyMatedpairs 				<< "\n";
	cout << "number of wrongly oriented pairs " 	<< wronglyOrientedPairs 			<< "\n";
	cout << "number of wrongly distance pairs " 	<< wronglyDistancePairs 			<< "\n";
	cout << "number of wrongly chromosomes pairs " 	<< matesDifferentChromosomesPairs 	<< "\n";
	cout << "number of Optical duplicates pairs "	<< matesOpticalDuplicates 			<< "\n";
	cout << "------------\n";
	cout << "read coverage across the genome " 		<< readCoverage << "\n";
	cout << "spanning coverage across the genome " 	<< spanningCoverage << "\n";

}





void findTranslocations(ofstream & outputFileDescriptor, samfile_t *fp, int32_t MinInsert, int32_t MaxInsert, uint64_t estimatedGenomeSize) {
	//Initialize bam entity
	bam1_t *b_read1 = bam_init1();
	bam1_t *b_read2 = bam_init1();

	bam_header_t* head = fp->header; // sam header
	for(int i=0; i< head->n_targets ; i++) {
		cout << i << " " << head->target_name[i] << "\n";
	}

	Translocations *Trans;
	Trans = new Translocations(contigsNumber);

	Trans->initTrans(fp);

	while (samread(fp, b_read1) >= 0) {
		const bam1_core_t *core_read1 = &b_read1->core;
		if(core_read1->flag&BAM_FMUNMAP) { // mate unmapped
			//DONOTHING as I have removed unmapped reads, if one read is without the mate it means that I need to read only one!!!!
		} else {
			samread(fp, b_read2);

			const bam1_core_t *core_read2 = &b_read2->core;
			if (core_read1 == NULL or core_read2 == NULL) {  //There is something wrong with the read/file
				printf("Input file is corrupt!");
				return ;
			}
			pairStatus currentPair = checkPair(core_read1,core_read2, MinInsert, MaxInsert );


			if(currentPair ==  pair_wrongChrs) {
				// possible trans-location event found
				uint32_t startRead_1 			= core_read1->pos; // start position on the chromosome
				uint32_t chromosomeRead_1		= core_read1->tid;
				uint32_t qualityAligRead_1		= core_read1->qual;

				uint32_t startRead_2 			= core_read2->pos; // start position on the chromosome
				uint32_t chromosomeRead_2		= core_read2->tid;
				uint32_t qualityAligRead_2		= core_read2->qual;

				if(qualityAligRead_1 >= 20 and qualityAligRead_2 >= 20) {
					if(chromosomeRead_1 < chromosomeRead_2) {
						Trans->insertConnection(chromosomeRead_1,startRead_1, chromosomeRead_2,startRead_2);
					} else {
						Trans->insertConnection(chromosomeRead_2,startRead_2, chromosomeRead_1,startRead_1);
					}
				}
			} else if (currentPair == pair_Deletion) {
				// possible deletion event found
				int32_t startRead_1 			= core_read1->pos; // start position on the chromosome
				int32_t chromosomeRead_1		= core_read1->tid;
				int32_t qualityAligRead_1		= core_read1->qual;

				int32_t startRead_2 			= core_read2->pos; // start position on the chromosome
				int32_t chromosomeRead_2		= core_read2->tid;
				int32_t qualityAligRead_2		= core_read2->qual;

				if(qualityAligRead_1 >= 20 and qualityAligRead_2 >= 20) {
					//if(chromosomeRead_1 == 13  ) { //and startRead_1 >= 16700000 and startRead_1 < 16850000
					//	int32_t difference = startRead_2 -  startRead_1;
					//	cout << chromosomeRead_1 << "\t" << startRead_1 << " --- " <<  chromosomeRead_2 << "\t" << startRead_2 << " --- " << difference <<"\n";
					//}
					if(startRead_1 < startRead_2) {
						Trans->insertConnection(chromosomeRead_1,startRead_1, chromosomeRead_2,startRead_2);
					} else {
						Trans->insertConnection(chromosomeRead_2,startRead_2, chromosomeRead_1,startRead_1);
					}
				}
			}
		}
	}


	uint32_t minimumNumberOfSupportingPairs = 10;
	readCoverage = 1.44958;
	float minCov = readCoverage/4;
	float maxCov = readCoverage*3;
	uint32_t windowSize = 4000;
	uint32_t windowStep = 1000;

	for (uint32_t i = 0; i<= Trans->chromosomesNum; i++) {
		for(uint32_t j = i +1; j<= Trans->chromosomesNum; j++) {
			Trans->findEvents(outputFileDescriptor, i,j, minimumNumberOfSupportingPairs, minCov, maxCov, windowSize, windowStep);
		}
	}


	ofstream outputDeletionFile;
	outputDeletionFile.open ("output_deletions.bed");
	minimumNumberOfSupportingPairs = ;
	minCov = 0;
	maxCov = 100;
	windowSize = 8000;
	windowStep = 1000;
	for (uint32_t i = 0; i<= Trans->chromosomesNum; i++) {
		Trans->findEvents(outputDeletionFile, i, i, minimumNumberOfSupportingPairs, minCov, maxCov, windowSize, windowStep);
	}
}



pairStatus checkPair(const bam1_core_t *core_read1, const bam1_core_t *core_read2, int32_t MinInsert, int32_t MaxInsert) {
	//uint32_t read1_startingPos = core_read1->mpos;
	if( (core_read1->flag&BAM_FDUP) || (core_read1->flag&BAM_FQCFAIL) || (core_read2->flag&BAM_FDUP) || (core_read2->flag&BAM_FQCFAIL)) {// both reads are mapped but at least one is an optical duplicat or fails vendor quality specifications
		return pair_wrongOptDup;
	} else if(core_read1->tid != core_read2->tid) { // read mapped on different chromosomes (the one I am interested in)
		return pair_wrongChrs;
	}
	// if it is neither an optical duplicate nor a pair aligned on different chromosomes
	int32_t read1_insertSize  = core_read1->isize;
	if(read1_insertSize < 0) {
		 read1_insertSize = -1*(read1_insertSize);
	}

	int32_t read2_insertSize  = core_read2->isize;
	if(read2_insertSize < 0) {
		read2_insertSize = -1*(read2_insertSize);
	}

	bool read1_forw; // true for forward and false otherwise
	if(!(core_read1->flag&BAM_FREVERSE)) {
		read1_forw = true;
	} else {
		read1_forw = false;
	}
	bool read1_isFirst;
	if(core_read1->flag&BAM_FREAD1) { // if this is the first read in the pair
		read1_isFirst = true;
	} else {
		read1_isFirst = false;
	}

	bool read2_forw; // true for forward and false otherwise
	if(!(core_read2->flag&BAM_FREVERSE)) {
		read2_forw = true;
	} else {
		read2_forw = false;
	}
	bool read2_isFirst;
	if(core_read2->flag&BAM_FREAD1) { // if this is the first read in the pair
		read2_isFirst = true;
	} else {
		read2_isFirst = false;
	}



	//if(read1_forw == read2_forw) {
	//	return pair_wrongOrientation;
	//} else
	if (read1_insertSize >= MinInsert and read1_insertSize <= MaxInsert) {
		return pair_Proper;
	} else if (read1_insertSize > MaxInsert) {
		return pair_Deletion;
	}
	return pair_Proper;

	/*

	if(read1_isFirst) {
		if(!read1_forw && read2_forw) {
			if(read1_insertSize >= MinInsert and read1_insertSize <= MaxInsert) {
				return pair_Proper;
			} else {
				if(read1_insertSize > MaxInsert) {
					return pair_Deletion;
				}
				return pair_wrongDistance;
			}
		} else {
			return pair_wrongOrientation;
		}
	} else if(read2_isFirst) {
		if(!read1_forw && read2_forw) {
			if(read2_insertSize >= MinInsert and read2_insertSize <= MaxInsert) {
				return pair_Proper;
			} else {
				if(read2_insertSize > MaxInsert) {
					return pair_Deletion;
				}
				return pair_wrongDistance;
			}
		} else {
			return pair_wrongOrientation;
		}
	}
	*/



}



