/*
   Francesco Vezzi
 */

#include "common.h"
#include "radix.h"
#include "samtools/sam.h"


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



void computeLibraryStats(samfile_t *fp, unsigned int minInsert, unsigned int maxInsert, uint64_t genomeLength);



int main(int argc, char *argv[]) {
	//MAIN VARIABLE


	string alignmentFile = "";
	int32_t MinInsert = 100;
	int32_t MaxInsert = 1000000;

	unsigned int WINDOW = 1000;
	uint64_t estimatedGenomeSize;



	string outputFile = "FRC.txt";
	string featureFile = "Features.txt";

	// PROCESS PARAMETERS
	stringstream ss;
	ss << package_description() << endl << endl << "Allowed options";
	po::options_description desc(ss.str().c_str());


	desc.add_options() ("help", "produce help message")
					("sam", po::value<string>(), "paired end alignment file (in sam or bam format). Expected orientation <- ->")
					("min-insert",  po::value<int>(), "paired reads minimum allowed insert size. Used in order to filter outliers. Insert size goeas from beginning of first read to end of second read")
					("max-insert",  po::value<int>(), "paired reads maximum allowed insert size. Used in order to filter outliers.")
					("window",  po::value<unsigned int>(), "window size for features computation")
					("genome-size", po::value<unsigned long int>(), "estimated genome size (if not supplied genome size is believed to be assembly length")
					("output",  po::value<string>(), "Header output file names (default FRC.txt and Features.txt)")
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

	if (vm.count("min-insert")) {
		MinInsert = vm["min-insert"].as<int>();
		if(MinInsert <= 0) {
			DEFAULT_CHANNEL << "minimum insert should be at least 1\n";
			DEFAULT_CHANNEL << desc << endl;
			exit(2);
		}
	}
	if (vm.count("pe-max-insert")) {
		MaxInsert = vm["pe-max-insert"].as<int>();
	}


	if(vm.count("sam")){
		cout << "sam file name is " << alignmentFile << endl;
		cout << "library min " << MinInsert << "\n";
		cout << "library max " << MaxInsert << "\n";
	}





	if (vm.count("window")) {
		WINDOW = vm["window"].as<unsigned int>();
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



	uint64_t genomeLength = 0;
	uint32_t contigsNumber = 0;
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
	cout << "assembly length " << genomeLength << "\n";
	cout << "estimated length " << estimatedGenomeSize << "\n";

	samclose(fp); // close the file

	//parse BAM


	fp = open_alignment_file(alignmentFile);



	EXIT_IF_NULL(fp);
	head = fp->header; // sam header


	computeLibraryStats(fp, 100, 5000, genomeLength ) ;

	return 0;


	int currentTid = -1;
	int reads = 0;
	bam1_t *b = bam_init1();

	// NOW PROCESS LIBRARIES

	fp = open_alignment_file(alignmentFile);
	EXIT_IF_NULL(fp);
	head = fp->header; // sam header
	while (samread(fp, b) >= 0) {
		//Get bam core.
		const bam1_core_t *core = &b->core;
		if (core == NULL) {
			printf("Input file is corrupt!");
			return -1;
		}
		++reads;

		// new contig
		if (is_mapped(core)) {
			if (core->tid != currentTid) { // another contig or simply the first one
				cout << "now porcessing contig " << core->tid << "\n";

				if(currentTid == -1) { // first read that I`m processing

				} else {

				}
			} else {
				//add information to current contig

			}
		}
	}

	samclose(fp); // close the file


}





void computeLibraryStats(samfile_t *fp, unsigned int minInsert, unsigned int maxInsert, uint64_t genomeLength) {


	//Initialize bam entity
	bam1_t *b = bam_init1();

	//All var declarations
	unsigned int contigs = 0; // number of contigs/scaffolds
	// total reads
	uint32_t reads = 0;
	uint64_t readsLength = 0;   // total length of reads
	uint32_t unmappedReads = 0;
	uint32_t mappedReads = 0;
	uint32_t zeroQualityReads = 0;
	uint32_t duplicates = 0;

	uint64_t contigSize = 0;
	uint64_t insertsLength = 0; // total inserts length
	float insertMean;
	float insertStd;

	// mated reads (not necessary correctly mated)
	uint32_t matedReads = 0;        // length of reads that align on a contig with the mate
	uint64_t matedReadsLength = 0;  // total length of mated reads
	// correctly aligned mates
	uint32_t correctlyMatedReads = 0; // total number of correctly mated reads
	uint64_t correctlyMatedReadsLength = 0; // length of correctly mated reads
	// wrongly oriented reads
	uint32_t wronglyOrientedReads = 0; // number of wrongly oriented reads
	uint64_t wronglyOrientedReadsLength = 0; // length of wrongly oriented reads
	// wrongly distance reads
	uint32_t wronglyDistanceReads       = 0; // number of reads at the wrong distance
	uint64_t wronglyDistanceReadsLength = 0;  // total length of reads placed in different contigs
	// singletons
	uint32_t singletonReads = 0; // number of singleton reads
	uint64_t singletonReadsLength = 0;     // total length of singleton reads
	// mates on different contigs
	uint32_t matedDifferentContig = 0; // number of contig placed in a different contig
	uint64_t matedDifferentContigLength = 0; // total number of reads placed in different contigs

	float C_A = 0; // total read coverage
	float S_A = 0; // total span coverage
	float C_M = 0; // coverage induced by correctly aligned pairs
	float C_W = 0; // coverage induced by wrongly mated pairs
	float C_S = 0; // coverage induced by singletons
	float C_C = 0; // coverage induced by reads with mate on a diferent contif

	// compute mean and std on the fly
	float Mk = 0;
	float Qk = 0;
	uint32_t counterK = 1;
	//Keep header for further reference
	bam_header_t* head = fp->header;
	int32_t currentTid = -1;
	int32_t iSize;

	while (samread(fp, b) >= 0) {
		//Get bam core.
		const bam1_core_t *core = &b->core;
		if (core == NULL) {  //There is something wrong with the read/file
			printf("Input file is corrupt!");
			return ;
		}
		++reads; // otherwise one more read is readed

		if (!is_mapped(core)) {
			++unmappedReads;
		} else {
			if (core->tid != currentTid) {
				//Get length of next section
				contigSize = head->target_len[core->tid];
				contigs++;
				if (contigSize < 1) {//We can't have such sizes! this can't be right
					fprintf(stderr,"%s has size %d, which can't be right!\nCheck bam header!",head->target_name[core->tid],contigSize);
				}
				currentTid = core->tid;
			}
			//&& !(core->flag&BAM_FSECONDARY)
			if(!(core->flag&BAM_FUNMAP) && !(core->flag&BAM_FDUP) && !(core->flag&BAM_FQCFAIL)) { // if read has been mapped and it is not a DUPLICATE or a SECONDARY alignment
				uint32_t* cigar = bam1_cigar(b);
				++mappedReads;
				uint32_t alignmentLength = bam_cigar2qlen(core,cigar);
				readsLength += alignmentLength;
				uint32_t startRead = core->pos; // start position on the contig
				uint32_t startPaired;
				//Now check if reads belong to a proper pair: both reads aligned on the same contig at the expected distance and orientation
				if ((core->flag&BAM_FREAD1) //First in pair
						&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
						&& (core->tid == core->mtid) /*Mate on the same chromosome*/
				) {
					//pair is mapped on the same contig and I'm looking the first pair
					startPaired = core->mpos;
					if(startRead < startPaired) {
						iSize = (startPaired + core->l_qseq -1) - startRead; // insert size, I consider both reads of the same length
						if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
							//here reads are correctly oriented
							if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
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
								correctlyMatedReads++;
								correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
							} else {
								wronglyDistanceReads++;
								wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
							}
						} else {
							//pair is wrongly oriented
							wronglyOrientedReads++;
							wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
						}
					} else {
						iSize = (startRead + alignmentLength - 1) - startPaired;
						if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
							//here reads are correctly oriented
							//here reads are correctly oriented
							if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
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
								correctlyMatedReads++;
								correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
							} else {
								wronglyDistanceReads++;
								wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
							}
						} else {
							//pair is wrongly oriented
							wronglyOrientedReads++;
							wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
						}
					}
				} else  if ((core->flag&BAM_FREAD2) //Second in pair
						&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
						&& (core->tid == core->mtid) /*Mate on the same chromosome*/
				)
					// if I'm considering the second read in a pair I must check it is is a correctly mated read and if this is the case update the right variables
				{
					startPaired = core->mpos;
					if(startRead > startPaired) {
						iSize = (startRead + alignmentLength -1) - startPaired;
						if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
							//here reads are correctly oriented
							if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert, no need to update insert coverage
								correctlyMatedReads++;
								correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
							} else {
								wronglyDistanceReads++;
								wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
							}
						} else {
							//pair is wrongly oriented
							wronglyOrientedReads++;
							wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
						}
					} else {
						iSize = (startPaired + core->l_qseq -1) - startRead;
						if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
							//here reads are correctly oriented
							if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert, no need to update insert coverage
								correctlyMatedReads++;
								correctlyMatedReadsLength +=  bam_cigar2qlen(core,cigar); // update number of correctly mapped and their length
							}else {
								wronglyDistanceReads++;
								wronglyDistanceReadsLength += bam_cigar2qlen(core,cigar);
							}
						} else {
							//pair is wrongly oriented
							wronglyOrientedReads++;
							wronglyOrientedReadsLength += bam_cigar2qlen(core,cigar);
						}
					}
				} else if (core->tid != core->mtid && !(core->flag&BAM_FMUNMAP)) {
					//Count inter-chrom pairs
					matedDifferentContig++;
					matedDifferentContigLength += bam_cigar2qlen(core,cigar);
				} else if(core->flag&BAM_FMUNMAP) {
					// if mate read is unmapped
					singletonReads++;
					singletonReadsLength += bam_cigar2qlen(core,cigar);
				}


				if (core->flag&BAM_FPROPER_PAIR) {
					//Is part of a proper pair
					matedReads ++; // increment number of mated reads
					matedReadsLength += bam_cigar2qlen(core,cigar); // add the length of the read aligne as proper mate (not necessary correctly mated)
				}

				if (core->flag&BAM_FDUP) {   //This is a duplicate. Don't count it!.
					++duplicates;
				}
			} else {
				++zeroQualityReads;

			}
		}
	}

	cout << "LIBRARY STATISTICS\n";
	cout << "\ttotal reads number " << reads << "\n";
	cout << "\ttotal mapped reads " << mappedReads << "\n";
	cout << "\ttotal unmapped reads " << unmappedReads << "\n";
	cout << "\tproper pairs " << matedReads << "\n";
	cout << "\tzero quality reads " << zeroQualityReads << "\n";
	cout << "\tcorrectly oriented " << correctlyMatedReads << "\n";
	cout << "\twrongly oriented " << wronglyOrientedReads << "\n";
	cout << "\twrongly distance " << wronglyDistanceReads << "\n";
	cout << "\twrongly contig " <<  matedDifferentContig << "\n";
	cout << "\tsingletons " << singletonReads << "\n";

	uint32_t total = correctlyMatedReads + wronglyOrientedReads + wronglyDistanceReads + matedDifferentContig + singletonReads;
	cout << "\ttotal " << total << "\n";
	cout << "\tCoverage statistics\n";


	cout << "\tC_A = " << C_A << endl;
	cout << "\tS_A = " << S_A << endl;
	cout << "\tC_M = " << C_M << endl;
	cout << "\tC_W = " << C_W << endl;
	cout << "\tC_S = " << C_S << endl;
	cout << "\tC_C = " << C_C << endl;
	cout << "\tMean Insert length = " << Mk << endl;
	cout << "\tStd Insert length = " << Qk << endl;
	cout << "----------\n";



}


