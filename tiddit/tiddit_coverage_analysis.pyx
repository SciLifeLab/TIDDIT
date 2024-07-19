import numpy
import re
cimport numpy
import pysam
import gzip

import tiddit.tiddit_coverage as tiddit_coverage

def determine_ploidy(dict coverage_data,contigs,dict library,int ploidy,str prefix,c, str reference_fasta,int bin_size,bam_header,gc):

	f=open( "{}.ploidies.tab".format(prefix),"w" )
	f.write("Chromosome\tPloidy\tPloidy_rounded\tMean_coverage\n")
	all_cov=[]
	for chromosome in coverage_data:
		tmp=[]
		for i in range(0,len(coverage_data[chromosome])):
			if coverage_data[chromosome][i] > 0 and gc[chromosome][i] != -1:

				tmp.append(coverage_data[chromosome][i])
				all_cov.append(coverage_data[chromosome][i])
		
		library[ "avg_coverage_{}".format(chromosome) ]=numpy.median(tmp)
		if numpy.isnan(library[ "avg_coverage_{}".format(chromosome) ]):
			library[ "avg_coverage_{}".format(chromosome) ]=0

	if not c:
		library["avg_coverage"]=numpy.median(all_cov)
	else:
		library["avg_coverage"]=c

	for chromosome in contigs:
		if not chromosome in coverage_data:
			continue

		avg_coverage_contig=library[ "avg_coverage_{}".format(chromosome) ]
		library["contig_ploidy_{}".format(chromosome)]=int(round(ploidy*avg_coverage_contig/library["avg_coverage"]))
		f.write("{}\t{}\t{}\t{}\n".format(chromosome,avg_coverage_contig/library["avg_coverage"]*ploidy,library["contig_ploidy_{}".format(chromosome)],avg_coverage_contig))
	

	f.close()
	return(library)

