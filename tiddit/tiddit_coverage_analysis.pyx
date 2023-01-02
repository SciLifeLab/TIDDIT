import numpy
import re
cimport numpy
import pysam
import gzip

import tiddit.tiddit_coverage as tiddit_coverage


def get_gc(str reference_fasta,bam_header,bin_size):

	gc=tiddit_coverage.create_coverage(bam_header,bin_size)[0]

	reference={}
	chromosomes=[]
	if not reference_fasta.endswith(".gz"):
		with open(reference_fasta, 'r') as f:
			sequence=f.read()
	else:
		with gzip.open(reference_fasta, 'r') as f:
			sequence=f.read()


	#split_reference=sequence.split(">")
	split_reference=re.split("\n>|^>", sequence)
	del sequence
	del split_reference[0]

	for chromosome in split_reference:
		content=chromosome.split("\n",1)
		reference[content[0].strip().split()[0]]=content[1].replace("\n","")

	N=set(["N","n"])
	GC=set(["G","g","C","c"])

	for chromosome in gc:
		for i in range(0,len(gc)):
			start=i*bin_size
			end=(i+1)*bin_size
			if end > len(reference[chromosome]):
				end=len(reference[chromosome])

			N_count=0
			GC_count=0
			for j in range(start,end):
				if reference[chromosome][j] in N:
					N_count +=1

			#masked bin
			if N_count > bin_size/2:
				gc[chromosome][i] = -1	

	return(gc)


def determine_ploidy(dict coverage_data,contigs,dict library,int ploidy,str prefix,c, str reference_fasta,int bin_size,bam_header):

	gc=get_gc(reference_fasta,bam_header,bin_size)

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

