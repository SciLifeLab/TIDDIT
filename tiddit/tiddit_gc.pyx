import pysam
import numpy
import math
from joblib import Parallel, delayed

def binned_gc(fasta_path,contig,bin_size,n_cutoff):
	fasta=pysam.FastaFile(fasta_path)
	contig_length=fasta.get_reference_length(contig)
	number_of_bins=int(math.ceil(contig_length/bin_size))
	
	contig_gc=numpy.zeros(number_of_bins,dtype=numpy.int8)

	next_start=0
	for bin in range(0,number_of_bins):
		slice=fasta.fetch(contig, next_start, next_start+bin_size)
		n=0
		gc=0
		number_of_chars=0

		for character in slice:
			number_of_chars += 1
			if character == "N" or character == "n":
				n+=1
			elif character == "C" or character == "c" or character == "G" or character == "g":
				gc+=1

		if n/bin_size > n_cutoff:
			contig_gc[bin]=-1
		else:
			result = round(100*gc/number_of_chars)
			contig_gc[bin] = result
			
		next_start+=bin_size
	return([contig,contig_gc])

def main(reference,contigs,threads,bin_size,n_cutoff):
	gc_list=Parallel(n_jobs=threads)( delayed(binned_gc)(reference,contig,bin_size,n_cutoff) for contig in contigs)

	gc_dictionary={}
	for gc in gc_list:
		gc_dictionary[gc[0]]=gc[1]

	return(gc_dictionary)


#contigs=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"]
