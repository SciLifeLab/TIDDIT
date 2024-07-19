import pysam
import numpy
import math
from joblib import Parallel, delayed

def binned_gc(fasta_path,contig,bin_size,n_cutoff):
	fasta=pysam.FastaFile(fasta_path)
	contig_length=fasta.get_reference_length(contig)
	elements=int(math.ceil(contig_length/bin_size))
	
	contig_gc=numpy.zeros(elements,dtype=numpy.int8)

	start=0
	for i in range(0,elements):
		slice=fasta.fetch(contig, start, start+bin_size)
		n=0
		gc=0

		for charachter in slice:
			if charachter == "N" or charachter == "n":
				n+=1
			elif charachter == "C" or charachter == "c" or charachter == "G" or charachter == "g":
				gc+=1

		if n/bin_size > n_cutoff:
			contig_gc[i]=-1

		else:
			contig_gc[i]=round(100*gc/elements)

		start+=bin_size

	return([contig,contig_gc])

def main(reference,contigs,threads,bin_size,n_cutoff):
	gc_list=Parallel(n_jobs=threads)( delayed(binned_gc)(reference,contig,bin_size,n_cutoff) for contig in contigs)

	gc_dictionary={}
	for gc in gc_list:
		gc_dictionary[gc[0]]=gc[1]

	return(gc_dictionary)


#contigs=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"]
