import sys
cimport numpy
import numpy
import math
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

def create_coverage(bam_header,bin_size):
	coverage_data={}
	end_bin_size={}

	for contig in bam_header["SQ"]:
		bins= int(math.ceil(contig["LN"]/float(bin_size)))
		coverage_data[ contig["SN"] ]=numpy.zeros(bins)
		end_bin_size[contig["SN"]]=contig["LN"]-math.floor( contig["LN"]/float(bin_size) )*bin_size

	return(coverage_data,end_bin_size)

def print_coverage(coverage_data,bam_header,bin_size,file_type,outfile):
	f=open(outfile,"w",buffering=819200)

	if file_type == "bed":
		f.write("#chromosome\tstart\tend\tcoverage\n")
	elif file_type == "wig":
		f.write("track type=wiggle_0 name=\"Coverage\" description=\"Per bin average coverage\"\n")

	for contig in bam_header["SQ"]:
		if file_type == "wig":
			f.write("fixedStep chrom={} start=1 step={}\n".format(contig["SN"],bin_size))

		for i in range(0,len(coverage_data[ contig["SN"] ]) ) :
			if file_type == "bed":
				bin_end=(i+1)*bin_size+1

				if i == len(coverage_data[ contig["SN"] ]) -1:
					bin_end=contig["LN"]

				f.write("{}\t{}\t{}\t{}\n".format(contig["SN"],1+i*bin_size,bin_end,coverage_data[ contig["SN"] ][i] ) )
			elif file_type == "wig":
				f.write("{}\n".format( coverage_data[ contig["SN"] ][i] ))

	f.close()

ctypedef numpy.double_t DTYPE_t
def update_coverage(read,int bin_size,numpy.ndarray[DTYPE_t, ndim=1] coverage_data,int min_q,int end_bin_size):

	cdef long ref_start=read.reference_start
	cdef long ref_end=read.reference_end

	cdef int first_bin=ref_start//bin_size
	cdef int end_bin=int(ref_end-1)//bin_size

	cdef int bases_first_bin

	if end_bin == first_bin:
		bases_first_bin=ref_end-ref_start
		coverage_data[first_bin]=float(bases_first_bin)/bin_size+coverage_data[first_bin]

		return(coverage_data)

	bases_first_bin=((first_bin+1)*bin_size)-ref_start
	coverage_data[first_bin]=float(bases_first_bin)/bin_size+coverage_data[first_bin]	
	
	cdef int bases_last_bin=(ref_end-1)-end_bin*bin_size

	if end_bin < len(coverage_data)-1:
		coverage_data[end_bin]+=float(bases_last_bin)/bin_size
	else:
		coverage_data[end_bin]+=float(bases_last_bin)/end_bin_size

	for i in range(first_bin+1,end_bin):
		coverage_data[i]+=1.0

	return(coverage_data)


#bam_file_name=sys.argv[1]

#samfile = pysam.AlignmentFile(bam_file_name, "r")
#bam_header=samfile.header
#bin_size=50
#file_type="bed"
#outfile=sys.argv[2]
#min_q=10

#coverage_data,end_bin_size=create_coverage(bam_header,bin_size)

#for read in samfile.fetch():
#	coverage_data=update_coverage(read,bin_size,coverage_data,min_q,end_bin_size)

#print_coverage(coverage_data,bam_header,bin_size,file_type,outfile)

#samfile.close()
