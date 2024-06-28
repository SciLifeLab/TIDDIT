import pysam
import numpy
import time

def statistics(bam_file_name,ref,min_mapq,max_ins_len,n_reads):
	library={}
	samfile = pysam.AlignmentFile(bam_file_name, "r",reference_filename=ref)

	insert_size=[]
	read_length=[]
	is_innie=0
	is_outtie=0

	n_sampled=0
	t=time.time()

	for read in samfile.fetch():

		read_length.append( read.query_length )	
		n_sampled+=1

		if n_sampled > n_reads:
			break

		if read.mate_is_unmapped:
			continue

		if read.is_reverse == read.mate_is_reverse:
			continue

		if read.next_reference_name != read.reference_name or (read.template_length > max_ins_len):
			continue

		if read.next_reference_start < read.reference_start:
			continue	

		if read.is_supplementary or read.is_secondary or read.is_duplicate or read.mapq < min_mapq:
			continue

		insert_size.append( read.template_length )

		if read.is_reverse and not read.mate_is_reverse:
			is_outtie+=1
		else:
			is_innie+=1


	samfile.close()

	library["avg_read_length"]=numpy.average(read_length)
	if len(insert_size):
		library["avg_insert_size"]=numpy.average(insert_size)
		library["std_insert_size"]=numpy.std(insert_size)
		library["percentile_insert_size"]=numpy.percentile(insert_size, 99.9)
	else:
		library["avg_insert_size"]=0
		library["std_insert_size"]=0
		library["percentile_insert_size"]=0

		

	print("LIBRARY STATISTICS")
	if is_innie > is_outtie:
		library["mp"]=False
		print("\tPair orientation = Forward-Reverse")
	else:
		print("\tPair orientation = Reverse-Forward")
		library["mp"]=True


	print("\tAverage Read length = {}".format(library["avg_read_length"]) )
	print("\tAverage insert size = {}".format(library["avg_insert_size"]) )
	print("\tStdev insert size = {}".format(library["std_insert_size"] ) )
	print("\t99.95 percentile insert size = {}".format( library["percentile_insert_size"]) )
	print("Calculated statistics in: " + str( t-time.time() ))
	print("")

	return(library)
