import pysam
import sys
import os
import itertools
import time

from . import tiddit_coverage

def find_SA_query_range(SA):
	a =pysam.AlignedSegment()
	a.reference_start=int( SA[1] )

	if SA[2] == "+":
		a.flag = 64
	else:
		a.flag = 80

	cdef list SA_cigar=[]
	SC = ["".join(x) for _, x in itertools.groupby(SA[3], key=str.isdigit)]		

	cdef dict s_to_op={"M":0,"S":4,"H":5,"D":2,"I":1}
	for i in range(0,int(len(SC)/2)):
		op=s_to_op[SC[i*2+1]]
		SA_cigar.append( (op,int(SC[i*2])) )

	a.cigar = tuple(SA_cigar)
	return(a)

def SA_analysis(read,min_q,splits,tag):
	#print(read.query_alignment_start,read.query_alignment_end,read.is_reverse,read.cigarstring)
	suplementary_alignments=read.get_tag(tag).rstrip(";").split(";")
	
	if len(suplementary_alignments) > 1:
		SA_lengths=[]
		ok_q=[]
		for i in range(0,len(suplementary_alignments)):
			SA_data=suplementary_alignments[0].split(",")
			if int(SA_data[4]) >= min_q:
				ok_q.append(i)
				supplementry_alignment=find_SA_query_range(SA_data)
				SA_lengths.append(supplementry_alignment.query_alignment_end-supplementry_alignment.query_alignment_start)


		longest_aln=0
		for i in range(0,len(ok_q)):
			if SA_lengths[i] > SA_lengths[longest_aln]:
				longest_aln=i

		#all alignments fail quality treshold
		if len(ok_q) == 0:
			return(splits)

		#only one SA pass mapping quality treshold		
		elif len(ok_q) == 1:
			suplementary_alignments[0]=suplementary_alignments[ ok_q[0] ]
		#many SA pass mapping quality treshold, pick the longest alignment.
		else:
			suplementary_alignments[0]=suplementary_alignments[ longest_aln ]

	SA_data=suplementary_alignments[0].split(",")
	SA_pos=int(SA_data[1])

	if int(SA_data[4]) < min_q:
		return(splits)

	if (read.query_alignment_start ) < (read.query_length - read.query_alignment_end):
		split_pos=read.reference_end+1
	else:
		split_pos=read.reference_start+1

	supplementry_alignment=find_SA_query_range(SA_data)
	SA_chr=SA_data[0]


	if (supplementry_alignment.query_alignment_start ) < (supplementry_alignment.query_length - read.query_alignment_end):
		SA_split_pos=supplementry_alignment.reference_end
	else:
		SA_split_pos=supplementry_alignment.reference_start


	if SA_chr < read.reference_name:
		chrA=SA_chr
		chrB=read.reference_name
		tmp=split_pos
		split_pos=SA_split_pos
		SA_split_pos=tmp

	else:
		chrA=read.reference_name
		chrB=SA_chr

		if chrA == chrB:
			if SA_split_pos < split_pos:
				tmp=split_pos
				split_pos=SA_split_pos
				SA_split_pos=tmp

	if not chrA in splits:
		splits[chrA]={}

	if not chrB in splits[chrA]:
		splits[chrA][chrB]={}

	if not read.query_name in splits[chrA][chrB]:
		splits[chrA][chrB][read.query_name]=[]

	if "-" == SA_data[2]:
		splits[chrA][chrB][read.query_name]+=[split_pos,read.is_reverse,SA_split_pos,True]
	else:
		splits[chrA][chrB][read.query_name]+=[split_pos,read.is_reverse,SA_split_pos,False]

	return(splits)

def main(str bam_file_name,str ref,str prefix,int min_q,int max_ins,str sample_id):

	samfile = pysam.AlignmentFile(bam_file_name, "r",reference_filename=ref)
	bam_header=samfile.header
	cdef int bin_size=50
	cdef str file_type="wig"
	cdef str outfile=prefix+".tiddit_coverage.wig"

	t_update=0
	t_split=0
	t_disc=0
	t_tot=0

	coverage_data,end_bin_size=tiddit_coverage.create_coverage(bam_header,bin_size)	

	cdef dict data={}
	cdef dict splits={}
	cdef dict clips={}

	for chrA in bam_header["SQ"]:

		clips[chrA["SN"]]=[]
		for chrB in bam_header["SQ"]:
			if chrA["SN"] <= chrB["SN"]:
				if not chrA["SN"] in data:
					data[chrA["SN"]] = {}
					splits[chrA["SN"]] = {}
				data[chrA["SN"]][chrB["SN"]]={}
				splits[chrA["SN"]][chrB["SN"]]={}


	chromosome_set=set([])
	chromosomes=[]

	clip_dist=100

	t_tot=time.time()
	f=open("{}_tiddit/clips_{}.fa".format(prefix,sample_id ),"w")
	for read in samfile.fetch(until_eof=True):

		if read.is_unmapped or read.is_duplicate:
			continue

		t=time.time()
		if read.mapq >= min_q:
			coverage_data[read.reference_name]=tiddit_coverage.update_coverage(read,bin_size,coverage_data[read.reference_name],min_q,end_bin_size[read.reference_name])
		t_update+=time.time()-t

		if not read.reference_name in chromosome_set:
			print("Collecting signals on contig: {}".format(read.reference_name))
			chromosome_set.add(read.reference_name)

			if len(chromosomes):
				for clip in clips[ chromosomes[-1] ]:
					f.write("".join( clip ))
				del clips[ chromosomes[-1] ]

			chromosomes.append(read.reference_name)


		t=time.time()

		if read.has_tag("SA") and not (read.is_supplementary or read.is_secondary) and read.mapq >= min_q:
			splits=SA_analysis(read,min_q,splits,"SA")

		if not (read.is_supplementary or read.is_secondary) and read.mapq > 1:
			if (read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 10) and (read.cigartuples[-1][0] == 0 and read.cigartuples[-1][1] > 30) and len(read.cigartuples) < 7:
				clips[read.reference_name].append([">{}|{}|{}\n".format(read.query_name,read.reference_name,read.reference_start+1),read.query_sequence+"\n"])

			elif read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 10 and (read.cigartuples[0][0] == 0 and read.cigartuples[0][1] > 30) and len(read.cigartuples) < 7:
				clips[read.reference_name].append([">{}|{}|{}\n".format(read.query_name,read.reference_name,read.reference_start+1),read.query_sequence+"\n"])

		t_split+=time.time()-t

		t=time.time()
		if ( abs(read.isize) > max_ins or read.next_reference_name != read.reference_name ) and read.mapq >= min_q and not (read.is_supplementary or read.is_secondary):
			if read.next_reference_name < read.reference_name:
				chrA=read.next_reference_name
				chrB=read.reference_name
			else:
				chrA=read.reference_name
				chrB=read.next_reference_name

			if not read.query_name in data[chrA][chrB]:
				data[chrA][chrB][read.query_name]=[]


			data[chrA][chrB][read.query_name].append([read.reference_start+1,read.reference_end+1,read.is_reverse,read.reference_name])
		t_disc+=time.time()-t

	f.close()

	print("total",time.time()-t_tot)
	print("coverage",t_update)
	print("split",t_split)
	print("disc",t_disc)

	#print("writing coverage wig")
	#tiddit_coverage.print_coverage(coverage_data,bam_header,bin_size,file_type,outfile)

	print("Writing signals to file")

	f=open("{}_tiddit/discordants_{}.tab".format(prefix,sample_id),"w")

	for chrA in data:
		for chrB in data[chrA]:

			for fragment in data[chrA][chrB]:
				if len(data[chrA][chrB][fragment]) < 2:
					continue

				if chrA == chrB:
					if data[chrA][chrB][fragment][1][-1] < data[chrA][chrB][fragment][0][-1]:
						out=data[chrA][chrB][fragment][1][0:-1]+data[chrA][chrB][fragment][0][0:-1]
					else:
						out=data[chrA][chrB][fragment][0][0:-1]+data[chrA][chrB][fragment][1][0:-1]
				else:
					if data[chrA][chrB][fragment][0][-1] == chrA:
						out=data[chrA][chrB][fragment][0][0:-1]+data[chrA][chrB][fragment][1][0:-1]
					else:
						out=data[chrA][chrB][fragment][1][0:-1]+data[chrA][chrB][fragment][0][0:-1]

				f.write("{}\t{}\t{}\t{}\n".format(fragment,chrA,chrB,"\t".join(map(str, out )))  )

	f.close()

	f=open("{}_tiddit/splits_{}.tab".format(prefix,sample_id),"w")

	for chrA in splits:
		for chrB in splits[chrA]:

			for fragment in splits[chrA][chrB]:
				f.write("{}\t{}\t{}\t{}\n".format(fragment,chrA,chrB,"\t".join(map(str, splits[chrA][chrB][fragment] )))  )

	f.close()

	f=open("{}_tiddit/clips_{}.fa".format(prefix,sample_id),"a")

	for chrA in clips:
		for clip in clips[chrA]:
			f.write("".join( clip ))
	f.close()


	return(coverage_data)
	#return(coverage_data,clips)
	#return(coverage_data,clips)
