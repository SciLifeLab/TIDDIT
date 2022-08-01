import pysam
import sys
import os
import itertools
import time
from joblib import Parallel, delayed

import tiddit.tiddit_coverage as tiddit_coverage

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

def SA_analysis(read,min_q,tag,reference_name):
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
			return()

		#only one SA pass mapping quality treshold		
		elif len(ok_q) == 1:
			suplementary_alignments[0]=suplementary_alignments[ ok_q[0] ]
		#many SA pass mapping quality treshold, pick the longest alignment.
		else:
			suplementary_alignments[0]=suplementary_alignments[ longest_aln ]

	SA_data=suplementary_alignments[0].split(",")
	SA_pos=int(SA_data[1])

	if int(SA_data[4]) < min_q:
		return()

	cdef long read_query_alignment_end=read.query_alignment_end

	if (read.query_alignment_start ) < (read.query_length - read_query_alignment_end):
		split_pos=read.reference_end+1
	else:
		split_pos=read.reference_start+1

	supplementry_alignment=find_SA_query_range(SA_data)
	SA_chr=SA_data[0]


	if (supplementry_alignment.query_alignment_start ) < (supplementry_alignment.query_length - read_query_alignment_end):
		SA_split_pos=supplementry_alignment.reference_end
	else:
		SA_split_pos=supplementry_alignment.reference_start


	if SA_chr < reference_name:
		chrA=SA_chr
		chrB=reference_name
		tmp=split_pos
		split_pos=SA_split_pos
		SA_split_pos=tmp

	else:
		chrA=reference_name
		chrB=SA_chr

		if chrA == chrB:
			if SA_split_pos < split_pos:
				tmp=split_pos
				split_pos=SA_split_pos
				SA_split_pos=tmp

	split=[]
	if "-" == SA_data[2]:
		split=[chrA,chrB,read.query_name,split_pos,read.is_reverse,SA_split_pos, True]
	else:
		split=[chrA,chrB,read.query_name,split_pos,read.is_reverse,SA_split_pos,False]
		#splits[chrA][chrB][read.query_name]+=[split_pos,read.is_reverse,SA_split_pos,False]

	return(split)

def worker(str chromosome, str bam_file_name,str ref,str prefix,int min_q,int max_ins,str sample_id, int bin_size):
	print("Collecting signals on contig: {}".format(chromosome))
	samfile = pysam.AlignmentFile(bam_file_name, "r",reference_filename=ref,index_filename="{}_tiddit/{}.csi".format(prefix,sample_id))
	bam_header=samfile.header
	coverage_data,end_bin_size=tiddit_coverage.create_coverage(bam_header,bin_size,chromosome)	

	clips=[]
	data=[]
	splits=[]

	clip_dist=100

	cdef long read_position
	cdef long read_end
	cdef int mapq

	for read in samfile.fetch(chromosome,until_eof=True):

		if read.is_unmapped or read.is_duplicate:
			continue

		read_chromosome=read.reference_name
		mate_chromosome=read.next_reference_name
		read_position=read.reference_start
		read_end=read.reference_end
		read_mapq=read.mapq
		read_supplementary=read.is_supplementary

		if read_mapq >= min_q:
			coverage_data=tiddit_coverage.update_coverage(read_position,read_end,bin_size,coverage_data,end_bin_size)

		if read_supplementary or read.is_secondary:
			continue

		if read_mapq > 1:
			cigar_tuple=read.cigartuples
			if (cigar_tuple[0][0] == 4 and cigar_tuple[0][1] > 10) and (cigar_tuple[-1][0] == 0 and cigar_tuple[-1][1] > 30) and len(cigar_tuple) < 7:
				clips.append([">{}|{}|{}\n".format(read.query_name,read_chromosome,read_position+1),read.query_sequence+"\n"])

			elif cigar_tuple[-1][0] == 4 and cigar_tuple[-1][1] > 10 and (cigar_tuple[0][0] == 0 and cigar_tuple[0][1] > 30) and len(cigar_tuple) < 7:
				clips.append([">{}|{}|{}\n".format(read.query_name,read_chromosome,read_position+1),read.query_sequence+"\n"])

		if read_mapq < min_q:
			continue

		if read.has_tag("SA"):
			split=SA_analysis(read,min_q,"SA",read_chromosome)
			if split:
				splits.append(split)

		if read.mate_is_unmapped:
			continue


		if ( abs(read.isize) > max_ins or mate_chromosome != read_chromosome ):
			read_query_name=read.query_name

			if mate_chromosome < read_chromosome:
				chrA=mate_chromosome
				chrB=read_chromosome
			else:
				chrA=read_chromosome
				chrB=mate_chromosome

			data.append([chrA,chrB,read_query_name,read_position+1,read_end+1,read.is_reverse,read_chromosome])

	f=open("{}_tiddit/clips/{}.fa".format(prefix,chromosome),"w")
	for clip in clips:
		f.write("".join(clip))
	f.close()

	return(chromosome,data,splits,coverage_data, "{}_tiddit/clips/{}.fa".format(prefix,chromosome) )

def main(str bam_file_name,str ref,str prefix,int min_q,int max_ins,str sample_id, int threads, int min_contig):

	samfile = pysam.AlignmentFile(bam_file_name, "r",reference_filename=ref,index_filename="{}_tiddit/{}.csi".format(prefix,sample_id))
	bam_header=samfile.header
	samfile.close()
	cdef int bin_size=50
	cdef str file_type="wig"
	cdef str outfile=prefix+".tiddit_coverage.wig"

	t_tot=0

	cdef dict data={}
	cdef dict splits={}
	cdef dict coverage_data={}
	cdef list clip_fasta=[]
	chromosomes=[]
	
	for chrA in bam_header["SQ"]:
		if chrA["LN"] < min_contig:
			continue

		chromosomes.append(chrA["SN"])
		data[chrA["SN"]]={}
		splits[chrA["SN"]]={}
		for chrB in bam_header["SQ"]:
			data[chrA["SN"]][chrB["SN"]]={}
			splits[chrA["SN"]][chrB["SN"]]={}

	t=time.time()
	res=Parallel(n_jobs=threads)( delayed(worker)(chromosome,bam_file_name,ref,prefix,min_q,max_ins,sample_id,bin_size) for chromosome in chromosomes )

	chromosomes=set(chromosomes)
	for i in range(0,len(res)):
		coverage_data[ res[i][0] ] = res[i][3]

		if not res[i][0] in chromosomes:
			continue

		for signal in res[i][1]:
			if not signal[0] in data:
				continue
	
			if not signal[2] in data[ signal[0] ][ signal[1] ]:
				data[ signal[0] ][signal[1]][signal[2]]=[]
			data[ signal[0] ][signal[1]][signal[2]].append(signal[3:])

		for signal in res[i][2]:
			if not signal[0] in splits:
				continue

			if not signal[2] in splits[ signal[0] ][ signal[1] ]:
				splits[ signal[0] ][signal[1]][signal[2]]=[]
			splits[ signal[0] ][signal[1]][signal[2]]+=signal[3:]

		clip_fasta.append(res[i][4])

	t_tot=time.time()-t

	print("total",t_tot)
	#print("coverage",t_update)
	#print("split",t_split)
	#print("disc",t_disc)

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

	f=open("{}_tiddit/clips_{}.fa".format(prefix,sample_id),"w")
	for clips in clip_fasta:
		for clip in open(clips):
			f.write(clip)
	f.close()

	return(coverage_data)
