import pysam
from subprocess import Popen, PIPE, DEVNULL
from joblib import Parallel, delayed

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment

import numpy
import math
import os

import time

import tiddit.DBSCAN as DBSCAN
import tiddit.silverfish as silverfish
import tiddit.tiddit_signal as tiddit_signal


def read_contigs(aligned_contigs,prefix,sample_id,min_size):
	samfile = pysam.AlignmentFile(aligned_contigs, "r")
	bam_header=samfile.header
	

	split_contigs={}
	for chrA in bam_header["SQ"]:
		for chrB in bam_header["SQ"]:
			if chrA["SN"] <= chrB["SN"]:

				if not chrA["SN"] in split_contigs:
					split_contigs[chrA["SN"]] = {}

				split_contigs[chrA["SN"]][chrB["SN"]]={}

	for read in samfile.fetch(until_eof=True):

		if read.is_unmapped:
			continue

		if read.has_tag("SA") and not (read.is_supplementary or read.is_secondary):
			split=tiddit_signal.SA_analysis(read,-2,"SA",read.reference_name)

			if split:
				if not split[2] in split_contigs[split[0]][split[1]]:
					split_contigs[split[0]][split[1]][split[2]]=[]
				split_contigs[split[0]][split[1]][split[2]]+=split[3:]

			
		elif read.has_tag("XA") and not (read.is_supplementary or read.is_secondary):
			XA=read.get_tag("XA")
			if XA.count(";") == 1:
				if ",-" in XA:
					XA=XA.replace(",-",",")
					xa_list=XA.split(",")
					xa_list.insert(2,"-")
					XA=",".join(xa_list)
				else:
					XA=XA.replace(",+",",")
					xa_list=XA.split(",")
					xa_list.insert(2,"+")
					XA=",".join(xa_list)

				read.set_tag("XA",XA)
				split=tiddit_signal.SA_analysis(read,-2,"XA",read.reference_name)

				if split:
					if not split[2] in split_contigs[split[0]][split[1]]:
						split_contigs[split[0]][split[1]][split[2]]=[]
					split_contigs[split[0]][split[1]][split[2]]+=split[3:]

		elif not (read.is_supplementary or read.is_secondary) and len(read.cigartuples) > 2:

			current_bp=read.reference_start
			for i in range(0,len(read.cigartuples)-1):
				if read.cigartuples[i][0] == 2 and read.cigartuples[i][1] > min_size:
					
					split_contigs[read.reference_name][read.reference_name]["{}_d_{}".format(read.query_name,i)]=[current_bp,read.is_reverse,current_bp+read.cigartuples[i][1],read.is_reverse,read.reference_start,current_bp,current_bp+read.cigartuples[i][1],read.reference_end]
				current_bp+=read.cigartuples[i][1]

	f=open("{}_tiddit/contigs_{}.tab".format(prefix,sample_id),"w")
	positions=set([])
	for chrA in split_contigs:
		for chrB in split_contigs[chrA]:
			for fragment in split_contigs[chrA][chrB]:

				p=(chrA,chrB,split_contigs[chrA][chrB][fragment][0],split_contigs[chrA][chrB][fragment][2])
				if p in positions:
					continue

				f.write("{}\t{}\t{}\t{}\n".format(fragment,chrA,chrB,"\t".join(map(str, split_contigs[chrA][chrB][fragment] )))  )
				positions.add(p)

	f.close()

def local_assembly(args,sample_id,prefix,regions,chr):

	if os.path.isfile("{}_tiddit/clips/clips.fa.assembly.{}.clean.mag".format(prefix,chr)):
		os.remove("{}_tiddit/clips/clips.fa.assembly.{}.clean.mag".format(prefix,chr))

	cdef AlignmentFile samfile  = AlignmentFile(args.bam, "r",reference_filename=args.ref,index_filename="{}_tiddit/{}.csi".format(args.o,sample_id))
	mag=open( "{}_tiddit/clips/clips.fa.assembly.{}.clean.mag".format(prefix,chr) ,"w")
	contig=1
	for region in regions[chr]:

		n_reads=0
		proper=0
		low_mapq=0

		if region[2]-region[1] > args.max_local_assembly_region:
			continue

		reads={}
		for read in samfile.fetch(region[0], region[1], region[2]):
			if read.is_supplementary or read.is_duplicate or read.is_secondary:
				continue
			n_reads+=1
			if read.mapq < 10:
				low_mapq+=1
			if read.is_proper_pair:
				proper+=1

			reads[str(n_reads)]=read.query_sequence
			if n_reads > 50000:
				break


		if n_reads > args.max_assembly_reads:
			continue

		if low_mapq/n_reads > 0.25 or proper/n_reads < 0.75:
			continue
	
		results=silverfish.main(reads,args.k,args.min_clip)
		del reads

		for result in results:
			if len(result) > args.min_contig_len:
				mag.write(f">{chr}_{region[1]}_{region[2]}_{contig}\n")
				mag.write(result+"\n")
			contig+=1

	mag.close()
	return( "{}_tiddit/clips/clips.fa.assembly.{}.clean.mag".format(prefix,chr) )

def main(prefix,sample_id,library,contigs,coverage_data,args):

	clips={}
	c=[]

	for line in open("{}_tiddit/clips_{}.fa".format(prefix,sample_id)):

		if line[0] == ">":
			c.append(line.strip())
			pos=int(line.strip().split("|")[-1])
			chr=line.strip().split("|")[-2]
			if not chr in clips:
				clips[chr]=[[],[]]

		else:
			c.append(line.strip())
			clips[chr][0].append( "\n".join(c) )
			clips[chr][1].append([pos,0])
			c=[]

	regions={}


	assembly_l=args.min_pts_clips
	
	for chr in clips:
		regions[chr]=[]


		l=assembly_l
		if library[ "avg_coverage_{}".format(chr) ]/library["avg_coverage"] > 5:
			l=args.l*int(round(library[ "avg_coverage_{}".format(chr) ]/library["avg_coverage"]/2.0))

		clusters,cluster_id = DBSCAN.x_coordinate_clustering(numpy.array(clips[chr][1]),50,l)
		cluster_stats={}

		for i in range(0,len(clusters)):
			if clusters[i] == -1:
				continue
			if not clusters[i] in cluster_stats:
				cluster_stats[clusters[i]]=[0,[]]
			cluster_stats[clusters[i]][0]+=1
			cluster_stats[clusters[i]][1].append( clips[chr][1][i][0] )

		for cluster in cluster_stats:

			if cluster_stats[cluster][0] < args.min_clip:
				continue

			clip_coverage=max(coverage_data[chr][ int(math.floor(min(cluster_stats[cluster][1])/50.0)):int(math.floor(max(cluster_stats[cluster][1])/50.0))+1 ])

			if clip_coverage/library[ "avg_coverage_{}".format(chr) ] > args.max_coverage:
				continue

			regions[chr].append([chr,min(cluster_stats[cluster][1] )-args.padding,max(cluster_stats[cluster][1])+args.padding])

			if regions[chr][-1][1] < 1:
				regions[chr][-1][1]=1

	del clips

	contigs=Parallel(n_jobs=args.threads,timeout=99999)( delayed(local_assembly)(args,sample_id,prefix,regions,chr) for chr in regions)

	mag=open(f"{prefix}_tiddit/clips.fa.assembly.clean.mag","w")
	for contig in contigs:
		for line in open(contig):
			mag.write(line.rstrip()+"\n")
	mag.close()

	os.system("{} mem -t {} -x intractg {} {}_tiddit/clips.fa.assembly.clean.mag  1> {}_tiddit/clips.sam 2> /dev/null".format(args.bwa,args.threads,args.ref,prefix,prefix))

	read_contigs("{}_tiddit/clips.sam".format(prefix) , prefix, sample_id, args.z)
