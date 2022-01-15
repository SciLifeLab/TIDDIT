import sys
import os
import argparse
import math
import numpy
import time
import DBSCAN

import tiddit_stats
import tiddit_signal
import tiddit_vcf_header
import tiddit_coverage_analysis
import tiddit_cluster
import tiddit_variant
import tiddit_contig_analysis

import pysam

version="3.0.0"
parser = argparse.ArgumentParser("""TIDDIT-{}""".format(version),add_help=False)
parser.add_argument("--sv"	 , help="call structural variation", required=False, action="store_true")
parser.add_argument("--cov"        , help="generate a coverage bed file", required=False, action="store_true")
args, unknown = parser.parse_known_args()
print(args)
print(args.sv)

if args.sv == True:

	parser = argparse.ArgumentParser("""TIDDIT --sv --bam inputfile [-o prefix] --ref ref.fasta""")
	parser.add_argument('--sv'       , help="call structural variation", required=False, action="store_true")
	parser.add_argument('--bam', type=str,required=True, help="coordinate sorted bam file(required)")
	parser.add_argument('-o', type=str,default="output", help="output prefix(default=output)")
	parser.add_argument('-i', type=int, help="paired reads maximum allowed insert size. Pairs aligning on the same chr at a distance higher than this are considered candidates for SV (default= 99.9th percentile of insert size)")
	parser.add_argument('-d', type=str,help="expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default: major orientation within the dataset")
	parser.add_argument('-p', type=int,default=3, help="Minimum number of supporting pairs in order to call a variant (default 3)")
	parser.add_argument('-r', type=int,default=3, help="Minimum number of supporting split reads to call a variant (default 3)")
	parser.add_argument('-q', type=int,default=5, help="Minimum mapping quality to consider an alignment (default 5)")
	parser.add_argument('-n', type=int,default=2, help="the ploidy of the organism,(default = 2)")
	parser.add_argument('-e', type=int, help="clustering distance  parameter, discordant pairs closer than this distance are considered to belong to the same variant(default = sqrt(insert-size*2)*12)")
	parser.add_argument('-c', type=float, help="average coverage, overwrites the estimated average coverage (useful for exome or panel data)")
	parser.add_argument('-l', type=int,default=3, help="min-pts parameter (default=3),must be set >= 2")
	parser.add_argument('-s', type=int,default=25000000, help="Number of reads to sample when computing library statistics(default=25000000)")
	parser.add_argument('-z', type=int,default=50, help="minimum variant size (default=50), variants smaller than this will not be printed ( z < 10 is not recomended)")
	parser.add_argument('--force_ploidy',action="store_true", help="force the ploidy to be set to -n across the entire genome (i.e skip coverage normalisation of chromosomes)")
	parser.add_argument('--n_mask',type=float,default=0.5, help="exclude regions from coverage calculation if they contain more than this fraction of N (default = 0.5)")
	parser.add_argument('--ref', type=str, help="reference fasta",required=True)
	parser.add_argument('--bwa', type=str,default="bwa", help="path to bwa executable file(default=bwa)")
	parser.add_argument('--fermi2', type=str,default="fermi2", help="path to fermi2 executable file (default=fermi2)")
	parser.add_argument('--ropebwt2', type=str , default="ropebwt2", help="path to ropebwt2 executable file (default=ropebwt2)")
	parser.add_argument('--p_ratio', type=float,default=0.1, help="minimum discordant pair/normal pair ratio at the breakpoint junction(default=0.1)")
	parser.add_argument('--r_ratio', type=float,default=0.1, help="minimum split read/coverage ratio at the breakpoint junction(default=0.1)")
	parser.add_argument('--max_coverage', type=float,default=4, help="filter call if X times higher than chromosome average coverage (default=4)")
	args= parser.parse_args()

	if args.l < 2:
		print ("error, too low --l value!")
		quit()

	if args.ref:
		if not os.path.isfile(args.ref):
			print ("error,  could not find the reference file")
			quit()

	if not (args.bam.endswith(".bam") or args.bam.endswith(".cram")):
		print ("error, the input file is not a bam file, make sure that the file extension is .bam or .cram")
		quit()

	if not os.path.isfile(args.bam):
		print ("error,  could not find the bam file")
		quit()

	bam_file_name=args.bam
	samfile = pysam.AlignmentFile(bam_file_name, "r")
	bam_header=samfile.header
	samfile.close()


	try:
		sample_id=header["RG"][0]["SM"]
	except:
		sample_id=bam_file_name.split("/")[-1].split(".")[0]


	samples=[sample_id]

	contigs=[]
	contig_number={}
	i=0
	for contig in bam_header["SQ"]:
		contigs.append(contig["SN"])
		contig_number[contig["SN"]]=i
		i+=0

	prefix=args.o
	try:
		os.mkdir( "{}_tiddit".format(prefix) )
		os.mkdir("{}_tiddit/clips".format(prefix) )
	except:
		print("Folder already exists")
	
	min_mapq=args.q
	max_ins_len=100000
	n_reads=args.s 

	library=tiddit_stats.statistics(bam_file_name,min_mapq,max_ins_len,n_reads)
	if args.i:
		max_ins_len=args.i
	else:
		max_ins_len=library["percentile_insert_size"]


	t=time.time()
	coverage_data=tiddit_signal.main(bam_file_name,prefix,min_mapq,max_ins_len,sample_id)
	print("extracted signals in")
	print(t-time.time())

	t=time.time()
	library=tiddit_coverage_analysis.determine_ploidy(coverage_data,contigs,library,args.n,prefix,args.c)
	print("calculated coverage in")
	print(time.time()-t)

	clustered_seq={}

	f=open("{}_tiddit/clips.fa".format(prefix),"w")	
	for chr in contigs:
		clips={}
		clips[chr]=[[],[]]
		c=[]
		for line in open("{}_tiddit/clips_{}_{}.fa".format(prefix,sample_id,chr)):
			if line[0] == ">":
				c.append(line.strip())
				pos=int(line.strip().split("|")[-1])
			else:
				c.append(line.strip())
				clips[chr][0].append( "\n".join(c) )
				clips[chr][1].append([pos,0])
				c=[]
				
		
		clusters,cluster_id = DBSCAN.x_coordinate_clustering(numpy.array(clips[chr][1]),50,args.l)
		cluster_stats={}

		for i in range(0,len(clusters)):
			if clusters[i] == -1:
				continue
			if not clusters[i] in cluster_stats:
				cluster_stats[clusters[i]]=[0,[]]
			cluster_stats[clusters[i]][0]+=1
			cluster_stats[clusters[i]][1].append( clips[chr][1][i][0] )

		for i in range(0,len(clusters)):
			if clusters[i] == -1:
				continue
			if cluster_stats[clusters[i]][0] < args.r:
				continue

			if cluster_stats[clusters[i]][0] > 2*library[ "avg_coverage_{}".format(chr) ]:
				continue

			clip_coverage=coverage_data[chr][ int(math.floor(clips[chr][1][i][0]/50.0)) ]
			if clip_coverage >  args.max_coverage/2*library[ "avg_coverage_{}".format(chr) ]:
				continue

			f.write( clips[chr][0][i].strip() +"\n")

	f.close()

	os.system("{} -dNCr {}_tiddit/clips.fa | {} assemble -l 81 - > {}_tiddit/clips.fa.assembly.mag".format(args.ropebwt2,args.fermi2,prefix,prefix))
	os.system("{} simplify -COS -d 0.8 {}_tiddit/clips.fa.assembly.mag > {}_tiddit/clips.fa.assembly.clean.mag".format(args.fermi2,prefix,prefix))
	os.system("{} mem -x intractg {} {}_tiddit/clips.fa.assembly.clean.mag  > {}_tiddit/clips.sam".format(args.bwa,args.ref,prefix,prefix))
	tiddit_contig_analysis.main("{}_tiddit/clips.sam".format(prefix) , prefix, sample_id, args.z )

	vcf_header=tiddit_vcf_header.main( bam_header,library,sample_id,version )
	
	if not args.e:
		args.e=int(library["avg_insert_size"]/2.0)

	t=time.time()
	sv_clusters=tiddit_cluster.main(prefix,contigs,samples,library["mp"],args.e,args.l,max_ins_len)
	print("generated clusters in")
	print(time.time()-t)
	
	f=open(prefix+".vcf","w")
	f.write(vcf_header+"\n")
	
	t=time.time()
	#print(sv_clusters)
	variants=tiddit_variant.main(bam_file_name,sv_clusters,args,library,min_mapq,samples,coverage_data,contig_number,max_ins_len)
	print("analyzed clusters in")
	print(time.time()-t)
	
	for variant in variants:
		f.write( "\t".join(variant)+"\n" ) 
	f.close()
	quit()

elif args.cov:
	pass

else:
	parser.print_help()

