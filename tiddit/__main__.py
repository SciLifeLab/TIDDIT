import sys
import json
import argparse
import time
import pysam
import os
import shutil
import glob

import tiddit.tiddit_stats as tiddit_stats
import tiddit.tiddit_signal as tiddit_signal
import tiddit.tiddit_vcf_header as tiddit_vcf_header
import tiddit.tiddit_coverage_analysis as tiddit_coverage_analysis
import tiddit.tiddit_coverage as tiddit_coverage
import tiddit.tiddit_cluster as tiddit_cluster
import tiddit.tiddit_variant as tiddit_variant
import tiddit.tiddit_contig_analysis as tiddit_contig_analysis

def main():
	version="3.4.0"
	parser = argparse.ArgumentParser("""tiddit-{}""".format(version),add_help=False)
	parser.add_argument("--sv"	 , help="call structural variation", required=False, action="store_true")
	parser.add_argument("--cov"        , help="generate a coverage bed file", required=False, action="store_true")
	args, unknown = parser.parse_known_args()

	if args.sv == True:

		parser = argparse.ArgumentParser("""tiddit --sv --bam inputfile [-o prefix] --ref ref.fasta""")
		parser.add_argument('--sv'       , help="call structural variation", required=False, action="store_true")
		parser.add_argument('--bam', type=str,required=True, help="coordinate sorted bam file(required)")
		parser.add_argument('-o', type=str,default="output", help="output prefix(default=output)")
		parser.add_argument('-i', type=int, help="paired reads maximum allowed insert size. Pairs aligning on the same chr at a distance higher than this are considered candidates for SV (default= 99.9th percentile of insert size)")
		parser.add_argument('-d', type=str,help="expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default: major orientation within the dataset")
		parser.add_argument('-p', type=int,default=3, help="Minimum number of supporting pairs in order to call a variant (default 3)")
		parser.add_argument('--threads', type=int,default=1, help="Number of threads (default=1)")
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
		parser.add_argument('--skip_assembly', action="store_true", help="Skip running local assembly, tiddit will perform worse, but wont require fermi2, bwa, ropebwt and bwa indexed ref")
		parser.add_argument('--p_ratio', type=float,default=0.1, help="minimum discordant pair/normal pair ratio at the breakpoint junction(default=0.1)")
		parser.add_argument('--r_ratio', type=float,default=0.1, help="minimum split read/coverage ratio at the breakpoint junction(default=0.1)")
		parser.add_argument('--max_coverage', type=float,default=4, help="filter call if X times higher than chromosome average coverage (default=4)")
		parser.add_argument('--min_contig', type=int,default=10000, help="Skip calling on small contigs (default < 10000 bp)")
		args= parser.parse_args()

		if args.l < 2:
			print ("error, too low --l value!")
			quit()

		if not args.skip_assembly:
			if not os.path.isfile(args.bwa) and not shutil.which(args.bwa):
				print("error, BWA executable missing, add BWA to path, or specify using --bwa")
				quit()

			if not os.path.isfile(args.fermi2) and not shutil.which(args.fermi2):
				print("error, fermi2 executable missing, add fermi2 to path, or specify using --fermi2")
				quit()

			if not os.path.isfile(args.ropebwt2) and not shutil.which(args.ropebwt2):
				print("error, ropebwt2 executable missing, add ropebwt2 to path, or specify using --ropebwt2")
				quit()

			if not glob.glob("{}*.bwt*".format(args.ref)):
				print ("error, The reference must be indexed using bwa index; run bwa index, or skip local assembly (--skip_assembly)")
				quit()


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
		samfile = pysam.AlignmentFile(bam_file_name, "r",reference_filename=args.ref)

		bam_header=samfile.header
		samfile.close()


		try:
			sample_id=bam_header["RG"][0]["SM"]
		except:
			sample_id=bam_file_name.split("/")[-1].split(".")[0]


		samples=[sample_id]

		contigs=[]
		contig_number={}
		contig_length={}
		i=0
		for contig in bam_header["SQ"]:
			contigs.append(contig["SN"])
			contig_number[contig["SN"]]=i
			contig_length[ contig["SN"] ]=contig["LN"]
			i+=1

		prefix=args.o
		try:
			os.mkdir( "{}_tiddit".format(prefix) )
			os.mkdir("{}_tiddit/clips".format(prefix) )
		except:
			print("Folder already exists")

		pysam.index("-c","-m","6","-@",str(args.threads),bam_file_name,"{}_tiddit/{}.csi".format(args.o,sample_id))

		min_mapq=args.q
		max_ins_len=100000
		n_reads=args.s

		library=tiddit_stats.statistics(bam_file_name,args.ref,min_mapq,max_ins_len,n_reads)
		if args.i:
			max_ins_len=args.i
		else:
			max_ins_len=library["percentile_insert_size"]


		t=time.time()
		coverage_data=tiddit_signal.main(bam_file_name,args.ref,prefix,min_mapq,max_ins_len,sample_id,args.threads,args.min_contig)
		print("extracted signals in:")
		print(t-time.time())

		t=time.time()
		library=tiddit_coverage_analysis.determine_ploidy(coverage_data,contigs,library,args.n,prefix,args.c,args.ref,50,bam_header)
		print("calculated coverage in:")
		print(time.time()-t)


		if not args.skip_assembly:

			t=time.time()
			tiddit_contig_analysis.main(prefix,sample_id,library,contigs,coverage_data,args)
			print("Clip read assembly in:")
			print(time.time()-t)

		vcf_header=tiddit_vcf_header.main( bam_header,library,sample_id,version )

		if not args.e:
			args.e=int(library["avg_insert_size"]/2.0)

		t=time.time()
		sv_clusters=tiddit_cluster.main(prefix,contigs,contig_length,samples,library["mp"],args.e,args.l,max_ins_len,args.min_contig,args.skip_assembly)

		print("generated clusters in")
		print(time.time()-t)

		f=open(prefix+".vcf","w")
		f.write(vcf_header+"\n")

		t=time.time()
		variants=tiddit_variant.main(bam_file_name,sv_clusters,args,library,min_mapq,samples,coverage_data,contig_number,max_ins_len)
		print("analyzed clusters in")
		print(time.time()-t)

		for chr in contigs:
			if not chr in variants:
				continue
			for variant in sorted(variants[chr], key=lambda x: x[0]):
				f.write( "\t".join(variant[1])+"\n" )
		f.close()
		quit()

	elif args.cov:
		parser = argparse.ArgumentParser("""tiddit --cov --bam inputfile [-o prefix]""")
		parser.add_argument('--cov'        , help="generate a coverage bed/wig file", required=False, action="store_true")
		parser.add_argument('--bam', type=str,required=True, help="coordinate sorted bam file(required)")
		parser.add_argument('-o', type=str,default="output", help="output prefix(default=output)")
		parser.add_argument('-z', type=int,default=500, help="use bins of specified size(default = 500bp) to measure the coverage of the entire bam file, set output to stdout to print to stdout")
		parser.add_argument('-w'        , help="generate wig instead of bed", required=False, action="store_true")
		parser.add_argument('-q'        , help="minimum mapping quality(default=20)", required=False, default=20)
		parser.add_argument('--ref', type=str, help="reference fasta, used for reading cram")
		args= parser.parse_args()

		if not os.path.isfile(args.bam):
			print ("error,  could not find the bam file")
			quit()

		samfile = pysam.AlignmentFile(args.bam, "r",reference_filename=args.ref)
		bam_header=samfile.header
		coverage_data,end_bin_size=tiddit_coverage.create_coverage(bam_header,args.z)
		n_reads=0
		for read in samfile.fetch(until_eof=True):

			if read.is_unmapped or read.is_duplicate:
				continue

			t=time.time()
			if read.mapq >= args.q:
				n_reads+=1

				read_position=read.reference_start
				read_end=read.reference_end
				read_reference_name=read.reference_name

				coverage_data[read_reference_name]=tiddit_coverage.update_coverage(read_position,read_end,args.z,coverage_data[read_reference_name],end_bin_size[read_reference_name])

		if args.w:
			tiddit_coverage.print_coverage(coverage_data,bam_header,args.z,"wig",args.o +".wig")
		else:
			tiddit_coverage.print_coverage(coverage_data,bam_header,args.z,"bed",args.o +".bed")

	else:
		parser.print_help()

if __name__ == '__main__':
	main()
