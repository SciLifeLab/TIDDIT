#!/usr/bin/env python
import argparse
import os
import sys
import time

wd=os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, '{}/src/'.format(wd))
import TIDDIT_calling

version = "2.7.1"
parser = argparse.ArgumentParser("""TIDDIT-{}""".format(version),add_help=False)
parser.add_argument('--sv'       , help="call structural variation", required=False, action="store_true")
parser.add_argument('--cov'        , help="generate a coverage bed file", required=False, action="store_true")
args, unknown = parser.parse_known_args()

if args.sv:
	parser = argparse.ArgumentParser("""TIDDIT --sv --bam inputfile [-o prefix] --ref ref.fasta""")
	parser.add_argument('--sv'       , help="call structural variation", required=False, action="store_true")
	parser.add_argument('--bam', type=str,required=True, help="coordinate sorted bam file(required)")
	parser.add_argument('-o', type=str,default="output", help="output prefix(default=output)")
	parser.add_argument('-i', type=int, help="paired reads maximum allowed insert size. Pairs aligning on the same chr at a distance higher than this are considered candidates for SV (default=3std + mean_insert_size)")
	parser.add_argument('-d', type=str,help="expected reads orientations, possible values \"innie\" (-> <-) or \"outtie\" (<- ->). Default: major orientation within the dataset")
	parser.add_argument('-p', type=int,default=3, help="Minimum number of supporting pairs in order to call a variation event (default 3)")
	parser.add_argument('-r', type=int,default=3, help="Minimum number of supporting split reads to call a small variant (default 3)")
	parser.add_argument('-q', type=int,default=10, help="Minimum mapping quality to consider an alignment (default 10)")
	parser.add_argument('-Q', type=int,default=20, help="Minimum regional mapping quality (default 20)")
	parser.add_argument('-n', type=int,default=2, help="the ploidy of the organism,(default = 2)")
	parser.add_argument('-e', type=int, help="clustering distance  parameter, discordant pairs closer than this distance are considered to belong to the same variant(default = sqrt(insert-size*2)*12)")
	parser.add_argument('-l', type=int,default=3, help="min-pts parameter (default=3),must be set > 2")
	parser.add_argument('-s', type=int,default=50000000, help="Number of reads to sample when computing library statistics(default=50000000)")
	parser.add_argument('-z', type=int,default=100, help="minimum variant size (default=100), variants smaller than this will not be printed( z < 10 is not recomended)")
	parser.add_argument('--force_ploidy',action="store_true", help="force the ploidy to be set to -n across the entire genome (i.e skip coverage normalisation of chromosomes)")
	parser.add_argument('--debug',action="store_true", help="rerun the tiddit clustering procedure")
	parser.add_argument('--n_mask',type=float,default=0.5, help="exclude regions from coverage calculation if they contain more than this fraction of N (default = 0.5)")
	parser.add_argument('--ref', type=str, help="reference fasta, used for GC correction")

	args= parser.parse_args()
	args.wd=os.path.dirname(os.path.realpath(__file__))
	if args.l < 2:
		print ("error, too low --l value!")
		quit()
	if args.ref:
		if not os.path.isfile(args.ref):
			print ("error,  could not find the reference file")
			quit()
	if not os.path.isfile(args.bam):
		print ("error,  could not find the bam file")
		quit()
	if not os.path.isfile("{}/bin/TIDDIT".format(args.wd)):
		print ("error,  could not find the TIDDIT executable file, try rerun the INSTALL.sh script")
		quit()

	command_str="{}/bin/TIDDIT --sv -b {} -o {} -p {} -r {} -q {} -n {} -s {}".format(args.wd,args.bam,args.o,args.p,args.r,args.q,args.n,args.s)
	if args.i:
		command_str += " -i {}".format(args.i)
	if args.d:
		command_str += " -d {}".format(args.d)

	if not args.debug:
		os.system(command_str)

		if args.ref:
			t=time.time()
			print ("Generating GC wig file")
			if args.ref.endswith(".gz"):
				os.system("zcat {} | {}/bin/TIDDIT --gc -z 50 -o {}".format(args.ref,args.wd,args.o))
			else:
				os.system("cat {} | {}/bin/TIDDIT --gc -z 50 -o {}".format(args.ref,args.wd,args.o))
			print ("Constructed GC wig in {} sec".format(time.time()-t))

	TIDDIT_calling.cluster(args)

elif args.cov:
	parser = argparse.ArgumentParser("""TIDDIT --cov --bam inputfile [-o prefix]""")
	parser.add_argument('--cov'        , help="generate a coverage bed/wig file", required=False, action="store_true")
	parser.add_argument('--bam', type=str,required=True, help="coordinate sorted bam file(required)")
	parser.add_argument('-o', type=str,default="output", help="output prefix(default=output)")
	parser.add_argument('-z', type=int,default=500, help="use bins of specified size(default = 500bp) to measure the coverage of the entire bam file, set output to stdout to print to stdout")
	parser.add_argument('-w'        , help="generate wig instead of bed", required=False, action="store_true")
	parser.add_argument('-u'        , help="skip per bin mapping quality", required=False, action="store_true")
	args= parser.parse_args()
	args.wd=os.path.dirname(os.path.realpath(__file__))
	command="{}/bin/TIDDIT --cov -b {} -o {} -z {}".format(args.wd,args.bam,args.o,args.z)
	if args.w:
		command += " -w"

	if args.u:
		command += " -u"

	if not os.path.isfile(args.bam):
		print ("error,  could not find the bam file")
		quit()
	if not os.path.isfile("{}/bin/TIDDIT".format(args.wd)):
		print ("error,  could not find the TIDDIT executable file, try rerun the INSTALL.sh script")
		quit()
	os.system(command)

else:
	parser.print_help()
