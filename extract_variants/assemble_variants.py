import argparse
import os



parser = argparse.ArgumentParser("""assemble the variants extracted through extract_variants.py, generates a sam file for each assembled variants""")
parser.add_argument('--vcf',type=str,required = True,help="the path to the TIDDIT vcf file")
parser.add_argument('--fa',type=str,required=True,help="reference genome fasta")
parser.add_argument('--t',type=int,default=8,help="number of threads(default = 8)")
parser.add_argument('--working_dir',required=True,type=str,default="",help="path to the folder containing the bam files of extracted variants")
args, unknown = parser.parse_known_args()

for line in open(args.vcf):
    if "#" == line[0]:
        continue
    if not "WINA" in line and not "WINB" in line:
        continue
    input_bam=os.path.join( args.working_dir,line.split("\t")[2]+".bam" )
    fastq=os.path.join( args.working_dir,line.split("\t")[2]+".fastq" )
    process="samtools bam2fq {} > {}".format(input_bam,fastq)
    os.system(process)
    
    variant_assembly=os.path.join( args.working_dir,line.split("\t")[2]+".fa" )
    process="mpirun -np {} ABYSS -k50 {} -o {}".format(args.t,fastq,variant_assembly)
    os.system(process)
    
    blat_results=os.path.join( args.working_dir,line.split("\t")[2]+".sam" )
    process="bwa bwasw -C {} {} -t {} > {}".format(args.fa,variant_assembly,args.t,blat_results)
    print process
    os.system(process)

