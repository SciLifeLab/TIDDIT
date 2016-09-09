import argparse
import subprocess
import os

parser = argparse.ArgumentParser("""assemble the variants extracted through extract_variants.py, generates a sam file for each assembled variants""")
parser.add_argument('--fa',type=str,required=True,help="reference genome fasta")
parser.add_argument('--t',type=int,default=8,help="number of threads(default = 8)")
parser.add_argument('--bam',required=True,type=str,help="path to the variant bam file")
args, unknown = parser.parse_known_args()

fastq=args.bam.replace(".bam",".fastq")
process=subprocess.check_output("samtools bam2fq {} > {}".format(args.bam,fastq), shell = True)

i=0
reads={}
read=[]
filtered_fastq=args.bam.replace(".bam","_filtered.fastq")
f=open(filtered_fastq,"w")
for line in open(fastq):
        if i == 4:
                bases= set(['A', 'C', 'G', 'T'])
                if len(read[1]) > 75 and bases | set(list(read[1])) == bases:
                        for string in read:
                                f.write(string+"\n")
                i = 0
                read=[]
        i +=1
	read.append( line.strip() )
f.close()

variant_assembly=args.bam.replace(".bam",".fa")
os.system( "mpirun -np {} ABYSS -k50 {} -o {}".format(args.t,filtered_fastq,variant_assembly))

