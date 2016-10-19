import argparse
import os
import sys
import random
import glob

parser = argparse.ArgumentParser("""extracts the alignments of structural variants called by TIDDIT""")
parser.add_argument('--fa',type=str,required = True,help="the fa file containing the contigs of the variant")
parser.add_argument('--TIDDIT',type=str,help="the path to TIDDIT(defaut = TIDDIT)", default="TIDDIT")
parser.add_argument('--bam',type=str,required = True,help="the path to the bam containing reads of the extracted region")
parser.add_argument('--prefix',type=str,default="assemlatron",help="path to the output dir(default is pwd)")
args, unknown = parser.parse_known_args()
contigs={}
current_header=""
contig_list=[]

for line in open(args.fa):
    if ">" in line:
        current_header=line[1:].strip().split()[0]
    else:
        contigs[current_header]={}
        
        contigs[current_header]["contig"]=line.strip()
        contig_list.append(line.strip())
os.system("bwa index {} ".format(args.fa))
os.system("samtools bam2fq  {} | bwa mem {} - | samtools view -bhS - | samtools sort - {}_contig_aln".format(args.bam,args.fa,args.prefix))
os.system("{} --cov --bin_size 1000000 -b {}_contig_aln.bam -o {}".format(args.TIDDIT, args.prefix , args.prefix))   

for line in open(args.prefix+".tab"):
    if line[0] == "#":
        continue
    content=line.split("\t")
    contigs[content[0]]["coverage"]=float(content[3])
    contigs[content[0]]["length"]=int(content[2])    
    
    
for contig in contigs:
    coverage=contigs[contig]["length"]
    
    p=0
    n=100
    for i in range(0,n):
        f=open(args.prefix+"_{}.fa".format(i),"w")   
        for test_contig in contigs:
            if test_contig == contig:
                continue
                
            f.write(">" + test_contig+"\n")
            f.write(contigs[test_contig]["contig"]+"\n")
        sequence=""
        while True:
            selected_contig=random.randint(0,len(contig_list)-1)

            selected_len=random.randint(5,40)
            selected_pos=random.randint(0,len(contig_list[selected_contig])-selected_len)
            sequence += contig_list[ selected_contig ][ selected_pos  : selected_pos+selected_len  ]
            if len(sequence) >= contigs[contig]["length"]:
                sequence=sequence[0:contigs[contig]["length"]]
                break
        f.write(">simulated_contig\n")
        f.write(sequence+"\n")
        f.close()
        os.system("bwa index {} ".format(args.prefix+"_{}.fa".format(i)))
        os.system("samtools bam2fq  {} | bwa mem {} - | samtools view -bhS - | samtools sort - {}_contig_aln_{}".format(args.bam, args.prefix+"_{}.fa".format(i) ,args.prefix,i))
        os.system("{} --cov --bin_size 1000000 -b {}_contig_aln_{}.bam -o {}_{}".format(args.TIDDIT, args.prefix,i,args.prefix,i))  
        
        for line in open("{}_{}.tab".format(args.prefix,i)):
            if line[0] == "#":
                continue
            content=line.split("\t")
            if "simulated_contig" in content[0]:
                if float(content[3]) >= contigs[contig]["coverage"]:
                    p += 1
        
        for f in glob.glob(args.prefix+"_{}*".format(i)):
            os.remove(f)
        for f in glob.glob(args.prefix+"_contig_aln_{}*".format(i)):
            os.remove(f)
    p=p/float(n)
    contigs[contig]["significant"]="True"
    contigs[contig]["p"]=p
    if p > 0.05:
        contigs[contig]["significant"]="False"
    
f=open(args.prefix+".fa","w")  
for contig in contigs:
    f.write(">" + contig+" length:{} significant:{} Pval:{} coverage:{}\n".format(contigs[contig]["length"],contigs[contig]["significant"],contigs[contig]["p"],contigs[contig]["coverage"]))
    f.write(contigs[contig]["contig"]+"\n")
f.close()
