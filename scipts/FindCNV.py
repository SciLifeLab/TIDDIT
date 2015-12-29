import sys, os, glob
import argparse
import math
import random
import itertools
import time

def mean_coverage(chromosome,start,end,bin_size,coverage_dict):
    start_bin=int(math.floor(start/float(bin_size)))
    end_bin=int(math.floor(end/float(bin_size)))
    
    cov = 0
    n=0;
    for i in range(start_bin,end_bin+1):
        cov += coverage_dict[chromosome][i];
        n+=1;
    cov=cov/float(n);
    return cov;

def chromosome_coverage(bin_size,coverage_dict):
    cov={}
    for chromosome in coverage_dict:
        cov[chromosome]=0;
        n=0;
        for bins in coverage_dict[chromosome]:
            cov[chromosome]+=bins;
            n+=1;
        cov[chromosome]=cov[chromosome]/float(n)
        #print(chromosome + "\t" + str(cov[chromosome]))
        
    return cov;
    
def perm_test(coverage,chromosome_cov,idx,n):
    p =0;
    for cov in idx:
        if(chromosome_cov[cov] <= coverage):
            p+=1
    p=p/float(n)
    return p
    
def perm_test_max(coverage,chromosome_cov,idx,n):
    p =0;
    for cov in idx:
        if(chromosome_cov[cov] >= coverage):
            p+=1
    p=p/float(n)
    return p


def main(args):
    bin_size=0;
    coverage={}
    
    #read the coverage file
    chromosomes=[]
    with open(args.tab) as fin:
        for line in fin:
            if(line[0]=="#"):
                pass
            else:
                content=line.strip().split("\t")
                if not content[0] == "":
                    if content[0] not in chromosomes:
                        chromosomes.append(content[0])
                    bin_size= int(content[2]) - int(content[1])
                    if not content[0] in coverage:
                        coverage[content[0]]=[round(float(content[3]),args.precision)]
                    else:
                        coverage[content[0]].append(round(float(content[3]),args.precision))
    #for chromosome in chromosomes:
    #    print chromosome
    
    #compute the coverage across each chromosome
    chr_cov=chromosome_coverage(bin_size,coverage);
    mean_chr_cov={}
    for chromosome in coverage:
      n=0;
      cov=0;
      for bin in coverage[chromosome]:
           cov+=bin
           n+=1
      mean_chr_cov[chromosome]=round(cov/float(n),args.precision)
      #print(chromosome + " " + str(mean_chr_cov[chromosome]))
       
       
    n=50000  
    P=0.025
    lower_lim={}
    upper_lim={}
    #use binary testing and permutations to decide the coverage limits to call a variant
    for chromosome in chromosomes:
      cov_set=set(coverage[chromosome])
      p=-1;
      N1=min(cov_set)
      
      chromosome_cov = [s for s in coverage[chromosome] if s > 0 and s  <  3*mean_chr_cov[chromosome] ]
      size = len(chromosome_cov)
      _random, _int = random.random, int
      idx=[_int(_random() * size) for i in itertools.repeat(None, n)]
      
      #print(chromosome)
      while p < P and cov_set:
        p = perm_test(N1,chromosome_cov,idx,n)
        cov_set=cov_set-set([N1])
        if len(cov_set) > 0:
            N1=min(cov_set)
      lower_lim[chromosome]=N1
      cov_set=set(coverage[chromosome])
      
      #print(N1)
      p=2;
      N1=min(cov_set)
      upper_lim[chromosome]=N1
      while p >= P and cov_set:
        p = perm_test_max(N1,chromosome_cov,idx,n)
        cov_set=cov_set-set([N1])
        try:
            N1=min(cov_set)
        except:
            pass
      upper_lim[chromosome]=N1
      #print(N1)
      
      #print(lower_lim[chromosome])
    #search for cnvs
    variant={};
    for chromosome in chromosomes:
        #print(chromosome)
        var=[]
        length=0;
        variant[chromosome]=[]
        for i in range(0,len(coverage[chromosome])):
            
            current_alt="DEL";
            if coverage[chromosome][i] > mean_chr_cov[chromosome]:
                current_alt="DUP"
                
            if coverage[chromosome][i] < lower_lim[chromosome] or coverage[chromosome][i] > upper_lim[chromosome]:               
                 if var:
                    #if a dup is close to a del
                    if not (current_alt == var[3]):
                        if length >= args.nbins:
                            var_list=var+[str(var_cov)]+[str(length)]
                            variant[chromosome].append(var_list)             
                        
                        var=[str(chromosome),str(i*bin_size),str((i*bin_size)+bin_size)]
                        alt=current_alt
                        length=1;
                        var_cov=coverage[chromosome][i]
                        var.append(current_alt)
                    #extend the old variant    
                    else:
                        var_cov=float(coverage[chromosome][i]+var_cov*length)/float(length+1)
                        var[2]=str((i*bin_size)+bin_size)
                        length+=1
                 #store a new variant
                 else:
                    var=[str(chromosome),str(i*bin_size),str((i*bin_size)+bin_size)]
                    alt=current_alt
                    length=1;
                    var_cov=coverage[chromosome][i]
                    var.append(current_alt)
                    
            #if the last bin of a  variant is found
            else: 
                if var:
                    if length >= args.nbins:
                        var_list=var+[str(var_cov)]+[str(length)]
                        variant[chromosome].append(var_list)
                    var = [];
                    var_cov=0;
                    length=0
                    
    #merge closely located bins              
    if(args.merge > 0):
        for chromosome in chromosomes:
            i=0;
            while i < len(variant[chromosome])-1:
                diff=(   float(variant[chromosome][i+1][1])-float(variant[chromosome][i][2])  ) /float(bin_size)
                if(diff <= args.merge and variant[chromosome][i][3] == variant[chromosome][i+1][3]):
                    variant[chromosome][i][2] = variant[chromosome][i+1][2]
                    variant[chromosome][i][4] =str( round( (float(variant[chromosome][i][4])*float(variant[chromosome][i][-1])+float(variant[chromosome][i+1][-1])*float(variant[chromosome][i+1][4]))/(float(variant[chromosome][i+1][-1])+float(variant[chromosome][i][-1])) ,args.precision) )
                    variant[chromosome][i][-1] = str( float(variant[chromosome][i+1][-1])+diff + float(variant[chromosome][i][-1]) )
                    del variant[chromosome][i+1]
                else:
                    i  += 1

    #print the data
    print("##fileformat=VCFv4.1")
    print("##source=FindCNV")
    print("##ALT=<ID=DEL,Description=\"Deletion>")
    print("##ALT=<ID=DUP,Description=\"Duplication\">")
    print("##INFO=<ID=RD,Number=1,Type=float,Description=\"The read depth of the variant\">")
    print("##INFO=<ID=RDR,Number=1,Type=float,Description=\"The read depth/chromosome read depth\">")
    print("##INFO=<ID=END,Number=1,Type=float,Description=\"The end position of the variant\">")
    print("##INFO=<ID=SVLEN,Number=1,Type=float,Description=\"The length of the variant\">")
    print("##INFO=<ID=SVTYPE,Number=1,Type=float,Description=\"The variant type\">")
    print("##INFO=<ID=BINS,Number=1,Type=float,Description=\"The number of bins used to call the variant\">")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    i=1;                 
    for chromosome in chromosomes:
        for line in variant[chromosome]:
            out_line=line[0]+"\t"+line[1]+"\t"+str(i)+"\t.\t"+"<"+line[3]+">"+"\t."+"\tPASS\t"+"END="+line[2]+";SVLEN="+str(int(float(line[-1])*bin_size))+";"
            out_line +="BINS="+line[-1]+";"+"RD="+line[4]+";RDR="+str(float(line[4])/float(mean_chr_cov[line[0]]))+";SVTYPE="+line[3]
            i += 1
            print(out_line)
parser = argparse.ArgumentParser("""
This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
and a variation file produced by FindTranslocations or CNVnator. The script will output the variation file sorting the variation
by number of occurences in the DB.
""")
parser.add_argument('--tab' , type=str, required = True, help="a tab file containing the binned coverage across the entire genomes(generated by FT or sambamba)")
parser.add_argument('--nbins' , type=int,default=2, help="the required numbe rof bins for calling a variant")
parser.add_argument('--merge' , type=int,default=2, help="merge bins of the same variant type, separated by n normal bins(default = 0)")
parser.add_argument('--precision' , type=int,default=0, help="round the coverage of each bin to reduce to noise,i.e the number of significant figures(default = 0)")
args = parser.parse_args()
#tic=time.clock()
main(args)
#print(time.clock()-tic)







