import sys, os, glob
import argparse
import readVCF
import random
from operator import itemgetter
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import copy
import subprocess

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'

#open the DB
def getDB(path):

    DB=[]
    variants=[]
    for line in open(path):
        data=line.split("\t");
        DB.append(data[0])
        variants.append(data[1]);


    return(DB,variants);

#get db and vcf via the path to their folders
def getDB_VCF(db_path,vcf_path):
    file_dict={}    
    files=glob.glob(os.path.join(db_path,"*.db"));
    #add all the .db files into the dictionary of path to the files
    for file in files:
        file_dict[file.split("/")[-1].strip(".db")]=[file,[]]
    files=glob.glob(os.path.join(vcf_path,"*.vcf"));
    #add the path to all vcf that have  db
    for file in files:
        for key in file_dict:
            if key in file:
                file_dict[key][1].append(file)
    #if a db file lacks a vcf file, remove it
    for file in file_dict:
        if( len(file_dict[file][1]) == 0):
            del file_dict[file];

    #add the paths to the DB and variant lists
    DB=[];variants=[];
    for file in file_dict:
        DB.append(file_dict[file][0]);variants.append(";".join(file_dict[file][1]));
    return(DB,variants);

#randomly picks N/(i+1) samples from the list
def getSample(DB,N):

    positions=range(0,len(DB))
    selection=random.sample(positions,N)

    tmpDB=[]
    for j in range(0,N):
        tmpDB.append(DB[selection[j]])
    return(tmpDB)
        
def main(dataBases,variations):
    frq={};
    command=["python","query_db.py","--variations",variations,"--files"]
    command=command+dataBases
    query=subprocess.check_output(command);
    query=query.strip().split("\n");
    for line in query:
        if not (line.startswith("#")):
            hits=int(line.split("OCC=")[-1].split(";")[0]);
            if hits in frq:
                frq[hits] +=1
            else:
                frq[hits] = 1
    return(frq);

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
    and a variation file produced by FindTranslocations or CNVnator. The script will output the variation file sorting the variation
    by number of occurences in the DB.
    """)
    parser.add_argument('--file'        , type=str,default =None, help="the inputfile, use either this option or the vcf,db options")
    parser.add_argument('--vcf'        , type=str, default =None, help="a path to the folder containing vcf files. the prefix of the vcf files is assumed do be the same as the prefix of the db files")
    parser.add_argument('--db'        , type=str, default =None,  help="a path to the folder containing db files. the prefix of the vcf files is assumed do be the same as the prefix of the vcf file")
    parser.add_argument("--n"   ,type=int, help="the number of samples, default = all samples")
    parser.add_argument("--N"   ,type=int, default=1, help="Itterations per sample,default = 1")
    parser.add_argument("--range"   ,type=int, default=None, help="histogram x axis range, default = number of samples")
    args = parser.parse_args()

    VAR=[];DATABASE=[]
    if(args.file and not args.vcf and not args.db):
        DATABASE,VAR=getDB(args.file)
    if(args.vcf and args.db):
        DATABASE,VAR=getDB_VCF(args.db,args.vcf)
    else:
        print("error");

    x=None
    if args.n:
        x=args.n;
    else:
        x=len(VAR)
    if not args.range:
        args.range=x
    
    frequencies={}
    for j in range(0,x):

       for k in range(0,args.N):
           tmpDATABASE=copy.copy(DATABASE)
           del tmpDATABASE[j]
           DB=getSample(tmpDATABASE,x-1)
           DB.append(DATABASE[j]);
           print("query " + str(j+1) + " , itteration " + str(k+1))
           vcf=VAR[j].rstrip()
           vcf=vcf.split(";");
           for vcfFile in vcf:
               #print(DB)
               V=main(DB,vcfFile)
               for element in V:
                   if element in frequencies:
                        frequencies[element] += V[element]
                   else:
                        frequencies[element] = V[element]
 
           #print(occurances)
    occurances =[]
    total_hits=0;
    for element in frequencies:
        total_hits += frequencies[element]
        occurances += [element]*frequencies[element]
    occurances = np.asarray(occurances)
    plt.hist(occurances, bins = range(0,x+2),normed=True)
    formatter = FuncFormatter(to_percent)
    #plt.xlim(0,20)
    plt.gca().yaxis.set_major_formatter(formatter)
    print("percentage of unique hits:" + str(frequencies[1]/float(total_hits)))
    plt.title("the percentage number of hits in a database of " + str(x) + " samples")
    plt.xlabel("Number of hits")
    plt.ylabel("Percentage")       
    plt.show()
            





