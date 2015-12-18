import sys
import os
import argparse
import math
import itertools
import subprocess
#accepts the path to a summary file as the first command line argument, and a prefix as the second, outputs a prefix.hst file in the same folder as the input summary
parser = argparse.ArgumentParser("""this script accepts the name of some callers, aswell as database files, returns output that may be used to construct a Venn-diagram
""")
parser.add_argument('--files'        ,required = True, type=str, nargs='*', help="the callers vcf output, if a caller returns multiple files, separate them with :")
parser.add_argument('--db'        ,required = True, type=str, nargs='*', help="the database files generated from the vcf, separate them with:")
parser.add_argument('--callers', required= True , type = str, nargs = '*',help='the name of the callers')
args = parser.parse_args()



tmp_databases=[];
for files in args.files:
        file_name=files.split(':')
        tmp_databases.append(file_name)

tmp_db=[]
for files in args.db:
        file_name=files.split(':')
        tmp_db.append(file_name)
databases={};i=0


for caller in args.callers:
        databases[caller]=[tmp_databases[i],tmp_db[i]]
        i += 1

del tmp_databases;del tmp_db
n=len(args.callers)

for i in range(0,n):
        permutations=list(itertools.combinations(args.callers, i+1));
        for combination in permutations:
                caller_db=[]
                for callers in combination:
                        caller_db += (databases[callers][1])
                for vcf in databases[combination[0]][0]:
                        command=["python", "query_db.py","--variations",vcf,"--overlap","0.9","--files"]
                        command +=caller_db
                        results=subprocess.check_output(command)
                        counter=0;
                        results=results.split("\n");
                        for line in results:
                                if(";OCC="+str(len(combination)) in line):
                                        counter +=1
                        print(" ".join(combination)+" "+str(counter))

