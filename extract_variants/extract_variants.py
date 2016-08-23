import pysam
import argparse
import os
import sys

def get_regions(vcf,padding):
    variants=[]
    for line in open(vcf):
        if line[0] == "#":
            continue
        if "WINA" in line and "WINB" in line:
            a=line.split("WINA=")[-1].split(";")[0].split(",")
            a[0]=int(a[0])-padding
            if a[0] < 0:
                a[0] =0
                
            a[1]=int(a[1])+padding
            if "CHRA=" in line:
               chra=line.split("CHRA=")[-1].split(";")[0]
            else:
               chra=line.split("\t")[0]
               chrb=chra
            b=line.split("WINB=")[-1].split(";")[0].split(",")
            b[0]=int(b[0])-padding
            if b[0] < 0:
                b[0] =0
            b[1]=int(b[1])+padding
            if "CHRB=" in line:
               chrb=line.split("CHRB=")[-1].split(";")[0]
            
            #the regions should not overlap
            if(chrb == chra and b[0] <= a[1]):
                b[0] = a[1]+1
            variant={"ID":line.split("\t")[2],"startA":a[0],"endA":a[1],"chrA":chra,"startB":b[0],"endB":b[1],"chrB":chrb}
            variants.append(variant)
    
    return variants
    
def get_tab(tab,padding):
    variants=[]
    i=0
    for line in open(tab):
        if line[0] == "#":
            continue
        content=line.strip().split("\t");
        if len(content) >= 6:
            a=[ content[1],content[2] ]
            a[0]=int(a[0])-padding
            if a[0] < 0:
                a[0] =0
                
            a[1]=int(a[1])+padding
            chra=content[0]
        
            b=[content[4],content[5]]
            b[0]=int(b[0])-padding
            if b[0] < 0:
                b[0] =0
            b[1]=int(b[1])+padding
            chrb=content[3]
            
            #the regions should not overlap
            if(chrb == chra and b[0] <= a[1]):
                b[0] = a[1]+1
            variant={"ID":"Region_"+str(i),"startA":a[0],"endA":a[1],"chrA":chra,"startB":b[0],"endB":b[1],"chrB":chrb}
            variants.append(variant)
            i +=1
        else:
            a=[ content[1],int(content[2]) ]
            a[0]=int(a[0])-padding
            if a[0] < 0:
                a[0] =0
                
            a[1] += padding
            region="Region_"+str(i)
            if len(content) > 3:
                region += "_" + content[3]
                
            variant={"ID":region,"startA":a[0],"endA":a[1],"chrA":content[0],"startB":-1,"endB":-1,"chrB":-1}
            variants.append(variant)
            i +=1
            
    return variants
def main(args):
    if not os.path.exists(args.working_dir) and not args.working_dir == "":
        os.makedirs(args.working_dir)
    alignmentFile = pysam.AlignmentFile(args.bam, "rb") 
    if args.vcf:
        variants=get_regions(args.vcf,args.padding)
    elif args.tab:
        variants=get_tab(args.tab,args.padding)
    else:
        print("Error: --vcf or --tab is required")
        sys.exit()
    #print the regions to a bam file
    print "please remember these numbers!"
    for variant in variants:
        output_dir=os.path.join(args.working_dir,variant["ID"])
        with pysam.AlignmentFile(output_dir+".bam", "wb", header=alignmentFile.header) as outf:
            if not variant["chrB"] == -1 and not variant["startB"] == -1:
                print("{}|{}:{}-{}|{}:{}-{}".format(variant["ID"],variant["chrA"], variant["startA"], variant["endA"],variant["chrB"], variant["startB"], variant["endB"]))
            else:
                print("{}|{}:{}-{}".format(variant["ID"],variant["chrA"], variant["startA"], variant["endA"]))
            iter = alignmentFile.fetch(variant["chrA"], variant["startA"], variant["endA"])
            for x in iter:
                outf.write(x)
                
            if not variant["chrB"] == -1 and not variant["startB"] == -1:
                iter = alignmentFile.fetch(variant["chrB"], variant["startB"], variant["endB"])
                for x in iter:
                    outf.write(x)
                    
    #sort the file   
    for variant in variants:
        try:  
            pysam.sort(output_dir+".bam","-o" + output_dir+".sorted.bam")
        except:
            pysam.sort(output_dir+".bam",output_dir+".sorted")
        os.rename(output_dir+".sorted.bam", output_dir+".bam")
    return()


parser = argparse.ArgumentParser("""extracts the alignments of structural variants called by TIDDIT""")
parser.add_argument('--vcf',type=str,help="the path to the TIDDIT vcf file")
parser.add_argument('--tab',type=str,help="the path to to a tab file(used instead of vcf)")
parser.add_argument('--bam',type=str,required = True,help="the path to the bam file")
parser.add_argument('--working_dir',type=str,default="",help="path to the output dir(default is pwd)")
parser.add_argument('--padding',type=int,default=150,help="extend the region by this number of bases(default is 150)")

args, unknown = parser.parse_known_args()

main(args)
