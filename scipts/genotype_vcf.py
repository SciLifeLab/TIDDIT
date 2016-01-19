import sys, os, glob
import argparse
import readVCF
import math

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

def main(args):
    bin_size=0;
    coverage={}
    #read the coverage file
    with open(args.tab) as fin:
        for line in fin:
            if(line[0]=="#"):
                pass
            else:
                content=line.strip().split("\t")
                bin_size= int(content[2]) - int(content[1])
                if not content[0].replace("chr","").replace("CHR","") in coverage:
                    coverage[content[0].replace("chr","").replace("CHR","")]=[float(content[3])]
                else:
                    coverage[content[0].replace("chr","").replace("CHR","")].append(float(content[3]))


    #compute the coverage across each chromosome
    chr_cov=chromosome_coverage(bin_size,coverage);

    #use the coverage to genotype the CNVS
    RD=0;
    PE=0;
    SR=0;
    GT=0;
    
    format_print=1;
    with open(args.vcf) as fin:
        for line in fin:
            if(line[0]=="#" and line[1] == "#" ):
                
                if ("#FORMAT=<ID=RD,Number=1,Type=Float,Description=\"The read depth fraction at the variant\">" in line):
                    RD = 1;
                elif("#FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"The number of discordant pairs\">" in line):
                    PE = 1;
                elif("#FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"The number of split reads\">" ):
                    SR = 1;
                elif("#FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" ):
                    GT = 1;
                print(line.strip());
                
            elif(line[0] == "#"):
                out_line=line.strip()+"\tFORMAT\t"+args.sample
                
                if format_print:
                    if not RD:
                        print("##FORMAT=<ID=RD,Number=1,Type=Float,Description=\"The read depth fraction at the variant\">");
                    if not PE:
                        print("##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"The number of discordant pairs\">" );
                    if not SR:
                        print("##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"The number of split reads\">");
                    if not GT:
                        print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
                    print(out_line)
                    format_print=0;
            else:

                RD="."
                PE="0";
                SR=0;                    
                GT="./."
                    
                chrA,startA,endA,chrB,startB,endB,event_type =readVCF.readVCFLine(None, line);
                out_line=line.strip();
                try:
                    chrA_cov=chr_cov[chrA];
                    #calculate the coverage across intra chromosomal variants
                    if(chrA == chrB):
                        var_cov=mean_coverage(chrA,startA,endA,bin_size,coverage)
                        RD=( var_cov/float( chr_cov[chrA] ) )
                        if("DUP" in event_type):
                            if(RD < 1.4)*chr_cov[chrA]:
                                GT="0/1"
                            
                        elif("DEL" in event_type):
                            if(RD > 0.1*chr_cov[chrA]):
                                GT="0/1"
                            else:
                                GT="1/1"
                except:
                    pass
                        
                RD=str(RD)
                #compute the number of discordant pairs
                if(";LTE=" in line):
                    LTE=line.split(";LTE=")[1];
                    LTE=LTE.split(";")[0]    
                    PE=LTE
                #generate the format data
                format_line=[GT,PE,RD,str(SR)]
                out_line=out_line+"\t"+"GT:PE:RD:SR" + "\t" + ":".join(format_line)
                print(out_line)
if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
    and a variation file produced by FindTranslocations or CNVnator. The script will output the variation file sorting the variation
    by number of occurences in the DB.
    """)
    parser.add_argument('--vcf', type=str, required = True, help="the structural variant vcf")
    parser.add_argument('--tab'        , type=str, required = True, help="a tab file containing the binned coverage across the entire genomes(generated by FT or sambamba)")
    parser.add_argument('--sample'        , type=str, required = True,help="the name of the sample")

    args = parser.parse_args()
    main(args)







