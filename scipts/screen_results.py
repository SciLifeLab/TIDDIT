import sys, os, glob
import argparse
from operator import itemgetter

def main(args):
    #print "memorizing masked reagions"
    toBeMaskedElements = {};
    for bed_file in [item for sublist in args.bed_files for item in sublist] :
        #sys.stdout.write("memorizing {} ...".format(bed_file))
        with open(bed_file) as fin:
            rows = ( line.rstrip().split('\t') for line in fin)
            for row in rows:
                if not row[0].startswith("#"):
                    try:
			 toBeMaskedElements[row[0]].append([row[1], row[2], row[3]])

                    except :
                         toBeMaskedElements[row[0]] = [[row[1], row[2], row[3]]]
            #sys.stdout.write("done\n")

    #print "sorting masked regions"
    for chr in toBeMaskedElements:
        sorted(toBeMaskedElements[chr], key=itemgetter(1))


    #print toBeMaskedElements

    #UncompressedVariations = [];
    #CompressedVariations = []
    #currentVariation     = UncompressedVariations[0]

    with open(args.variations) as fin:

        UncompressedVariations = ( line.rstrip().split('\t') for line in fin)
        noFeatureTag=1;
        infoFound=0;
        for row in UncompressedVariations:
            #the first charachter is # in the metadata, otherwise it is c,C or a number
            if(row[0][0] != "#"):
                INFO=row[7].split(";");
	        chr_1=INFO[1].split("=")[1];
	        posA=INFO[2].split("=")[1];
                posA=posA.split(",");
	        chr_1_start=int(posA[0]);
	        chr_1_end=int(posA[1]);

	        chr_2=INFO[3].split("=")[1];
	        posB=INFO[4].split("=")[1];
	        posB=posB.split(",");
	        chr_2_start=int(posB[0]);
	        chr_2_end=int(posB[1]);


                outrow="\t".join(row);
                outRow=outrow.replace("\n","");
                sys.stdout.write(outRow)
            
                ## Chategorise the two breackpoints
                for j in range(2):
                        if j == 0:
                                variation_chr   = chr_1
                                variation_start = chr_1_start
                                variation_end   = chr_1_end
                        else:
                                variation_chr   = chr_2
                                variation_start = chr_2_start
                                variation_end   = chr_2_end
                                stopSearching = 0

                        if variation_chr in toBeMaskedElements:
                                categorized   = 0
                                i = 0
                                max_Size = 0
                                max_Size_type = ""

                                while i < len(toBeMaskedElements[variation_chr]) :
                                        repeat_start = int(toBeMaskedElements[variation_chr][i][0])
                                        repeat_end   = int(toBeMaskedElements[variation_chr][i][1])
                                        repeat_type  = toBeMaskedElements[variation_chr][i][2]
                                        i += 1
                                        overlap = 0
                                        if variation_start > repeat_end:
                                                stopSearching = 1
                                        elif variation_start <= repeat_start and variation_end >= repeat_end and variation_start < repeat_end:
                                                overlap = repeat_end - repeat_start + 1
                                        elif variation_start <= repeat_start and variation_end <= repeat_end and variation_start < repeat_end:
                                                overlap = variation_end - repeat_start + 1
                                        elif variation_start >= repeat_start and variation_end >= repeat_end and variation_start < repeat_end:
                                                overlap = repeat_end - variation_start + 1
                                        elif variation_start >= repeat_start and variation_end <= repeat_end and variation_start < repeat_end:
                                                overlap = variation_end - variation_start + 1

                                        if overlap > max_Size:
                                                max_Size       = overlap
                                                categorized = 1
                                                max_Size_type  = "{}".format(repeat_type)
                    
                                if categorized == 0:
			                if( j == 0):
                      	                        sys.stdout.write(";FEATURE=NoCategory")
			                else:
			                        sys.stdout.write(",NoCategory")
                                else:
			                if(j == 0):
                     	                        sys.stdout.write(";FEATURE={}".format(max_Size_type))
			                else:
			                        sys.stdout.write(",{}".format(max_Size_type))

    
                        else:
	                        if(j == 0):
                    	                sys.stdout.write(";FEATURE=NotPresentInRef")
		                else:
			                sys.stdout.write(",NotPresentInRef")
                sys.stdout.write("\n")
            else:
                lookForFilter=row[0].split("=");
                #the last infotag will be the Feature tag
                if(lookForFilter[0] !="##INFO" and noFeatureTag and infoFound==1):
                        
                        sys.stdout.write("##INFO=<ID=FEATURE,Number=2,Type=String,Description=\"The features of regions A and B\">\n");
                        sys.stdout.write(row[0]+"\n");
                        infoFound=0;noFeatureTag=0;
                elif(lookForFilter[0] == "##INFO"):
                        sys.stdout.write(row[0]+"\n");
                        infoFound=1;
                        #there should only be one feature tag per vf file
                        if(row[0] == "##INFO=<ID=FEATURE,Number=2,Type=String,Description=\"The features of regions A and B\">"):
                                noFeatureTag=0;
                else:
                        sys.stdout.write("\t".join(row)+"\n");
            

    return 0






if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This script takes as input one or more bed files containing genomics areas with specif features (bed entry must look like [chr start end feature]
    and a variation file produced by FindTranslocations. For each line of the variation file the script tries to match the current feature to one of the
    entries in the bed files. It adds two columns to the variation file.
    """)
    parser.add_argument('--bed-files', type=str, required=True, action='append', nargs='+', help="bed files containing the features of interest [chr start end feature]")
    parser.add_argument('--variations', help="vcf file produced by FindTransloactions containing variations (either inter of intra chromosomal)", type=str)
    args = parser.parse_args()

    main(args)



