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

  	#Define fileformat and source
    	sys.stdout.write("##fileformat=VCFv4.1\n");
    	sys.stdout.write("##source=FindTranslocations\n");
    	#define the alowed events
    	sys.stdout.write("##ALT=<ID=DEL,Description=\"Deletion\">\n");
    	sys.stdout.write("##ALT=<ID=DUP,Description=\"Duplication\">\n");
    	sys.stdout.write("##ALT=<ID=BND,Description=\"Duplication\">\n");
    	#Define the info fields
    	sys.stdout.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
    	sys.stdout.write("##INFO=<ID=LFW,Number=1,Type=Integer,Description=\"Links from window\">\n");
    	sys.stdout.write("##INFO=<ID=LCB,Number=1,Type=Integer,Description=\"Links to chromosome B\">\n");
    	sys.stdout.write("##INFO=<ID=LTE,Number=1,Type=Integer,Description=\"Links to event\"\b>\n");
    	sys.stdout.write("##INFO=<ID=COVA,Number=1,Type=Integer,Description=\"Coverage on window A\">\n");
    	sys.stdout.write("##INFO=<ID=COVB,Number=1,Type=Integer,Description=\"Coverage on window B\">\n");
    	sys.stdout.write("##INFO=<ID=OA,Number=1,Type=Integer,Description=\"Orientation of the reads in window A\">\n");
    	sys.stdout.write("##INFO=<ID=OB,Number=1,Type=Integer,Description=\"Orientation of the mates in window B\">\n");
    	sys.stdout.write("##INFO=<ID=CHRA,Number=1,Type=String,Description=\"The chromosome of window A\">\n");
        sys.stdout.write("##INFO=<ID=CHRB,Number=1,Type=String,Description=\"The chromosome of window B\">\n");
        sys.stdout.write("##INFO=<ID=WINA,Number=2,Type=Integer,Description=\"start and stop positon of window A\">\n");
    	sys.stdout.write("##INFO=<ID=WINB,Number=2,Type=Integer,Description=\"start and stop position of window B\">\n");
    	sys.stdout.write("##INFO=<ID=EL,Number=1,Type=Integer,Description=\"Expected links to window B\">\n");
    	sys.stdout.write("##INFO=<ID=RATIO,Number=1,Type=Integer,Description=\"The number of links divided by the expected number of links\">\n");
    	sys.stdout.write("##INFO=<ID=ED,Number=1,Type=Integer,Description=\"The average estimated distance between paired ends within the window\">\n");
        sys.stdout.write("##INFO=<ID=FEATURE,Number=2,Type=String,Description=\"The features of regions A and B\">\n");
        sys.stdout.write("##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n");
    	#set filters
    	sys.stdout.write("##FILTER=<ID=BelowExpectedLinks,Description=\"The number of links between A and B is less than 40\% of the expected value\">\n");
    	sys.stdout.write("##FILTER=<ID=FewLinks,Description=\"Fewer than 40% of the links in window A link to chromosome B\">\n");
        sys.stdout.write("##FILTER=<ID=UnexpectedDistance,Description=\"The average paired reads distance is deviating\">\n");
    	sys.stdout.write("##FILTER=<ID=UnexpectedCoverage,Description=\"The coverage of the window on chromosome B or A is higher than 10*average coverage\">\n");
        #Header
    	sys.stdout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        UncompressedVariations = ( line.rstrip().split('\t') for line in fin if not line.startswith("#"))
        
        for row in UncompressedVariations:
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
                            categorized = 1
                            overlap = repeat_end - repeat_start + 1
                        elif variation_start <= repeat_start and variation_end <= repeat_end and variation_start < repeat_end:
                            categorized = 1
                            overlap = variation_end - repeat_start + 1
                        elif variation_start >= repeat_start and variation_end >= repeat_end and variation_start < repeat_end:
                            categorized = 1
                            overlap = repeat_end - variation_start + 1
                        elif variation_start >= repeat_start and variation_end <= repeat_end and variation_start < repeat_end:
                            categorized = 1
                            overlap = variation_end - variation_start + 1
                        
                        if overlap > max_Size:
                            max_Size       = overlap
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



