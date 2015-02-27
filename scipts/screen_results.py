import sys, os, glob
import argparse
from operator import itemgetter

def main(args):
    #print "memorizing masked reagions"
    toBeMaskedElements = {};
    for bed_file in [item for sublist in args.bed_files for item in sublist] :
        #sys.stdout.write("memorizing {} ...".format(bed_file))
        with open(bed_file) as fin:
            rows = ( line.rstrip().split('\t') for line in fin if not line.startswith("chrA") )
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
        # 1        2          3      4       5        6            7              8              9              10              11             12
        #chrA    startOnA  endOnA   chrB  startOnB  endOnB  LinksFromWindow   LinksToChrB   LinksToEvent    CoverageOnChrA  OrientationA    OrientationB
        sys.stdout.write("chrA\tstartOnA\tendOnA\tchrB\tstartOnB\tendOnB\tLinksFromWindow\tLinksToChrB\LinksToEvent\tCoverageOnChrA\tOrientationA\tOrientationB\n")
        UncompressedVariations = ( line.rstrip().split('\t') for line in fin if not line.startswith("chrA"))
        
        for row in UncompressedVariations:
            chr_1       = row[0]
            chr_1_start = int(row[1])
            chr_1_end   = int(row[2])

            chr_2       = row[3]
            chr_2_start = int(row[4])
            chr_2_end   = int(row[5])

            LinksFromWindow     = row[6]
            LinksToChrB         = row[7]
            LinksInCurrentEvent = row[8]
            Coverage            = row[9]
            OrientationA        = row[10]
            OrientationB        = row[11]
           
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chr_1,chr_1_start,
                 chr_1_end, chr_2,  chr_2_start, chr_2_end, LinksFromWindow,  LinksToChrB, LinksInCurrentEvent,
                 Coverage , OrientationA , OrientationB))
            
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
                        sys.stdout.write("\tNoCaterogy")
                    else:
                        sys.stdout.write("\t{}".format(max_Size_type))
    
                else:
                    sys.stdout.write("\tNoPresentInRef")
            sys.stdout.write("\n")



            

    return 0






if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This script takes as input one or more bed files containing genomics areas with specif features (bed entry must look like [chr start end feature]
    and a variation file produced by FindTranslocations. For each line of the variation file the script tries to match the current feature to one of the
    entries in the bed files. It adds two columns to the variation file.
    """)
    parser.add_argument('--bed-files', type=str, required=True, action='append', nargs='+', help="bed files containing the features of interest [chr start end feature]")
    parser.add_argument('--variations', help="tabular file produced by FindTransloactions containing variations (either inter of intra chromosomal)", type=str)
    args = parser.parse_args()

    main(args)



