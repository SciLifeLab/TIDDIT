import sys, os, glob
import argparse
import collections


def main(args):
    #load the DB
    allVariations = {};
    with open(args.db) as fDB:
        Variations = ( line.rstrip().split('\t') for line in fDB if not line.startswith("#"))
        for variation in Variations:
            if variation[0] in allVariations:
                if variation[1] in allVariations[variation[0]]:
                    allVariations[variation[0]][variation[1]].append([variation[2], variation[3], variation[4], variation[5], variation[6]])
                else:
                    allVariations[variation[0]][variation[1]] = [[variation[2], variation[3], variation[4], variation[5], variation[6]]]
            else:
                allVariations[variation[0]] = {}
                allVariations[variation[0]][variation[1]] = [[variation[2], variation[3], variation[4], variation[5], variation[6]]]


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

    variationsByOcc = {}
    with open(args.variations) as fin:
        Variations = [ line.rstrip().split('\t') for line in fin if not line.startswith("#")]
	INFO=[infoCol[7].split(";") for infoCol in Variations]
        for i in range(0,len(Variations)):
		
		INFO=Variations[i][7].split(";");
		chrA=INFO[1].split("=")[1];
		posA=INFO[2].split("=")[1];
		posA=posA.split(",");
		startA=posA[0];
		endA=posA[1];

		chrB=INFO[3].split("=")[1];
		posB=INFO[4].split("=")[1];
		posB=posB.split(",");
		startB=posB[0];
		endB=posB[1];

                current_variation = [chrA, int(startA), int(endA), chrB, int(startB), int(endB)]
                hit = isVariationInDB(allVariations, current_variation)
                hit_string = ";".join(INFO)
		otherFields=Variations[i][0:7];
		otherFields="\t".join(otherFields);
                if hit == None:
                    if 0 in variationsByOcc:
                        variationsByOcc[0].append(otherFields+"\t"+"{};OCC=0".format(hit_string))
                    else:
                        variationsByOcc[0] = [otherFields+"\t"+"{};OCC=0".format(hit_string)]
                    #print "0\t{}".format(hit_string)
                else:
                    if int(hit[4]) in variationsByOcc:
                        variationsByOcc[int(hit[4])].append(otherFields+"\t"+"{};OCC={}".format(hit_string,hit[4]))
                    else:
                        variationsByOcc[int(hit[4])] = [otherFields+"\t"+"{};OCC={}".format(hit_string,hit[4])]
                    #print "{}\t{}".format(hit[4], hit_string)

    #now print them in order of occurence
    for key in sorted(variationsByOcc):
        for var in variationsByOcc[key]:
            print "{}".format(var)





def isVariationInDB(allVariations, variation):
    """
isVariationInDB requires a dictionaty like this
chrA
    chrB
        var1: [chrA, charA_start, chrA_end, chrB, chrB_start, chrB_end]
        var2: [chrA, charA_start, chrA_end, chrB, chrB_start, chrB_end]
while variation is an array of the form
[chrA, charA_start, chrA_end, chrB, chrB_start, chrB_end]

the function checks if there is an hit and returns it
    """

    chrA          = variation[0]
    chrA_start    = int(variation[1])
    chrA_end      = int(variation[2])
    chrB          = variation[3]
    chrB_start    = int(variation[4])
    chrB_end      = int(variation[5])
    
    hit = None
    if chrA in allVariations:
        # now look if chrB is here
        if chrB in allVariations[chrA]:
            # check if this variation is already present
            variationsBetweenChrAChrB = allVariations[chrA][chrB]
            hit = None
            for event in variationsBetweenChrAChrB:
                hit_tmp = isSameVariation(event, [chrA_start, chrA_end, chrB_start, chrB_end])
                if hit_tmp != None:
                    #if hit != None:
                    #    print "attention multiple hits possible case not conisedere so fa...."
                    hit = hit_tmp
    return hit



    
def isSameVariation(event, variation): #event is in the DB, variation is the new variation I want to insert
    event_chrA_start    = int(event[0])
    event_chrA_end      = int(event[1])
    event_chrB_start    = int(event[2])
    event_chrB_end      = int(event[3])

    variation_chrA_start      = int(variation[0])
    variation_chrA_end        = int(variation[1])
    variation_chrB_start      = int(variation[2])
    variation_chrB_end        = int(variation[3])

    found               = 0
    
    for j in range(2):
        if j == 0:
            variation_start    = variation_chrA_start
            variation_end      = variation_chrA_end
            event_start  = event_chrA_start
            event_end    = event_chrA_end
        else:
            variation_start    = variation_chrB_start
            variation_end      = variation_chrB_end
            event_start  = event_chrB_start
            event_end    = event_chrB_end

        overlap     = 0
        
        if variation_start <= event_start and variation_end >= event_end and variation_start < event_end:
            overlap = event_end - event_start + 1
        elif variation_start <= event_start and variation_end <= event_end and variation_start < event_end:
            overlap = variation_end - event_start + 1
        elif variation_start >= event_start and variation_end >= event_end and variation_start < event_end:
            overlap = event_end - variation_start + 1
        elif variation_start >= event_start and variation_end <= event_end and variation_start < event_end:
            overlap = variation_end - variation_start + 1

        if overlap  > 0:
            found += 1

    if found == 2:
        return event
    else:
        return None




if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
    and a variation file produced by FindTranslocations. The script will output the variation file sorting the variation
    by number of occurences in the DB.
    """)
    parser.add_argument('--variations', type=str, required=True, help="vcf file containing variations")
    parser.add_argument('--db'        , type=str,                help="File containing the database to be queryed")
    args = parser.parse_args()

    main(args)







