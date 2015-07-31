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

    variationsByOcc = {}
    with open(args.variations) as fin:
        Variations = [ line.rstrip().split('\t') for line in fin]
	#INFO=[infoCol[7].split(";") for infoCol in Variations]
        noOCCTag=1;
        infoFound=0;
        for i in range(0,len(Variations)):
                #if we are not in the metadata
		if(Variations[i][0][0] != "#"):
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
                #print the metadata and add the number of occurances
                else:
                        lookForFilter=Variations[i][0].split("=");
                        #the last infotag will be the Feature tag
                        if(lookForFilter[0] !="##INFO" and noOCCTag and infoFound==1):
                        
                                sys.stdout.write("##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n");
                                sys.stdout.write(Variations[i][0]+"\n");
                                infoFound=0;noFeatureTag=0;

                        elif(lookForFilter[0] == "##INFO"):
                                sys.stdout.write(Variations[i][0]+"\n");
                                infoFound=1;
                                #there should only be one feature tag per vf file
                                if(Variations[i][0] == "##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">"):
                                        noOCCTag=0;

                        else:
                                sys.stdout.write("\t".join(Variations[i])+"\n");

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







