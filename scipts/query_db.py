import sys, os, glob
import argparse
import collections
import readVCF
from operator import itemgetter


def main(args):
    #start by loading the variations
    variations = args.variations
    queries = []
    with open(variations) as fin:
        noOCCTag=1;
        infoFound=0;
        for line in fin:
            #process and output the metadata
            if line.startswith("#") or line.startswith("="):
                #find the output source(cnvnator or Findtranslocations)
                meta_line=line.replace("#","");
                content=meta_line.split("=");
                if(content[0] == "source"):
                    outputSource=content[1].rstrip().split()[0];

                lookForFilter=meta_line.split("=");
                #the last infotag will be the Feature tag
                if(lookForFilter[0] != "##INFO" and noOCCTag and infoFound==1):
                    sys.stdout.write("##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n");
                    sys.stdout.write(line);
                    infoFound=0;noFeatureTag=0;
                elif(lookForFilter[0] == "##INFO"):
                    sys.stdout.write(line);
                    infoFound=1;
                    #there should only be one feature tag per vf file
                    if(line == "##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">"):
                        noOCCTag=0
                else:
                    sys.stdout.write(line)
            else:
                #in this case I need to store a query
                chrA,startA,endA,chrB,startB,endB =readVCF.readVCFLine(outputSource, line);
                current_variation = [chrA, int(startA), int(endA), chrB, int(startB), int(endB), 0, line] # plus a counter and the variatio
                queries.append(current_variation)
    # at this point queries contains an entry for each variation
    #now query each sample.db present in the given folder and store the occurences
    for sample_db in glob.glob("{}/*.db".format(os.path.abspath(args.db))):
        allVariations = {}
        with open(sample_db) as fDB:
            for line in fDB:
                db_entry = (line.rstrip().split('\t'))
                if db_entry[0] in allVariations:
                    if db_entry[1] in allVariations[db_entry[0]]:
                        allVariations[db_entry[0]][db_entry[1]].append([db_entry[2], db_entry[3], db_entry[4], db_entry[5], db_entry[6]])
                    else:
                            allVariations[db_entry[0]][db_entry[1]] = [[db_entry[2], db_entry[3], db_entry[4], db_entry[5], db_entry[6]]]
                else:
                    allVariations[db_entry[0]] = {}
                    allVariations[db_entry[0]][db_entry[1]] = [[db_entry[2], db_entry[3], db_entry[4], db_entry[5], db_entry[6]]]
    
        for query in queries:
            hit = isVariationInDB(allVariations, query)
            if hit is not None:
                query[6] += 1 # found hit

    for query in sorted(queries, key=itemgetter(6)):
        vcf_entry = query[7].rstrip()
        sys.stdout.write("{};OCC={}\n".format(vcf_entry, query[6]))

    """
    variationsByOcc = {}
    with open(args.variations) as fin:
        Variations = [ line.rstrip().split('\t') for line in fin]
	#INFO=[infoCol[7].split(";") for infoCol in Variations]
        noOCCTag=1;
        infoFound=0;
        for i in range(0,len(Variations)):
                #if we are not in the metadata
                if(Variations[i][0][0] != "#"):
                        chrA,startA,endA,chrB,startB,endB =readVCF.readVCFLine(outputSource,"\t".join(Variations[i]));
                        
                        current_variation = [chrA, int(startA), int(endA), chrB, int(startB), int(endB)]
                        hit = isVariationInDB(allVariations, current_variation)
                        hit_string =Variations[i][7]
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


                        else:


    #now print them in order of occurence
    for key in sorted(variationsByOcc):
        for var in variationsByOcc[key]:
            print "{}".format(var)

    """



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
    
    new_event = []
    found     = 0
    
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
        
        if variation_start < event_start and variation_end > event_end:
            overlap = event_end - event_start + 1
            new_event.append(variation_start)
            new_event.append(variation_end)
        #event         ---------------------
        #variaton        ---------------      take into account special cases
        elif variation_start >= event_start and variation_end <= event_end:
            overlap = variation_end - variation_start + 1
            new_event.append(event_start)
            new_event.append(event_end)
        #event           ---------------------
        #variaton     ------------------
        elif variation_start < event_start and variation_end < event_end and variation_end >= event_start: #variation_end > event_start
            overlap = variation_end - event_start + 1
            new_event.append(variation_start)
            new_event.append(event_end)
        #event         ---------------------
        #variaton           ----------------------
        elif variation_start >= event_start and variation_end > event_end and variation_start <= event_end:
            overlap = event_end - variation_start + 1
            new_event.append(event_start)
            new_event.append(variation_end)

        if overlap  > 0:
            found += 1

    if found == 2:
        return new_event
    else:
        return None




if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
    and a variation file produced by FindTranslocations or CNVnator. The script will output the variation file sorting the variation
    by number of occurences in the DB.
    """)
    parser.add_argument('--variations', type=str, required=True, help="vcf file containing variations")
    parser.add_argument('--db'        , type=str,                help="path to DB (a folder containing samples .db files")
    args = parser.parse_args()

    main(args)







