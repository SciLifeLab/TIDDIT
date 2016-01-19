import sys, os, glob
import argparse
import collections
import readVCF
from operator import itemgetter
import threading


def main(args):
    #start by loading the variations
    variations = args.variations
    ratio= args.overlap
    queries = []
    with open(variations) as fin:
        noOCCTag=1;
        infoFound=0;
        outputSource="Not_specified"
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
                if(lookForFilter[0] != "INFO" and noOCCTag and infoFound==1):
                    sys.stdout.write("##INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">\n");
                    sys.stdout.write("##INFO=<ID=FRQ,Number=1,Type=Integer,Description=\"The frequency of the event in the database\">\n");
                    sys.stdout.write(line);
                    infoFound=0;noFeatureTag=0;
                elif(lookForFilter[0] == "INFO"):
                    sys.stdout.write(line);
                    infoFound=1;
                    #there should only be one feature tag per vf file
                    if(line == "INFO=<ID=OCC,Number=1,Type=Integer,Description=\"The number of occurances of the event in the database\">"):
                        noOCCTag=0
                else:
                    sys.stdout.write(line)
            else:
                #in this case I need to store a query
                chrA,startA,endA,chrB,startB,endB,event_type =readVCF.readVCFLine(outputSource, line);
                current_variation = [chrA.replace("chr",""), int(startA), int(endA), chrB.replace("chr",""), int(startB), int(endB),event_type, 0, line] # plus a counter and the variation
                queries.append(current_variation)
    # at this point queries contains an entry for each variation
    #now query each sample.db present in the given folder and store the occurences
    if(args.db):
        dataBases=glob.glob("{}/*.db".format(os.path.abspath(args.db)))
    else:
        dataBases=args.files

    for sample_db in dataBases:
        allVariations = {}
        with open(sample_db) as fDB:
            for line in fDB:
                db_entry = (line.rstrip().split('\t'))
                db_entry[0] = db_entry[0].replace("chr","")
                db_entry[1] = db_entry[1].replace("chr","")
                if db_entry[0] in allVariations:
                        if db_entry[1] in allVariations[db_entry[0]]:
                            allVariations[db_entry[0]][db_entry[1]].append([db_entry[2], db_entry[3], db_entry[4], db_entry[5], db_entry[6],db_entry[7]])
                        else:
                                allVariations[db_entry[0]][db_entry[1]] = [[db_entry[2], db_entry[3], db_entry[4], db_entry[5], db_entry[6],db_entry[7]]]
                else:
                        allVariations[db_entry[0]] = {}
                        allVariations[db_entry[0]][db_entry[1]] = [[db_entry[2], db_entry[3], db_entry[4], db_entry[5], db_entry[6],db_entry[7]]]

                    
        for query in queries:
            hit = isVariationInDB(allVariations, query,ratio,args)
            if hit:
                query[7] += 1 # found hit

    for query in sorted(queries, key=itemgetter(7)):
        vcf_entry = query[8].rstrip()
        sys.stdout.write("{};OCC={};FRQ={}\n".format(vcf_entry, query[7],(query[7]/float(len(dataBases))) ) )




def isVariationInDB(allVariations, variation,ratio,args):
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
    variation_type=variation[6]

    if chrA in allVariations:
        # now look if chrB is here
        if chrB in allVariations[chrA]:
            # check if this variation is already present
            for event in allVariations[chrA][chrB]:
                #check if the variant type of the events is the same
                if(event[4] == variation_type or args.no_var):
                    #check so that the events are not disjunct on chromosome A
                    event[0]=int(event[0])
                    event[1]=int(event[1])
                    event[2]=int(event[2])
                    event[3]=int(event[3])
                    if not (chrA == chrB):
                        hit_tmp=precise_overlap(event, [chrA_start, chrA_end, chrB_start, chrB_end],args.bnd_distance)
                        if hit_tmp != None:
                            return(hit_tmp)
                    elif(event[0] <= chrA_end and event[1] >= chrA_start):
                        hit_tmp = isSameVariation(event, [chrA_start, chrA_end, chrB_start, chrB_end],ratio)
                        if hit_tmp != None:
                            return(hit_tmp)
    return None


#check the "overlap" of translocations
def precise_overlap(db_var,query_var,distance):    
    Adist=abs(db_var[0]-query_var[0]);
    Bdist=abs(db_var[2]-query_var[2]);
    if(Adist <= distance and Bdist <= distance ):
        return(True)


def isSameVariation(event, variation,ratio): #event is in the DB, variation is the new variation I want to insert
    event_chrA_start    = event[0]
    event_chrA_end      = event[1]
    event_chrB_start    = event[2]
    event_chrB_end      = event[3]
    variation_chrA_start      = variation[0]
    variation_chrA_end        = variation[1]
    variation_chrB_start      = variation[2]
    variation_chrB_end        = variation[3]

    variation_start    = variation_chrA_start
    variation_end      = variation_chrA_end
    event_start  = event_chrA_start
    event_end    = event_chrA_end

      
    overlapRatio=overlap_ratio(variation_start,variation_end,event_start,event_end)

    if overlapRatio  >= ratio:
        return(True);
    return None

#compute the total area spanned by the two events(overlaping events), calculate the intersect of the two events, return the ratio of the length of these two regions
def overlap_ratio(variation_start,variation_end,event_start,event_end):
    if(variation_start < event_start):
        region_start=variation_start;
        overlap_start=event_start;
    else:
        region_start=event_start;
        overlap_start=variation_start;

    if(variation_end > event_end):
        region_end=variation_end
        overlap_end=event_end

    else:
        region_end=event_end
        overlap_end=variation_end
    try:
        event_ratio=float(overlap_end-overlap_start+1)/float(region_end-region_start+1)
    except:
        event_ratio=0;
    return(event_ratio)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts wants as input a databse containing varaition and counts [chrA chrB startA endA startB endB occurences]
    and a variation file produced by FindTranslocations or CNVnator. The script will output the variation file sorting the variation
    by number of occurences in the DB.
    """)
    parser.add_argument('--variations', type=str, required = True, help="vcf file containing variations")
    parser.add_argument('--files'        , type=str, nargs='*', help="the paths to the db files are given as a command line arguments")
    parser.add_argument('--db'        , type=str,                help="path to DB (a folder containing samples .db files")
    parser.add_argument('--bnd_distance'        , type=int,default= 10000,help="the maximum distance between two similar precise breakpoints")
    parser.add_argument('--overlap', type=float, default = 0.6,help="the overlap required to merge two events(0 means anything that touches will be merged, 1 means that two events must be identical to be merged), default = 0.6")
    parser.add_argument('--no_var',help="count overlaping variants of different type as hits in the db", required=False, action="store_true")
    args = parser.parse_args()
    if not (args.files or args.db):
        print("a DB input method must be selected(file or db arguments)")
        quit()
    elif (args.files and args.db):
        print("only one DB input method may be chosen");
        quit()
    main(args)







