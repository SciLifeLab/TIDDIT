import sys, os, glob
import argparse
from operator import itemgetter

def main(args):
    allVariations = {};
    for variation_file in [item for sublist in args.variations for item in sublist] :
        with open(variation_file) as fin:
            #memorize all variations
            Variations = [ line.rstrip().split('\t') for line in fin if not line.startswith("chrA")]
            #now I need to collapse similar variations
            # 1        2          3      4       5        6            7              8              9              10              11             12
            #chrA    startOnA  endOnA   chrB  startOnB  endOnB  LinksFromWindow   LinksToChrB   LinksToEvent    CoverageOnChrA  OrientationA    OrientationB
            collapsedVariations = {}
            for variation in Variations:
                tollerance = 25000
                startA = int(variation[1]) - tollerance
                if startA < 0:
                    startA = 0
                startB = int(variation[4]) - tollerance
                if startB < 0:
                    startB = 0
                current_variation = [variation[0], startA , int(variation[2]), variation[3], startB, int(variation[5])]
                collapsedVariations = populate_DB(collapsedVariations, current_variation, True, 0)


            #now populate the DB
            for chrA in collapsedVariations:
                for chrB in collapsedVariations[chrA] :
                    for collapsedVariation in collapsedVariations[chrA][chrB]:
                        current_variation = [chrA, collapsedVariation[0],  collapsedVariation[1], chrB, collapsedVariation[2], collapsedVariation[3]]
                        allVariations = populate_DB(allVariations, current_variation, False, 0)
    

    for chrA in allVariations:
        for chrB in allVariations[chrA] :
            for event in allVariations[chrA][chrB]:
                print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrA, chrB,  event[0],  event[1],  event[2], event[3], event[4])






def populate_DB(allVariations, variation, local, tollerance):
    """
populate_DB requires a dictionaty like this
chrA
    chrB
        var1: [chrA, charA_start, chrA_end, chrB, chrB_start, chrB_end]
        var2: [chrA, charA_start, chrA_end, chrB, chrB_start, chrB_end]
while variation is an array of the form
[chrA, charA_start, chrA_end, chrB, chrB_start, chrB_end]

the function checks if this is a new element or if it is already present. 
In this case the boundaries are updtated.
If local is set to true it means we are simply collapsing a variation file,
if local is set to false it means we are building the reaDB  
    """

    chrA          = variation[0]
    chrA_start    = int(variation[1])
    chrA_end      = int(variation[2])
    chrB          = variation[3]
    chrB_start    = int(variation[4])
    chrB_end      = int(variation[5])
    
    
    if chrA in allVariations:
        # now look if chrB is here
        if chrB in allVariations[chrA]:
            # check if this variation is already present
            variationsBetweenChrAChrB = allVariations[chrA][chrB]
            found = False
            for event in variationsBetweenChrAChrB:
                new_event = mergeIfSimilar(event, [chrA_start, chrA_end, chrB_start, chrB_end])
                if  event[4] != new_event[4]:
                    event[0] = new_event[0]
                    event[1] = new_event[1]
                    event[2] = new_event[2]
                    event[3] = new_event[3]
                    if local != True: # do not change multiplicity while collapsing a set of variations
                        event[4] = new_event[4]
                    found     = True
            #otherwise add new event
            if found == False:
                startA = chrA_start - tollerance
                if startA < 0:
                    startA = 0
                startB = chrB_start - tollerance
                if startB < 0:
                    startB = 0
                allVariations[chrA][chrB].append([ startA ,chrA_end + tollerance,  startB,  chrB_end + tollerance,  1])
        else:
            startA = chrA_start - tollerance
            if startA < 0:
                startA = 0
            startB = chrB_start - tollerance
            if startB < 0:
                startB = 0
            allVariations[chrA][chrB] =  [[ startA ,chrA_end + tollerance,  startB,  chrB_end + tollerance,  1]]
        
    else:
        allVariations[chrA]       = {}
        startA = chrA_start - tollerance
        if startA < 0:
            startA = 0
        startB = chrB_start - tollerance
        if startB < 0:
            startB = 0
        allVariations[chrA][chrB] =  [[ startA ,chrA_end + tollerance,  startB,  chrB_end + tollerance,  1]]

    return allVariations
    
    
def mergeIfSimilar(event, variation): #event is in the DB, variation is the new variation I want to insert
    event_chrA_start    = int(event[0])
    event_chrA_end      = int(event[1])
    event_chrB_start    = int(event[2])
    event_chrB_end      = int(event[3])

    variation_chrA_start      = int(variation[0])
    variation_chrA_end        = int(variation[1])
    variation_chrB_start      = int(variation[2])
    variation_chrB_end        = int(variation[3])

    new_event = []

    totalOverlapSize    = 0
    totalEventSize      = 0
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

        event_size  = event_end - event_start + 1
        totalOverlapSize += event_size
        overlap     = 0
        
        
        if variation_start <= event_start and variation_end >= event_end and variation_start < event_end:
            overlap = event_end - event_start + 1
            new_event.append(variation_start)
            new_event.append(variation_end)
        elif variation_start <= event_start and variation_end <= event_end and variation_start < event_end:
            overlap = variation_end - event_start + 1
            new_event.append(variation_start)
            new_event.append(event_end)
        elif variation_start >= event_start and variation_end >= event_end and variation_start < event_end:
            overlap = event_end - variation_start + 1
            new_event.append(event_start)
            new_event.append(variation_end)
        elif variation_start >= event_start and variation_end <= event_end and variation_start < event_end:
            overlap = variation_end - variation_start + 1
            new_event.append(event_start)
            new_event.append(event_end)
        
        totalOverlapSize += overlap
        if overlap  > 0:
            found += 1

    if found == 2:
        new_event.append(event[4] +1)
        return new_event
    else:
        return event



if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts takes as input a set of variation files generated by FindTranslocations and build a databases of common varaitions it 
    first collapse similar variations belonging to the same sample and then it pushes new varaition in the DB or it updates the number
    of times a variation has been seen """)
    parser.add_argument('--variations', help="tab file containing variations", type=str,  required=True, action='append', nargs='+')
    args = parser.parse_args()

    main(args)




