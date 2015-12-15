import sys, os, glob
import argparse
from operator import itemgetter
import readVCF
import subprocess

def main(args):
    allVariations       = {}
    tollerance    = args.tollerance
    fixed =0
    if args.fixed:
        fixed =1
    for variation_file in [item for sublist in args.variations for item in sublist] :
        outputSource=None;
        collapsedVariations = {} # this will contain the SVs of this file, but collapsing those that are close
        with open(variation_file) as fin:
            #memorize all variations
            for line in fin:
                if line.startswith("#") or line.startswith("="):
                        #find the output source(cnvnator or Findtranslocations
                        line=line.replace("#","");
                        content=line.split("=");
                        if(content[0] == "source"):
                                outputSource=content[1].rstrip().split()[0];
                        continue
                
                chrA,startA,endA,chrB,startB,endB,event_type =readVCF.readVCFLine(outputSource,line);

                startA = startA - tollerance
                if startA < 0:
                        startA = 0
                startB = startB - tollerance
                if startB < 0:
                        startB = 0

                current_variation = [chrA, startA , endA + tollerance, chrB, startB, endB + tollerance,event_type]
                collapsedVariations = populate_DB(collapsedVariations, current_variation, True, 0, fixed )

        ##collapse again in order to avoid problems with areas that have become too close one to onther
        elemnets_before = 0
        for chrA in collapsedVariations:
            for chrB in collapsedVariations[chrA] :
                for collapsedVariation in collapsedVariations[chrA][chrB]:
                    elemnets_before += 1
        elemnets_after  = 0

        while elemnets_before != elemnets_after:
            collapsedVariationsFinal = {}
            elemnets_before = 0
            elemnets_after  = 0
            for chrA in collapsedVariations:
                for chrB in collapsedVariations[chrA] :
                    for collapsedVariation in collapsedVariations[chrA][chrB]:
                        current_variation = [chrA, collapsedVariation[0],  collapsedVariation[1], chrB, collapsedVariation[2], collapsedVariation[3],collapsedVariation[4]]
                        collapsedVariationsFinal = populate_DB(collapsedVariationsFinal, current_variation, True, 0,fixed)
                        elemnets_before += 1
            for chrA in collapsedVariationsFinal:
                for chrB in collapsedVariationsFinal[chrA] :
                    for collapsedVariation in collapsedVariationsFinal[chrA][chrB]:
                        elemnets_after += 1
            collapsedVariations.clear()
            collapsedVariations = collapsedVariationsFinal


        #now populate the DB
        for chrA in collapsedVariations:
            for chrB in collapsedVariations[chrA] :
                for collapsedVariation in collapsedVariations[chrA][chrB]:
                    current_variation = [chrA, collapsedVariation[0],  collapsedVariation[1], chrB, collapsedVariation[2], collapsedVariation[3],collapsedVariation[4]]
                    allVariations = populate_DB(allVariations, current_variation, False, 0,fixed)
                        
        
    

    for chrA in allVariations:
        for chrB in allVariations[chrA] :
            for event in allVariations[chrA][chrB]:
                print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrA, chrB,  event[0],  event[1],  event[2], event[3], event[4], event[5])



def populate_DB(allVariations, variation, local, tollerance,fixed):
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
    event_type    = variation[6];
    
    
    if chrA in allVariations:
        # now look if chrB is here
        if chrB in allVariations[chrA]:
            if not fixed:
                # check if this variation is already present
                variationsBetweenChrAChrB = allVariations[chrA][chrB]
                found = False
                for event in variationsBetweenChrAChrB:
                    new_event = mergeIfSimilar(event, [chrA_start, chrA_end, chrB_start, chrB_end,event_type])
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
                    allVariations[chrA][chrB].append([ startA ,chrA_end + tollerance,  startB,  chrB_end + tollerance,event_type,  1])
            else:
                allVariations[chrA][chrB].append([ chrA_start ,chrA_end, chrB_start,  chrB_end,event_type,  1])
        else:
            startA = chrA_start - tollerance
            if startA < 0:
                startA = 0
            startB = chrB_start - tollerance
            if startB < 0:
                startB = 0
            allVariations[chrA][chrB] =  [[ startA ,chrA_end + tollerance,  startB,  chrB_end + tollerance,event_type,  1]]
        
    else:
        allVariations[chrA]       = {}
        startA = chrA_start - tollerance
        if startA < 0:
            startA = 0
        startB = chrB_start - tollerance
        if startB < 0:
            startB = 0
        allVariations[chrA][chrB] =  [[ startA ,chrA_end + tollerance,  startB,  chrB_end + tollerance,event_type,  1]]

    return allVariations
    
    
def mergeIfSimilar(event, variation): #event is in the DB, variation is the new variation I want to insert
    event_chrA_start    = int(event[0])
    event_chrA_end      = int(event[1])
    event_chrB_start    = int(event[2])
    event_chrB_end      = int(event[3])
    event_type = event[4];

    variation_chrA_start      = int(variation[0])
    variation_chrA_end        = int(variation[1])
    variation_chrB_start      = int(variation[2])
    variation_chrB_end        = int(variation[3])
    variation_type=variation[4];

    new_event = []

    totalEventSize      = 0
    found               = 0
    
    for j in range(2):
        if j == 0:
            variation_start    = variation_chrA_start
            variation_end      = variation_chrA_end
            event_start        = event_chrA_start
            event_end          = event_chrA_end
        else:
            variation_start    = variation_chrB_start
            variation_end      = variation_chrB_end
            event_start        = event_chrB_start
            event_end          = event_chrB_end

        overlap     = 0
        #event         ---------------------
        #variaton   ------------------------------
        #do not merge if the events are of different type
        if(event_type == variation_type):
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
        new_event.append(event[5] +1)
        return new_event
    else:
        return event



if __name__ == '__main__':
    parser = argparse.ArgumentParser("""
    This scripts takes as input any number of vcf files and generates a structural variant database""")
    parser.add_argument('--variations', help="VCF file input", type=str, action='append', nargs='+')
    parser.add_argument('--folder', help="folder containing vcf file, one db file will be generated per prefix", type=str)
    parser.add_argument('--merge_FT', help="merge all inter and intra chromosomal FT db files within teh specified folder", required=False, type=str)
    parser.add_argument('--tollerance', help="expand variations right and left in order to merge similar ones", type=float,  required=False, default=0)
    parser.add_argument('--fixed', help="no expansion or merging of the variants", required=False, action="store_true")
    args = parser.parse_args()
    if(args.variations):
        main(args)
    elif(args.folder):
        vcf_folder = glob.glob(os.path.join(args.folder,"*.vcf"));
        for vcf in vcf_folder:
            #scriptception, a script withing a script calling a script
            command=["python","build_db.py","--variations",vcf]
            if args.fixed:
                command += ["--fixed"]
            if args.tollerance:
                command += ["--tollerance",args.tollerance]
            db=subprocess.check_output(command);
            file_name=vcf.strip(".vcf")+".db"
            f=open(file_name,"w")
            f.write(db)
            f.close();
    elif(args.merge_FT):
        inter_db=glob.glob(os.path.join(args.merge_FT,"*_inter_chr_events.db"))
        intra_db=glob.glob(os.path.join(args.merge_FT,"*_intra_chr_events.db"))

        sample_dict={}
        for db in inter_db:
            print(db)
            sample_name=db.replace("_inter_chr_events.db","")
            sample_dict[sample_name]=[db]
            print(sample_name)
        for db in intra_db:
            sample_name=db.replace("_intra_chr_events.db","")
            sample_dict[sample_name].append(db)

        for sample in sample_dict:
            f=open(sample+".db","w")
            for db in sample_dict[sample]:
                 with open(db) as fin:
                    for line in fin:
                        f.write(line)
                 os.remove(db)
            f.close()             
              
    else:
        print("Error: select a vcf using the variations command or select a folder containing vcf using the folder command")




