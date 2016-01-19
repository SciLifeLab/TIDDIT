import sys,re

def readVCFLine(source,line):
    variation   = line.rstrip().split("\t")
    event_type=""
    chrA=variation[0];
    startA=int(variation[1]);
    startB=0;

    description ={}
    INFO=variation[7].split(";");
    for tag in INFO:
        tag=tag.split("=")
        if(len(tag) > 1):
            description[tag[0]]=tag[1];
    #Delly translocations
    if("TRA" in variation[4]):
        endA=startA

        event_type="BND"
        chrB=description["CHR2"]
        startB=int(description["END"]);
        endB=int(description["END"]);       

    #intrachromosomal variant
    elif(not  "]" in variation[4] and not "[" in variation[4]):
        chrB=chrA;
                
        startB=startA;
        endA=int(description["END"]);
        endB=endA
        #sometimes the fermikit intra chromosomal events are inverted i.e the end pos is a lower position than the start pos
        if(endB <startB):
            tmp=endB;
            endB=startB;
            startB=tmp;

            tmp = endA;
            endA=startA;
            startA=tmp;

        event_type=variation[4].strip("<").rstrip(">");
    #if the variant is given as a breakpoint, it is stored as a precise variant in the db
    else:
        B=variation[4];
        endA=startA;

        B=re.split("[],[]",B);
        for string in B:
            if string.count(":"):
                lst=string.split(":");
                chrB=lst[0]
                startB=int(lst[1]);
                endB=startB

        event_type="BND"
        if(chrA == chrB):
            endA=endB
            startB=startA;
                
    return( chrA.replace("chr","").replace("CHR",""), startA,endA , chrB.replace("chr","").replace("CHR",""), startB, endB, event_type);
