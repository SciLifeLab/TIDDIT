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

    if(not  "]" in variation[4] and not "[" in variation[4]):
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

    else:
        B=variation[4];
        #Create a small interval around the break points
        endA=startA+1000;
        if(startA-1000 > 0):
            startA=startA-1000;
        else:
            startA=1;

        B=re.split("[],[]",B);
        for string in B:
            if string.count(":"):
                lst=string.split(":");
                chrB=lst[0]
                startB=int(lst[1]);
                endB=startB+1000;
            if(startB -1000 > 0):
                startB=startB-1000;
            else:
                startB=1

        event_type="BND"
                
    return( chrA.replace("chr","").strip("CHR",""), startA,endA , chrB.replace("chr","").strip("CHR",""), startB, endB, event_type);
