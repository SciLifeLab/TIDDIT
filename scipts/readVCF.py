import sys,re

def readVCFLine(source,line):

        variation   = line.rstrip().split("\t")
        if(source == "CNVnator"):
                #print(variation)
                chrA=variation[0];
                chrB=chrA;

                startA=int(variation[1]);
                startB=startA;

                endA=int(variation[7].split(";")[0].split("=")[1]);
                endB=endA;

                event_type=variation[4].strip("<").rstrip(">");

        elif(source == "FindTranslocations"):

                description = dict(item.split("=") for item in variation[7].split(";"))
                if(not "]" in variation[4] and not "[" in variation[4]):
                    #if the event is an intrachromosomal event
                    chrA    = description["CHRA"]
                    chrB=chrA;
                    startA=int(variation[1]);
                    startB=startA;

                    endA=int(description["END"]);
                    endB=endA        

                    event_type=description["SVTYPE"];

                else:

                    
                    #now I need to collapse similar variations
                    chrA    = description["CHRA"]
                    startA  = int(description["WINA"].split(",")[0])
                    endA    = int(description["WINA"].split(",")[1])
                    chrB    = description["CHRB"]
                    startB  = int(description["WINB"].split(",")[0])
                    endB    = int(description["WINB"].split(",")[1])
                    event_type=description["SVTYPE"]

        #if the source is fermikit
        elif(source == "htsbox-abreak-r303"):

            chrA=variation[0];
            startA=int(variation[1]);

            description = dict(item.split("=") for item in variation[7].split(";"))
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
                event_type=description["SVTYPE"];

            else:
                B=variation[4];
                #fermikit assigns translocations as precise events, thus we need to create a small interval around the positions of the translocations or we will be unlikely to overlap similar fermikit events
                endA=startA+500;
                if(startA-500 > 0):
                    startA=startA-500;
                else:
                    startA=1;

                B=re.split("[],[]",B);
                for string in B:
                    if string.count(":"):
                        lst=string.split(":");
                        chrB=lst[0]
                        startB=int(lst[1]);
                        endB=startB+500;
                        if(startB -500 > 0):
                            startB=startB-500;
                        else:
                            startB=1

                event_type=description["SVTYPE"]
                
                
                
        else:
                print("error, non supported vcf source");
                sys.exit();
        return(chrA, startA,endA,chrB, startB, endB, event_type);
