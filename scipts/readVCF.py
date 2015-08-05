import sys

def readVCFLine(source,line):
        if(source == "CNVnator"):
                variation   = line.rstrip().split("\t")
                chrA=variation[0];
                chrB=chrA;

                startA=int(variation[1]);
                startB=startA;

                endA=int(variation[7].split(";")[0].split("=")[1]);
                endB=endA;

        elif(source == "FindTranslocations"):

                variation   = line.rstrip().split("\t")[7]
                description = dict(item.split("=") for item in variation.split(";"))
                #now I need to collapse similar variations
                chrA    = description["CHRA"]
                startA  = int(description["WINA"].split(",")[0])
                endA    = int(description["WINA"].split(",")[1])
                chrB    = description["CHRB"]
                startB  = int(description["WINB"].split(",")[0])
                endB    = int(description["WINB"].split(",")[1])
                
        else:
                print("error, non supported vcf source");
                sys.exit();
        return(chrA, startA,endA,chrB, startB, endB);
