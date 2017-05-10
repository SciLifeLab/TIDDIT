import sys

read_lim=int(sys.argv[2])
pair_lim=int(sys.argv[3])

for line in open(sys.argv[1]):
    if line[0] == "#":
        print line.strip()
        continue
    content=line.strip().split("\t")
    format_fields=content[9:]
    print_line = False
    for sample in format_fields:
        reads=int( sample.split(":")[-1] )
        pairs=int( sample.split(":")[-2] )
        if not pairs:
            if reads >= read_lim:
                print_line=True
        else:
            if reads+pairs >= pair_lim:
                print_line=True
 
    if print_line:
        print line.strip()
