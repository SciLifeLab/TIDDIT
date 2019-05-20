import itertools
import sqlite3
import math
import numpy

#The functions of this script are used to read the input signals, and to create the database of signals

#extract the vcf header of the signals tab file
def find_contigs(header):
	chromosomes=[]
	library_stats={}
	chromosome_len={}
	for line in header.split("\n"):
		if "##contig=<ID=" in line:
			chromosomes.append(line.split("##contig=<ID=")[-1].split(",")[0])
			chromosome_len[chromosomes[-1]]=int(line.split("length=")[-1].split(">")[0])
		if "##LibraryStats=" in line:
			stats=line.strip().split(" ")
			library_stats={"Coverage":float(stats[1].split("=")[-1]),"ReadLength":int(stats[2].split("=")[-1]),"MeanInsertSize":int(stats[3].split("=")[-1]),"STDInsertSize":int(stats[4].split("=")[-1]),"Orientation":stats[5].split("=")[-1]}
	return (chromosomes,library_stats,chromosome_len)

def sample(args,coverage_data,span_data):
	n=2000
	bridges=[]
	bridges_reads=[]
	bin_size=50

	for chromosome in coverage_data:

		if bin_size*len(coverage_data[chromosome][:,0]) < 200000:
			continue
		index=numpy.random.randint(low=1000, high=len(coverage_data[chromosome][:,0])-1000,size=n)
		for i in range(0,n):
			
			coverage=coverage_data[chromosome][index[i],0]
			if coverage == 0:
				continue

			bridging=span_data[chromosome][index[i],0]
			bridging_reads=span_data[chromosome][index[i],1]
			bridges_reads.append(bridging_reads/coverage)
			bridges.append(bridging/coverage)

	percentile_list=numpy.arange(0,101,1)
	
	if len(bridges):
		percentiles_disc=numpy.percentile(bridges,percentile_list)
		percentiles_splits=numpy.percentile(bridges_reads,percentile_list)
	else:
		print ("Warning - to few reads in the bam file, skipping permutation tests")
		percentiles_disc=[]
		percentiles_splits=[]

	#for i in range(0,len(percentiles_disc)):
	#	print "{} {} {}".format(i,percentiles_disc[i], percentiles_splits[i])

	return (percentiles_disc,percentiles_splits)

#Read the signals tab file
def signals(args,coverage_data):
	signal_data=[]
	header=""
	first_signal=True
	read_list={}
	n_reads=0
	for line in open(args.o+".signals.tab"):
		if line[0] == "#":
			header += line
			continue

		if first_signal:
			chromosomes,library_stats,chromosome_len=find_contigs(header)
			first_signal=False

		content=line.strip().split()

		read_name=content[0]
		if read_name in read_list:
			name=read_list[read_name]
		else:
			name=n_reads
			read_list[read_name]=n_reads
			n_reads+=1

		chrA=content[1]
		posA=int(content[2])
		if content[3] == "+":
			orientationA=1
		else:
			orientationA=0

		cigarA=content[4]
		if ("S" in cigarA or "H" in cigarA) and not cigarA == "NA":
			deletions,insertions,length,clip_after,aligned_range=read_cigar(cigarA)
			if clip_after:
				posA+=length
		else:
			if "M" in cigarA:
				deletions,insertions,length,clip_after,aligned_range=read_cigar(cigarA)
			else:
				length=library_stats["ReadLength"]

			if library_stats["Orientation"] == "innie":			
				if orientationA:
					posA+=length-1
			else:
				if not orientationA:
					posA+=length-1

		if posA >= chromosome_len[chrA]:
			posA=chromosome_len[chrA]-1

		qualA=int(content[5])

		chrB=content[6]
		posB=int(content[7])
		if content[8] == "+":
			orientationB=1
		else:
			orientationB=0
		cigarB=content[9]
		qualB=int(content[10])

		if ("S" in cigarB or "H" in cigarB) and not cigarB == "NA":
			deletions,insertions,length,clip_after,aligned_range=read_cigar(cigarB)
			if clip_after:
				posB+=length
		else:
			if "M" in cigarB:
				deletions,insertions,length,clip_after,aligned_range=read_cigar(cigarB)
			else:
				length=library_stats["ReadLength"]

			if library_stats["Orientation"] == "innie":			
				if orientationB:
					posB+=length-1
			else:
				if not orientationB:
					posB+=length-1

		if posB >= chromosome_len[chrB]:
			posB=chromosome_len[chrB]-1

		resolution=int(content[-1])
		

		if chrA > chrB or (posB < posA and chrA == chrB):
			signal_data.append([chrB,chrA,posB,posA,orientationB,orientationA,qualB,qualA,cigarB,cigarA,resolution,name])
		else:
			signal_data.append([chrA,chrB,posA,posB,orientationA,orientationB,qualA,qualB,cigarA,cigarB,resolution,name])
		

		if len (signal_data) > 1000000:
			args.c.executemany('INSERT INTO TIDDITcall VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',signal_data)  
			signal_data=[]

		idx_b=int(math.floor(posB/50.0))
		if idx_b >= len(coverage_data[chrB]):
			idx_b =len(coverage_data[chrB])-1
		coverage_data[chrB][idx_b,2]+=1
                
		idx_a=int(math.floor(posA/50.0))
		if idx_a >= len(coverage_data[chrA]):
			idx_a=len(coverage_data[chrA])-1
		coverage_data[chrA][idx_a,2]+=1


	if len(signal_data):
		args.c.executemany('INSERT INTO TIDDITcall VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',signal_data)  

	#If no signals were found, parse and return the header
	if first_signal:
		chromosomes,library_stats,chromosome_len=find_contigs(header)

	return(header,chromosomes,library_stats)

#analyse the cigar operation of the read
def read_cigar(cigar):
	deletions=0
	insertions=0
	SC = ["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
	length=0
	first=True
	clip_after=True

	aligned_range=[]
	current_pos=1
	for i in range(0,len(SC)/2):
		if first and SC[i*2+1] == "M":
			first = False
		elif first and SC[i*2+1] == "S":
			first = False
			clip_after=False

		if SC[i*2+1] == "M":
			length += int( SC[i*2] )
			bases=list(range(0,int( SC[i*2] )))
			for j in range(0,len(bases)):
				bases[j] += current_pos

			aligned_range += bases
			current_pos += int( SC[i*2] )
		elif SC[i*2+1] == "I":
			insertions+=1
			bases=list(range(0,int( SC[i*2] )))
			for j in range(0,len(bases)):
				bases[j] += current_pos
			aligned_range += bases

			current_pos += int( SC[i*2] )
		elif SC[i*2+1] == "D":
			length += int( SC[i*2] )
			deletions +=1
		else:
			current_pos += int( SC[i*2] )

	return deletions,insertions,length,clip_after,aligned_range
