import numpy
import math
import copy
import DBSCAN
import gzip
import sys
import sqlite3
import itertools

from scipy.stats import norm

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
			bases=range(0,int( SC[i*2] ))
			for j in range(0,len(bases)):
				bases[j] += current_pos

			aligned_range += bases
			current_pos += int( SC[i*2] )
		elif SC[i*2+1] == "I":
			insertions+=1
			length += int( SC[i*2] )
			bases=range(0,int( SC[i*2] ))
			for j in range(0,len(bases)):
				bases[j] += current_pos
			aligned_range += bases

			current_pos += int( SC[i*2] )
		elif SC[i*2+1] == "D":
			deletions +=1
		else:
			current_pos += int( SC[i*2] )

	return deletions,insertions,length,clip_after,aligned_range
def coverage(args):
	coverage_data={}
	for line in open(args.o+".tab"):
		if line[0] == "#":
			continue
		content=line.strip().split()
		if not content[0] in coverage_data:
			coverage_data[content[0]]=[]
		coverage_data[content[0]].append([ float(content[3]),float(content[4]),0 ])

	for chromosome in coverage_data:
		coverage_data[chromosome]=numpy.array(coverage_data[chromosome])

	return(coverage_data)

def signals(args,coverage_data):
	signal_data=[]
	header=""
	for line in open(args.o+".signals.tab"):
		if line[0] == "#":
			header += line
			continue
		content=line.strip().split()

		chrA=content[0]
		posA=int(content[1])
		if content[2] == "+":
			orientationA=1
		else:
			orientationA=0
		cigarA=content[3]
		if ("S" in cigarA or "H" in cigarA) and not cigarA == "NA":
			deletions,insertions,length,clip_after,aligned_range=read_cigar(cigarA)
			if clip_after:
				posA+=length

		qualA=int(content[4])

		chrB=content[5]
		posB=int(content[6])
		if content[7] == "+":
			orientationB=1
		else:
			orientationB=0
		cigarB=content[8]
		qualB=int(content[9])

		if ("S" in cigarB or "H" in cigarB) and not cigarB == "NA":
			deletions,insertions,length,clip_after,aligned_range=read_cigar(cigarB)
			if clip_after:
				posB+=length

		resolution=int(content[-1])
		

		if chrA > chrB or (posB < posA and chrA == chrB):
			#signal_data[chrB][chrA].append([posB,posA,orientationB,qualB,orientationA,qualA,resolution])
			signal_data.append([chrB,chrA,posB,posA,orientationB,orientationA,qualB,qualA,cigarB,cigarA,resolution])
		else:
			#signal_data[chrA][chrB].append([posA,posB,orientationA,qualA,orientationB,qualB,resolution])
			signal_data.append([chrA,chrB,posA,posB,orientationA,orientationB,qualA,qualB,cigarA,cigarB,resolution])
		

		if len (signal_data) > 1000000:
			args.c.executemany('INSERT INTO TIDDITcall VALUES (?,?,?,?,?,?,?,?,?,?,?)',signal_data)  
			signal_data=[]
		coverage_data[chrB][int(math.floor(posB/100)),2]+=1
		coverage_data[chrA][int(math.floor(posA/100)),2]+=1

	if len(signal_data):
		args.c.executemany('INSERT INTO TIDDITcall VALUES (?,?,?,?,?,?,?,?,?,?,?)',signal_data)  
	return(header)

def find_contigs(header):
	chromosomes=[]
	library_stats={}
	for line in header.split("\n"):
		if "##contig=<ID=" in line:
			chromosomes.append(line.split("##contig=<ID=")[-1].split(",")[0])
		if "##LibraryStats=" in line:
			stats=line.strip().split(" ")
			library_stats={"Coverage":float(stats[1].split("=")[-1]),"ReadLength":int(stats[2].split("=")[-1]),"MeanInsertSize":int(stats[3].split("=")[-1]),"STDInsertSize":int(stats[4].split("=")[-1]),"Orientation":stats[5].split("=")[-1]}
	return (chromosomes,library_stats)

def analyse_pos(candidate_signals,discordants,library_stats,args):
	analysed_signals={}

	max_a=max(candidate_signals[:,0:1])[0]
	min_a=min(candidate_signals[:,0:1])[0]

	max_b=max(candidate_signals[:,1:2])[0]
	min_b=min(candidate_signals[:,1:2])[0]

	
	FF=0
	FR=0
	RF=0
	RR=0

	splitsINV=0
	splits=0

	discs=0

	avgQb=[]
	avgQa=[]

	for i in range(0,len(candidate_signals)):
		if candidate_signals[i][-1] == 1:
			if(candidate_signals[i][2] != candidate_signals[i][4]):
				splitsINV +=1
		else:
			if(candidate_signals[i][2] and not candidate_signals[i][4]):
				FR+=1
			elif( not candidate_signals[i][2] and candidate_signals[i][4]):
				RF+=1
			elif( not candidate_signals[i][2] and not candidate_signals[i][4]):
				RR+=1
			else:
				FF+=1

		if candidate_signals[i][-1] == 1:
			splits +=1
		else:
			discs+=1

		if candidate_signals[i][3] != -1:
			avgQa.append(candidate_signals[i][3])
		if candidate_signals[i][5] != -1:
			avgQb.append(candidate_signals[i][5])

	if avgQa:
		avgQa=numpy.mean(avgQa)
	else:
		avgQa=-1

	if avgQb:
		avgQb=numpy.mean(avgQb)
	else:
		avgQb=-1

	major_orientation=max([FF,FR,RF,RR])
	if not discs:
				posA=max_a
				posB=min_b		
	else:
		if library_stats["Orientation"] == "innie":
			if major_orientation == FF:
				posA=max_a
				posB=max_b
			elif major_orientation == FR:
				posA=max_a
				posB=min_b
			elif major_orientation == RF:
				posA=min_a
				posB=max_b
			else:
				posA=min_a
				posB=min_b
		else:
			if major_orientation == FF:
				posA=min_a
				posB=min_b
			elif major_orientation == FR:
				posA=min_a
				posB=max_b
			elif major_orientation == RF:
				posA=max_a
				posB=min_b
			else:
				posA=max_a
				posB=max_b
	analysed_signals={"posA":posA,"posB":posB,"min_A":min_a,"max_A":max_a,"min_B":min_b,"max_B":max_b,"QA":avgQa,"QB":avgQb,"FF":FF,"FR":FR,"RF":RF,"RR":RR,"splitsINV":splitsINV,"splits":splits,"discs":discs,"signals":candidate_signals}
	return(analysed_signals)


def generate_clusters(chrA,chrB,coordinates,library_stats,args):
	candidates=[]
	coordinates=coordinates[numpy.lexsort((coordinates[:,1],coordinates[:,0]))]
	db=DBSCAN.main(coordinates[:,0:2],args.e,args.l)
	unique_labels = set(db)

	for var in unique_labels:
		if var == -1:
			continue
		class_member_mask = (db == var)
		candidate_signals =coordinates[class_member_mask]
		resolution=candidate_signals[:,-1]
		discordants=True
		if len(set(resolution)) == 1 and max(resolution) == 1:
			disordants=False
		if discordants and len(candidate_signals) >= args.p:
			candidates.append( analyse_pos(candidate_signals,discordants,library_stats,args) )
		elif not discordants and len(candidate_signals) >= args.r and chrA == chrB:
			candidates.append( analyse_pos(candidate_signals,discordants,library_stats,args) )

	return(candidates)

def retrieve_coverage(chromosome,start,end,coverage_data):
	start_index=int(math.floor(start/100.0))
	end_index=int(math.floor(end/100.0))+1
	coverage=numpy.average(coverage_data[chromosome][start_index:end_index,0])
	max_coverage=max(coverage_data[chromosome][start_index:end_index,0])
	quality=numpy.average(coverage_data[chromosome][start_index:end_index,1])
	return (coverage,max_coverage,int(round(quality)))


def retrieve_discs(chromosome,start,end,coverage_data):
	discs=0
	start_index=int(math.floor(start/100.0))
	end_index=int(math.floor(end/100.0))+1
	discs=sum(coverage_data[chromosome][start_index:end_index,2])
	return(discs)



def Part(a, b, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength ):

	readfrequency = 2 * readLength / coverage;
	expr1 = (min([sizeA, sizeB]) - (readLength - 0)) / readfrequency * norm.cdf(a, 0, 1);
	expr2 = -(- 0 ) / readfrequency * norm.cdf(b, 0, 1);
	expr3 = (b * insert_stddev) / readfrequency * (norm.cdf(b, 0, 1) - norm.cdf(a, 0, 1));
	expr4 = (insert_stddev / readfrequency) * (norm.pdf(b, 0, 1) - norm.pdf(a, 0, 1));
	value = expr1 + expr2 + expr3 + expr4
	return(value)

def expected_links(coverageA,coverageB,sizeA,sizeB,gap,insert_mean,insert_stddev,readLength):
	coverage=numpy.average([coverageA,coverageB])

	b1 = (sizeA + sizeB + gap - insert_mean) / insert_stddev
	a1 = (max([sizeA, sizeB]) + gap + readLength  - insert_mean) / insert_stddev
	b2 = (min([sizeA, sizeB]) + gap + readLength  - insert_mean) / insert_stddev
	a2 = (gap + 2 * readLength - insert_mean) / insert_stddev

	e = Part(a1, b1, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength ) - Part(a2, b2, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength)

	return(e)


def fetch_variant_type(chrA,chrB,candidate,args,library_stats):
	variant_type="SVTYPE=BND"
	var="N[{}:{}[".format(chrB,candidate["posB"])
	if not library_stats["ploidies"][chrA] and chrA == chrB:
		variant_type="SVTYPE=DUP"
		var="<DUP>"		

	if chrA == chrB and library_stats["ploidies"][chrA]:
		if candidate["discs"]:
			if candidate["FF"] + candidate["RR"] > candidate["RF"] + candidate["FR"]:
				variant_type="SVTYPE=INV"
				var="<INV>"
			elif library_stats["Orientation"] == "innie":
				if candidate["covM"]/library_stats["chr_cov"][chrA] > (args.n+0.5)/args.n:
					variant_type="SVTYPE=DUP"
					var="<DUP>"
					if candidate["RF"] > candidate["FR"]: 
						variant_type="SVTYPE=TDUP"
						var="<TDUP>"
				elif candidate["covM"]/library_stats["chr_cov"][chrA] < (args.n-0.5)/args.n:
					variant_type="SVTYPE=DEL"
					var="<DEL>"	

			else:
				if candidate["covM"]/library_stats["chr_cov"][chrA] > (args.n+0.5)/args.n:
					variant_type="SVTYPE=DUP"
					var="<DUP>"
					if candidate["RF"] < candidate["FR"]: 
						variant_type="SVTYPE=TDUP"
						var="<TDUP>"

				elif candidate["covM"]/library_stats["chr_cov"][chrA] < (args.n-0.5)/args.n:
					variant_type="SVTYPE=DEL"
					var="<DEL>"
		else:
			if candidate["splitsINV"] > candidate["splits"]-candidate["splitsINV"]:
				variant_type="SVTYPE=INV"
				var="<INV>"
			elif candidate["covM"]/library_stats["chr_cov"][chrA] > (args.n+0.5)/args.n:
					variant_type="SVTYPE=DUP"
					var="<DUP>"
			elif candidate["covM"]/library_stats["chr_cov"][chrA] < (args.n-0.5)/args.n:
					variant_type="SVTYPE=DEL"
					var="<DEL>"		

	if candidate["discs"]:
		if candidate["e2"]*1.6 <= (candidate["discs"]):
			GT="1/1"
		else:
			GT="0/1"
	else:
		if candidate["e2"]*1.6 <= (candidate["splits"]):
			GT="1/1"
		else:
			GT="0/1"

	if "DUP" in var or var == "<DEL>":
		if var == "DEL" and candidate["covM"]/library_stats["chr_cov"][chrA] < 0.1:
			GT="1/1"
		elif var == "DUP" and candidate["covM"]/library_stats["chr_cov"][chrA] > 1.8: 
			GT="1/1"
 
	return(var,variant_type,GT)

def fetch_filter(chrA,chrB,candidate,args,library_stats):
	filt="PASS"
	#fewer links than expected
	if chrA == chrB and (abs(candidate["max_A"]-candidate["min_B"]) < args.z):
		filt="MinSize"
	if candidate["discs"]:
		if candidate["e1"]*0.4 >= candidate["discs"]:
			filt = "BelowExpectedLinks"
	else:
		if candidate["e1"]*0.4 >= candidate["splits"]:
			filt = "BelowExpectedLinks"
	if library_stats["ploidies"][chrA] == 0:
		return("Ploidy")
	if candidate["MaxcovA"] >= library_stats["chr_cov"][chrA]*(library_stats["ploidies"][chrA]+2) or candidate["MaxcovB"] >= library_stats["chr_cov"][chrA]*(library_stats["ploidies"][chrA]+2):
		filt = "UnexpectedCoverage"
	elif candidate["discsA"] > (candidate["discs"]+candidate["splits"])*(1+library_stats["ploidies"][chrA]) or candidate["discsB"] > (candidate["discs"]+candidate["splits"])*(1+library_stats["ploidies"][chrA]):
		filt= "FewLinks"
	elif chrA == chrB and candidate["max_A"] > candidate["min_B"]:
		filt = "Smear"
	elif chrA != chrB or (abs(candidate["posA"]-candidate["posB"]) > library_stats["MeanInsertSize"]+3*library_stats["STDInsertSize"] ):
		if not candidate["discs"] or candidate["discs"]*4 <  candidate["splits"] or candidate["discs"] <= args.p/2:
			filt= "SplitsVSDiscs"
	elif candidate["QRA"] < args.Q or candidate["QRB"] < args.Q:
		filt = "RegionalQ"

	return(filt)

def redefine_inv(vcf_line,signals,library_stats,args):
	analysed_signals=analyse_pos(signals,True,library_stats,args)
	vcf_line[1]=analysed_signals["posA"]
	vcf_line[7]=vcf_line[7].split(";END=")[0]+";END="+str(analysed_signals["posB"]) +";SVLEN=" + vcf_line[7].split(";SVLEN=")[-1]
	vcf_line[7]=vcf_line[7].split("SVLEN=")[0]+"SVLEN="+str(analysed_signals["posB"]-analysed_signals["posA"]+1) +";COVM=" + vcf_line[7].split(";COVM=")[-1]

	vcf_line[7]=vcf_line[7].split("CIPOS=")[0]+"CIPOS={},{}".format(analysed_signals["min_A"]-analysed_signals["posA"],analysed_signals["max_A"]-analysed_signals["posA"]) + ";CIEND" + vcf_line[7].split(";CIEND")[-1]
	vcf_line[7]=vcf_line[7].split("CIEND=")[0]+"CIEND={},{}".format(analysed_signals["min_B"]-analysed_signals["posB"],analysed_signals["max_B"]-analysed_signals["posB"]) + ";END=" + vcf_line[7].split(";END=")[-1]
	return(vcf_line)

def inv_recluster(vcf_line,candidate,library_stats,args):
	RR=[]
	FF=[]

	for signal in candidate["signals"]:
		if( not signal[2] and not signal[4] ):
			RR.append(signal)
		else:
			FF.append(signal)

	vcf_line_rr=copy.copy(vcf_line)
	vcf_line_ff=copy.copy(vcf_line)
	vcf_line_rr=redefine_inv(vcf_line_rr,numpy.array(RR),library_stats,args)
	vcf_line_ff=redefine_inv(vcf_line_ff,numpy.array(FF),library_stats,args)
	var_id=list(vcf_line_ff[2])
	var_id[-1]="2"
	vcf_line_ff[2]="".join(var_id)

	return([vcf_line_rr,vcf_line_ff])

def generate_vcf_line(chrA,chrB,n,candidate,args,library_stats):
	vcf_line=[]
	if chrA == chrB and  candidate["posA"] > candidate["posB"]:
		candidate["posA"]=candidate["min_A"]
		candidate["posB"]=candidate["max_B"]

	vcf_line.append(chrA)
	vcf_line.append(candidate["posA"])
	vcf_line.append("SV_{}_1".format(n))
	vcf_line.append("N")
	var,variant_type,GT=fetch_variant_type(chrA,chrB,candidate,args,library_stats)
	vcf_line.append(var)
	vcf_line.append(".")
	vcf_line.append(fetch_filter(chrA,chrB,candidate,args,library_stats))
	INFO="{};CIPOS={},{};CIEND={},{}".format(variant_type,candidate["min_A"]-candidate["posA"],candidate["max_A"]-candidate["posA"],candidate["min_B"]-candidate["posB"],candidate["max_B"]-candidate["posB"])
	if chrA == chrB:
		if not "BND" in variant_type:
			INFO += ";END={};SVLEN={}".format(candidate["posB"],abs(candidate["posB"]-candidate["posA"]+1))
		INFO += ";COVM={}".format(candidate["covM"])
	stats=";COVA={};COVB={};LFA={};LFB={};LTE={};E1={};E2={};OR={},{},{},{};ORSR={},{};QUALA={};QUALB={}".format(candidate["covA"],candidate["covB"],int(round(candidate["discsA"])),int(round(candidate["discsB"])),candidate["discs"]+candidate["splits"],candidate["e1"],candidate["e2"],candidate["FF"],candidate["RR"] ,candidate["RF"],candidate["FR"],candidate["splitsINV"],candidate["splits"]-candidate["splitsINV"],candidate["QRA"],candidate["QRB"])

	INFO+=stats

	vcf_line.append(INFO)
	FORMAT_FORMAT="GT:CN:DV:RV"
	vcf_line.append(FORMAT_FORMAT)
	CN="."
	if "DEL" in var or "DUP" in var:
		CN=int(round(candidate["covM"]/(library_stats["Coverage"]/args.n)))
	if "DEL" in var:
		CN=library_stats["ploidies"][chrA]-CN
		if CN < 0:
			CN=library_stats["ploidies"][chrA]

	FORMAT_STR="{}:{}:{}:{}".format(GT,CN,candidate["discs"],candidate["splits"])
	vcf_line.append(FORMAT_STR)
	if not "BND" in variant_type and not "INV" in variant_type:
		return([vcf_line])
	elif "INV" in variant_type:
		if candidate["FF"] and candidate["RR"] and candidate["RR"] > args.p/2 and candidate["FF"] > args.p/2 and candidate["discs"]:
			return(inv_recluster(vcf_line,candidate,library_stats,args))
		else:
			return([vcf_line])
	else:
		vcf_line_a=copy.copy(vcf_line)
		vcf_line_b=copy.copy(vcf_line)
		inverted=False
		before=True
		if candidate["posA"] == candidate["max_A"]:
			before=False

		if candidate["FF"] + candidate["RR"]  > candidate["RF"] + candidate["FR"]:
			inverted=True

		if not inverted and not before:
			alt_str_a="N[{}:{}[".format(chrB,candidate["posB"])
			alt_str_b="]{}:{}]N".format(chrA,candidate["posA"])
		elif not inverted and before:
			alt_str_a="]{}:{}]N".format(chrB,candidate["posB"])
			alt_str_b="N[{}:{}[".format(chrA,candidate["posA"])
		elif inverted and  not before:
			alt_str_a="N]{}:{}]".format(chrB,candidate["posB"])
			alt_str_b="[{}:{}[N".format(chrA,candidate["posA"])
		else:
			alt_str_a="[{}:{}[N".format(chrB,candidate["posB"])
			alt_str_b="N]{}:{}]".format(chrA,candidate["posA"])

		vcf_line_a[4]=alt_str_a
		vcf_line_b[4]=alt_str_b
		vcf_line_b[2]="SV_{}_2".format(n)
		vcf_line_b[0]=chrB
		vcf_line_b[1]=candidate["posB"]
		return([vcf_line_a,vcf_line_b])

def determine_ploidy(args,chromosomes,coverage_data,Ncontent,sequence_length,library_stats):

	library_stats["chr_cov"]={}

	normalising_chromosomes={}
	ploidies={}
	avg_coverage=[]
	coverage_norm=0
	if not args.force_ploidy:
		try:
			for chromosome in args.s.split(","):
				ploidies[chromosome]=args.n
				N_count=Ncontent[chromosome]
				chromosomal_average=numpy.median(coverage_data[chromosome][numpy.where(N_count > 0),0])
				avg_coverage.append( chromosomal_average )
				library_stats["chr_cov"][chromosome]=chromosomal_average
		except:
			print "error: reference mismatch!"
			print "make sure that the contigs of the bam file and the reference match"
			print "also make sure that the chromosomes suplied through the -s parameter match the reference"
			print "you may use --force_ploidy to skip the per chromosome ploidy estimation"
			quit()

		coverage_norm=numpy.median(avg_coverage)

	chromosomal_average=0
	outfile=open(args.o+".ploidy.tab", 'w')
	outfile.write("Contig\tploidy_rounded\tploidy_raw\tmedian_coverage\n")
	for chromosome in chromosomes:
		if not chromosome in ploidies:
   
			chromosome_length=sequence_length[chromosome]
			N_count=Ncontent[chromosome]
			chromosomal_average=numpy.median(coverage_data[chromosome][numpy.where(N_count > 0),0])
			if not args.force_ploidy:
				try:
					ploidies[chromosome]=int(round((chromosomal_average)/coverage_norm*args.n))
				except:
					ploidies[chromosome]=args.n
			else:
				ploidies[chromosome]=args.n  
			library_stats["chr_cov"][chromosome]=chromosomal_average
		
		outfile.write("{}\t{}\t{}\t{}\n".format(chromosome,ploidies[chromosome],round( library_stats["chr_cov"][chromosome]/coverage_norm*args.n,2),library_stats["chr_cov"][chromosome]))

	outfile.close()
	return(ploidies,library_stats)

def retrieve_N_content(args):
	if not args.ref.endswith(".gz"):
		with open(args.ref, 'r+') as f:
			sequence=f.read().split(">")
	else:
		with gzip.open(args.ref, 'r+') as f:
			sequence=f.read().split(">")

	del sequence[0]
	Ncontent={}
	sequence_length={}
	for chromosome in sequence:
		content=chromosome.split("\n",1)
		sequence=content[1].replace("\n","")
		contig=content[0].split()[0]
		regions=[sequence[i:i+100] for i in range(0, len(sequence), 100)]
		Ncontent[contig]=[]
		for region in regions:
			if region.upper().count("N")/100.0 > args.n_mask:
				#print region.upper().count("N")/100.0
				Ncontent[contig].append(0)
			else:  
				Ncontent[contig].append(1)
		sequence_length[contig]=len(sequence)		
		Ncontent[contig]=numpy.array(Ncontent[contig])
		

	return(Ncontent,sequence_length)

def cluster(args):

	Ncontent,sequence_length=retrieve_N_content(args)
	coverage_data=coverage(args)

	conn = sqlite3.connect(args.o+".db")
	args.c = conn.cursor()

	tableListQuery = "SELECT name FROM sqlite_master WHERE type=\'table\'"
	args.c.execute(tableListQuery)
	tables = map(lambda t: t[0], args.c.fetchall())
	if "TIDDITcall" in tables:
		args.c.execute("DROP TABLE TIDDITcall")
	A="CREATE TABLE TIDDITcall (chrA TEXT,chrB TEXT,posA INT,posB INT,forwardA INT,forwardB INT,qualA INT, qualB INT,cigarA TEXT,cigarB TEXT, resolution INT)"
	args.c.execute(A)
	header=signals(args,coverage_data)
	A="CREATE INDEX CHR ON TIDDITcall (chrA, chrB)"
	args.c.execute(A)

	chromosomes,library_stats=find_contigs(header)
	ploidies,library_stats=determine_ploidy(args,chromosomes,coverage_data,Ncontent,sequence_length,library_stats)
	library_stats["ploidies"]=ploidies

	if not args.e:
		args.e=int(math.sqrt(library_stats["MeanInsertSize"]*2)*12)
	n=1
	print "clustering signals on chromosome:"
	calls={}
	for chrA in chromosomes:
		calls[chrA] =[]
		print "{}".format(chrA)
		for chrB in chromosomes:
			signal_data=numpy.array([ [hit[0],hit[1],hit[2],hit[3],hit[4],hit[5],hit[6]] for hit in args.c.execute('SELECT posA,posB,forwardA,qualA,forwardB,qualB,resolution FROM TIDDITcall WHERE chrA == \'{}\' AND chrB == \'{}\''.format(chrA,chrB)).fetchall()])

			if not len(signal_data):
				continue

			candidates=generate_clusters(chrA,chrB,signal_data,library_stats,args)

			for i in range(0,len(candidates)):
				candidates[i]["covA"],candidates[i]["MaxcovA"],candidates[i]["QRA"]=retrieve_coverage(chrA,candidates[i]["min_A"],candidates[i]["max_A"],coverage_data)
				candidates[i]["covB"],candidates[i]["MaxcovB"],candidates[i]["QRB"]=retrieve_coverage(chrB,candidates[i]["min_B"],candidates[i]["max_B"],coverage_data)
				if chrA == chrB:
					if candidates[i]["posB"] > candidates[i]["posA"]:
						candidates[i]["covM"],candidates[i]["MaxcovM"],candidates[i]["QRM"]=retrieve_coverage(chrB,candidates[i]["posA"],candidates[i]["posB"],coverage_data)
					else:
						candidates[i]["covM"],candidates[i]["MaxcovM"],candidates[i]["QRM"]=retrieve_coverage(chrB,candidates[i]["posB"],candidates[i]["posA"],coverage_data)
				else:
					candidates[i]["covM"]=0
					candidates[i]["QRM"]=0
				candidates[i]["discsA"]=retrieve_discs(chrA,candidates[i]["min_A"],candidates[i]["max_A"],coverage_data)
				candidates[i]["discsB"]=retrieve_discs(chrB,candidates[i]["min_B"],candidates[i]["max_B"],coverage_data)
				sizeA=candidates[i]["max_A"]-candidates[i]["min_A"]+library_stats["ReadLength"]
				sizeB=candidates[i]["max_B"]-candidates[i]["min_B"]+library_stats["ReadLength"]
				gap=[]
				for signal in candidates[i]["signals"]:
					gap.append(abs(sizeA-(signal[0]-candidates[i]["min_A"])+(signal[1]-candidates[i]["min_B"])))
				gap=numpy.average(gap)
				if gap > sizeA:
					gap = gap - sizeA
				else:
					gap=1
				if ploidies[chrA]:
					coverageA=candidates[i]["covA"]/ploidies[chrA]
				else:
					coverageA=candidates[i]["covA"]
				if ploidies[chrB]:
					coverageB=candidates[i]["covB"]/ploidies[chrB]
				else:
					coverageB=candidates[i]["covB"]

				if candidates[i]["discs"]:
					candidates[i]["e1"]=int(round(expected_links(coverageA,coverageB,sizeA,sizeB,gap,library_stats["MeanInsertSize"],library_stats["STDInsertSize"],library_stats["ReadLength"])))
					candidates[i]["e2"]=int(round(expected_links(coverageA,coverageB,library_stats["MeanInsertSize"],library_stats["MeanInsertSize"],0,library_stats["MeanInsertSize"],library_stats["STDInsertSize"],library_stats["ReadLength"])))
				else:
					candidates[i]["e1"]=min([coverageA,coverageB])
					candidates[i]["e2"]=min([coverageA,coverageB])
				vcf_line=generate_vcf_line(chrA,chrB,n,candidates[i],args,library_stats)

				if len(vcf_line) == 1:
					calls[chrA].append(vcf_line[0])
				else:
					calls[chrA].append(vcf_line[0])
					if not chrB in calls:
						calls[chrB]=[]
					calls[chrB].append(vcf_line[1])
				n+=1

	outfile=open(args.o+".vcf", 'w')
	new_header=[]
	for line in header.strip().split("\n"):
		if "##TIDDITcmd" in line:
			new_header.append("##TIDDITcmd=\"{}\"\n".format(" ".join(sys.argv)))
			continue
		new_header.append(line+"\n")
	header="".join(new_header)
	outfile.write(header)
	for chromosome in chromosomes:
		for call in sorted(calls[chromosome],key=lambda x: x[1]):
			output=call
			output[1]=str(output[1])
			output="\t".join(output)+"\n"
			outfile.write(output)

	return()

