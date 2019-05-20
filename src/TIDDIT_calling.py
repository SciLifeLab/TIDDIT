import numpy
import math
import copy
import sys
import sqlite3
import os
import time

import DBSCAN
import TIDDIT_coverage
import TIDDIT_filtering
import TIDDIT_signals

def retrieve_coverage(chromosome,start,end,coverage_data):
	start_index=int(math.floor(start/50.0))
	end_index=int(math.floor(end/50.0))+1
	regional_coverage=coverage_data[chromosome][start_index:end_index,0]
	coverage=numpy.average(regional_coverage[numpy.where(regional_coverage > -1)])
	max_coverage=max(regional_coverage[numpy.where(regional_coverage > -1)])
	quality=numpy.average(coverage_data[chromosome][start_index:end_index,1])
	return (coverage,max_coverage,int(round(quality)))


#Find the number of discordant pairs within a region
def retrieve_discs(chromosome,start,end,coverage_data):
	discs=0
	start_index=int(math.floor(start/50.0))
	end_index=int(math.floor(end/50.0))+1
	discs=sum(coverage_data[chromosome][start_index:end_index,2])
	return(discs)

def count_ref(chromosome,pos,span_data):
	index=int(math.floor(pos/50.0))
	dr=span_data[chromosome][index,0]
	rr=span_data[chromosome][index,1]
	return(dr,rr)

#split inversions into two separate calls
def redefine_inv(vcf_line,signals,library_stats,args):
	analysed_signals=DBSCAN.analyse_pos(signals,True,library_stats,args)
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

def generate_vcf_line(chrA,chrB,n,candidate,args,library_stats,percentiles,span_data):
	vcf_line=[]
	if chrA == chrB and  candidate["posA"] > candidate["posB"]:
		candidate["posA"]=candidate["min_A"]
		candidate["posB"]=candidate["max_B"]

	if candidate["discs"] > candidate["splits"]:
		pairs=candidate["discs"]
	else:
		pairs=candidate["splits"]

	percentile_list=numpy.arange(0,101,1)
	p=0
	if candidate["covA"] < candidate["covB"]:
		coverage=candidate["covA"]
		ploidy=library_stats["ploidies"][chrA]
	else:
		coverage=candidate["covB"]
		ploidy=library_stats["ploidies"][chrB]
	i=0
	perc=0
	p=0
	for percentile in percentiles:
		if pairs*ploidy > percentile*coverage:
			p=percentile_list[i]
			perc=percentile
		i+=1
	if p:
		qual=str(p)
	elif not p and not len(percentiles):
		qual="."
	else:
		qual="0"


	a_dr,a_rr=count_ref(chrA,candidate["posA"],span_data)
	b_dr,b_rr=count_ref(chrB,candidate["posB"],span_data)

	split_ratio=0
	if candidate["splits"]:
		split_ratio=candidate["splits"]/float(candidate["splits"]+min([a_rr,b_rr]))
	disc_ratio=0
	if candidate["discs"]:
		disc_ratio=candidate["discs"]/float(candidate["discs"]+min([a_dr,b_dr]))

	candidate["ratio"]=1
	if candidate["splits"] and candidate["splits"] > candidate["discs"]:
		candidate["ratio"]=split_ratio
	elif candidate["discs"]:
		candidate["ratio"]=disc_ratio


	vcf_line.append(chrA)
	vcf_line.append(candidate["posA"])
	vcf_line.append("SV_{}_1".format(n))
	vcf_line.append("N")
	var,variant_type,GT=TIDDIT_filtering.fetch_variant_type(chrA,chrB,candidate,args,library_stats,disc_ratio,split_ratio)
	vcf_line.append(var)
	vcf_line.append(str(qual))
	vcf_line.append(TIDDIT_filtering.fetch_filter(chrA,chrB,candidate,args,library_stats))
	INFO="{};CIPOS={},{};CIEND={},{}".format(variant_type,candidate["min_A"]-candidate["posA"],candidate["max_A"]-candidate["posA"],candidate["min_B"]-candidate["posB"],candidate["max_B"]-candidate["posB"])
	if chrA == chrB:
		if not "BND" in variant_type:
			INFO += ";END={};SVLEN={}".format(candidate["posB"],abs(candidate["posB"]-candidate["posA"]+1))
		INFO += ";COVM={}".format(candidate["covM"])
	stats=";COVA={};COVB={};LFA={};LFB={};LTE={};OR={},{},{},{};ORSR={},{};QUALA={};QUALB={}".format(candidate["covA"],candidate["covB"],int(round(candidate["discsA"])),int(round(candidate["discsB"])),candidate["discs"]+candidate["splits"],candidate["FF"],candidate["RR"] ,candidate["RF"],candidate["FR"],candidate["splitsINV"],candidate["splits"]-candidate["splitsINV"],candidate["QRA"],candidate["QRB"])

	INFO+=stats

	vcf_line.append(INFO)
	FORMAT_FORMAT="GT:CN:DV:RV:DR:RR"
	vcf_line.append(FORMAT_FORMAT)
	CN="."
	if "DEL" in var or "DUP" in var:
		CN=int(round(candidate["covM"]/(library_stats["Coverage"]/args.n)))
	if "DEL" in var:
		CN=library_stats["ploidies"][chrA]-CN
		if CN < 0:
			CN=library_stats["ploidies"][chrA]

	FORMAT_STR="{}:{}:{}:{}:{},{}:{},{}".format(GT,CN,candidate["discs"],candidate["splits"],a_dr,b_dr,a_rr,b_rr)
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

#main function
def cluster(args):
	start_time=time.time()

	if args.ref:
		Ncontent=TIDDIT_coverage.retrieve_N_content(args)
	else:
		Ncontent=[]

	coverage_data,span_data=TIDDIT_coverage.coverage(args)


	conn = sqlite3.connect(args.o+".db")
	args.c = conn.cursor()

	tableListQuery = "SELECT name FROM sqlite_master WHERE type=\'table\'"
	args.c.execute(tableListQuery)
	tables = map(lambda t: t[0], args.c.fetchall())
	if "TIDDITcall" in tables:
		args.c.execute("DROP TABLE TIDDITcall")

	args.c.execute("CREATE TABLE TIDDITcall (chrA TEXT,chrB TEXT,posA INT,posB INT,forwardA INT,forwardB INT,qualA INT, qualB INT,cigarA TEXT,cigarB TEXT, resolution INT,name INT)")
	header,chromosomes,library_stats=TIDDIT_signals.signals(args,coverage_data)
	args.c.execute("CREATE INDEX CHR ON TIDDITcall (chrA, chrB)")

	ploidies,library_stats,coverage_data=TIDDIT_coverage.determine_ploidy(args,chromosomes,coverage_data,Ncontent,library_stats)
	library_stats["ploidies"]=ploidies

	percentiles_disc,percentiles_splits=TIDDIT_signals.sample(args,coverage_data,span_data)

	if not args.e:
		args.e=int(math.sqrt(library_stats["MeanInsertSize"]*2)*12)
	n=1
	print ("clustering signals on chromosome:")
	calls={}
	for chrA in chromosomes:
		calls[chrA] =[]
		print ("{}".format(chrA))
		for chrB in chromosomes:
			signal_data=numpy.array([ [hit[0],hit[1],hit[2],hit[3],hit[4],hit[5],hit[6],hit[7]] for hit in args.c.execute('SELECT posA,posB,forwardA,qualA,forwardB,qualB,resolution,name FROM TIDDITcall WHERE chrA == \'{}\' AND chrB == \'{}\''.format(chrA,chrB)).fetchall()])

			if not len(signal_data):
				continue

			candidates=DBSCAN.generate_clusters(chrA,chrB,signal_data,library_stats,args)

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
				if sizeA > library_stats["MeanInsertSize"] or sizeB > library_stats["MeanInsertSize"]:
					gap=1
				elif gap > sizeA:
					gap = gap - sizeA
				else:
					gap=1
				if ploidies[chrA]:
					coverageA=candidates[i]["covA"]/ploidies[chrA]
				else:
					coverageA=candidates[i]["covA"]/float(args.n)
				if ploidies[chrB]:
					coverageB=candidates[i]["covB"]/ploidies[chrB]
				else:
					coverageB=candidates[i]["covB"]/float(args.n)

				if candidates[i]["discs"] > candidates[i]["splits"]:
					vcf_line=generate_vcf_line(chrA,chrB,n,candidates[i],args,library_stats,percentiles_disc,span_data)
				else:
					vcf_line=generate_vcf_line(chrA,chrB,n,candidates[i],args,library_stats,percentiles_splits,span_data)

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
			if "MinSize" == output[6]:
				continue

			output="\t".join(output)+"\n"
			outfile.write(output)

	print ("variant clustering completed in {}".format(time.time()-start_time))
	print ("Work complete!")
	os.remove("{}.db".format(args.o))
	return()

