import time
import math
import numpy
from joblib import Parallel, delayed

import pysam
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment

def percentile(a, q):
	size = len(a)
	percentiles=[]

	sorted_a=sorted(a)

	for v in q:
		if not size:
			percentiles.append(0)
		else:
			percentiles.append( sorted_a[ int(math.ceil((size * v) / 100.0)) - 1 ] )

	return(percentiles)


def scoring(scoring_dict,percentiles):
	score=[0]
	if scoring_dict["n_contigs"]:
		score.append(50)

	if scoring_dict["n_discordants"]:
		score.append(0)
		for p in percentiles["FA"]:
			if scoring_dict["n_discordants"]/(scoring_dict["refFA"]+scoring_dict["n_discordants"]) >= p:
				score[-1]+=5

		score.append(0)
		for p in percentiles["FB"]:
			if scoring_dict["n_discordants"]/(scoring_dict["refFB"]+scoring_dict["n_discordants"]) >= p:
				score[-1]+=5


	if scoring_dict["n_splits"]:
		score.append(0)
		for p in percentiles["RA"]:
			if scoring_dict["n_splits"]/(scoring_dict["refRA"]+scoring_dict["n_splits"]) >= p:
				score[-1]+=5

		score.append(0)
		for p in percentiles["RB"]:
			if scoring_dict["n_splits"]/(scoring_dict["refRB"]+scoring_dict["n_splits"]) >= p:
				score[-1]+=5

	return(max(score))

def get_region(AlignmentFile samfile,str chr,int start,int end,int bp,int min_q,int max_ins, contig_number):

	cdef int low_q=0
	cdef int n_reads=0
	cdef long bases=0
	cdef int n_discs=0
	cdef int n_splits=0

	cdef int crossing_r=0
	cdef int crossing_f=0

	bam_header=samfile.header
	contig_length=bam_header["SQ"][contig_number[chr]]["LN"]

	q_start=start
	q_end=end+max_ins

	if q_end > contig_length:
		q_end=contig_length

	if q_start >= q_end:
		q_start=q_end-10

	cdef long read_reference_start
	cdef long read_reference_end

	cdef long r_start
	cdef long r_end

	cdef AlignedSegment read

	for read in samfile.fetch(chr, q_start, q_end):
		if read.is_unmapped:
			continue

		read_reference_start=read.reference_start

		if not read.mate_is_unmapped:
			if read.next_reference_start > end and read_reference_start > end:
				continue
		else:
			if read_reference_start > end:
				continue

		if read.is_duplicate:
			continue

		if not (read_reference_start > end):
			n_reads+=1
			if read.mapq < min_q:
				low_q+=1

		if read.mapq < min_q:
			continue

		read_reference_end=read.reference_end
		read_reference_name=read.reference_name
		read_next_reference_name=read.next_reference_name

		r_start=read_reference_start
		r_end=read_reference_end

		if read_reference_start < bp-20 and r_end > bp+20:
			crossing_r+=1

		mate_bp_read= (read.next_reference_start < bp-50 and r_end > bp+50)
		discordant= ( abs(read.isize) > max_ins or read_next_reference_name != read_reference_name )

		if mate_bp_read and not discordant:
			crossing_f+=1

		if read_reference_end < start:
			continue
		elif read_reference_start > end:
			continue

		if read_reference_start < start:
			r_start=start

		if read_reference_end > end:
			r_end=end

		bases+=r_end-r_start+1

		if read.has_tag("SA"):
			n_splits+=1

		if discordant:
			n_discs+=1

	coverage= bases/(end-start+1)

	if n_reads > 0:
		frac_low_q=low_q/float(n_reads)
	else:
		frac_low_q=0

	return(coverage,frac_low_q,n_discs,n_splits,crossing_f,crossing_r)

def find_sv_type(chrA,chrB,inverted,non_inverted,args,sample_data,samples,library):
	if chrA != chrB:
		return("BND",".")

	p=library["contig_ploidy_{}".format(chrA)]
	for sample in samples:
		if library[ "avg_coverage_{}".format(chrA) ] != 0:
			cn=int(round(sample_data[sample]["covM"]*p/library[ "avg_coverage_{}".format(chrA) ]))
		else:
			cn=int(round(sample_data[sample]["covM"]*args.n/library[ "avg_coverage" ]))


	#mitochondria or similar
	if p > args.n*10:
		if cn > p*1.05:
			if inverted:
				return("DUP:INV",cn)
			else:
				return("DUP:TANDEM",cn)
		elif cn < p*0.95:
			return("DEL",cn)
		elif inverted > non_inverted:
			return("INV",cn)
		else:
			return("BND",cn)
	else:
		if cn > p:
			if inverted:
				return("DUP:INV",cn)
			else:
				return("DUP:TANDEM",cn)

		if inverted > non_inverted:
			return("INV",cn)
		elif cn < p:
			return("DEL",cn)
		else:
			return("BND",cn)

def sv_filter(sample_data,args,chrA,chrB,posA,posB,max_ins_len,n_discordants,n_splits,library,n_discs_regionA,n_discs_regionB,n_splits_regionA,n_splits_regionB,n_contigs):
	filt="PASS"
	for sample in sample_data:

		#filter SV in high coverage regions
		if sample_data[sample]["covA"] > args.max_coverage*library[ "avg_coverage_{}".format(chrA) ]:
			return("UnexpectedCoverage")
		if sample_data[sample]["covB"] > args.max_coverage*library[ "avg_coverage_{}".format(chrB) ]:
			return("UnexpectedCoverage")
		if sample_data[sample]["covM"] > args.max_coverage*library[ "avg_coverage_{}".format(chrB) ]:
			return("UnexpectedCoverage")

		#if abs(posA-posB) > max_ins_len:
		if not n_contigs:
			if (n_discordants < args.p_ratio * sample_data[sample]["refFA"] or n_discordants < args.p_ratio * sample_data[sample]["refFB"]) and (n_splits < args.r_ratio * sample_data[sample]["refRA"] or n_splits < args.r_ratio * sample_data[sample]["refRB"]):
				return("BelowExpectedLinks")

			if (n_discordants < args.p_ratio * sample_data[sample]["covA"] or n_discordants < args.p_ratio * sample_data[sample]["covB"]) and (n_splits < args.r_ratio * sample_data[sample]["covA"] or n_splits < args.r_ratio * sample_data[sample]["covB"]):
				return("BelowExpectedLinks")

			if n_discordants > n_splits:
				if n_discordants < 0.25*n_discs_regionA or n_discordants < 0.25*n_discs_regionB:
					return("FewLinks")

			else:
				if n_splits < 0.25*n_splits_regionA or n_splits < 0.25*n_splits_regionB:
					return("FewLinks")

			if sample_data[sample]["QA"] > 0.2 or sample_data[sample]["QB"] > 0.2:
				return("RegionalQ")

			if n_discordants == 0 and (chrA != chrB):
				return("SplitsVSDiscs")

		#interchromsomal translocation, supported only by contigs
		if n_contigs and (chrA != chrB):
			if n_discordants < args.p:
				return("BelowExpectedLinks")

		#large variant, supported only by contigs but not discordant pairs
		elif n_contigs and (chrA == chrB and max_ins_len*3 < abs(posB-posA) ):
			if n_discordants < args.p:
				return("BelowExpectedLinks")

	return(filt)

def define_variant(str chrA, str bam_file_name,dict sv_clusters,args,dict library,int min_mapq,samples,dict coverage_data,contig_number,max_ins_len,contig_seqs,gc):
	cdef AlignmentFile samfile  = AlignmentFile(bam_file_name, "r",reference_filename=args.ref,index_filename="{}_tiddit/{}.csi".format(args.o,samples[0]))
	variants=[]

	var_n=0
	for chrB in sv_clusters[chrA]:

		for cluster in sv_clusters[chrA][chrB]:

			n_discordants=sv_clusters[chrA][chrB][cluster]["N_discordants"]
			n_splits=sv_clusters[chrA][chrB][cluster]["N_splits"]
			n_contigs=sv_clusters[chrA][chrB][cluster]["N_contigs"]

			if (n_discordants < args.p and n_splits < args.r) and not n_contigs:
					continue

			posA=sv_clusters[chrA][chrB][cluster]["posA"]
			posB=sv_clusters[chrA][chrB][cluster]["posB"]

			if chrA == chrB and posA > posB:
				posT=posA
				posA=posB
				posB=posT

			if chrA == chrB and abs(posA-posB) < args.z:
				continue

			s=int(math.floor(sv_clusters[chrA][chrB][cluster]["startA"]/50.0))
			e=int(math.floor(sv_clusters[chrA][chrB][cluster]["endA"]/50.0))+1
			avg_a=numpy.average(coverage_data[chrA][s:e])

			if avg_a > args.max_coverage*library[ "avg_coverage_{}".format(chrA) ]:
				continue
			elif (args.max_coverage*n_discordants/avg_a < args.p_ratio/2 and args.max_coverage*n_splits/avg_a < args.r_ratio/2) and not n_contigs:
				continue

			s=int(math.floor(sv_clusters[chrA][chrB][cluster]["startB"]/50.0))
			e=int(math.floor(sv_clusters[chrA][chrB][cluster]["endB"]/50.0))+1

			avg_b=numpy.average(coverage_data[chrB][s:e])

			if avg_b == 0:
				continue
			elif avg_b > args.max_coverage*library[ "avg_coverage_{}".format(chrB) ]:
				continue
			elif (args.max_coverage*n_discordants/avg_b < args.p_ratio/2 and args.max_coverage*n_splits/avg_b < args.r_ratio/2) and not n_contigs:
				continue

			var_n+=1
			sample_data={}
			for sample in samples:

				coverageA,frac_low_qA,n_discsA,n_splitsA,crossing_f_A,crossing_r_A=get_region(samfile,chrA,sv_clusters[chrA][chrB][cluster]["startA"],sv_clusters[chrA][chrB][cluster]["endA"],posA,min_mapq,max_ins_len,contig_number)
				coverageB,frac_low_qB,n_discsB,n_splitsB,crossing_f_B,crossing_r_B=get_region(samfile,chrB,sv_clusters[chrA][chrB][cluster]["startB"],sv_clusters[chrA][chrB][cluster]["endB"],posB,min_mapq,max_ins_len,contig_number)

				sample_data[sample]={}
				sample_data[sample]={"covA":coverageA,"QA":frac_low_qA,"discA":n_discsA,"splitA":n_splitsA,"refRA":crossing_r_A,"refFA":crossing_f_A}
				sample_data[sample].update({"covB":coverageB,"QB":frac_low_qB,"discB":n_discsB,"splitB":n_splitsB,"refRB":crossing_r_B,"refFB":crossing_f_B})

				if chrA != chrB:
					sample_data[sample]["covM"]=0
				elif abs(posB - posA) < 1000:
					if posA < posB:
						coverageM,_,_,_,_,_=get_region(samfile,chrA,posA,posB,posA,min_mapq,max_ins_len,contig_number)
					else:
						coverageM,_,_,_,_,_=get_region(samfile,chrA,posB,posA,posB,min_mapq,max_ins_len,contig_number)

					sample_data[sample]["covM"]=coverageM
				else:
					s=int(math.floor(posA/50.0))
					e=int(math.floor(posB/50.0))+1
					coverage_between=coverage_data[chrA][s:e]
					gc_between=gc[chrA][s:e]
					coverage_between=coverage_between[ gc_between > -1 ]
					if len(coverage_between) > 4:
						sample_data[sample]["covM"]=numpy.average(coverage_between)
					else:
						sample_data[sample]["covM"]=library[ "avg_coverage_{}".format(chrA) ]
						

			inverted=0
			non_inverted=0
			for i in range(0,len(sv_clusters[chrA][chrB][cluster]["positions_A"]["orientation_discordants"]) ):
				if sv_clusters[chrA][chrB][cluster]["positions_A"]["orientation_discordants"][i] == sv_clusters[chrA][chrB][cluster]["positions_B"]["orientation_discordants"][i]:
					inverted+=1
				else:
					non_inverted+=1

			for i in range(0,len(sv_clusters[chrA][chrB][cluster]["positions_A"]["orientation_splits"]) ):
				if not sv_clusters[chrA][chrB][cluster]["positions_A"]["orientation_splits"][i] == sv_clusters[chrA][chrB][cluster]["positions_B"]["orientation_splits"][i]:
					inverted+=1
				else:
					non_inverted+=1

			for i in range(0,len(sv_clusters[chrA][chrB][cluster]["positions_A"]["orientation_contigs"]) ):
				if not sv_clusters[chrA][chrB][cluster]["positions_A"]["orientation_contigs"][i] == sv_clusters[chrA][chrB][cluster]["positions_B"]["orientation_contigs"][i]:
					inverted+=1
				else:
					non_inverted+=1

			svtype,cn=find_sv_type(chrA,chrB,inverted,non_inverted,args,sample_data,samples,library)

			filt=sv_filter(sample_data,args,chrA,chrB,posA,posB,max_ins_len,n_discordants,n_splits,library,sample_data[sample]["discA"],sample_data[sample]["discB"],sample_data[sample]["splitA"],sample_data[sample]["splitB"],n_contigs)
			format_col="GT:CN:COV:DV:RV:LQ:RR:DR"

			#configure filters for CNV based on Read depth
			for sample in samples:

				covA=sample_data[sample]["covA"]
				covM=sample_data[sample]["covM"]
				covB=sample_data[sample]["covB"]

				if "DEL" in svtype:
					#homozygout del based on coverage
					if cn == 0:
						filt="PASS"


					#normal coverage on the flanking regions, abnormal inbetween
					if covA > covM*(cn+0.9) and covB > covM*(cn+0.9):
						filt="PASS"

				#too few reads, but clear DR signal
				elif "DUP" in svtype and filt == "BelowExpectedLinks":
					filt="PASS"
				scoring_dict={"n_contigs":n_contigs, "n_discordants":n_discordants,"n_splits":n_splits,"covA":covA,"covM":covM,"covB":covB,"refRA":sample_data[sample]["refRA"],"refRB":sample_data[sample]["refRB"],"refFA":sample_data[sample]["refFA"],"refFB":sample_data[sample]["refFB"]}

			if svtype != "BND":
				info=["SVTYPE={}".format(svtype),"SVLEN={}".format(posB-posA),"END={}".format(posB)]
				alt="<{}>".format(svtype)

				info+=["REGIONA={},{}".format(sv_clusters[chrA][chrB][cluster]["startA"],sv_clusters[chrA][chrB][cluster]["endA"])]
				info+=["REGIONB={},{}".format(sv_clusters[chrA][chrB][cluster]["startB"],sv_clusters[chrA][chrB][cluster]["endB"])]
				info+=["LFA={},{}".format(sample_data[sample]["discA"],sample_data[sample]["splitA"])]
				info+=["LFB={},{}".format(sample_data[sample]["discB"],sample_data[sample]["splitB"])]
				info+=["LTE={},{}".format(n_discordants,n_splits)]

				if n_contigs:
					for c in sv_clusters[chrA][chrB][cluster]["contigs"]:
						if "_d_" in c:
							c=c.split("_d_")[0]

						ctgs=[ contig_seqs[c] ]
					info+=["CTG={}".format("|".join(ctgs) )]

				else:
					info+=["CTG=."]



				info=";".join(info)
				variant=[chrA,str(posA),"SV_{}_1".format(var_n),"N",alt,".",filt,info,format_col]
				for sample in samples:
					GT="./."

					if len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) >= args.r or len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) >= args.p:
						GT="0/1"
					if sample_data[sample]["refRB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) or sample_data[sample]["refRA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]):
						GT="1/1"
					if sample_data[sample]["refFB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) or sample_data[sample]["refFA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]):
						GT="1/1"
					if n_contigs and (not len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) and not len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample])):
						if sample_data[sample]["covB"]:
							if sample_data[sample]["refRB"]/sample_data[sample]["covB"] < 0.2:
								GT="1/1"
							else:
								GT="0/1"
						else:
							GT="1/1"

						if sample_data[sample]["covA"]:
							if sample_data[sample]["refRA"]/sample_data[sample]["covA"] < 0.2:
								GT="1/1"
							else:
								GT="0/1"
						else:
							GT="1/1"


					if "DEL" in alt:
						if cn == 0:
							GT = "1/1"
						else:
							GT= "0/1"
					elif "DUP" in alt:
						if cn >= 2*library["contig_ploidy_{}".format(chrA)]:
							GT= "1/1"
						else:
							GT="0/1"

					variant.append( "{}:{}:{},{},{}:{}:{}:{},{}:{},{}:{},{}".format(GT,cn,sample_data[sample]["covA"],sample_data[sample]["covM"],sample_data[sample]["covB"],n_discordants,n_splits,sample_data[sample]["QA"],sample_data[sample]["QB"],sample_data[sample]["refRA"],sample_data[sample]["refRB"],sample_data[sample]["refFA"],sample_data[sample]["refFB"]) )
				variants.append([chrA,posA,variant,scoring_dict])
			else:
				info=["SVTYPE=BND".format(svtype)]
				inverted=False
				before=True

				if posA == sv_clusters[chrA][chrB][cluster]["endA"]:
					before=False

				if inverted  > non_inverted:
					inverted=True

				if not inverted and not before:
					alt_str_a="N[{}:{}[".format(chrB,posB)
					alt_str_b="]{}:{}]N".format(chrA,posA)
				elif not inverted and before:
					alt_str_a="]{}:{}]N".format(chrB,posB)
					alt_str_b="N[{}:{}[".format(chrA,posA)
				elif inverted and  not before:
					alt_str_a="N]{}:{}]".format(chrB,posB)
					alt_str_b="[{}:{}[N".format(chrA,posA)
				else:
					alt_str_a="[{}:{}[N".format(chrB,posB)
					alt_str_b="N]{}:{}]".format(chrA,posA)

				info+=["REGIONA={},{}".format(sv_clusters[chrA][chrB][cluster]["startA"],sv_clusters[chrA][chrB][cluster]["endA"])]
				info+=["REGIONB={},{}".format(sv_clusters[chrA][chrB][cluster]["startB"],sv_clusters[chrA][chrB][cluster]["endB"])]
				info+=["LFA={},{}".format(sample_data[sample]["discA"],sample_data[sample]["splitA"])]
				info+=["LFB={},{}".format(sample_data[sample]["discA"],sample_data[sample]["splitA"])]
				info+=["LTE={},{}".format(n_discordants,n_splits)]

				if n_contigs:
					for c in sv_clusters[chrA][chrB][cluster]["contigs"]:
						if "_d_" in c:
							c=c.split("_d_")[0]
 
						ctgs=[ contig_seqs[c] ]
					info+=["CTG={}".format("|".join(ctgs) )]

				else:
					info+=["CTG=."]



				info=";".join(info)
				variant=[chrA,str(posA),"SV_{}_1".format(var_n),"N",alt_str_a,".",filt,info,format_col]
				for sample in samples:
					GT="./."
					if len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) >= args.r or len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) >= args.p:
						GT="0/1"
					if sample_data[sample]["refRB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) or sample_data[sample]["refRA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]):
						GT="1/1"
					if sample_data[sample]["refFB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) or sample_data[sample]["refFA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]):
						GT="1/1"
					if n_contigs and (not len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) and not len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample])):
						if sample_data[sample]["covB"]:
							if sample_data[sample]["refRB"]/sample_data[sample]["covB"] < 0.2:
								GT="1/1"
							else:
								GT="0/1"
						else:
								GT="1/1"

						if sample_data[sample]["covA"]:
							if sample_data[sample]["refRA"]/sample_data[sample]["covA"] < 0.2:
								GT="1/1"
							else:
								GT="0/1"


						else:
								GT="1/1"



					variant.append( "{}:{}:{},{},{}:{}:{}:{},{}:{},{}:{},{}".format(GT,cn,sample_data[sample]["covA"],sample_data[sample]["covM"],sample_data[sample]["covB"],n_discordants,n_splits,sample_data[sample]["QA"],sample_data[sample]["QB"],sample_data[sample]["refRA"],sample_data[sample]["refRB"],sample_data[sample]["refFA"],sample_data[sample]["refFB"]) )
				variants.append([chrA,posA,variant,scoring_dict])


				variant=[chrB,str(posB),"SV_{}_2".format(var_n),"N",alt_str_b,".",filt,info,format_col]
				for sample in samples:
					GT="./."
					if len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) >= args.r or len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) >= args.p:
						GT="0/1"
					if sample_data[sample]["refRB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) or sample_data[sample]["refRA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]):
						GT="1/1"
					if sample_data[sample]["refFB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) or sample_data[sample]["refFA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]):
						GT="1/1"
					if n_contigs and (not len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) and not len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample])):
						if sample_data[sample]["covB"]:
							if sample_data[sample]["refRB"]/sample_data[sample]["covB"] < 0.2:
								GT="1/1"
							else:
								GT="0/1"
						else:
								GT="1/1"

						if sample_data[sample]["covA"]:
							if sample_data[sample]["refRA"]/sample_data[sample]["covA"] < 0.2:
								GT="1/1"
							else:
								GT="0/1"
						else:
								GT="1/1"




					variant.append( "{}:{}:{},{},{}:{}:{}:{},{}:{},{}:{},{}".format(GT,cn,sample_data[sample]["covA"],sample_data[sample]["covM"],sample_data[sample]["covB"],n_discordants,n_splits,sample_data[sample]["QA"],sample_data[sample]["QB"],sample_data[sample]["refRA"],sample_data[sample]["refRB"],sample_data[sample]["refFA"],sample_data[sample]["refFB"]) )
				variants.append([chrB,posB,variant, scoring_dict ])

	samfile.close()
	return(variants)

def main(str bam_file_name,dict sv_clusters,args,dict library,int min_mapq,samples,dict coverage_data,contig_number,max_ins_len,gc):
	contig_seqs={}
	new_seq=False
	if not args.skip_assembly:
		for line in open("{}_tiddit/clips.fa.assembly.clean.mag".format(args.o)):

			if line[0] == ">":
				name=line[1:].rstrip()
			else:
				contig_seqs[name]=line.strip("\n")

			#if not new_seq and line[0] == "@" and "\t" in line:
			#	name=line.split("\t")[0][1:]
			#	new_seq=True

			#elif new_seq:
			#	contig_seqs[name]=line.strip("\n")
			#	new_seq=False


	variants={}
	for chrA in sv_clusters:
		variants[chrA]=[]
		for chrB in sv_clusters[chrA]:
			variants[chrB]=[]

	variants_list=Parallel(n_jobs=args.threads,prefer="threads")( delayed(define_variant)(chrA,bam_file_name,sv_clusters,args,library,min_mapq,samples,coverage_data,contig_number,max_ins_len,contig_seqs,gc) for chrA in sv_clusters)

	ratios={"fragments_A":[],"fragments_B":[],"reads_A":[],"reads_B":[]}
	for v in variants_list:
		for variant in v:
			if variant[3]["n_discordants"]:
				ratios["fragments_A"].append(variant[3]["n_discordants"]/(variant[3]["refFA"]+variant[3]["n_discordants"]) )
				ratios["fragments_B"].append(variant[3]["n_discordants"]/(variant[3]["refFB"]+variant[3]["n_discordants"]) )

			if variant[3]["n_splits"]:
				ratios["reads_A"].append(variant[3]["n_splits"]/(variant[3]["refRA"]+variant[3]["n_splits"]) )
				ratios["reads_B"].append(variant[3]["n_splits"]/(variant[3]["refRB"]+variant[3]["n_splits"]) )


	p=[1,5,10,20,30,40,50,60,70,75,80,85,90,95,97.5,99]

	percentiles={"FA":percentile(ratios["fragments_A"],p),"FB":percentile(ratios["fragments_B"],p),"RA":percentile(ratios["reads_A"],p),"RB":percentile(ratios["reads_B"],p)}

	for v in variants_list:
		for variant in v:	
			score=scoring(variant[3],percentiles)
			variant[2][5]=str(score)
			variants[ variant[0] ].append( [ variant[1],variant[2] ] )

	return(variants)
