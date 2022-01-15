import time
import math
import numpy

#from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam import AlignmentFile, AlignedSegment

def get_region(samfile,str chr,int start,int end,int bp,int min_q,int max_ins, contig_number):

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

	for read in samfile.fetch(chr, q_start, q_end):
		if read.is_unmapped:
			continue
		if not read.mate_is_unmapped:
			if read.next_reference_start > end and read.reference_start > end:
				continue		
		else:
			if read.reference_start > end:
				continue


		if read.is_duplicate:
			continue
		
		if not (read.reference_start > end):
			n_reads+=1
			if read.mapq < min_q:
				low_q+=1

		if read.mapq < min_q:
			continue

		r_start=read.reference_start
		r_end=read.reference_end

		if r_start < bp-10 and r_end > bp:
			crossing_r+=1

		mate_bp_read= (read.next_reference_start < bp and r_end > bp)
		discordant= ( abs(read.isize) > max_ins or read.next_reference_name != read.reference_name )

		if mate_bp_read and not discordant:
			crossing_f+=1

		if read.reference_end < start:
			continue
		elif read.reference_start > end:
			continue

		if read.reference_start < start:
			r_start=start
		
		if read.reference_end > end:
			r_end=end

		bases+=r_end-r_start+1

		if read.has_tag("SA"):
			n_splits+=1

		if (abs(read.isize) > max_ins) or (read.next_reference_name != read.reference_name) :
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
		cn=int(round(sample_data[sample]["covM"]*p/library[ "avg_coverage_{}".format(chrA) ]))

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
		elif cn < p:
			return("DEL",cn)
		
		elif inverted > non_inverted:
			return("INV",cn)
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

		if n_splits or n_discordants:
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

	return(filt)

def main(str bam_file_name,dict sv_clusters,args,dict library,int min_mapq,samples,dict coverage_data,contig_number,max_ins_len):

	variants=[]
	#cdef AlignmentFile samfile  = AlignmentFile(bam_file_name, "r")
	samfile  = AlignmentFile(bam_file_name, "r")

	var_n=0
	for chrA in sv_clusters:
		for chrB in sv_clusters[chrA]:
			for cluster in sv_clusters[chrA][chrB]:
				n_discordants=sv_clusters[chrA][chrB][cluster]["N_discordants"]
				n_splits=sv_clusters[chrA][chrB][cluster]["N_splits"]
				n_contigs=sv_clusters[chrA][chrB][cluster]["N_contigs"]

				if (n_discordants < args.p and n_splits < args.r) and not n_contigs:
					continue

				posA=sv_clusters[chrA][chrB][cluster]["posA"]
				posB=sv_clusters[chrA][chrB][cluster]["posB"]

				if chrA == chrB and abs(posA-posB) < args.z:
					continue

				s=int(math.floor(sv_clusters[chrA][chrB][cluster]["startA"]/50.0))
				e=int(math.floor(sv_clusters[chrA][chrB][cluster]["endA"]/50.0))+1
				avg_a=numpy.average(coverage_data[chrA][s:e])

				if avg_a > args.max_coverage*library[ "avg_coverage_{}".format(chrA) ]:
					continue
				elif (args.max_coverage*n_discordants/avg_a < args.p_ratio/2 and args.max_coverage*n_splits/avg_a < args.r_ratio/2) and not n_contigs:
					continue

				avg_b=numpy.average(coverage_data[chrA][s:e])
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
					sample_data[sample].update({"covB":coverageB,"QB":frac_low_qA,"discB":n_discsB,"splitB":n_splitsB,"refRB":crossing_r_B,"refFB":crossing_f_B})

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
						sample_data[sample]["covM"]=numpy.average(coverage_data[chrA][s:e] )
				
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
				format_col="GT:CN:DP:DR:SR:LQ:RR:RD"

				#configure filters for CNV based on Read depth
				for sample in samples:
					if "DEL" in svtype:
						#homozygout del based on coverage 
						if cn == 0:
							filt="PASS"

						covA=sample_data[sample]["covA"]
						covM=sample_data[sample]["covM"]
						covB=sample_data[sample]["covB"]

						#normal coverage on the flanking regions, abnormal inbetween
						if covA > covM*(cn+0.9) and covB > covM*(cn+0.9):
							filt="PASS"

					#too few reads, but clear RD signal
					elif "DUP" in svtype and filt == "BelowExpectedLinks":
						filt="PASS"



				if svtype != "BND":
					info=["SVTYPE={}".format(svtype),"SVLEN={}".format(posB-posA),"END={}".format(posB)]
					alt="<{}>".format(svtype)

					info+=["REGIONA={},{}".format(sv_clusters[chrA][chrB][cluster]["startA"],sv_clusters[chrA][chrB][cluster]["endA"])]
					info+=["REGIONB={},{}".format(sv_clusters[chrA][chrB][cluster]["startB"],sv_clusters[chrA][chrB][cluster]["endB"])]
					info+=["LFA={},{}".format(sample_data[sample]["discA"],sample_data[sample]["splitA"])]
					info+=["LFB={},{}".format(sample_data[sample]["discB"],sample_data[sample]["splitB"])]
					info+=["LTE={},{}".format(n_discordants,n_splits)]

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

						if "DEL" in alt:
							if cn == 0:
								GT = "1/1"
						elif "DUP" in alt:
							if cn >= 2*library["contig_ploidy_{}".format(chrA)]:
								GT= "1/1"
							
						variant.append( "{}:{}:{},{},{}:{}:{}:{},{}:{},{}:{},{}".format(GT,cn,sample_data[sample]["covA"],sample_data[sample]["covM"],sample_data[sample]["covB"],n_discordants,n_splits,sample_data[sample]["QA"],sample_data[sample]["QB"],sample_data[sample]["refRA"],sample_data[sample]["refRB"],sample_data[sample]["refFA"],sample_data[sample]["refFB"]) )
					variants.append(variant)				

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

						variant.append( "{}:{}:{},{},{}:{}:{}:{},{}:{},{}:{},{}".format(GT,cn,sample_data[sample]["covA"],sample_data[sample]["covM"],sample_data[sample]["covB"],n_discordants,n_splits,sample_data[sample]["QA"],sample_data[sample]["QB"],sample_data[sample]["refRA"],sample_data[sample]["refRB"],sample_data[sample]["refFA"],sample_data[sample]["refFB"]) )
					variants.append(variant)				

					variant=[chrB,str(posB),"SV_{}_2".format(var_n),"N",alt_str_b,".",filt,info,format_col]
					for sample in samples:
						if len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) >= args.r or len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) >= args.p:
							GT="0/1"
						if sample_data[sample]["refRB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]) or sample_data[sample]["refRA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_splits"][sample]):
							GT="1/1"
						if sample_data[sample]["refFB"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]) or sample_data[sample]["refFA"] < 0.1*len(sv_clusters[chrA][chrB][cluster]["sample_discordants"][sample]):
							GT="1/1"

						GT="./."
						variant.append( "{}:{}:{},{},{}:{}:{}:{},{}:{},{}:{},{}".format(GT,cn,sample_data[sample]["covA"],sample_data[sample]["covM"],sample_data[sample]["covB"],n_discordants,n_splits,sample_data[sample]["QA"],sample_data[sample]["QB"],sample_data[sample]["refRA"],sample_data[sample]["refRB"],sample_data[sample]["refFA"],sample_data[sample]["refFB"]) )
					variants.append(variant)				

	samfile.close()
	return(variants)
