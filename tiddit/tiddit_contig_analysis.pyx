import pysam
import tiddit_signal

def main(aligned_contigs,prefix,sample_id,min_size):
	samfile = pysam.AlignmentFile(aligned_contigs, "r")
	bam_header=samfile.header
	

	split_contigs={}
	for chrA in bam_header["SQ"]:
		for chrB in bam_header["SQ"]:
			if chrA["SN"] <= chrB["SN"]:

				if not chrA["SN"] in split_contigs:
					split_contigs[chrA["SN"]] = {}

				split_contigs[chrA["SN"]][chrB["SN"]]={}

	for read in samfile.fetch(until_eof=True):

		if read.is_unmapped:
			continue

		if read.has_tag("SA") and not (read.is_supplementary or read.is_secondary):
			split_contigs=tiddit_signal.SA_analysis(read,-2,split_contigs,"SA")
		#elif read.has_tag("XA") and not (read.is_supplementary or read.is_secondary):
		#	split_contigs=tiddit_signal.SA_analysis(read,-2,split_contigs,"XA")
		elif not (read.is_supplementary or read.is_secondary) and len(read.cigartuples) > 2:

			current_bp=read.reference_start
			for i in range(0,len(read.cigartuples)-1):
				if read.cigartuples[i][0] == 2 and read.cigartuples[i][1] > min_size:
					split_contigs[read.reference_name][read.reference_name]["{}_d_{}".format(read.query_name,i)]=[current_bp,read.is_reverse,current_bp+read.cigartuples[i][1],read.is_reverse]
				current_bp+=read.cigartuples[i][1]

	for chrA in split_contigs:
		for chrB in split_contigs[chrA]:
			f=open("{}_tiddit/contigs_{}_{}_{}.tab".format(prefix,sample_id,chrA,chrB),"w")

			for fragment in split_contigs[chrA][chrB]:

				f.write("{}\t{}\n".format(fragment,"\t".join(map(str, split_contigs[chrA][chrB][fragment] )))  )

			f.close()

#main("NA12878_default_tiddit/clips.sam","hej_ctg","test")



