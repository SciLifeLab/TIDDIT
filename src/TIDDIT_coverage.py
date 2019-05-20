import numpy

#The functions in this script are used to read and proccess the coverage tab file, as well as to perform GC correction

#read the coverage tab file
def coverage(args):
	chromosome_size={}
	for line in open(args.o+".wig"):
		if "chrom=" in line:
			chromosome=line.split("chrom=")[-1].split()[0]
			continue
		elif "quality" in line:
			break
		elif "descrip" in line:
			continue

		if not chromosome in chromosome_size:
			chromosome_size[chromosome]=0
		chromosome_size[chromosome]+=1
	
	coverage_data={}
	for chromosome in chromosome_size:
		coverage_data[chromosome]=numpy.zeros( (chromosome_size[chromosome],3), dtype=numpy.float32 )

	span_coverage={}
	for chromosome in chromosome_size:
		span_coverage[chromosome]=numpy.zeros( (chromosome_size[chromosome],2), dtype=numpy.uint32 )


	mapq=False
	spanning_pair=False
	spanning_read=False
	for line in open(args.o+".wig"):
		if "name=\"Coverage\"" in line:
			continue
		elif "chrom=" in line:
			chromosome=line.split("chrom=")[-1].split()[0]
			chromosome_index=0
			continue
		elif "quality" in line:
			mapq=True
			continue
		elif "SpanPairs" in line:
			spanning_pair=True
			continue
		elif "SpanReads" in line:
			spanning_read=True
			continue

		content=line.strip()	
		if mapq and not spanning_read and not spanning_pair:
			coverage_data[chromosome][chromosome_index][1]=float(content)
		elif mapq and spanning_pair and not spanning_read:
			span_coverage[chromosome][chromosome_index][0]=int(content)
		elif mapq and spanning_pair and  spanning_read:
			span_coverage[chromosome][chromosome_index][1]=int(content)
		else:
			coverage_data[chromosome][chromosome_index][0]=float(content)

		chromosome_index+=1

	return(coverage_data,span_coverage)

#estimate the ploidy of each chromosome
def determine_ploidy(args,chromosomes,coverage_data,Ncontent,library_stats):
	library_stats["chr_cov"]={}
	ploidies={}
	avg_coverage=[]

	cov=numpy.array([],dtype=numpy.float32)
	for chromosome in chromosomes:		
		try:
			if args.ref:
				N_count=Ncontent[chromosome]
				chr_cov=coverage_data[chromosome][numpy.where( (N_count > 0) & (coverage_data[chromosome][:,1] > args.Q)  ),0][0]
			else:
				chr_cov=coverage_data[chromosome][numpy.where( coverage_data[chromosome][:,1] > args.Q  ),0][0]

			if len(chr_cov):
				chromosomal_average=numpy.median(chr_cov)
				cov = numpy.append(cov,chr_cov)
			else:
				chromosomal_average=0
			library_stats["chr_cov"][chromosome]=chromosomal_average

		except:
			print ("error: reference mismatch!")
			print ("make sure that the contigs of the bam file and the reference match")
			quit()

	if len(cov):
		coverage_norm=numpy.median(cov)
	else:
		coverage_norm=1

	if args.ref:
		coverage_data=gc_norm(args,coverage_norm,chromosomes,coverage_data,Ncontent)

	chromosomal_average=0
	outfile=open(args.o+".ploidy.tab", 'w')
	outfile.write("Contig\tploidy_rounded\tploidy_raw\tmedian_coverage\n")
	for chromosome in chromosomes:
	
		if args.ref:
			N_count=Ncontent[chromosome]
			cov=coverage_data[chromosome][numpy.where( (N_count > -1) & ( (coverage_data[chromosome][:,1] > args.Q) | (coverage_data[chromosome][:,1] == 0) ) ),0]
		else:
			cov=coverage_data[chromosome][numpy.where( (coverage_data[chromosome][:,1] > args.Q) | (coverage_data[chromosome][:,1] == 0) ),0]

		chromosomal_average=numpy.median(cov)
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
	return(ploidies,library_stats,coverage_data)

#normalise the coverage based on GC content
def gc_norm(args,median_coverage,normalising_chromosomes,coverage_data,Ncontent):
	gc_vals=range(0,101)
	gc_dict={}

	for gc in gc_vals:
		tmp=[]
		gc_dict[gc]=median_coverage
		for chromosome in normalising_chromosomes:
			selected= Ncontent[chromosome] == gc
			selected_coverage=coverage_data[chromosome][selected]
			tmp+=list(selected_coverage[ numpy.where( selected_coverage[:,1] > args.q)[0] ,0])
		if len(tmp):
			gc_dict[gc]=numpy.median(tmp)
		if 0 == gc_dict[gc]:
			gc_dict[gc]=median_coverage

	for chromosome in coverage_data:
		for i in range(0,len(coverage_data[chromosome])):
			if not Ncontent[chromosome][i] == -1:
				coverage_data[chromosome][i,0]=median_coverage*coverage_data[chromosome][i,0]/gc_dict[ Ncontent[chromosome][i] ]

	return(coverage_data)

#load the GC tab, and find which bins contain too many N
def retrieve_N_content(args):

	chromosome_size={}
	for line in open(args.o+".gc.wig"):
		if "Per bin GC values" in line:
			continue
		elif "chrom=" in line:
			chromosome=line.split("chrom=")[-1].split()[0]
			continue

		elif "Per bin fraction of N" in line:
			break
		if not chromosome in chromosome_size:
			chromosome_size[chromosome]=0
		chromosome_size[chromosome]+=1
	
	Ncontent={}
	for chromosome in chromosome_size:
		Ncontent[chromosome]=numpy.zeros( chromosome_size[chromosome],dtype=numpy.int8 )
	read_n=False
	for line in open(args.o+".gc.wig"):
		if "Per bin GC values" in line:
			continue
		elif "chrom=" in line:
			chromosome=line.split("chrom=")[-1].split()[0]
			chromosome_index=0
			continue
		elif "Per bin fraction of N" in line:
			read_n=True
			continue

		content=line.strip()
	
		if read_n:
			n=float(content[-1])
			if n > args.n_mask:
				Ncontent[chromosome][chromosome_index]=-1
		else:
			gc=int(round(float(content)*100))
			Ncontent[chromosome][chromosome_index]=gc 

		chromosome_index+=1
	return(Ncontent)
