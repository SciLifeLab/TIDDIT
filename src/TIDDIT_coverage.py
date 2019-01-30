import numpy

#The functions in this script are used to read and proccess the coverage tab file, as well as to perform GC correction

#read the coverage tab file
def coverage(args):
	chromosome_size={}
	for line in open(args.o+".tab"):
		if line[0] == "#":
			continue
		content=line.strip().split()
		if not content[0] in chromosome_size:
			chromosome_size[content[0]]=0
		chromosome_size[content[0]]+=1
	
	coverage_data={}
	chromosome_index={}
	for chromosome in chromosome_size:
		coverage_data[chromosome]=numpy.zeros( (chromosome_size[chromosome],3) )
		chromosome_index[chromosome]=0

	for line in open(args.o+".tab"):
		if line[0] == "#":
			continue
		content=line.strip().split()
		
		coverage_data[content[0]][chromosome_index[content[0]]][0]=float(content[3])
		coverage_data[content[0]][chromosome_index[content[0]]][1]=float(content[4])
		chromosome_index[content[0]]+=1

	return(coverage_data)

#estimate the ploidy of each chromosome
def determine_ploidy(args,chromosomes,coverage_data,Ncontent,library_stats):
	library_stats["chr_cov"]={}
	ploidies={}
	avg_coverage=[]
	cov=[]
	for chromosome in chromosomes:		
		try:
			if args.ref:
				N_count=Ncontent[chromosome]
				chr_cov=coverage_data[chromosome][numpy.where( (N_count > 0) & (coverage_data[chromosome][:,1] > args.Q)  ),0][0]
			else:
				chr_cov=coverage_data[chromosome][numpy.where( coverage_data[chromosome][:,1] > args.Q  ),0][0]

			if len(chr_cov):
				chromosomal_average=numpy.median(chr_cov)
				cov+= list(chr_cov)
			else:
				chromosomal_average=0
			library_stats["chr_cov"][chromosome]=chromosomal_average

		except:
			print "error: reference mismatch!"
			print "make sure that the contigs of the bam file and the reference match"
			quit()

	cov=numpy.array(cov)
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
	gc_vals=[0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 , 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3 , 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4 , 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5 , 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6 , 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7 , 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8 , 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9 , 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1]
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
	for line in open(args.o+".gc.tab"):
		if line[0] == "#":
			continue
		content=line.strip().split()
		if not content[0] in chromosome_size:
			chromosome_size[content[0]]=0
		chromosome_size[content[0]]+=1
	
	Ncontent={}
	chromosome_index={}
	for chromosome in chromosome_size:
		Ncontent[chromosome]=numpy.zeros( chromosome_size[chromosome] )
		chromosome_index[chromosome]=0

	for line in open(args.o+".gc.tab"):
		if line[0] == "#":
			continue

		content=line.strip().split()
		contig=content[0]
		gc=round(float(content[-2]),2)
		n=float(content[-1])
		if n > args.n_mask:
			Ncontent[contig][chromosome_index[contig]]=-1
		else:  
			Ncontent[contig][chromosome_index[contig]]=gc 
		chromosome_index[contig]+=1

	return(Ncontent)
