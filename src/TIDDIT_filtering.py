import numpy

#compute the filters
def fetch_filter(chrA,chrB,candidate,args,library_stats):
	filt="PASS"

	#Less than the expected number of signals
	if candidate["ratio"] <= 0.2 and candidate["discs"] > candidate["splits"]:
		filt = "BelowExpectedLinks"
	elif candidate["ratio"] <= 0.1 and candidate["discs"] < candidate["splits"]:
		filt = "BelowExpectedLinks"

	#The ploidy of this contig is 0, hence there shoud be no variant here
	if library_stats["ploidies"][chrA] == 0 or library_stats["ploidies"][chrB] == 0:
		return("Ploidy")

	#coverage is too high
	if candidate["covA"] >= library_stats["chr_cov"][chrA]*(library_stats["ploidies"][chrA]*4+library_stats["ploidies"][chrA]) or candidate["covB"] >= library_stats["chr_cov"][chrB]*(library_stats["ploidies"][chrB]*4+library_stats["ploidies"][chrA]):
		filt = "UnexpectedCoverage"
	elif candidate["discsA"] > (candidate["discs"]+candidate["splits"])*(library_stats["ploidies"][chrA]*2) or candidate["discsB"] > (candidate["discs"]+candidate["splits"])*(library_stats["ploidies"][chrA]*2):
		filt= "FewLinks"
	elif chrA == chrB and candidate["max_A"] > candidate["min_B"]:
		filt = "Smear"
	elif chrA != chrB or (abs(candidate["posA"]-candidate["posB"]) > library_stats["MeanInsertSize"]+3*library_stats["STDInsertSize"] ):
		if not candidate["discs"] or candidate["discs"]*4 <  candidate["splits"] or candidate["discs"] <= args.p/2:
			filt= "SplitsVSDiscs"
	#fewer links than expected
	if chrA == chrB and (abs(candidate["max_A"]-candidate["min_B"]) < args.z):
		filt="MinSize"

	return(filt)

#determine the variant type
def fetch_variant_type(chrA,chrB,candidate,args,library_stats,disc_ratio,split_ratio):
	variant_type="SVTYPE=BND"
	var="N[{}:{}[".format(chrB,candidate["posB"])
	if not library_stats["ploidies"][chrA] and chrA == chrB:
		variant_type="SVTYPE=DUP"
		var="<DUP>"		

	if chrA == chrB and library_stats["ploidies"][chrA]:
		ploidy=library_stats["ploidies"][chrA]
		if ploidy > 10:
			if candidate["discs"] and abs(candidate["covM"]/library_stats["chr_cov"][chrA]-1) < 0.05:
				if candidate["FF"] + candidate["RR"] > candidate["RF"] + candidate["FR"]:
					variant_type="SVTYPE=INV"
					var="<INV>"
			elif not candidate["discs"] and abs(candidate["covM"]/library_stats["chr_cov"][chrA]-1) < 0.05:
				if candidate["splitsINV"] > candidate["splits"]-candidate["splitsINV"]:
					variant_type="SVTYPE=INV"
					var="<INV>"
			elif candidate["covM"]/library_stats["chr_cov"][chrA]-1 > 0.05:
				variant_type="SVTYPE=DUP"
				var="<DUP>"
			elif candidate["covM"]/library_stats["chr_cov"][chrA]-1 < -0.05:
				variant_type="SVTYPE=DEL"
				var="<DEL>"		

		elif candidate["discs"]:
			if candidate["FF"] + candidate["RR"] > candidate["RF"] + candidate["FR"]:
				variant_type="SVTYPE=INV"
				var="<INV>"
			elif library_stats["Orientation"] == "innie":
				if candidate["covM"]/library_stats["chr_cov"][chrA] > (ploidy+0.5)/float(ploidy):
					variant_type="SVTYPE=DUP"
					var="<DUP>"
					if candidate["RF"] > candidate["FR"]: 
						variant_type="SVTYPE=TDUP"
						var="<TDUP>"
				elif candidate["covM"]/library_stats["chr_cov"][chrA] < (ploidy-0.5)/float(ploidy):
					variant_type="SVTYPE=DEL"
					var="<DEL>"	

			else:
				if candidate["covM"]/library_stats["chr_cov"][chrA] > (ploidy+0.5)/float(ploidy):
					variant_type="SVTYPE=DUP"
					var="<DUP>"
					if candidate["RF"] < candidate["FR"]: 
						variant_type="SVTYPE=TDUP"
						var="<TDUP>"

				elif candidate["covM"]/library_stats["chr_cov"][chrA] < (ploidy-0.5)/float(ploidy):
					variant_type="SVTYPE=DEL"
					var="<DEL>"
		else:
			if candidate["splitsINV"] > candidate["splits"]-candidate["splitsINV"]:
				variant_type="SVTYPE=INV"
				var="<INV>"
			elif candidate["covM"]/library_stats["chr_cov"][chrA] >(ploidy+0.5)/float(ploidy):
					variant_type="SVTYPE=DUP"
					var="<DUP>"
			elif candidate["covM"]/library_stats["chr_cov"][chrA] < (ploidy-0.5)/float(ploidy):
					variant_type="SVTYPE=DEL"
					var="<DEL>"		


	if candidate["discs"] > candidate["splits"]:

		if disc_ratio >= 0.8:
			GT="1/1"
		else:
			GT="0/1"
	else:
		if split_ratio >= 0.8:
			GT="1/1"
		else:
			GT="0/1"

	if ("DUP" in var or var == "<DEL>") and library_stats["chr_cov"][chrA] != 0:
		if var == "<DEL>" and candidate["covM"]/library_stats["chr_cov"][chrA] < 0.2:
			GT="1/1"
		elif "DUP" in var and candidate["covM"]/library_stats["chr_cov"][chrA] > 1.8: 
			GT="1/1"
		else:
			gt="0/1"

	return(var,variant_type,GT)
