import numpy
from scipy.stats import norm

#This script contains functions for applying and computing the filters and statistics
def expected_links(coverageA,coverageB,sizeA,sizeB,gap,insert_mean,insert_stddev,readLength):
	coverage=numpy.average([coverageA,coverageB])

	b1 = (sizeA + sizeB + gap - insert_mean) / float(insert_stddev)
	a1 = (max([sizeA, sizeB]) + gap + readLength  - insert_mean) / float(insert_stddev)
	b2 = (min([sizeA, sizeB]) + gap + readLength  - insert_mean) / float(insert_stddev)
	a2 = (gap + 2 * readLength - insert_mean) / float(insert_stddev)

	e = Part(a1, b1, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength ) - Part(a2, b2, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength)
	return(e)

#compute the filters
def fetch_filter(chrA,chrB,candidate,args,library_stats):
	filt="PASS"

	#Less than the expected number of signals
	if candidate["discs"] and not ( candidate["splits"] and abs(candidate["posA"]-candidate["posB"]) < 3*library_stats["STDInsertSize"] and chrA == chrB ):
		if candidate["e1"]*0.6 >= candidate["discs"]+candidate["splits"]:
			filt = "BelowExpectedLinks"
		elif candidate["e2"]*library_stats["ReadLength"]/float(library_stats["MeanInsertSize"]) >= candidate["discs"]+candidate["splits"]:
			filt = "BelowExpectedLinks"
	else:
		if candidate["e1"]*0.4 >= (candidate["splits"]+candidate["discs"]):
			filt = "BelowExpectedLinks"
	#The ploidy of this contig is 0, hence there shoud be no variant here
	if library_stats["ploidies"][chrA] == 0 or library_stats["ploidies"][chrB] == 0:
		return("Ploidy")

	#coverage is too high
	if candidate["MaxcovA"] >= library_stats["chr_cov"][chrA]*(library_stats["ploidies"][chrA]+2) or candidate["MaxcovB"] >= library_stats["chr_cov"][chrB]*(library_stats["ploidies"][chrB]+2):
		filt = "UnexpectedCoverage"
	elif candidate["discsA"] > (candidate["discs"]+candidate["splits"])*(library_stats["ploidies"][chrA]+1) or candidate["discsB"] > (candidate["discs"]+candidate["splits"])*(library_stats["ploidies"][chrA]+1):
		filt= "FewLinks"
	elif chrA == chrB and candidate["max_A"] > candidate["min_B"]:
		filt = "Smear"
	elif chrA != chrB or (abs(candidate["posA"]-candidate["posB"]) > library_stats["MeanInsertSize"]+3*library_stats["STDInsertSize"] ):
		if not candidate["discs"] or candidate["discs"]*4 <  candidate["splits"] or candidate["discs"] <= args.p/2:
			filt= "SplitsVSDiscs"
	elif candidate["QRA"] < args.Q or candidate["QRB"] < args.Q:
		filt = "RegionalQ"

	#fewer links than expected
	if chrA == chrB and (abs(candidate["max_A"]-candidate["min_B"]) < args.z):
		filt="MinSize"

	return(filt)

#determine the variant type
def fetch_variant_type(chrA,chrB,candidate,args,library_stats):
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

def Part(a, b, sizeA, sizeB, gap, insert_mean, insert_stddev, coverage, readLength ):

	readfrequency = 2 * readLength / coverage;
	expr1 = (min([sizeA, sizeB]) - (readLength - 0)) / readfrequency * norm.cdf(a, 0, 1);
	expr2 = -(- 0 ) / readfrequency * norm.cdf(b, 0, 1);
	expr3 = (b * insert_stddev) / readfrequency * (norm.cdf(b, 0, 1) - norm.cdf(a, 0, 1));
	expr4 = (insert_stddev / readfrequency) * (norm.pdf(b, 0, 1) - norm.pdf(a, 0, 1));
	value = expr1 + expr2 + expr3 + expr4
	return(value)


