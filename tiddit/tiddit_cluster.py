import sys
import os
import DBSCAN
import numpy
import statistics 
from statistics import mode

def find_discordant_pos(fragment,is_mp):
	if is_mp:
		if fragment[3] == "False" and fragment[6] == "True":
			posA=fragment[2]
			posB=fragment[4]
		elif fragment[3] == "False" and fragment[6] == "False":
			posA=fragment[1]
			posB=fragment[4]
		elif fragment[3] == "True" and fragment[6] == "True":		
			posA=fragment[2]
			posB=fragment[5]
		else:
			posA=fragment[1]
			posB=fragment[5]

	else:
		if fragment[3] == "False" and fragment[6] == "True":
			posA=fragment[2]
			posB=fragment[4]
		elif fragment[3] == "False" and fragment[6] == "False":
			posA=fragment[2]
			posB=fragment[5]
		elif fragment[3] == "True" and fragment[6] == "True":
			posA=fragment[1]
			posB=fragment[4]

		else:
			posA=fragment[1]
			posB=fragment[5]

	return(posA,posB)

def main(prefix,chromosomes,samples,is_mp,epsilon,m,max_ins_len):

	discordants={}
	splits={}
	candidates={}
	for chrA in chromosomes:

		discordants[chrA]={}
		splits[chrA]={}
		candidates[chrA] = {}

		for chrB in chromosomes:
			if chrA > chrB:
				continue

			candidates[chrA][chrB]={}
			discordants[chrA][chrB]=[]
			splits[chrA][chrB]=[]

			positions=[]
			i=0
			contigs=set([])
			for sample in samples:
				discordants_path="{}_tiddit/discordants_{}_{}_{}.tab".format(prefix,sample,chrA,chrB)
				splits_path="{}_tiddit/splits_{}_{}_{}.tab".format(prefix,sample,chrA,chrB)
				contigs_path="{}_tiddit/contigs_{}_{}_{}.tab".format(prefix,sample,chrA,chrB)

				if not  os.path.exists(discordants_path):
					continue
				if not  os.path.exists(splits_path):
					continue

				for line in open(discordants_path):
					content=line.rstrip().split("\t")
					posA,posB=find_discordant_pos(content,is_mp)
					discordants[chrA][chrB].append([content[0],sample,"D",posA,content[3],posB,content[6],i])
					positions.append([int(posA),int(posB),i])
					i+=1

				for line in open(splits_path):
					content=line.rstrip().split("\t")
					posA=content[1]
					posB=content[3]

					discordants[chrA][chrB].append([content[0],sample,"S",posA,content[2],posB,content[4],i])
					positions.append([int(posA),int(posB),i])
					i+=1
				for line in open(contigs_path):
					content=line.rstrip().split("\t")
					posA=content[1]
					posB=content[3]

					discordants[chrA][chrB].append([content[0],sample,"A",posA,content[2],posB,content[4],i])
					positions.append([int(posA),int(posB),i])
					contigs.add(i)
					i+=1
				
			positions=numpy.array(sorted(positions,key=lambda l:l[0]))
			#print(positions[:,[0,1]])
			
			clusters=DBSCAN.main(positions,epsilon,m)
			positions=list(positions)
			cluster_pos=[]
			for i in range(0,len(positions)):
				cluster_pos.append(list(positions[i])+[clusters[i]] )
				#positions[i].append(clusters[i])

			cluster_pos= sorted(cluster_pos,key=lambda l:l[2])
			n_ctg_clusters=0
			for i in range(0,len(cluster_pos)):
				candidate=int(cluster_pos[i][-1])
				if candidate == -1 and not (chrA == chrB and discordants[chrA][chrB][i][2] == "A" and ( int(discordants[chrA][chrB][i][5])- int(discordants[chrA][chrB][i][3]) ) < max_ins_len*2):
					continue
				elif candidate == -1 and discordants[chrA][chrB][i][2] == "A":
					candidate=len(cluster_pos)+n_ctg_clusters
					n_ctg_clusters+=1

				if not candidate in candidates[chrA][chrB]:
					candidates[chrA][chrB][candidate]={}
					candidates[chrA][chrB][candidate]["signal_type"]={}
					candidates[chrA][chrB][candidate]["samples"]=set([])
					candidates[chrA][chrB][candidate]["sample_discordants"]={}
					candidates[chrA][chrB][candidate]["sample_splits"]={}
					candidates[chrA][chrB][candidate]["sample_contigs"]={}


					candidates[chrA][chrB][candidate]["N_discordants"]=0
					candidates[chrA][chrB][candidate]["discordants"]=set([])
					candidates[chrA][chrB][candidate]["N_splits"]=0
					candidates[chrA][chrB][candidate]["splits"]=set([])
					candidates[chrA][chrB][candidate]["N_contigs"]=0
					candidates[chrA][chrB][candidate]["contigs"]=set([])


					candidates[chrA][chrB][candidate]["n_signals"]=0

					candidates[chrA][chrB][candidate]["posA"]=0
					candidates[chrA][chrB][candidate]["positions_A"]={}
					candidates[chrA][chrB][candidate]["positions_A"]["contigs"]=[]
					candidates[chrA][chrB][candidate]["positions_A"]["splits"]=[]
					candidates[chrA][chrB][candidate]["positions_A"]["discordants"]=[]
					candidates[chrA][chrB][candidate]["positions_A"]["orientation_contigs"]=[]
					candidates[chrA][chrB][candidate]["positions_A"]["orientation_splits"]=[]
					candidates[chrA][chrB][candidate]["positions_A"]["orientation_discordants"]=[]
					candidates[chrA][chrB][candidate]["start_A"]=0
					candidates[chrA][chrB][candidate]["end_A"]=0

					candidates[chrA][chrB][candidate]["posB"]=0
					candidates[chrA][chrB][candidate]["positions_B"]={}
					candidates[chrA][chrB][candidate]["positions_B"]["contigs"]=[]
					candidates[chrA][chrB][candidate]["positions_B"]["splits"]=[]
					candidates[chrA][chrB][candidate]["positions_B"]["discordants"]=[]
					candidates[chrA][chrB][candidate]["positions_B"]["orientation_contigs"]=[]
					candidates[chrA][chrB][candidate]["positions_B"]["orientation_splits"]=[]
					candidates[chrA][chrB][candidate]["positions_B"]["orientation_discordants"]=[]

					candidates[chrA][chrB][candidate]["start_B"]=0
					candidates[chrA][chrB][candidate]["end_B"]=0

				if not discordants[chrA][chrB][i][1] in candidates[chrA][chrB][candidate]["samples"]:
					candidates[chrA][chrB][candidate]["sample_discordants"][discordants[chrA][chrB][i][1]]=set([])
					candidates[chrA][chrB][candidate]["sample_splits"][discordants[chrA][chrB][i][1]]=set([])
					candidates[chrA][chrB][candidate]["sample_contigs"][discordants[chrA][chrB][i][1]]=set([])

				candidates[chrA][chrB][candidate]["samples"].add(discordants[chrA][chrB][i][1])
				
				if discordants[chrA][chrB][i][2] == "D":
					candidates[chrA][chrB][candidate]["discordants"].add(discordants[chrA][chrB][i][0])
					candidates[chrA][chrB][candidate]["positions_A"]["discordants"].append(int(discordants[chrA][chrB][i][3]))
					candidates[chrA][chrB][candidate]["positions_A"]["orientation_discordants"].append(discordants[chrA][chrB][i][4])

					candidates[chrA][chrB][candidate]["positions_B"]["discordants"].append(int(discordants[chrA][chrB][i][5]))
					candidates[chrA][chrB][candidate]["positions_B"]["orientation_discordants"].append(discordants[chrA][chrB][i][6])
					candidates[chrA][chrB][candidate]["sample_discordants"][discordants[chrA][chrB][i][1]].add(discordants[chrA][chrB][i][0])

				elif discordants[chrA][chrB][i][2] == "S":
					candidates[chrA][chrB][candidate]["splits"].add(discordants[chrA][chrB][i][0])
					candidates[chrA][chrB][candidate]["positions_A"]["splits"].append(int(discordants[chrA][chrB][i][3]))
					candidates[chrA][chrB][candidate]["positions_A"]["orientation_splits"].append(discordants[chrA][chrB][i][4])

					candidates[chrA][chrB][candidate]["positions_B"]["splits"].append(int(discordants[chrA][chrB][i][5]))
					candidates[chrA][chrB][candidate]["positions_B"]["orientation_splits"].append(discordants[chrA][chrB][i][6])
					candidates[chrA][chrB][candidate]["sample_splits"][discordants[chrA][chrB][i][1]].add(discordants[chrA][chrB][i][0])
				else:
					candidates[chrA][chrB][candidate]["contigs"].add(discordants[chrA][chrB][i][0])
					candidates[chrA][chrB][candidate]["positions_A"]["contigs"].append(int(discordants[chrA][chrB][i][3]))
					candidates[chrA][chrB][candidate]["positions_A"]["orientation_contigs"].append(discordants[chrA][chrB][i][4])

					candidates[chrA][chrB][candidate]["positions_B"]["contigs"].append(int(discordants[chrA][chrB][i][5]))
					candidates[chrA][chrB][candidate]["positions_B"]["orientation_contigs"].append(discordants[chrA][chrB][i][6])
					candidates[chrA][chrB][candidate]["sample_contigs"][discordants[chrA][chrB][i][1]].add(discordants[chrA][chrB][i][0])



	for chrA in candidates:
		for chrB in candidates[chrA]:
			for candidate in candidates[chrA][chrB]:
				candidates[chrA][chrB][candidate]["N_discordants"]=len(candidates[chrA][chrB][candidate]["discordants"])
				candidates[chrA][chrB][candidate]["N_splits"]=len(candidates[chrA][chrB][candidate]["splits"])
				candidates[chrA][chrB][candidate]["N_contigs"]=len(candidates[chrA][chrB][candidate]["contigs"])

				if candidates[chrA][chrB][candidate]["N_contigs"]:
					candidates[chrA][chrB][candidate]["posA"]=mode(candidates[chrA][chrB][candidate]["positions_A"]["contigs"])
					candidates[chrA][chrB][candidate]["posB"]=mode(candidates[chrA][chrB][candidate]["positions_B"]["contigs"])
				elif candidates[chrA][chrB][candidate]["N_splits"]:
					candidates[chrA][chrB][candidate]["posA"]=mode(candidates[chrA][chrB][candidate]["positions_A"]["splits"])
					candidates[chrA][chrB][candidate]["posB"]=mode(candidates[chrA][chrB][candidate]["positions_B"]["splits"])


				else:
					reverse_A = candidates[chrA][chrB][candidate]["positions_A"]["orientation_discordants"].count("True")
					forward_A = candidates[chrA][chrB][candidate]["positions_A"]["orientation_discordants"].count("False")

					reverse_B = candidates[chrA][chrB][candidate]["positions_B"]["orientation_discordants"].count("True") 
					forward_B = candidates[chrA][chrB][candidate]["positions_B"]["orientation_discordants"].count("False")

					if ( reverse_A >= 5*forward_A or reverse_A*5 <= forward_A ) and ( reverse_B >= 5*forward_B or reverse_B*5 <= forward_B ):
						A_reverse=False
						if reverse_A > forward_A:
							A_reverse=True

						B_reverse=False
						if reverse_B > forward_B:
							B_reverse=True

						if is_mp:
							if A_reverse and not B_reverse:
								candidates[chrA][chrB][candidate]["posA"]=max(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=min(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

							elif not A_reverse and B_reverse:
								candidates[chrA][chrB][candidate]["posA"]=min(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=max(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

							elif A_reverse and B_reverse:
								candidates[chrA][chrB][candidate]["posA"]=max(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=max(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

							else:
								candidates[chrA][chrB][candidate]["posA"]=min(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=min(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

						else:
							if not A_reverse and B_reverse:
								candidates[chrA][chrB][candidate]["posA"]=max(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=min(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

							elif A_reverse and not B_reverse:
								candidates[chrA][chrB][candidate]["posA"]=min(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=max(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

							elif not A_reverse and not B_reverse:
								candidates[chrA][chrB][candidate]["posA"]=max(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=max(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

							else:
								candidates[chrA][chrB][candidate]["posA"]=min(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
								candidates[chrA][chrB][candidate]["posB"]=min(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

					else:
						candidates[chrA][chrB][candidate]["posA"]=mode(candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
						candidates[chrA][chrB][candidate]["posB"]=mode(candidates[chrA][chrB][candidate]["positions_B"]["discordants"])



					candidates[chrA][chrB][candidate]["posA"]
					candidates[chrA][chrB][candidate]["posB"]

				candidates[chrA][chrB][candidate]["startB"]=min(candidates[chrA][chrB][candidate]["positions_B"]["contigs"]+candidates[chrA][chrB][candidate]["positions_B"]["splits"]+candidates[chrA][chrB][candidate]["positions_B"]["discordants"])
				candidates[chrA][chrB][candidate]["endB"]=max(candidates[chrA][chrB][candidate]["positions_B"]["contigs"]+candidates[chrA][chrB][candidate]["positions_B"]["splits"]+candidates[chrA][chrB][candidate]["positions_B"]["discordants"])

				candidates[chrA][chrB][candidate]["startA"]=min(candidates[chrA][chrB][candidate]["positions_A"]["contigs"]+candidates[chrA][chrB][candidate]["positions_A"]["splits"]+candidates[chrA][chrB][candidate]["positions_A"]["discordants"])
				candidates[chrA][chrB][candidate]["endA"]=max(candidates[chrA][chrB][candidate]["positions_A"]["contigs"]+candidates[chrA][chrB][candidate]["positions_A"]["splits"]+candidates[chrA][chrB][candidate]["positions_A"]["discordants"])

	return(candidates)

#chromosomes=["1","2","3","4","5"]
#samples=["SweGen0001"]
#prefix=sys.argv[1]
#hej=main(prefix,chromosomes,samples,False,150,2)
##print(hej["3"]["5"])
#for entry in hej["3"]["5"]:
#	print(hej["3"]["5"][entry])
#
