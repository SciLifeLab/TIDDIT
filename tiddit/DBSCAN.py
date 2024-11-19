import numpy
#The functions of this script perform DBSCAN and returns the variatn clusters

#Call the cluser algorithm, the statistics function, and returns the final cluster
def generate_clusters(chrA,chrB,coordinates,library_stats,args):
	candidates=[]
	coordinates=coordinates[numpy.lexsort((coordinates[:,1],coordinates[:,0]))]
	min_pts=args.l
	if chrA == chrB and library_stats["ploidies"][chrA] > args.n*2:
		min_pts=int(round(args.l/float(args.n)*library_stats["ploidies"][chrA]))

	db=main(coordinates[:,0:2],args.e,min_pts)
	unique_labels = set(db)

	for var in unique_labels:
		if var == -1:
			continue
		class_member_mask = (db == var)
		candidate_signals =coordinates[class_member_mask]
		resolution=candidate_signals[:,-2]
		support=len(set(candidate_signals[:,-1]))
		discordants=True
		if len(set(resolution)) == 1 and max(resolution) == 1:
			disordants=False

		if discordants and support >= args.p:
			candidates.append( analyse_pos(candidate_signals,discordants,library_stats,args) )
		elif not discordants and support >= args.r and chrA == chrB:
			candidates.append( analyse_pos(candidate_signals,discordants,library_stats,args) )

	return(candidates)

def x_coordinate_clustering(data,epsilon,m):
	clusters=numpy.zeros(len(data))
	for i in range(0,len(clusters)):
		clusters[i]=-1
	cluster_id=-1
	cluster=False

	for i in range(0,len(data)-m+1):

		distances=[]
		current=data[i,:]
		points=data[i+1:i+m+1,:]

		#print points
		distances=[]
		for point in points:
			distances.append(abs(point[0]-current[0]))

		if max(distances) < epsilon:
			#add to the cluster
			if cluster:
				clusters[i+m-1]=cluster_id
				#define a new cluster
			else:
				cluster_id+=1
				cluster=True
				for j in range(i,i+m):
					clusters[j]=cluster_id
		else:
			cluster=False

	return(clusters,cluster_id)

def y_coordinate_clustering(data,epsilon,m,cluster_id,clusters):

	cluster_id_list=set(clusters)
	for cluster in cluster_id_list:
		if cluster == -1:
			continue
		class_member_mask = (clusters == cluster)
		indexes=numpy.where(class_member_mask)[0]
		signals=data[class_member_mask]

		y_coordinates=[]


		for i  in range(0,len(signals)):
			y_coordinates.append([signals[i][1],indexes[i]])
		y_coordinates.sort(key=lambda x:x[0])
		
		sub_clusters=numpy.zeros(len(indexes))
		for i in range(0,len(sub_clusters)):
			sub_clusters[i]=-1

		active_cluster=False
		sub_cluster_id=0
		y_coordinates=numpy.array(y_coordinates)
		for i in range(0,len(y_coordinates)-m+1):
			distances=[]
			current=y_coordinates[i,:]
			next=y_coordinates[i+1:i+m,:]

			distances=[]
			for pos in next:
				distances.append(abs(pos[0]-current[0]))	

			if max(distances) < epsilon:
				#add to the cluster
				if active_cluster:
					sub_clusters[i+m-1]=sub_cluster_id
					#define a new cluster
				else:
					sub_cluster_id+=1
					active_cluster=True
					for j in range(i,i+m):
						sub_clusters[j]=sub_cluster_id
			else:
				active_cluster=False

		for i in range(0,len(sub_clusters)):
			if sub_clusters[i] == 1:
				clusters[ y_coordinates[i][1] ]	= cluster

			elif sub_clusters[i] > -1:
				clusters[ y_coordinates[i][1] ]=sub_clusters[i] +cluster_id-1
			elif sub_clusters[i] == -1:
				clusters[ y_coordinates[i][1] ] = -1

		if sub_cluster_id > 1:
			cluster_id += sub_cluster_id-1
	return(clusters,cluster_id)

def main(data,epsilon,m):
	clusters,cluster_id=x_coordinate_clustering(data,epsilon,m)
	clusters,cluster_id=y_coordinate_clustering(data,epsilon,m,cluster_id,clusters)

	return(clusters)


#hej=numpy.array([[1,2],[1,2],[1,2],[10,11]])
#print(main(hej,0.1,2))
#print(hej)
