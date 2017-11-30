import numpy
#supports up to 4D data

def x_coordinate_clustering(data,epsilon,m):
	clusters=numpy.zeros(len(data))
	for i in range(0,len(clusters)):
		clusters[i]=-1
	cluster_id=-1
	cluster=False

	for i in range(0,len(data)-m+1):

		distances=[]
		current=data[i,:]
		points=data[i+1:i+m,:]
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


data=numpy.array([[1,1],[10481823,10483880],[10481947,10483785],[10481947,1],[10481947,1],[10481947,1],[10482033,10483984],[10482079,10483801],[10482111,10483972],[10482121,10483788],[10482125,10483769],[10482126,10484204],[10482163,10483811],[10482177,10483909],[10482186,10483906],[10482191,10483836],[10482202,10484150],[10482262,10483947],[10482285,10483797],[10482342,10483968],[10483770,10482390],[10482390,10483814],[10483770,10482405],[10483769,10482428],[10483770,10482405],[10483770,10482405],[10483770,1],[10483770,1],[10483770,1],[10483770,1]])


#data=data[numpy.lexsort((data[:,1],data[:,0]))]
#print data
#data=numpy.array([[1,2],[1,3],[1,4],[1,5],[5,100],[5,100],[10,2],[10,3],[10,4],[10,5],[20,100]])
#print main(data,398,3)

