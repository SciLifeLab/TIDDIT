import time
import sys
import graphlib
import copy

def build_kmer_hist(seq,kmer_hist,k):
	seq_len=len(seq)
	kmers=[]
	for i in range(0,seq_len-k+1):

		kmer=seq[i:i+k]
		if len(kmer) < k:
			break
		kmers.append(kmer)

		if not kmer in kmer_hist:
			kmer_hist[kmer]=0

		kmer_hist[kmer]+=1

	return(kmers,kmer_hist)	

def find_chain(graph,start,ends):
	chain=[start]
	current_node=start
	if start in ends:
		return(chain)
		
	while True:
		current_node=graph.sucessors[current_node]

		if not current_node or len(current_node) > 1 or current_node == start or current_node in ends:
			return(chain)
		current_node=list(current_node)[0]
		chain.append(current_node)
		if current_node in ends:
			return(chain)

def drop_kmers(graph,min_support):
	kmers=list(graph.kmers.keys())
	for kmer in kmers:
		if len(graph.kmers[kmer]) < min_support:
			graph.delete_kmer(kmer)
	return(graph)

def trim_edges(graph,min_weight):
	edge_list=list(graph.vertice_set)
	for edge in edge_list:
		if len(graph.vertices[edge[0]][edge[1]]) < min_weight:
			graph.delete_vertice(edge[0],edge[1])
	return(graph)

def remove_tips(graph,min_tip_length):
	branch_start=graph.in_branch_points
	branch_end=graph.out_branch_points
	starting_point=graph.starting_points

	switches=branch_end.union(branch_start)
	for start in starting_point.union(branch_start,branch_end):	
		chains=[]
		for node in graph.sucessors[start]:
			chains.append([start]+find_chain(graph,node,switches))

		for chain in chains:
			if len(chain) < 20 and chain[-1] in graph.end_points:
				for node in chain:
					graph.delete_kmer(node)

	return(graph)


def chain_typer(chain,graph):

	if chain[0] in graph.starting_points:
		return("starting_point")
	elif chain[-1] in graph.end_points:
		return("end_point")

	elif chain[0] in graph.in_branch_points:
		if chain[-1] in graph.out_branch_points:
			return("in_out")
		elif chain[-1] in graph.in_branch_points:
			return("in_in")

	elif chain[0] in graph.out_branch_points:
		if chain[-1] in graph.out_branch_points:
			return("out_out")

		elif chain[-1] in graph.in_branch_points:
			return("out_in")

	return("unknown")

def forward_scaffold(scaffold,chains,graph,chain_numbers):
	results=[]
	
	for i in range(0,len(chains)):
		if i in chain_numbers:
			continue

		if chains[i][0][0] == chains[scaffold][0][-1]:
			r=forward_scaffold(i,chains,graph,set([i]) | chain_numbers )

			for j in range(0,len(r)):
				results.append( [ chains[scaffold][0]+r[j][0][1:],r[j][1] | set([scaffold]) ] )
			
	if not results:
		results=[ [chains[scaffold][0], set([scaffold]) ] ]

	return(results)

def backward_scaffold(scaffold,chains,graph,chain_numbers):
	results=[]
	
	for i in range(0,len(chains)):
		if i in chain_numbers:
			continue

		if chains[i][0][-1] == chains[scaffold][0][0]:
			r=backward_scaffold(i,chains,graph,set([i]) | chain_numbers )

			for j in range(0,len(r)):
				results.append( [ r[j][0]+chains[scaffold][0][1:],r[j][1] | set([scaffold]) ] )
			
	if not results:
		results=[ [chains[scaffold][0], set([scaffold]) ] ]

	return(results)

def main(reads,k,min_support):

	time_all=time.time()

	kmers={}
	time_kmerize=time.time()
	graph = graphlib.graph()

	kmer_hist={}
	for read in reads:
		if len(reads[read]) < k:
			continue

		read_kmers,kmer_hist=build_kmer_hist(reads[read],kmer_hist,k)
		kmers[read]=read_kmers

	for read in kmers:
		if len(reads[read]) < k+1:
			continue


		for i in range(1,len(kmers[read])):

			if kmer_hist[kmers[read][i-1]] < min_support and kmer_hist[kmers[read][i]] < min_support:
				continue

			if kmer_hist[kmers[read][i]] < min_support and kmer_hist[kmers[read][i-1]] >= min_support:
				graph.add_kmer(kmers[read][i-1],read)

			elif kmer_hist[kmers[read][i]] >= min_support and kmer_hist[kmers[read][i-1]] < min_support:
				graph.add_kmer(kmers[read][i],read)

			if kmer_hist[kmers[read][i]] >= min_support and kmer_hist[kmers[read][i]] >= min_support:
				graph.add_vertice(kmers[read][i-1],kmers[read][i],read)		

	if not reads:
		return([])
	
	graph=drop_kmers(graph,min_support)
	graph=trim_edges(graph,min_support)
	graph=remove_tips(graph,10)

	branch_start=graph.in_branch_points
	branch_end=graph.out_branch_points
	starting_point=graph.starting_points

	chains=[]
	switches=branch_end.union(branch_start)
	for start in starting_point.union(branch_start,branch_end):	
		for node in graph.sucessors[start]:

			chain=[start]+find_chain(graph,node,switches)
			chain_type=chain_typer(chain,graph)
			chains.append([chain,chain_type])

	scaffolds=[]
	for i in range(0,len(chains)):
		chain=chains[i][0]
		start=chain[0]
		end=chain[-1]

		scaffold=[]	

		if chains[i][1] == "end_point":
			results=backward_scaffold(i,chains,graph,set([i]) )
			scaffolds+=results

		elif chains[i] == "start_point":
			results=forward_scaffold(i,chains,graph,set([i]) )
			scaffolds+=results
		else:
			forward=forward_scaffold(i,chains,graph,set([i]) )
			for forward_result in forward:
				backward_result=backward_scaffold(i,chains,graph,forward_result[1] )
				for result in backward_result:
					scaffolds.append([  result[0]+forward_result[0][len(chains[i][0])-1:],forward_result[1] | result[1]])

	results=[]
	for i in range(0,len(scaffolds)):
	
		skip=False
		for j in range(0,len(scaffolds)):
			if j ==i or j < i:
				continue
			if not len(scaffolds[i][-1]-scaffolds[j][-1]):
				skip=True

		if skip:
			continue

		scaffolded_chains=list(map(str,scaffolds[i][-1]))

		out=[]
		for j in range(1,len(scaffolds[i][0])):
			out.append(scaffolds[i][0][j][-1])

		results.append(scaffolds[i][0][0]+"".join(out))
	return(results)

min_branch_length=2
min_overlap=0.2
max_overlap=0.8
