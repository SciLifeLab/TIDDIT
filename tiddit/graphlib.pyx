class graph:

	def __init__(self):
		self.predecessors={}
		self.sucessors={}
		self.kmers={}
		self.vertices={}
		self.vertice_set=set([])
		self.in_branch_points=set([])
		self.out_branch_points=set([])
		self.starting_points=set([])
		self.end_points=set([])

	#add node to the graph
	def add_kmer(self,kmer,read):

		if not kmer in self.kmers:
			self.kmers[kmer]=set([])
			self.predecessors[kmer]=set([])
			self.sucessors[kmer]=set([])
			self.starting_points.add(kmer)
			self.end_points.add(kmer)

		self.kmers[kmer].add(read)
		
	#add vertices between nodes
	def add_vertice(self,kmer1,kmer2,read):

		self.add_kmer(kmer1,read)
		self.add_kmer(kmer2,read)
		
		if not kmer1 in self.vertices:
			self.vertices[kmer1]={}

		if kmer1 in self.end_points:
			self.end_points.remove(kmer1)

		if not kmer2 in self.vertices[kmer1]:
			self.vertices[kmer1][kmer2]=set([])

		if kmer2 in self.starting_points:
			self.starting_points.remove(kmer2)

		self.vertices[kmer1][kmer2].add(read)
			
		self.vertice_set.add((kmer1,kmer2))

		self.sucessors[kmer1].add(kmer2)
		if len(self.sucessors[kmer1]) > 1:
			self.in_branch_points.add(kmer1)

		self.predecessors[kmer2].add(kmer1)
		if len(self.predecessors[kmer2]) > 1:
			self.out_branch_points.add(kmer2)

	def delete_vertice(self,kmer1,kmer2):
	
		if kmer1 in self.vertices:
			if kmer2 in self.vertices[kmer1]:
				del self.vertices[kmer1][kmer2]

			if len(self.vertices[kmer1]) == 0:
				del self.vertices[kmer1]

		if kmer1 in self.in_branch_points:
			if len(self.sucessors[kmer1]) < 3:
				self.in_branch_points.remove(kmer1)

		if kmer1 in self.sucessors:
			if kmer2 in self.sucessors[kmer1]:
				self.sucessors[kmer1].remove(kmer2)

			if not self.sucessors[kmer1]:
				self.end_points.add(kmer1)

		if kmer2 in self.out_branch_points:
			if len(self.predecessors[kmer2]) < 3:
				self.out_branch_points.remove(kmer2)

		if kmer1 in self.predecessors[kmer2]:
			self.predecessors[kmer2].remove(kmer1)

		if not len(self.predecessors[kmer2]):
			self.starting_points.add(kmer2)

		if (kmer1,kmer2) in self.vertice_set:
			self.vertice_set.remove((kmer1,kmer2))

	def delete_kmer(self,kmer):
		if kmer in self.kmers:
			del self.kmers[kmer]

		sucessors=set([])
		for sucessor in self.sucessors[kmer]:
			sucessors.add(sucessor)

		predecessors=set([])
		for predecessor in self.predecessors[kmer]:
			predecessors.add(predecessor)

		for predecessor in predecessors:
			self.delete_vertice(predecessor,kmer)

		for sucessor in sucessors:
			self.delete_vertice(kmer,sucessor)
