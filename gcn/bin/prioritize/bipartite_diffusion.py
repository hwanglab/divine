#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''

import numpy as np
import os, dill
from scipy.sparse import coo_matrix
from gcn.lib.utils.lib_utils import py_struct, swap_keyVal, normalize_dic

def build_bipartitie_graph(kegg_genes_fn):
	fp = open(kegg_genes_fn, 'r')
	heads = fp.next()[1:-1].split('\t')
	hsa_to_genes = [[],[],[]]
	hsa_idx = 0
	gene_idx = 0
	Hsa = {}
	Gene = {}
	for i in fp:
		hsa, desc, genes = i[:-1].split('\t')
		if hsa not in Hsa:
			Hsa[hsa] = hsa_idx
			hsa_idx += 1
		for gene in genes.split(','):
			gene = gene.strip()
			if gene:
				if gene not in Gene:
					Gene[gene] = gene_idx
					gene_idx += 1
				hsa_to_genes[0].append(Gene[gene])
				hsa_to_genes[1].append(Hsa[hsa])
				hsa_to_genes[2].append(1.)
	fp.close()

	return py_struct(heads = heads,
									 rows = Gene,
									 cols = Hsa,
									 rcw = hsa_to_genes,
									 prior = {})

class bipartitie:
	def __init__(self,pyEdge):
		'''
		:param pyEdge:
		.rcw[0]: row index
		.rcw[1]: column index
		.rcw[2]: edge value
		.rows: dict{node_name:row_index)
		.cols: dict{node_name:col_index)
		.initial_values
		'''
		self.A = None
		self.dimR = 0 #column
		self.dimC = 0 #row
		self.priorR = None
		self.load_graph(pyEdge)

	def load_graph(self,pyEdge):
		#build a sparse matrix
		self.A = coo_matrix((pyEdge.rcw[2], (pyEdge.rcw[0], pyEdge.rcw[1])), \
									 dtype=np.float, shape=(len(pyEdge.rows), len(pyEdge.cols)))

		self.A = self.A.tocsr()
		self.dimR, self.dimC = self.A.get_shape()
		self.normalize()
		self.pyEdge = pyEdge

	def normalize(self):

		print('normalize the matrix (graph laplcian) ...')
		#to get column sum
		col_sum_sqrt = np.sqrt(self.A.sum(axis=0).transpose())

		#to get row sum
		row_sum_sqrt = np.sqrt(self.A.sum(axis=1))

		for i in range(self.dimR):
			for j in range(self.dimC):
				if self.A[i,j]!= 0.:
					self.A[i,j] = self.A[i,j]/(row_sum_sqrt[i]*col_sum_sqrt[j])

		print('Done.')

	def set_prior(self, genes, prior_vals):
		print('loading prior/activation values into the bipartite graph ...')
		priorR = np.zeros(shape=(len(self.pyEdge.rows), 1), dtype=float)

		for gene, prior_val in zip(genes,prior_vals):
			if gene in self.pyEdge.rows:
				self.pyEdge.prior[gene] = prior_val
				priorR[self.pyEdge.rows[gene]] = prior_val

		if np.sum(priorR) > 0.:
			self.priorR = priorR / np.sum(priorR)
		else:
			print('[Warning]:all activiation values are 0!')
		del priorR
		print('Done')

	def net_diff(self,a=0.75,maxiter=50):
		Rv = self.priorR
		initialRv = np.copy(Rv)
		Cv = np.zeros(shape=(self.dimC, 1), dtype=float)
		if np.sum(self.priorR) > 0.:

			initialCv = np.copy(Cv)
			Atr = self.A.transpose()
			eps2 = 1e-7
			print("net diffusion start to run at epsilon[%g] or maxiter[%d]..."%(eps2,maxiter))
			for k in range(maxiter):
				oldRv = np.copy(Rv)
				oldCv = np.copy(Cv)
				Rv = a * self.A.dot(oldCv) + (1 - a) * initialRv
				Cv =  a * Atr.dot(oldRv) + (1 - a) * initialCv

				if (max(abs(Rv - oldRv))[0]<eps2 and max(abs(Cv - oldCv))[0]<eps2):
					break
			print("Done.")
		return Rv,Cv

	def print_rows(self,Rv,topK=0):
		print('collecting topK[%d] nodes in biparitite[rows]'%topK)
		ridx2key = swap_keyVal(self.pyEdge.rows)

		Rv2 = np.zeros(shape=(self.dimR, 2), dtype=float)
		Rv2[:, 0] = range(self.dimR) #store index
		Rv2[:, 1] = Rv.transpose()  # store row difussion scores from bipartitie
		Rv2 = Rv2[(-Rv2[:, 1]).argsort()] #descending order sort

		if topK > self.dimR or topK==0:
			topK = self.dimR

		outputs = []
		for k in range(topK):
			if Rv2[k,0] in ridx2key:
				outputs.append([int(Rv2[k,0]),ridx2key[Rv2[k,0]],Rv2[k,1]])
		print('Done.')
		return outputs

	def get_final_scores(self, Rv, topK=0, nodes2sel=[]):

		dNodeScores = {}
		
		if np.sum(Rv) > 0.:
			print('collecting topK[%d] nodes in biparitite[rows]'%topK)
			ridx2key = swap_keyVal(self.pyEdge.rows)

			Rv2 = np.zeros(shape=(self.dimR, 2), dtype=float)
			Rv2[:, 0] = range(self.dimR) #store index
			Rv2[:, 1] = Rv.transpose()  # store row difussion scores from bipartitie
			Rv2 = Rv2[(-Rv2[:, 1]).argsort()] #descending order sort

			if topK > self.dimR or topK==0:
				topK = self.dimR

			k = 0
			v = 0
			while(k < self.dimR and v < topK):
				if Rv2[k,0] in ridx2key:
					node = ridx2key[Rv2[k, 0]]
					if nodes2sel:
						if node in nodes2sel:
							dNodeScores[node] = Rv2[k,1]
							v += 1
					else:
						dNodeScores[node] = Rv2[k, 1]
						v += 1
				k += 1
			dNodeScores = normalize_dic(dNodeScores,mode='max')
			print('Done.')
		return dNodeScores

def run_bp(genes,prior_vals,nodes2sel,kegg_genes_fn):
	pyv = '%s.dill'%kegg_genes_fn
	if os.path.exists(pyv):
		fp = open(pyv, 'rb')
		cB,pyEdge = dill.load(fp)
	else:
		pyEdge = build_bipartitie_graph(kegg_genes_fn)
		cB = bipartitie(pyEdge)
		fp = open(pyv, 'wb')
		dill.dump([cB, pyEdge], fp)
	fp.close()

	cB.set_prior(genes,prior_vals)
	Rv,Cv = cB.net_diff()
	return cB.get_final_scores(Rv,topK=25,nodes2sel=nodes2sel)

if __name__ == '__main__':
	genes = ['HLCS','PRMT1','C3','C3AR1']
	prior_vals=[0.25,0.15,0.30,0.3]
	fn="/home/hongc2/mount/vm_shared/projects/divine-0.1.1/gcndata/kegg_pathway/kegg_genes.tsv"
	run_bp(genes,prior_vals,{},kegg_genes_fn=fn)