#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''

import os, re, dill
import numpy as np
from scipy.sparse import coo_matrix
#from sklearn.preprocessing import normalize
from gcn.lib.utils import lib_utils
from gcn.lib.databases import kegg_pathway
import divine_inc

def extract_ensembl_protein(protein):
	mObj=re.search(r'\d+\.(\w+)',protein)
	return mObj.group(1)

def cal_array_distance(npVec1, npVec2):
	return sum(abs(npVec1-npVec2))[0]

class HeatDiff:
	def __init__(self,cDivine):
		self.dv = cDivine
		self.min_edge_weight = 350
		self.nNodes = 0
		self.Prots = []
		self.dProt2Idx = {}
		self.dPPI = {}
		self.dGene2Prot = {}
		self.dProt2Gene = {}
		self.dangledGenes = []
		self.harmonic_sc = []
		self.harmonic_dng_sc = []
		self.ppi = [[], [], []]  # from protein, to protein, link weight
		self.Y = None
		self.Y0 = None

	def get_sparse_elements(self):
		'''
		to store ppi network
		input: self.dProt2Gene, dGenes(whether the gene is in ppi or not)- protein-gene relation; proteinLinkFile- ppi link
		output: update self.dProt2Gene, dGenes when add_dangled is enabled. Store ppi and Prots
		'''
		#read string DB and assign an integer to each protein symbol
		fp = lib_utils.open2(self.dv.entries['string_link'],'r')

		linked = [-1,-1]

		self.nNodes = 0
		self.Prots = []
		self.dProt2Idx = {}

		lib_utils.msgout('notice', 'preparing a genetic network matrix. Please, be patient ...', 'pagerank|heat_diffusion')
		#store col,row,weight from ppi file

		fp.next()
		for i in fp:
			#print '%s'%i #debug
			linked[0],linked[1],weight=i.rstrip().split()
			weight = float(weight)
			if weight < self.min_edge_weight: continue

			for c in range(2):
				protein=extract_ensembl_protein(linked[c])

				#to register a protein node
				if not protein in self.dProt2Idx:
					self.dProt2Idx[protein]=self.nNodes

					# item index corresponds to a node number of the protein
					self.Prots.append(protein)
					self.nNodes+=1

				self.ppi[c].append(self.dProt2Idx[protein])
			self.ppi[2].append(weight)
		fp.close()

	def gen_adj_matrix(self,logger=None,reuse=True):

		if not self.ppi:
			raise RuntimeError('edge info is not available yet. run get_sparse_elements() first to load ppi edge info ...')

		dill_fn = self.dv.entries['string_link']+'.dill'
		if reuse and os.path.exists(dill_fn):
			msg = "loading adjacent matrix computed previously and stored in [%s]"%dill_fn
			lib_utils.msgout('notice',msg)
			if logger: logger.info(msg)

			with open(dill_fn, 'rb') as in_strm:
				self.A = dill.load(in_strm)
		else:
			self.A = coo_matrix((self.ppi[2], (self.ppi[0], self.ppi[1])), \
										 dtype=np.float, shape=(self.nNodes, self.nNodes))

			job_name = 'gen_adj_matrix'
			# convert to csr_matrix for faster/reliable matrix operation
			msg = 'reformatting the genetic network matrix.'
			lib_utils.msgout('notice', msg, job_name)
			if logger: logger.info(msg)
			self.A = self.A.tocsr()

			# normalize PPI matrix
			msg = 'normalizing (graph laplacian) the genetic network matrix. (it will take 4 hours!)'
			lib_utils.msgout('notice', msg, job_name)
			if logger: logger.info(msg)

			self.A = normalize_glap(self.A)
			#self.A = normalize(self.A, norm='l1', axis=0)
			with open(dill_fn, 'wb') as out_strm:
				dill.dump(self.A, out_strm)

	def conv_ppi_to_ragged_array(self):
		'''
		to convert ppi[[protein1_idx],[protein2_idx],[wt]]
		into {protein1_idx:[protein2_idx,...]} dic format
		'''
		self.dPPI = {}
		N = len(self.ppi[0])

		for n in range(N):
			p1,p2,wt = [self.ppi[0][n],self.ppi[1][n],self.ppi[2][n]]
			if not p1 in self.dPPI:
				self.dPPI[p1]=[]
			self.dPPI[p1].append(p2)

			if not p2 in self.dPPI:
				self.dPPI[p2]=[]
			self.dPPI[p2].append(p1)
		self.ppi = None

	def protein_to_gene(self):
		'''
		to create a map from protein to ref gene symbol
		'''
		dProt2Gene = {}

		fp = lib_utils.open2(self.dv.entries['esp_to_gene'],'r')
		fp.next()
		for i in fp:
			gene,protein = i[:-1].split('\t')
			if gene and protein:
				if protein not in dProt2Gene:
					dProt2Gene[protein]=gene

		fp.close()
		return dProt2Gene

	def gene_to_protein(self):
		'''
		to create a map from ref gene symbol to protein
		'''
		dGene2Prot = {}

		fp = lib_utils.open2(self.dv.entries['esp_to_gene'],'r')
		fp.next()
		for i in fp:
			gene,protein = i[:-1].split('\t')
			if gene and protein:
				if gene not in dGene2Prot:
					dGene2Prot[gene]=protein

		fp.close()
		return dGene2Prot

	def assign_node_prior(self,verbose=False):

		self.Y = np.zeros(shape=(self.nNodes,1),dtype=float)
		
		#to get the genes in a dic whose dmg_sc > 0.
		perturbedGeneOnLnk = lib_utils.list_to_dic(\
													list(self.dv.gene_dmg.keys()),False)
		
		#assign the combined dmg score to PPI
		for n,protein in enumerate(self.Prots): #scanning all ppi
			if protein in self.dProt2Gene:
				gene=self.dProt2Gene[protein] #protein-> gene
				if gene in self.dv.gene_dmg:
					#to assign combined dmg score to node
					self.Y[n] = self.dv.gene_dmg[gene][0]
					
					#to remember this gene is assigned
					perturbedGeneOnLnk[gene] = True
			elif verbose:
				print 'protein[%s] is not found from a gene-product mapping database!'%protein

		#count the number of perturbed genes not on the network (i.e. dangled nodes)
		nNodes0 = len([1 for isLnk in perturbedGeneOnLnk.values() if isLnk is False])
		self.Y0 = np.zeros(shape=(nNodes0,1),dtype=float)

		#now, remember which dmg gene is not assigned 
		self.dangledGenes = []
		n = 0
		for gene,isLnk in perturbedGeneOnLnk.iteritems():
			#print gene
			if not isLnk:
				self.dangledGenes.append(gene)
				self.Y0[n] = self.dv.gene_dmg[gene][0]
				n += 1
		
		denom = np.sum(self.Y) + np.sum(self.Y0)

		#normalize in a range of 0 to 1.
		self.Y /= denom
		self.Y0 /= denom

	def heat_diffusion_core(self,gamma=2.,M=100,alpha=0.9,\
													maxIter=150,logger=None):
		job_name = 'pagerank|heat_diffusion_core'

		N=len(self.Y)
		s = np.zeros(shape=(N,1))

		N0=len(self.Y0)
		s0 = np.zeros(shape=(N0,1))

		epsilon = 1e-4
		iter = 0

		msg = 'running heat diffusion on [%dx%d, gamma=%g, alpha=%g, max_iter=%d, M=%d]. Please, be patient ...' % (N,N,gamma,alpha,maxIter,M)
		lib_utils.msgout('notice', msg, job_name)
		if logger: logger.info(msg)

		e = 1.
		while (e>=epsilon and iter<maxIter):
			#heat diffusion
			s_new  = (1.-gamma/M)*s  + (gamma/M)*(alpha*self.A.dot(s)+(1.-alpha)*self.Y)
			s0_new = (1.-gamma/M)*s0 + (gamma/M)*(1.-alpha)*self.Y0

			#normalize
			denom = np.sum(s_new) + np.sum(s0_new)
			s_new = s_new/denom
			s0_new = s0_new/denom

			e = cal_array_distance(s_new,s)+cal_array_distance(s0_new,s0)

			s = np.copy(s_new)
			s0 = np.copy(s0_new)

			iter+=1

		msg = 'done. [iteration:%d/%d,e:%g]'%(iter,maxIter,e)
		lib_utils.msgout('notice', msg, job_name)
		if logger: logger.info(msg)
		
		return s,s0

	def convert_node2gene(self):
		'''
		for only gt_dmg genes, print out gene, harmonic score, and seed score 
		'''
		
		rank_fn_tmp = '%s.tmp'%self.dv.gene_rank_fn
		fp2=lib_utils.open2(rank_fn_tmp,'w')
		fp2.write('#gene\tpredicted_score\tseed_score\tgt_dmg_score\tpheno_score\tcontain_known_pathogenic\n')
		genes_printed = {}
		#browsing each node in the whole (original) ppi network
		for n,protein in enumerate(self.Prots):
			seed_score = 0.
			gene = protein
			
			#check if this node (restart value) was assigned previously
			if protein in self.dProt2Gene:
				gene = self.dProt2Gene[protein]

				if gene in self.dv.gene_dmg:
					seed_score = self.dv.gene_dmg[gene][0]
					
			#to get harmonic score and save into dv.gene_dmg
			pred_score = 0.
			if self.harmonic_sc[n][0]>0.:
				pred_score = self.harmonic_sc[n][0]
				if gene in self.dv.gene_dmg:
					self.dv.gene_dmg[gene][1] = pred_score
			
			#NOTE that print only a gene having at one mutation
			if (not self.dv.gt_dmg) or \
					(gene in self.dv.gt_dmg and self.dv.gt_dmg[gene].score>0.):
				
				pheno_sc = 0.
				if gene in self.dv.pheno_dmg:
					pheno_sc = self.dv.pheno_dmg[gene].score
				
				if self.dv.vknown:
					if gene in self.dv.vknown_genes: is_vknown = 'Y'
					else: is_vknown = 'N'
				else: is_vknown = 'NA'

				if gene in genes_printed:
					gene2 = '%s|%s'%(gene,protein)
				else:
					gene2 = gene
					genes_printed[gene]=True

				fp2.write('%s\t%g\t%g\t%g\t%g\t%s\n'%\
								(gene2,pred_score,seed_score,\
								self.dv.gt_dmg[gene].score,pheno_sc,is_vknown))

		#repeat the same procedure to dangled nodes
		for n,gene in enumerate(self.dangledGenes):
			
			self.dv.gene_dmg[gene][1] = self.harmonic_dng_sc[n][0]

			if (not self.dv.gt_dmg) or \
					(gene in self.dv.gt_dmg and self.dv.gt_dmg[gene].score>0.):
				
				pheno_sc = 0.
				if gene in self.dv.pheno_dmg:
					pheno_sc = self.dv.pheno_dmg[gene].score
				
				if self.dv.vknown:
					if gene in self.dv.vknown_genes: is_vknown = 'Y'
					else: is_vknown = 'N'
				else: is_vknown = 'NA'
					
				fp2.write('%s\t%g\t%g\t%s\t%g\t%s\n'%\
					(gene,self.dv.gene_dmg[gene][1],self.dv.gene_dmg[gene][0],\
					self.dv.gt_dmg[gene].score,pheno_sc,is_vknown))

		fp2.close()

		#sort by score
		lib_utils.sort_tsv_by_col2(\
			rank_fn_tmp, [2], ['gr'], False, self.dv.gene_rank_fn)
		os.unlink(rank_fn_tmp)
		
	def create_disease_rank_tab(self):
		
		fpw = open(self.dv.disease_rank_fn,'w')

		headStr = """
			disease_ID
			disease_description
			inheritance
			assoc_pheno_genes(^:mutated,*:known_pathogenic)
			num_of_assoc_pheno_genes
			num_of_gt_dmg_genes
			pheno_match_score
			avg_combined_dmg_score
			max_combined_dmg_score
			avg_harmonic_score
			max_harmonic_score
			external_genes_of_interest(kegg-ppi_or_GO_enriched[harmonic_score])
			PPI-KEGG_pathway_desc
			"""

		headCols = headStr.split()

		cell_delim = ';'
		fpw.write('#%s\n'%lib_utils.joined(headCols,'\t'))
		
		#
		cKegg = kegg_pathway.Kegg(hsa_fn=self.dv.entries['kegg_hsa'])
		cKegg.get_hsa()
		
		#annotate kegg_pathway to disease
		self.dv.omim.to_kegg_hsa(cKegg.cHsa)
		
		#browsing whole known disease entries whose HPO sim score with the patient > 0.
		for cD in self.dv.omim.cDis.itervalues():
			if cD.pheno_score == 0.: continue

			to_print = []
			to_print.append(cD.id) #disID
			to_print.append(cD.desc) #disDesc
			to_print.append(divine_inc.inheritStr[cD.inherit]) #inheritance

			Genes = [[],[]]
			max_rw_score = [0., 0.]
			sum_rw_score = [0., 0.]
			cnt_gene_dmg = 0
			gt2_dmg = None
			
			for gene in cD.genes:#for each gene assoc with the disease
				#split into two groups (one having gt_dmg, or else), and collect max & sum act score
				if gene in self.dv.gt_dmg:
					if self.dv.vknown and gene in self.dv.vknown_genes:
						Genes[0].append('%s*'%gene)
					else:
						Genes[0].append('%s^'%gene)
						
					if self.dv.gene_dmg[gene][0] > max_rw_score[0]:
						max_rw_score[0] = self.dv.gene_dmg[gene][0]
					sum_rw_score[0] += self.dv.gene_dmg[gene][0]
				else:
					Genes[1].append(gene)
				
				#to collect max & sum on harmonic scores
				if gene in self.dv.gene_dmg:
					if self.dv.gene_dmg[gene][1] > max_rw_score[1]:
						max_rw_score[1] = self.dv.gene_dmg[gene][1]
					sum_rw_score[1] += self.dv.gene_dmg[gene][1]
					cnt_gene_dmg += 1
			
			#bring KEGG genes (PPI) interacted with non-mutated phenotype genes
			goi,hsa_desc = self.external_goi(\
							Genes[1],cD.kegg_hsa,cKegg.cHsa)
			
			#bring GO enriched genes
			for gene2 in cD.enriched_genes:
				geneStr2 = 'go(%s:%s' % (cD.enriched_genes[gene2],gene2)
				if self.dv.vknown and (gene2 in self.dv.vknown_genes):
					geneStr2 = geneStr2 + '*'
				else:
					geneStr2 = geneStr2 + '^'
				
				if gene2 in self.dv.gene_dmg:
					goi.append('%s[%g])' % \
						(geneStr2,self.dv.gene_dmg[gene2][1]))

			to_print.append(cell_delim.join(Genes[0]+Genes[1])) #assoc_pheno_genes
			G = len(cD.genes)
			G_mt = len(Genes[0])
			to_print.append(G) #num_of_assoc_pheno_genes
			to_print.append(G_mt) #num_of_gt_dmg_genes
			to_print.append(cD.pheno_score) #pheno_match_score
			if G_mt > 0:
				to_print.append(sum_rw_score[0]/G_mt) #avg_combined_dmg_score
			else:
				to_print.append(0.)
			to_print.append(max_rw_score[0]) #max_combined_dmg_score
			
			if cnt_gene_dmg > 0:
				to_print.append(sum_rw_score[1]/cnt_gene_dmg) #avg_harmonic_score
				to_print.append(max_rw_score[1]) #max_harmonic_score
			else:
				to_print.append(0.)
				to_print.append(0.)
			to_print.append(cell_delim.join(goi)) #partner_in_protein_network_of_interest
			if hsa_desc:
				to_print.append(cell_delim.join(hsa_desc)) #kegg-pathway desc if exist
			else:
				to_print.append('NA') #kegg-pathway desc if exist
			fpw.write('%s\n'%(lib_utils.joined(to_print,'\t')))

		fpw.close()

	def external_goi(self,genes,hsas,cHsa):

		genes2 = []
		hsa_descs = []
		for gene in genes:
			#gene1 -> protein1 -> pIdx1
			if gene in self.dGene2Prot:
				prot = self.dGene2Prot[gene]
	
				if prot in self.dProt2Idx:
	
					for pIdx2 in self.dPPI[self.dProt2Idx[prot]]:
	
						#pIdx2 -> protein2 -> gene2
						if self.Prots[pIdx2] in self.dProt2Gene:
							gene2 = self.dProt2Gene[self.Prots[pIdx2]]
						else:
							gene2 = self.Prots[pIdx2]
	
						if gene2 in self.dv.gt_dmg:
							for hsa in hsas:
								if gene2 in cHsa[hsa].genes:
									geneStr2 = 'kegg_ppi(%s:%s' % (gene,gene2)
									if self.dv.vknown and (gene2 in self.dv.vknown_genes):
										geneStr2 = geneStr2 +'*'
									else:
										geneStr2 = geneStr2 +'^'
										
									genes2.append('%s[%g])' % \
										(geneStr2,self.dv.gene_dmg[gene2][1]))
									
									hsa_descs.append(cHsa[hsa].desc)

		return list(set(genes2)),list(set(hsa_descs))

def normalize_glap(A):

	dimR,dimC = A.get_shape()
	#to get column sum
	col_sum_sqrt = np.sqrt(A.sum(axis=0).transpose())

	#to get row sum
	row_sum_sqrt = np.sqrt(A.sum(axis=1))

	for i in range(dimR):
		for j in range(dimC):
			if A[i,j]!= 0.:
				A[i,j] = A[i,j]/(row_sum_sqrt[i]*col_sum_sqrt[j])

	return A

def run_heatdiffusion(cDivine,logger):

	cRW = HeatDiff(cDivine)

	job_name = "pagerank.run_heatdiffusion"

	msg = 'transferring gene product to a matrix [%s;%s]'%\
		(cRW.dv.entries['string_link'],cRW.dv.entries['esp_to_gene'])
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	
	#store original ppi file into self.Prots, self.dProt2Idx, self.nNodes
	cRW.get_sparse_elements()

	msg = 'creating a genetic network matrix[%d x %d].'%\
		(cRW.nNodes,cRW.nNodes)
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)

	cRW.gen_adj_matrix()

	#convert ppi to ragged array
	cRW.conv_ppi_to_ragged_array()

	cRW.dProt2Gene = cRW.protein_to_gene()
	cRW.dGene2Prot = cRW.gene_to_protein()

	#to get a map between protein and gene
	msg = 'prepping mapping btn gene and protein.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)

	msg = 'assign an activation value from cDivine into the protein network.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	cRW.assign_node_prior() #self.dangledGene

	#core heat diffusion in recursion
	msg = 'running heat diffusion on genetic networks labeled by perturbed genes.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	
	cRW.harmonic_sc,cRW.harmonic_dng_sc = \
		cRW.heat_diffusion_core(logger=logger)

	cRW.A = None; cRW.Y = None; cRW.Y0 = None

	#annotate gene to node
	msg = 'reporting ranked genes.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	cRW.convert_node2gene()

	#annotation disease ranking
	if cDivine.hpo_query:
		msg = 'reporting a disease rank.'
		lib_utils.msgout('notice',msg,job_name);logger.info(msg)
		
		cRW.create_disease_rank_tab()

	msg = 'done.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
