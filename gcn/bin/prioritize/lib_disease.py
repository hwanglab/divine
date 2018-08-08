#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''
import pandas as pd
from gcn.lib.utils import lib_utils
import divine_inc

class Disease:
	
	def __init__(self,disId):
		self.id = disId
		self.desc = None
		self.pheno_score = 0.
		self.genes = [] #its value represents gt_dmg_score
		self.enriched_genes = {}
		self.inherit = divine_inc.INH_UNKNOWN
		self.hpos = []
		self.kegg_hsa = []
	
class Omim:
	
	def __init__(self,omim_disgenet_fn,desc_fn):
		self.disease_fn = omim_disgenet_fn
		self.desc_fn = desc_fn
		self.cDis = {}
	
	def to_kegg_hsa(self,cHsa,min_cvd_ratio=0.9):
		
		for disId,cD in self.cDis.iteritems():
			for hsaId,cH in cHsa.iteritems():
				shared_genes = lib_utils.intersect(cD.genes,cH.genes)
				G = len(cD.genes)
				if G > 0:
					cvd_ratio = 1.*len(shared_genes)/G
					if cvd_ratio > min_cvd_ratio:
						cD.kegg_hsa.append(hsaId)
			self.cDis[disId] = cD

	def get_disId_sorted_by_pheno_score(self,order='desc'):

		disIds = []
		pheno_scores = []
		for disId,cD in self.cDis.iteritems():
			disIds.append(disId)
			pheno_scores.append(cD.pheno_score)

		disPhenoScDf = pd.DataFrame({'idx': range(len(disIds)),
																 'disId': disIds,
																 'pheno_score': pheno_scores})

		disPhenoScDf.sort_values(by=["pheno_score"],
														 ascending=[False],
														 inplace=True)
		return disPhenoScDf

	def load(self,logger):
		
		#append disease description
		msg = "appending a description into disease rankikng file ..."
		lib_utils.msgout('notice',msg);logger.info(msg)
		self.cDis = {}
		
		#OMIM:614652\tPDSS2\t57107\tHP:0002133\tStatus epilepticus
		fp = open(self.disease_fn,'r')
		for i in fp:
			if i[0] == '#':continue
			disId,gene,_,hpo,disDesc = i[:-1].split('\t')
			if disId not in self.cDis:
				self.cDis[disId] = Disease(disId)
			
			if gene not in self.cDis[disId].genes:
				self.cDis[disId].genes.append(gene)
				
			if hpo not in self.cDis[disId].hpos:
				self.cDis[disId].hpos.append(hpo)
				
		fp.close()
	
		#append disease description
		msg = "appending a disease description ..."
		lib_utils.msgout('notice',msg);logger.info(msg)
		
		fp = open(self.desc_fn,'r')
		DisDesc = {}
		for i in fp:
			j=i.split('\t')
			disId = '%s:%s'%(j[0],j[1])
			if disId in self.cDis:
				if self.cDis[disId].inherit == divine_inc.INH_UNKNOWN:
					if j[4] == 'HP:0000006':
						self.cDis[disId].inherit = divine_inc.INH_DOMINANT #autosomal dominant inheri
					elif j[4] == 'HP:0000007':
						self.cDis[disId].inherit = divine_inc.INH_RECESSIVE #autosomal recessive inheri
				
				if not self.cDis[disId].desc:
					self.cDis[disId].desc = j[2]
		fp.close()
		
		msg = "done."
		lib_utils.msgout('notice',msg);logger.info(msg)
				