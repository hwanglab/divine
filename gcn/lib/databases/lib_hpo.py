"""
.. module:: annotpipe
    :platform: Unix, Windows, MacOSX
    :synopsis: A wraper to call VARANT

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is a wrapper to call VARANT. The inputs are -
1. Unannotated VCF file path
2. Path to create annotated vcf file (Option)
"""
import os
from collections import namedtuple

hpo_flat = namedtuple('hpo_flat',['diseaseId','geneName','geneId','hpoId','termDesc'])

class Disease:
	def __init__(self):
		self.desc = None
		self.hpos = []
		self.genes = []
		
class Hpo:
	def __init__(self,db_dir=os.path.join(os.environ.get("DIVINE"),'gcndata','hpo')):
		self.db_dir = db_dir
		self.flat_fn = os.path.join(db_dir,'ALL_SOURCES_TYPICAL_FEATURES_diseases_to_genes_to_phenotypes.txt')
		if not os.path.exists(self.flat_fn):
			raise IOError('[%s] does not exists'%self.flat_fn)
		self.disease2hpo_fn = os.path.join(db_dir,'phenotype_annotation.tab')
		if not os.path.exists(self.disease2hpo_fn):
			raise IOError('[%s] does not exists'%self.disease2hpo_fn)

		self.diseases = {}
		
	def get_disease_to_hpos(self):
		'''
		update self.diseases.hpos/desc
		'''
		#open disease2geen fn and store source:disease_id key
		fp = open(self.disease2hpo_fn,'r')
		for i in fp:
			itms = i.split('\t')
			disId = '%s:%s'%(itms[0],itms[1])
			if disId not in self.diseases:
				self.diseases[disId] = Disease()
			self.diseases[disId].desc = itms[2]
			self.diseases[disId].hpos.append(itms[4])
		fp.close()

	def get_disease_to_genes(self):
		
		for r in self._iterate_fn():
			if r.diseaseId in self.diseases:
				if r.geneName not in self.diseases[r.diseaseId].genes:
					self.diseases[r.diseaseId].genes.append(r.geneName)
			else:
				raise RuntimeError('disease [%s] does not appear in [%s]'%(r.diseaseId,self.disease2hpo_fn))
	
	def get_diseases_no_gene_assoc(self):
		
		#get all disease ids
		self.get_disease_to_hpos()
		
		self.get_disease_to_genes()
		
		diseases_no_gene = []
		#get disease ids having gene association
		for disId,cDis in self.diseases.iteritems():
			if not cDis.genes:
				diseases_no_gene.append(disId)
		return list(set(diseases_no_gene))
	
	def _iterate_fn(self):
		
		fp = open(self.flat_fn,'r')
		fp.next()
		for i in fp:
			j = i.rstrip().split('\t')
			hflat = hpo_flat(j[0],j[1],j[2],j[3],j[4])
			yield hflat
		fp.close()
	
