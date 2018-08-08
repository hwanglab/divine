"""
.. module:: annotpipe
    :platform: Unix, Windows, MacOSX
    :synopsis: A wraper to call VARANT

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is a wrapper to call VARANT. The inputs are -
1. Unannotated VCF file path
2. Path to create annotated vcf file (Option)
"""
import os, re
from collections import namedtuple

disgenet = namedtuple('disgenet',['diseaseId','geneId','score','geneName','description','diseaseName','sourceId','NofPmids','NofSnps'])
		
class Disgenet:
	def __init__(self,db_dir=os.path.join(os.environ.get("DIVINE"),'gcndata','disgenet')):
		
		db_dir = db_dir
		self.disease_id_maps = []
		self.disease_id_maps.append(os.path.join(db_dir,'ls-umls2omim-exactMatch.ttl'))
		self.disease_id_maps.append(os.path.join(db_dir,'ls-umls2hp-orphanet-exactMatch.ttl'))
		self.disease_id_maps.append(os.path.join(db_dir,'ls-umls2hp-decipher-exactMatch.ttl'))
		self.map2umls = {}
		
		for fn in self.disease_id_maps:
			if not os.path.exists(fn):
				raise IOError('[%s] does not exists'%fn)
		
		self.all_gene_disease_fn = os.path.join(db_dir,'curated_gene_disease_associations.tsv')
	
	def map_to_umls_core(self,tag,disease_source,map_fn,map2umls):
		fp = open(map_fn,'r')
		for i in fp:
			i=i.strip()
			
			if not i: continue
			if i[0]=='@': continue
			
			if 'umls' not in i: continue
			
			mObj = re.search(r'umls\/id\/(\w+)>',i)
			if mObj:
				umls_id = mObj.group(1)
			else:
				raise RuntimeError('error in parsing [%s] for umls_id in [%s]'%(i,map_fn))

			i=fp.next().strip()
			itms=i.split(',')
			for itm in itms:
				itm = itm.strip()
				mObj = re.search(r'\/%s\/(\d+)>'%tag,itm)
				if mObj:
					disease_id = mObj.group(1)
				else:
					raise RuntimeError('error in parsing [%s] for %s in [%s]'%(itm,tag,map_fn))
				dkey = '%s:%s'%(disease_source,disease_id)
				if dkey not in map2umls:
					map2umls[dkey] = []
				map2umls[dkey].append(umls_id)
		fp.close()
		return map2umls
	
	def map_to_umls(self):
		tag_prefix = ['omim','orphanet','syndrome']
		disease_sources = ['OMIM','ORPHANET','DECIPHER']

		for tag,disease_source,map_fn in zip(tag_prefix,disease_sources,self.disease_id_maps):
			self.map2umls = self.map_to_umls_core(tag,disease_source,map_fn,self.map2umls)
	
	def get_umls_ids(self,disease_ids):
		self.map_to_umls() #update map2umls

		umls2did = {}
		for did in disease_ids:
			if did in self.map2umls:
				for umls in self.map2umls[did]:
					umls2did[umls]=did
		return umls2did
	
	def _iterate_fn(self):

		fp = open(self.all_gene_disease_fn,'r')
		fp.next()
		for i in fp:
			j=i.rstrip().split('\t')
			dnet = disgenet(j[0].split(':')[1],j[1],float(j[2]),j[3],j[4],j[5],j[6],j[7],j[8])
			yield dnet
			
		fp.close()
		
	def get_genes_by_umls(self,umls2did,min_score=0.10):
		disGenes = {}
		for dnet in self._iterate_fn():
			if dnet.score<min_score: continue
			if dnet.diseaseId in umls2did:
				if umls2did[dnet.diseaseId] not in disGenes:
					disGenes[umls2did[dnet.diseaseId]] = []
				disGenes[umls2did[dnet.diseaseId]].append(dnet.geneName)
		
		return disGenes
