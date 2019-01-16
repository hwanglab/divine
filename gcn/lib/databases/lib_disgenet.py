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
from gcn.lib.utils.lib_utils import py_struct, check_if_file_valid, open2
import string

disgenet = namedtuple('disgenet',['diseaseId','geneId','score','geneName','description','diseaseName','sourceId','NofPmids','NofSnps'])

class Disgenet5:
	def __init__(self,db_dir=os.path.join(os.environ.get("DIVINE"),'gcndata','disgenet')):
		self.umls_to_hpo_fpath = os.path.join(db_dir,'latest','ls-umls2hpo.ttl')

		if not check_if_file_valid(self.umls_to_hpo_fpath):
			raise IOError('check if the file [%s] exists'%self.umls_to_hpo_fpath)

		self.gene_to_disease_fpath = py_struct(all=os.path.join(db_dir,'latest','all_gene_disease_associations.tsv.gz'),
		                                 curated=os.path.join(db_dir,'latest','all_gene_disease_associations.tsv.gz'))

		if not check_if_file_valid(self.gene_to_disease_fpath.all):
			raise IOError('check if the file [%s] exists'%self.gene_to_disease_fpath.all)

		if not check_if_file_valid(self.gene_to_disease_fpath.curated):
			raise IOError('check if the file [%s] exists'%self.gene_to_disease_fpath.curated)

		self.umls2hpo = {}
		self.umls2genes = {}
		self.umls2desc={}

	def parse_umls_to_hpo_ttl(self):

		self.umls2hpo = {}
		umls_id = None
		print "start to parse [%s] to get a umls to HPOs map"%self.umls_to_hpo_fpath
		fp = open2(self.umls_to_hpo_fpath, 'r')
		
		for i in fp:
			i = i.strip()

			if not i: continue
			if i[0] == '@': continue

			if i.startswith('<http://linkedlifedata.com'):
				mObj = re.search(r'umls\/id\/(\w+)>', i)
				if mObj:
					umls_id = mObj.group(1)
					if umls_id not in self.umls2hpo:
						self.umls2hpo[umls_id] = []
				else:
					raise RuntimeError('error in parsing [%s] for umls_id in [%s]' % (i, self.umls_to_hpo_fpath))
			elif 'hpo:' in i:
				mList = re.findall(r'hpo:(\w+)', i)
				if mList:
					for hpo in mList:
						if umls_id:
							hpo = string.replace(hpo.strip(), '_', ':')
							self.umls2hpo[umls_id].append(hpo)
						else:
							raise RuntimeError('an associated umls_id is not defined for HPO[%s]'%mList)
				else:
					raise RuntimeError('error in parsing [%s] for HPOs in [%s]' % (i, self.umls_to_hpo_fpath))
		
		fp.close()

	def get_genes_by_umls(self, min_score=0.0, which_file='curated'):

		self.umls2genes = {}
		self.umls2desc = {}
		gene2disease = namedtuple('gene2disease', ['geneSymbol', 'diseaseId', 'diseaseName', 'score'])

		if which_file == 'all':
			fp = open2(self.gene_to_disease_fpath.all, 'r')
			print("header of the file[%s]"%self.gene_to_disease_fpath.all)
			print fp.next()
		else:
			fp = open2(self.gene_to_disease_fpath.curated, 'r')
			print("header of the file[%s]"%self.gene_to_disease_fpath.all)

		print fp.next()
		for i in fp:
			j = i.split('\t')
			dnet = gene2disease(j[1],j[2],j[3],j[4])
			if float(dnet.score) < min_score: continue

			if dnet.diseaseId not in self.umls2genes:
				self.umls2genes[dnet.diseaseId] = []
			self.umls2genes[dnet.diseaseId].append(dnet.geneSymbol)

			if dnet.diseaseId not in self.umls2desc:
				self.umls2desc[dnet.diseaseId] = dnet.diseaseName

		fp.close()


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
