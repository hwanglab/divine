
import os, sys, argparse, re, math, datetime, time, logging
import lib_utils

class Disease:
	def __init__(self,disId):
		self.id = disId
		self.desc = None
		self.pheno_score = 0.
		self.genes = {}
		self.gdmg = {}
		self.enriched_genes = {}
		self.combined = {}
		self.inherit = 0
		self.hpos = []
		
class HpoSim:
	def __init__(self):
		#divine_cfg = read_divine_config()
		
		hpo_base_dir='gcndata/hpo'
		
		divine_cfg = {}
		divine_cfg['hposim'] = 'python_libs/pkg/hpo_similarity/tests/cj_hpo_similarity.py'
		divine_cfg['hpo_obo'] = '%s/hp.obo' % hpo_base_dir
		divine_cfg['ext_disease_to_gene'] = '%s/ALL_SOURCES_disgnet_TYPICAL_FEATURES_diseases_to_genes_to_phenotypes.txt' % hpo_base_dir
		
		self.hposim = divine_cfg['hposim']
		self.hpo_obo = divine_cfg['hpo_obo']
		self.disease_to_gene = divine_cfg['ext_disease_to_gene']

	def load_disease_pheno(self):

		fp = open(self.disease_to_gene,'r')
		disIds = {}
		
		for i in fp:
			if i[0] == '#':continue
			disId,gene,_,hpo,hpo_desc = i[:-1].split('\t')
			disId = disId.replace(":", "_")
			if disId not in disIds:
				disIds[disId] = Disease(disId)
			
			if gene not in disIds[disId].genes:
				disIds[disId].genes[gene] = 0.
			if hpo not in disIds[disId].hpos:
				disIds[disId].hpos.append(hpo)
		fp.close()
		
		return disIds
	
	def disease_hpo_to_file(self,out_dir):
		disease_hpos = []
		
		for disId,cDis in self.load_disease_pheno().iteritems():
			
			fnw = '%s/%s.hpo'%(out_dir,disId)
			fpw = open(fnw,'w')
			fpw.write('#%s\n'%disId)
			fpw.write('%s'%('\n'.join(cDis.hpos)))
			fpw.close()
			disease_hpos.append([disId,fnw])
			
		return disease_hpos

	def run(self,hpo_query,out_fn,log_dir=None,logger=None):
		job_name = 'hpo_to_diseases'
		# prepare output file
	
		msg = 'matching query phenotypes to diseases in semantic HPO ontology[%s;%s]'%(job_name,out_fn)
		lib_utils.msgout('notice',msg)
		if logger:
			logger.info(msg)

		# run hpo similarity
		cmd = ["python", self.hposim, \
					"-q", hpo_query, \
					"-b", self.hpo_obo, \
					"-f", self.disease_to_gene, \
					"-s", "funSimAvg",\
					"--normalize", \
					"-o", out_fn]

		lib_utils.runcmd2(cmd)
		msg = 'done. [%s on %s]' % (job_name,hpo_query)
		lib_utils.msgout('notice',msg)
