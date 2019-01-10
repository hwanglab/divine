#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com|hongc2@ccf.org
'''
import os, sys, argparse, math, datetime, time, logging
import pandas as pd
import numpy as np
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from ConfigParser import SafeConfigParser
from gcn.bin.prioritize.bipartite_diffusion import run_bp
from gcn.bin.prioritize import damaging_model, pagerank, divine_inc, lib_disease
import dill
import gcn.lib.io.vcf as vcf
from gcn.lib.databases.refgene import Refgene
from gcn.lib.utils import lib_utils, lib_ped
from gcn.lib.varann.vartype.varant.annotator import get_min_maf
from gcn.lib.varann.vartype.varant import varant_parser as vp
from gcn.lib.varann.vartype.varant import annotateRegion, vcf_mask
from gcn.lib.io.vcfutils import normalize_variant, get_var_type
from gcn.lib.databases import geneontology
from gcn.etc import fileconfig

VERSION = '0.1.2'
author_email = 'hongc2@ccf.org'

class PhenoGene:
	def __init__(self):
		self.score = 0.
		self.disId = None
		self.num_assoc_genes = 0

class SnvGene:
	def __init__(self):
		self.score = 0.
		self.pheno_score = 0.
		self.zyg_cnt = [0, 0] #HET,HOM
		self.scores = [[],[]]
		
class Divine:
	'''
	collect program configuration, input parameters, and computational resource
	'''
	def __init__(self, uargs):
		#transferring user input arguments to class member variables
		
		self.to_delete_fns = []
		self.exp_tag = uargs.exp_tag
		self.vknown = uargs.vknown
		self.cadd = uargs.cadd
		self.top_k_disease = uargs.top_k_disease
		
		self.excl_non_coding = False
		self.sparser = SafeConfigParser()
		
		self.omim = None
		
		self.pheno_dmg = {}
		self.gt_dmg = {}
		self.gene_dmg = {}
		self.vknown_genes = {}
		
		lib_utils.msgout('notice','initializing Divine ...','Divine')
		
		divine_root_dir = os.environ.get("DIVINE")
		if not divine_root_dir:
			raise EnvironmentError("set DIVINE variable properly!")
		
		config_fn = os.path.join(divine_root_dir,'gcn','config','divine.conf')

		if not lib_utils.check_if_file_valid(config_fn):
			raise IOError("check if the configuration file[%s] is valid!" % config_fn)
		
		self.config_fn = config_fn
		self.entries = {'divine_root':divine_root_dir}
		self._set_args(uargs)

		self.hpo_query = uargs.hpo_query
		if self.hpo_query is None:
			self.hpo2disease_fn = None
			self.pheno_dmg_fn = None
			self.disease_rank_fn = None
		else:
			self.hpo2disease_fn = self._assign_out_fn('hpo_to_diseases','tsv')
			self.pheno_dmg_fn = self._assign_out_fn('pheno_gene_rank','tsv')
			self.disease_rank_fn = self._assign_out_fn('diseases_rank','tsv')

		self.gene_rank_fn = self._assign_out_fn('gene_rank', 'tsv')
		self.vcf = uargs.vcf
		self.ped = None
		self.proband_id = None
		self.genotype = True
		
		if self.vcf:
			self.is_family_vcf = False
			if uargs.ped:
				self.is_family_vcf = True
				if uargs.proband_id:
					proband_idx = lib_ped.check_consistency_ped_vcf(\
															self.vcf,uargs.ped,uargs.proband_id)
					self.ped = uargs.ped
					self.proband_id = uargs.proband_id
				else:
					msg = "A family file [%s] was provided but you didn't provide a proband ID to examine. Specify the probrand ID available in the VCF [%s] using an option -p."\
						%(uargs.ped,self.vcf)
					print(msg)
					raise RuntimeError(msg)

			else:
				#get sample_ids contained into VCF file
				v = vcf.VCFParser(self.vcf)
				if len(v.samples) > 1:
					raise RuntimeError('VCF file [%s] contains more than two samples. Let me know which sample is a proband to diagnose!'%self.vcf)
				elif len(v.samples) == 1:
					#search sample_id and create a temp ped for the proband
					self.ped = os.path.join(self.out_dir,'proband_tmp.ped')
					self.proband_id = lib_ped.create_proband_ped(self.vcf,self.ped)
					self.to_delete_fns.append(self.ped)
				else:
					self.genotype = False
		
		self.xls = None
		self.hgmd = uargs.hgmd
		self.cosmic = uargs.cosmic
		self.dblink = uargs.dblink
		
		# damage factor w.r.t the location of variant within the transcript
		self.dm = damaging_model.DmgCoeff(\
			uargs.indel_fidel,uargs.go_seed_k,self.logger)
		
		if uargs.ref_exon_only==1:
			msg = 'VCF is going to be masked by RefGene coding region'
			lib_utils.msgout('notice',msg);self.logger.info(msg)

		self.ref_exon_only = uargs.ref_exon_only

		lib_utils.msgout('notice','done. initialization')
	
	def _assign_out_fn(self,fbase,fext='tsv'):
		if self.exp_tag:
			fn = os.path.join(self.out_dir,'%s_%s.%s'%(fbase,self.exp_tag,fext))
		else:
			fn = os.path.join(self.out_dir,'%s.%s'%(fbase,fext))
		return fn
	
	def _set_args(self, uargs):
		'''
		-objective: checking input parameters, reading config, and storing user command line
		-input: uargs (args from main())
		-output: class initialization 
		'''
		job_name = '_set_args'
		lib_utils.msgout('notice','storing input condition ...',job_name)
		
		if not uargs.hpo_query and not uargs.vcf:
			raise RuntimeError('either VCF (-v) or query phenotype (-q) file should be provided!')

		# check sanity of the input files
		if uargs.hpo_query:
			if lib_utils.check_if_file_valid(uargs.hpo_query):
				self.hpo_query = uargs.hpo_query
			else:
				raise IOError('check if [%s] is valid' % uargs.hpo_query)
		
		if uargs.vcf:
			if lib_utils.check_if_file_valid(uargs.vcf):
				self.vcf = uargs.vcf
			else:
				raise IOError('check if [%s] is valid' % uargs.vcf)
		
		if uargs.capkit:
			if uargs.capkit in ['SureSelect_V6', 'SeqCapEZ_Exome']:
				self.capkit = uargs.capkit
			else:
				raise RuntimeError("revise capture kit symbol[%s]" % uargs.capkit)
		else:
			self.capkit = None
			
		# check input condition
		if uargs.out_dir is None:
			if self.vcf:
				uargs.out_dir = os.path.join(os.path.dirname(self.vcf), 'divine')
			else:
				uargs.out_dir = os.path.join(os.path.dirname(self.hpo_query), 'divine')
		
		#create the output directory user specifies
		if uargs.out_dir.endswith('/'):
			uargs.out_dir = uargs.out_dir[:-1]
		self.out_dir = uargs.out_dir
		lib_utils.ensure_dir(self.out_dir)
		
		#prepare output file name
		self.log_dir = os.path.join(self.out_dir, 'logs')
		lib_utils.ensure_dir(self.log_dir)
		
		msg = 'prepared log directory[%s]  ...'%self.log_dir
		lib_utils.msgout('notice',msg,job_name)
				
		#prepare loggig handler
		ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')

		FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
		logging.basicConfig(filename=os.path.join(self.log_dir, 'divine_%s.log' % ts),\
											 filemode="w", level=logging.DEBUG, format=FORMAT)
		
		# ------------------------
		self.logger = logging.getLogger('divine')
		# ------------------------
		self.logger.info(msg)
		
		#read configuration file containing 3rd parties s/w path and database locations
		self._read_config(uargs.vcf_filter_cfg)

		#record user command line
		self.record_commandline()
		msg = '<divine> initialization completed [%s]'%job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

	def _set_config(self, section, entry):
		'''
		objective: read item in the configuration file
		input:
			-section: part/section in the configuration file
			-entry: item to search in the section
		output:
			-entries: member var updated
		'''
		try:
			self.entries[entry] = self.sparser.get(section, entry)
			if section in ['program_paths', 'database', 'config']:
				if not self.entries[entry].startswith('/'):
					if self.entries['divine_root'] is None:
						raise ValueError('[divine_root] should be defined first!')
					self.entries[entry] = '%s/%s'%(self.entries['divine_root'],self.entries[entry])
					
					if section in ['program_paths','database']:
						if not os.path.exists(self.entries[entry]):
							raise IOError('check if the path/file [%s] exists'%self.entries[entry])
		except:
			raise IOError('check if [%s] exists in %s[%s]' % (entry, self.config_fn, section))

	def _read_config(self,vcf_filter_cfg=None):
		'''
		objective: read configuration file
		'''
		job_name = '_read_config'
		msg = 'reading configuration file [%s;%s] ...'%(job_name,self.config_fn)
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		self.sparser.read(self.config_fn)

		self._set_config('program_paths', 'varant')
		self._set_config('program_paths', 'hposim')
		self._set_config('program_paths', 'vcf2xls')
		
		self._set_config('config', 'temp_dir')
		if not vcf_filter_cfg:
			self._set_config('config', 'vcf_filter_conf')
		else:
			if os.path.exists(vcf_filter_cfg):
				self.entries['vcf_filter_conf'] = vcf_filter_cfg
			else:
				raise RuntimeError('check if the file [%s] is valid'%vcf_filter_cfg)

		self._set_config('database', 'ext_disease_to_gene')

		self._set_config('database', 'disease_desc')

		self._set_config('database', 'hpo_obo')

		self._set_config('database', 'beta_fit')
		self._set_config('database', 'string_link')
		
		'''
		to access to UCSC mysql database(hg19)
		select e2g.value, gtp.protein from ensGtp as gtp
		inner join ensemblToGeneName as e2g on e2g.name=gtp.transcript;
		'''
		self._set_config('database', 'esp_to_gene')
		self._set_config('database', 'kegg_hsa')

		# check if the file or directory all exists before long journey!
		for key, path2 in self.entries.iteritems():
			if not lib_utils.check_if_file_valid(path2):
				raise IOError('check [%s = %s] in the file [%s]' %\
										(key, path2, self.config_fn))

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		return self.entries
	
	def record_commandline(self):
		'''
		objective: record the divine run condition into logger
		'''
		import socket
		job_name = 'record_commandline'
		msg='capturing user command line [%s] ...'%job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		try:
			host_name = socket.gethostname()
		except:
			host_name = 'N/A'
		self.logger.info('host:%s'%host_name)
		
		try:
			user = os.environ.get('USER')
		except:
			user = 'N/A'
		self.logger.info('user:%s'%user)
		
		try:
			pwd = os.environ.get('PWD')
		except:
			pwd = 'N/A'
		self.logger.info('pwd:%s'%pwd)

		self.logger.info('cmd:%s'%(' '.join(sys.argv)))
		self.logger.info("divine configuration file:%s" % self.config_fn)
		
		self.logger.info('exclude_non_coding:%s'%self.excl_non_coding)
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
	
	def hpo_to_diseases(self,top_k_disease=0):
		'''
		objective: match HPO IDs from a given patint phenotype to known disease database
		input: hpo_query, hpo database
		method: hposim (funSimMax)
		output: phenotype similarity between patient and known diseases, store the HPO similarity into pheno_dmg  
		'''
		job_name = 'hpo_to_diseases'

		msg = 'matching query phenotypes to diseases in semantic HPO ontology[%s;%s]'%(job_name,self.hpo2disease_fn)
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		# run hpo similarity
		cmd = ["python", self.entries['hposim'], \
					"-q", self.hpo_query, \
					"-b", self.entries['hpo_obo'], \
					"-f", self.entries['ext_disease_to_gene'], \
					"--normalize", \
					"-o", self.hpo2disease_fn]

		lib_utils.runcmd2(cmd,self.log_dir,self.logger,job_name)

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

	def annotate_comphet_inherit(self, reuse=False):
		if not self.ped: return
		msg = "for a multi-sample VCF containing parent genotypes, append inheritance model"
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

		# to build database for transcript exon regions if not exist
		refgene_tx_fn = fileconfig.FILECONFIG['REFGENE']
		refgene_tx_fn_dir = os.path.dirname(os.path.abspath(refgene_tx_fn))
		genmod_db_dir = os.path.join(refgene_tx_fn_dir, 'genmod_db')
		if not os.path.exists(genmod_db_dir):
			os.makedirs(genmod_db_dir)

		if not os.path.exists(os.path.join(genmod_db_dir, 'genes.db')) or \
				not os.path.exists(os.path.join(genmod_db_dir, 'exons.db')):
			cmd = ["genmod", "build", \
						 "-t", "gene_pred", \
						 "--splice_padding", "2", \
						 "-o", genmod_db_dir, \
						 refgene_tx_fn]

			job_name = "annotate_inheritance.genmod_build_db"
			lib_utils.runcmd2(cmd, self.log_dir, self.logger, job_name)

		vcf_genmod_out = lib_utils.file_tag2(self.vcf, 'genmod', '')
		cmd = ["genmod", "annotate", \
					 "-r", \
					 self.vcf, \
					 "|", \
					 "genmod", "models", \
					 "-", \
					 "--family_file", self.ped, \
					 "-o", vcf_genmod_out]

		job_name = "annotate_inheritance.genmod_model"
		msg = "annotating inheritance model into VCF ..."
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

		if not reuse or not lib_utils.check_if_file_valid(vcf_genmod_out):
			lib_utils.runcmd2(cmd, self.log_dir, self.logger, job_name)

		msg = "Done."
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)
		self.vcf = vcf_genmod_out

	def _store_hposim_outfn(self,hpo2disease_fn,top_k_disease,mutated_genes=[]):

		fp = open(hpo2disease_fn,'r') #assuming there is no duplicated entry; also sorted by desc score
		k = 0
		self.pheno_dmg = {} #initialize
		for i in fp:
			if top_k_disease>0 and k>top_k_disease:
				break
			if i[0]=='#':continue
			disId,genes,score = i.strip().split('\t')

			score = float(score)
			Genes = genes.split(',')
			effGenes = lib_utils.intersect(mutated_genes, Genes)
			scale = np.log(len(effGenes) + 3)

			if score>0. and score > self.omim.cDis[disId].pheno_score:
				self.omim.cDis[disId].pheno_score = score
				# considering the number of genes associated with the disease, adjust the score to assign to each gene; keep only max pheno-match disease
				score /= scale
				for gene in Genes:
					#add one disease into each gene
					if gene not in self.pheno_dmg:
						self.pheno_dmg[gene] = PhenoGene()

					if score > self.pheno_dmg[gene].score:
						self.pheno_dmg[gene].score = score
						self.pheno_dmg[gene].disId = disId

			k += 1
		fp.close()

	def vannotate(self,reuse=False):
		'''
		objective: run varant (GCN) annotator
		input: self.vcf
		output: annotated vcf
		'''
		job_name = 'vannotate'
		msg = 'annotating VCF file[%s;%s] ...'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		# prepare output file
		varant_vcf = os.path.join(self.out_dir,'divine.vcf')
		
		# if necessary, masking the raw vcf file
		coding_vcf = None
		if self.ref_exon_only>0:
			if not lib_utils.check_if_file_valid(varant_vcf) or not reuse:
				cRef = annotateRegion.RefGeneUcscTB(work_dir=self.out_dir,logger=self.logger)
				coding_bed_fn = cRef.create_bed(ext_bp=20,reuse=False)
				
				msg = 'extracting variants in coding region from [%s] @ %s ...'%(self.vcf,job_name)
				lib_utils.msgout('notice',msg);self.logger.info(msg)
				
				coding_vcf = os.path.join(self.out_dir,'refgene_e20.vcf')
				self.vcf = vcf_mask.by_bed(self.vcf,coding_bed_fn,coding_vcf,logger=self.logger)
				
				msg = 'done.@ %s'%job_name
				lib_utils.msgout('notice',msg);self.logger.info(msg)

		if not lib_utils.check_if_file_valid(varant_vcf) or not reuse:
			self.logger.info('annotating [%s,%s] ...' % (job_name, self.vcf))
			
			cmd = ["python", self.entries['varant'], \
						"-i", self.vcf, \
						"-o", varant_vcf, \
						"-l", self.log_dir]
			if self.capkit:
				cmd.extend(["-c", self.capkit, "-e", "180"])

			if self.hgmd>0:
				cmd.extend(["--hgmd"])
			
			if self.cosmic>0:
				cmd.extend(["--cosmic"])

			if self.dblink>0:
				cmd.extend(["--dblink"])

			lib_utils.runcmd2(cmd,self.log_dir,self.logger,job_name)
			
		self.vcf = varant_vcf
		
		if coding_vcf:#cleanup intermediary file
			os.unlink(coding_vcf)
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
	
	def vfilter(self):
		'''
		objective:apply a standard filter to VCF file and classify variants
		input: annotated vcf from varant (GCN) annotator
		output: filtered vcf
		'''
		job_name = 'vfilter'
		msg = 'filtering the annotated VCF [%s;%s] ...'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		filtered_vcf = self._assign_out_fn(job_name,'vcf')

		msg='applying a standard filter/class tagging [%s]' % self.vcf
		lib_utils.msgout('notice',msg,job_name);self.logger.info(msg)
		
		gcn_filter = os.path.join(self.entries['divine_root'], \
														'gcn', 'lib', 'utils', 'filter_cj.py')
		
		cmd = ["python", gcn_filter, \
					"-i", self.vcf, \
					"-o", filtered_vcf]
		
		if not self.genotype:
			cmd.append("--no_genotype")
			
		filter_conf = self.entries['vcf_filter_conf']
		cmd.extend(["-f", filter_conf])
		
		self.logger.info('filter config [%s] is applied' % filter_conf)
		lib_utils.runcmd2(cmd,self.log_dir,self.logger,job_name)
		self.vcf = filtered_vcf
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

	def gather_pdomain_scores(self, vcfParser):

		msg = 'gathering pathogenic variant density in domains ...'
		lib_utils.msgout('notice', msg)
		self.logger.info(msg)

		pdomains = lib_utils.py_struct(ridx=[],
																	 denoms=[],
																	 benign_dens=[],
																	 vus_dens=[],
																	 patho_dens=[])

		ridx = 0
		for rec in vcfParser:
			vcfParser.parseinfo(rec)
			# to collect pdomain info
			if rec.info.PATHO_DOMAIN:
				pdoms = [float(pdom) for pdom in rec.info.PATHO_DOMAIN.split(',')]
				pdomains.ridx.append(ridx)
				pdomains.denoms.append(pdoms[0])
				pdomains.benign_dens.append(pdoms[1])
				pdomains.vus_dens.append(pdoms[2])
				pdomains.patho_dens.append(pdoms[3])

			ridx += 1

		pdomains = pd.DataFrame({'ridx': pdomains.ridx,
														 'denoms': pdomains.denoms,
														 'benign_dens': pdomains.benign_dens,
														 'vus_dens': pdomains.vus_dens,
														 'patho_dens': pdomains.patho_dens,
														 'phat_lo':None,
														 'patho_dens_p':None})

		phat = pdomains.patho_dens / (pdomains.benign_dens + pdomains.vus_dens + pdomains.patho_dens)

		tgt_z = damaging_model.get_z(confidence=0.75)

		pdomains['phat_lo'] = map(lambda x1,x2: damaging_model.ci_lower_bound(x1, x2, z=tgt_z), phat, pdomains.denoms)

		tgt_pctile = 50
		pdensl = np.log10(pdomains.patho_dens+1e-12)

		tgt_pctile_sc = np.percentile(pdensl, tgt_pctile)

		y = (pdensl >= tgt_pctile_sc).astype(np.float)
		X = pdensl[:, np.newaxis]

		model2 = LogisticRegression().fit(X, y)
		pdomains['patho_dens_p'] = model2.predict_proba(X)[:, 1]

		pdomains_default = lib_utils.py_struct(phat_lo=np.percentile(pdomains['phat_lo'], 15),
																			 patho_dens_p=np.percentile(pdomains['patho_dens_p'], 15))

		return pdomains, pdomains_default

	def _extract_mutation_info(self,beta_fits):
		'''
		objective: to extract (gene_qsymbol,mutation_type,variant_class_tag,transcript_length,insillico_prediction_score,MAF_significance_offset,zygosity) from annotated/filtered VCF file
			to transfer genmod information to class_tag and also get rid of some redundancy in VCF info
		'''
		job_name = '_extract_mutation_info'
		msg='collecting variant information and class label to determine genetic damage [%s;%s]...'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		mutation_info = []
		if self.proband_id:
			rewrite_vcf = True
			v = vcf.VCFParser(self.vcf,sampleids=[self.proband_id])
			pdom,pdom0 = self.gather_pdomain_scores(v)
			v.stream.close()
			vcf_tmp = self.vcf+'.tmp'
			ostream = open(vcf_tmp, 'w')
			rmInfo = ['Exonic','Annotation','Compounds']
			v = vcf.VCFParser(self.vcf)
			v.writeheader(ostream,to_del_info = rmInfo)
		else:
			rewrite_vcf = False
			v = vcf.VCFParser(self.vcf)
			pdom,pdom0 = self.gather_pdomain_scores(v)
			v.stream.close()
			v = vcf.VCFParser(self.vcf)

		msg = 'Importing max transcript length for each gene ...'
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

		refgene = Refgene()
		cds_lens = refgene.get_max_cds_length()
		tx_lens = {}
		for gene, cds_len in cds_lens.iteritems():
			tx_lens[gene] = int(cds_len/3.)

		ridx = 0
		for rec in v:
			
			v.parseinfo(rec)
			
			#to remove redundant gene symbols annotated by genmod but add transcript version
			if rewrite_vcf:
				for rkey in rmInfo:
					v.delete_info(rec, rkey)

				if rec.info.GeneticModels:
					genmod_tag = lib_ped.parse_genmod_inherit_model(\
												rec.info.GeneticModels[0].split(':')[1])
					rec.info.CLASS_TAG += genmod_tag
					
				v.write(ostream, rec)
			
			varlist = normalize_variant(rec.chrom, rec.pos, rec.ref, rec.alt)
			mut_type = varlist[0][-1]

			if ':' in rec.id[0]:
				mut_type = 'mnp'
			
			# collect conservation prediction score (CADD and GERP++) 
			cadd_aa = './.'
			px_cadd = None
			if rec.info.CADD_raw:
				# to get CADD_raw (average)
				px_cadd, cadd_aa = vcf.get_CADD_scores(mut_type, rec.info.CADD_aa, rec.info.CADD_raw, beta_fits)
			
			# to get GERP++ score
			px_gerp = None
			if rec.info.GerpConserve:
				px_gerp = vcf.get_GERP_scores(mut_type, cadd_aa, rec.info.GerpRSScore, beta_fits)
			
			# which score can be chosen
			px = 0.5
			if self.cadd>0 and px_cadd is not None:
				px = px_cadd
			elif px_gerp is not None:
				px = px_gerp

			vpop = vp.parse(rec.info)
			genes = []
			
			# to get MAF in the order of ExAC, ESP, and 1K
			if rec.info.EXACDB:
				maf = get_min_maf(rec.info.EXACAF[0])
			elif rec.info.ESPDB:
				maf = get_min_maf(rec.info.ESPAF[0])
			elif rec.info.KGDB:
				maf = get_min_maf(rec.info.KGAF[0])
			else:
				maf = 0.
			
			# to compute a significance of MAF
			maf_offset = self.dm.get_maf_xoffset(maf)

			#pdom.iloc[ridx]==ridx
			pdom_idx = pdom.index[pdom.ridx == ridx].tolist()
			if pdom_idx:
				patho_p = pdom.phat_lo[pdom_idx[0]]
				patho_pden = pdom.patho_dens_p[pdom_idx[0]]
			else:
				# assign a default pathogenic domain value (15% quantile value)
				patho_p = pdom0.phat_lo
				patho_pden = pdom0.patho_dens_p

			vartype = get_var_type(rec.ref,rec.alt)

			# to get transcript length
			for altnum, val in vpop.items():
				# for each gene involved with the variant
				for gene, gd in val.items():
					protein_len = self.dm.avg_protein_len
					if gene in tx_lens:
						protein_len = tx_lens[gene]

					# store a set of essential annotation to be used for genetic damage
					if gene not in genes:
						mutation_info.append([gene, vartype, rec.info.CLASS_TAG, protein_len, px, maf_offset, patho_p, patho_pden])
						genes.append(gene)

			ridx += 1
			
		# done reading filterd VCF file
		if rewrite_vcf:
			v.stream.close()
			ostream.close()
			os.rename(vcf_tmp,self.vcf)
			
		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)

		return mutation_info
	
	def get_kth_score(self,dmg,topR):
		msg = 'getting [%d]-th top pheno_score...'%topR
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		scores = []
		for scDid in dmg.itervalues():
			scores.append(scDid.score)
		
		if topR<1.:
			s1 = round(topR*len(scores))
		else:
			s1 = topR
		
		scores.sort(reverse=True)
		
		msg = 'selected pheno score:%g'%scores[s1]
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		return scores[s1]

	def _predict_gt_dmg(self, mutation_info):
		#
		job_name = '_predict_gt_dmg'
		msg = 'combining impact by variant location in tx and conservation pred score [%s] ...' % job_name
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

		pheno_sc_dom = 0.
		if self.hpo_query:
			# to get k-th score
			pheno_sc_dom = self.get_kth_score( \
				self.pheno_dmg, self.dm.conf_pheno_rank_dom)

		# estimate genetic damage scores per gene
		for gene, vtype, tag, protein_len, px, maf_offset, pdom, pdom_dens in mutation_info:

			vreg_dmg = 0.
			# filter out frequent or known-benign
			if '1' in tag or '2' in tag or 'b' in tag: continue

			# collect inheritance model of disease associated with the pheno gene
			dgInherit = divine_inc.INH_RECESSIVE
			if self.hpo_query and gene in self.pheno_dmg:
				disId = self.pheno_dmg[gene].disId
				if disId:
					dgInherit = self.omim.cDis[disId].inherit
				else:
					dgInherit = divine_inc.INH_UNKNOWN
				if dgInherit == divine_inc.INH_DOMINANT:
					if self.pheno_dmg[gene].score < pheno_sc_dom:
						dgInherit = divine_inc.INH_UNKNOWN
				elif dgInherit == divine_inc.INH_UNKNOWN:
					if self.pheno_dmg[gene].score < pheno_sc_dom:
						dgInherit = divine_inc.INH_RECESSIVE

			# to locate variant location
			if 'n' in tag:
				if self.excl_non_coding:
					continue
				else:
					vreg_dmg += self.dm.vncoding
			elif 'i' in tag:
				vreg_dmg += self.dm.vintronic
			elif 'e' in tag:
				vreg_dmg += self.dm.vexonic
			elif 'S' in tag:
				vreg_dmg += self.dm.vsplice
			elif 's' in tag:
				vreg_dmg += self.dm.vsplice_syn
			elif 'w' in tag:
				vreg_dmg += self.dm.warning

			clipatho = False
			if self.vknown:
				if 'c' in tag:  # previously known pathogenic (from ClinVar or HGMD)?
					clipatho = True
					if gene in self.vknown_genes: continue
					vreg_dmg += self.dm.cexonic
					self.vknown_genes[gene] = True
				elif 'I' in tag:  # was it in intronic/intergenic/non-coding?
					vreg_dmg += self.dm.cintronic
				elif 'g' in tag:  # pathogenic gene?
					vreg_dmg += self.dm.cgene

			vreg_dmg += maf_offset

			if gene not in self.gt_dmg:
				self.gt_dmg[gene] = SnvGene()

			if 'H' in tag:
				zygosity = divine_inc.ZYG_HOM
			else:
				zygosity = divine_inc.ZYG_HET

			self.gt_dmg[gene].zyg_cnt[zygosity] += 1

			if clipatho:
				if zygosity == divine_inc.ZYG_HOM:
					dmg_allele_wt = self.dm.hom_denovo_snp
				elif 'h' in tag:
					dmg_allele_wt = self.dm.het2_denovo_snp
				else:
					dmg_allele_wt = self.dm.het_denovo_snp
			elif '3' in tag:  # rare
				if vtype != 'mismatch':
					if zygosity == divine_inc.ZYG_HOM or \
							dgInherit == divine_inc.INH_DOMINANT:
						dmg_allele_wt = self.dm.hom_rare_indel
					elif 'h' in tag or dgInherit == divine_inc.INH_UNKNOWN:  # is it compound het?
						dmg_allele_wt = self.dm.het2_rare_indel
					else:
						dmg_allele_wt = self.dm.het_rare_indel
				else:
					if zygosity == divine_inc.ZYG_HOM or \
							dgInherit == divine_inc.INH_DOMINANT:
						dmg_allele_wt = self.dm.hom_rare
					elif 'h' in tag or dgInherit == divine_inc.INH_UNKNOWN:  # is it compound het?
						dmg_allele_wt = self.dm.het2_rare
					else:
						dmg_allele_wt = self.dm.het_rare
			elif '4' in tag:  # de-novo
				if vtype != 'mismatch':
					if zygosity == divine_inc.ZYG_HOM or \
							dgInherit == divine_inc.INH_DOMINANT:
						dmg_allele_wt = self.dm.hom_denovo_indel
					elif 'h' in tag or dgInherit == divine_inc.INH_UNKNOWN:
						dmg_allele_wt = self.dm.het2_denovo_indel
					else:
						dmg_allele_wt = self.dm.het_denovo_indel
				else:
					if zygosity == divine_inc.ZYG_HOM or \
							dgInherit == divine_inc.INH_DOMINANT:
						dmg_allele_wt = self.dm.hom_denovo_snp
					elif 'h' in tag or dgInherit == divine_inc.INH_UNKNOWN:
						dmg_allele_wt = self.dm.het2_denovo_snp
					else:
						dmg_allele_wt = self.dm.het_denovo_snp

			if vtype in ['complex', 'mnp']:
				dmg_allele_wt *= self.dm.cmp_mnp_penalty

			score = ((1. - self.dm.prwt) * (vreg_dmg+pdom_dens) + self.dm.prwt * (px+pdom)) * dmg_allele_wt / np.log10(protein_len)

			if self.gt_dmg[gene].zyg_cnt[zygosity] > 1:
				if self.gt_dmg[gene].scores[zygosity]:
					m, mscore = lib_utils.argmin(self.gt_dmg[gene].scores[zygosity])
					if score > mscore:
						self.gt_dmg[gene].scores[zygosity][m] = score
				else:
					self.gt_dmg[gene].scores[zygosity].append(score)
			else:
				self.gt_dmg[gene].scores[zygosity].append(score)

		for gene, snvGene in self.gt_dmg.iteritems():
			scores = snvGene.scores[0] + snvGene.scores[1]
			scores.sort(reverse=True)
			self.gt_dmg[gene].score = sum(scores[:2])

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

	def normalize2(self):

		if self.vcf:
			self.norm_genetic_dmg()

		# to normalize phenogene scores
		if self.hpo_query:
			self.norm_pheno_dmg()

	def preprocess_dmg_scores(self):
		'''
		-objective:
		-output: dictionary {gene:genetic damaged score}
		'''
		job_name = 'preprocess_dmg_scores'
		gdmg = []

		if self.vcf:
			msg='start to predict genetic damage score from variants in the provided VCF [%s]' % (job_name)
			lib_utils.msgout('notice',msg);self.logger.info(msg)

			msg = 'loading training model of CADD/GERP w.r.t AA change...'
			lib_utils.msgout('notice',msg);self.logger.info(msg)
			try:
				beta_fit_dill = self.entries['beta_fit']
				msg='loading beta fit cdf[%s] for conservation score w.r.t. AA'%beta_fit_dill
				lib_utils.msgout('notice',msg); self.logger.info(msg)
				fp = open(beta_fit_dill, 'rb')
				beta_fits = dill.load(fp)
				fp.close()
			except:
				beta_fits = [None, None, None]

			# to extract some info from annotated/filterd VCF to evaluate the genetic mutation damage
			# [gene, indel, class_tag, protein_len, in-sillico pred score, maf_offset, zygosity]
			mutation_info = self._extract_mutation_info(beta_fits)

			# to get a gene list having genetic mutations

			for minfo in mutation_info:
				if minfo[0] not in gdmg:
					gdmg.append(minfo[0])
				gdmg = list(set(gdmg))

		if self.hpo2disease_fn:
			self._store_hposim_outfn(self.hpo2disease_fn, self.top_k_disease, gdmg)

		# to enrich phenogenes (update self.pheno_dmg)
		if self.hpo_query and self.dm.go_seed_k>0 and gdmg:
			self.enrich_pheno_genes(gdmg)

		if self.vcf:
			# combine variant location and conservation pred dmg
			self._predict_gt_dmg(mutation_info)
		elif self.hpo_query:
			for gene in self.pheno_dmg.iterkeys():
				if gene not in self.gt_dmg:
					self.gt_dmg[gene] = SnvGene()
				self.gt_dmg[gene].score = self.pheno_dmg[gene].score

		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)

	def get_seed_genes(self, top_k):
		'''
		to collect genes associated a disease whose matching score to HPO is relatively high
		'''
		job_name = 'get_seed_genes'
		msg = 'collecting genes associated with top[%d] diseases showing high HPO similarity' % top_k
		lib_utils.msgout('notice',msg); self.logger.info(msg)

		msg = 'Loading adjusted phenotype-to-disease matching scores [omim.cDis.pheno_score] ...'
		lib_utils.msgout('notice', msg); self.logger.info(msg)

		genes = []
		scores = []
		pheno_scores = []
		k = 1
		disPhenoScDf=self.omim.get_disId_sorted_by_pheno_score(order='desc')
		#TOOD: for loop with pd.df upto top_k
		for i, row in disPhenoScDf.iterrows():
			if row.disId in self.omim.cDis:
				cD = self.omim.cDis[row.disId]
				registered = False
				for gene in cD.genes:
					if gene in self.pheno_dmg:
						genes.append(gene)
						scores.append(self.pheno_dmg[gene].score)
						registered = True
				if registered:
					pheno_scores.append(cD.pheno_score)
					k += 1
				if k > top_k:
					break


		msg = 'total [%d] genes are chosen for GO seeds in top[%d] diseases'%(len(genes), top_k)
		msg += ', done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)

		return genes, scores, sum(pheno_scores)/k

	def enrich_pheno_genes(self, ggenes):
		'''
		Objective:Gene-ontology enrichment (select private members of purturbed gene that highly matched with phenotypic-scored genes and assign predicted phenotypic score instead of assigning de-novo prior)
		Input:
			-pheno_dmg = {gene1:0.2,gene2:0.9,...} #e.g. phenotype score
			-genetic_dmg = {gene2:0.4,gene3:0.3,...} #e.g. genetic score
		'''

		job_name = 'enrich_pheno_genes'
		msg = 'enriching perturbed genes with both GO semantic similarity and KEGG pathways [%s] ...' % job_name
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

		# collect genes from both phenotype and genotype perturbation
		pgenes = list(self.pheno_dmg.keys())  # assuming that it's score >0
		P = len(pgenes)

		msg = 'total phenotypic genes before enrichment:%d' % P
		lib_utils.msgout('notice', msg, job_name);
		self.logger.info(msg)

		msg = 'total perturbed genes:%d' % len(ggenes)
		lib_utils.msgout('notice', msg, job_name);
		self.logger.info(msg)

		# draw a venn diagram and get genes not reported by phenotypes among genetic dmg genes
		priv_ggenes = lib_utils.difference(pgenes, ggenes)
		msg = 'the number of genes not associated with the given phenotypes:%d' % len(priv_ggenes)
		lib_utils.msgout('notice', msg, job_name);
		self.logger.info(msg)

		# to collect genes highly matched to do GO enrichment
		# Gene-ontology enrichment (select private members of purturbed gene that highly matched with phenotypic-scored genes and assign predicted phenotypic score instead of assigning de-novo prior)
		seed_pheno_genes, seed_scores, _ = \
			self.get_seed_genes(self.dm.go_seed_k)

		# query high-scored phenotype genes against private genetic-perturbed genes and bring high-matched ones
		msg = 'Using [%d] seed genes to enrich [%d] genetic variant genes with GO similarity ...' % (len(seed_pheno_genes),len(priv_ggenes))
		lib_utils.msgout('notice', msg, job_name);
		self.logger.info(msg)
		go = geneontology.Geneontology()
		goSimScores = go.get_funsim(seed_pheno_genes, priv_ggenes, min_score=self.dm.gosim_min)

		msg = 'Using [%d] seed genes to enrich [%d] genetic variant genes with KEGG similarity ...' % (len(seed_pheno_genes), len(priv_ggenes))
		lib_utils.msgout('notice', msg, job_name);
		self.logger.info(msg)

		# updating the original phenotype damage score
		# weighting enriched phenotype matching score to the gene not reported in the original phenotypes

		delta_pheno = {}
		if goSimScores:
			msg = 'Collecting [%d] GO enriched genes, enrichment_penality_ratio [%g] ...' % (len(goSimScores),self.dm.go_penalty)
			lib_utils.msgout('notice', msg, job_name);
			self.logger.info(msg)

			for pair, go_sc in goSimScores.iteritems():
				# search for a gene matched to seed pheno gene
				if pair[0] in priv_ggenes:
					new_gene = pair[0]
					seed_gene = pair[1]
				else:
					new_gene = pair[1]
					seed_gene = pair[0]

				score2 = go_sc * self.dm.go_penalty * self.pheno_dmg[seed_gene].score

				if score2 > 0.:
					# register enriched genes
					if new_gene not in delta_pheno:

						delta_pheno[new_gene] = lib_utils.py_struct(go=[0., None, None],
																												kegg=[0.,0.],
																												score=0.)

						delta_pheno[new_gene].go[2] = self.pheno_dmg[seed_gene].disId

					if score2 > delta_pheno[new_gene].go[0]: #keep only max score
						delta_pheno[new_gene].go[0] = score2
						delta_pheno[new_gene].go[1] = seed_gene
						delta_pheno[new_gene].go[2] = self.pheno_dmg[seed_gene].disId
						delta_pheno[new_gene].score = delta_pheno[new_gene].go[0]

			msg = 'Genes enriched by GO similarity:[%s]' % lib_utils.joined(delta_pheno.keys(), ',')
			lib_utils.msgout('notice', msg)
			self.logger.info(msg)

		seed_pheno_genes, seed_scores, mean_seed_score = \
			self.get_seed_genes(self.dm.go_seed_k * 4)  # update self.go_seeds

		msg = 'Using [%d] seed genes to enrich [%d] genetic variant genes with KEGG pathway genes ...' % (
		len(seed_pheno_genes), len(priv_ggenes))
		lib_utils.msgout('notice', msg, job_name);
		self.logger.info(msg)

		# query seed_pheno_genes to KEGG matrix and normalize the matched genes and ranking!
		keggEnriched = run_bp(seed_pheno_genes, seed_scores, priv_ggenes, kegg_genes_fn=self.entries['kegg_hsa'])

		if keggEnriched:
			msg = 'Collecting [%d] KEGG enriched genes with mean seed score [%g]...' % (len(keggEnriched),mean_seed_score)
			lib_utils.msgout('notice', msg, job_name);
			self.logger.info(msg)

			for kgene, kscore in keggEnriched.iteritems():
				# search for a gene matched to seed pheno gene
				score2 = kscore * mean_seed_score

				if score2 > 0.:
					# register enriched genes
					if kgene not in delta_pheno:

						delta_pheno[kgene] = lib_utils.py_struct(go=[0., None, None],
																										 kegg=[0],
																										 score=0.)
					if score2 > delta_pheno[kgene].kegg[0]: #keep only max score and sum two enriched scores
						delta_pheno[kgene].kegg[0] = score2
						delta_pheno[kgene].score = delta_pheno[kgene].go[0] + delta_pheno[kgene].kegg[0]

			msg = 'Genes enriched by KEGG bipartite network difussion:[%s]' % lib_utils.joined(keggEnriched.keys(), ',')
			lib_utils.msgout('notice', msg)
			self.logger.info(msg)

		max_score = -1.
		max_seed_gene = None
		msg = 'Total [%d] mutated genes that did not have any phenotype score previously are enriched. Assigning a new phenotype score to each enriched gene ...' % len(delta_pheno)
		lib_utils.msgout('notice', msg, job_name)
		self.logger.info(msg)
		if delta_pheno:
			for gene, deltaP in delta_pheno.iteritems():
				if deltaP.score > max_score:
					max_score = deltaP.score
					if deltaP.go[1]:
						max_seed_gene = deltaP.go[1]

		if max_score > 0:
			if max_seed_gene:
				max_enriched_score = self.pheno_dmg[max_seed_gene].score
			else:
				max_seed_gene = self.get_max_pheno_dmg()
				max_enriched_score = self.pheno_dmg[max_seed_gene].score
			max_scaled = max_enriched_score * self.dm.go_penalty * 2

			for ngene,deltaP in delta_pheno.iteritems():
				self.pheno_dmg[ngene] = PhenoGene()
				self.pheno_dmg[ngene].score = delta_pheno[ngene].score*max_scaled/max_score
				self.pheno_dmg[ngene].disId = deltaP.go[2]

				if deltaP.go[2]:
					self.omim.cDis[deltaP.go[2]].enriched_genes[ngene] = None
					if deltaP.go[1]:
						self.omim.cDis[deltaP.go[2]].enriched_genes[ngene] = deltaP.go[1]

		msg = 'max scaled phenotype score[%g], raw max enriched score[%g]' % (max_scaled,max_score)
		lib_utils.msgout('notice', msg, job_name)
		self.logger.info(msg)

	def get_max_pheno_dmg(self):
		gene_max = None
		max_score = -1.
		for gene,pdmg in self.pheno_dmg.iteritems():
			if pdmg.score > max_score:
				gene_max = gene
				max_score = pdmg.score
		return gene_max

	def simple_bayesian_pred(self, pdmg, gdmg):

		if pdmg > 0. and gdmg > 0.:
			comb_score = [pdmg * gdmg / (pdmg * gdmg \
																						+ (1. - pdmg) * (1. - gdmg)), 0.]
		else:
			comb_score = [0., 0.]
		return comb_score

	def combine_pheno_gt_dmg(self):

		job_name = 'combine_pheno_gt_dmg'

		msg = 'combining both phenotypes[%s] and geneotype[%s] damage scores ... [%s]' % \
					(self.hpo_query, self.vcf, job_name)
		lib_utils.msgout('notice', msg)
		self.logger.info(msg)

		L = len(self.gt_dmg.keys())
		msg = "total number of genes to investigate [%d]" % L
		lib_utils.msgout('notice', msg)
		# to prepare final gene-level dmg score self.gene_dmg
		if L==0:
			msg = 'combine_phenotype_gt_dmg() should not be called when neither VCF nor HPO query is given!'
			lib_utils.msgout('error',msg)
			raise RuntimeError(msg)
		elif not self.vcf:
			gdmg0 = 1. / L
			for gene in self.gt_dmg.iterkeys():
				pdmg = (1. - self.dm.ptwt) * self.gt_dmg[gene].score
				gdmg = gdmg0
				self.gene_dmg[gene] = self.simple_bayesian_pred(pdmg, gdmg)
		elif not self.hpo_query:
			pdmg0 = 1. / L
			for gene in self.gt_dmg.iterkeys():
				gdmg = (1. - self.dm.ptwt) * self.gt_dmg[gene].score
				pdmg = pdmg0
				self.gene_dmg[gene] = self.simple_bayesian_pred(pdmg, gdmg)
		else:
			self.logger.info(msg)
			for gene in self.gt_dmg.iterkeys():
				gdmg = (1. - self.dm.ptwt) * self.gt_dmg[gene].score
				pdmg = self.dm.ptwt * self.gt_dmg[gene].pheno_score
				self.gene_dmg[gene] = self.simple_bayesian_pred(pdmg, gdmg)

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice', msg)
		self.logger.info(msg)

		return self.gene_dmg

	def norm_pheno_dmg(self):

		msg = 'normalizing phenogenes by sum ...'
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

		gt_dmg = pd.DataFrame({
			'gene':self.gt_dmg.keys(),
			'score':[gt.score for gt in self.gt_dmg.itervalues()]}
		)

		pn_dmg = pd.DataFrame({
			'gene':self.pheno_dmg.keys(),
			'pheno_score':[pn.score for pn in self.pheno_dmg.itervalues()]}
		)

		gt_dmg = pd.merge(gt_dmg, pn_dmg, how='left', on='gene')
		gt_dmg.loc[gt_dmg.pheno_score.isna(), 'pheno_score'] = \
			gt_dmg.pheno_score.min() * self.dm.min_dmg_prior

		gt_dmg.pheno_score /= gt_dmg.pheno_score.sum()

		for r in gt_dmg.itertuples():
			self.gt_dmg[r.gene].pheno_score = r.pheno_score
			if r.gene in self.pheno_dmg:
				self.pheno_dmg[r.gene].score = r.pheno_score

		for pgene in self.pheno_dmg:
			if not any(gt_dmg.gene == pgene):
				self.pheno_dmg[pgene] = None

		msg += ', done.'
		lib_utils.msgout('notice', msg);
		self.logger.info(msg)

	def norm_genetic_dmg(self):

		msg = 'normalizing genetic_dmg by sum...'
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		gt_dmg_min = 1.
		denom = 0.
		
		for cSnvGene in self.gt_dmg.itervalues():

			if cSnvGene.score < gt_dmg_min:
				gt_dmg_min = cSnvGene.score
			denom += cSnvGene.score
		
		msg = '# of mutated genes:%d'%len(self.gt_dmg.keys())
		msg += ', denom for normalization:%g'%denom
				
		for gene in self.gt_dmg.iterkeys():
			self.gt_dmg[gene].score /= denom
			
		gt_dmg_min /= denom

		msg += ', done.'
		lib_utils.msgout('notice',msg); self.logger.info(msg)

	def run_vcf2xls(self):
		job_name = 'run_vcf2xls'
		msg = 'converting vcf file to excel file [%s] ...'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)

		if not os.path.exists(self.gene_rank_fn):
			msg = "check if gene rank file [%s] exist"%self.gene_rank_fn
			print(msg); self.logger(msg); RuntimeError(msg)
		
		rank_fn_tmp = self.gene_rank_fn + '.tmp'
			
		cmd = ["cut","-f1,2",self.gene_rank_fn,"|","grep","-v","'#'",">",rank_fn_tmp]
		lib_utils.runcmd2(cmd,self.log_dir,self.logger,"extract_pred_rank")
		
		self.xls = self._assign_out_fn('divine','xls')
		
		cmd = ["python", self.entries['vcf2xls'], \
					"-i", self.vcf, \
					"-o", self.xls, \
					"-l", self.log_dir, \
					"-g", rank_fn_tmp, \
					"-k", self.vknown]
		
		lib_utils.runcmd2(cmd,self.log_dir,self.logger,job_name)
		
		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		os.unlink(rank_fn_tmp)
	
	def rank_pheno_gene(self):
		job_name = 'rank_pheno_gene'
		
		msg = 'selecting genes matched by patient phenotypes ... [%s;%s]'%(job_name,self.hpo_query)
		lib_utils.msgout('notice',msg); self.logger.info(msg)
	
		tmp_fn = '%s.tmp' % self.gene_rank_fn
		fp2=open(tmp_fn,'w')
		fp2.write('#gene\tphenotypic_score\n')
		for gene,cPhenoGene in self.pheno_dmg.iteritems():
			fp2.write('%s\t%g\n'%(gene,cPhenoGene.score))
		fp2.close()
		
		lib_utils.sort_tsv_by_col2(tmp_fn,[2],['gr'],False,self.gene_rank_fn)
		msg = 'done. [%s]'%job_name
		os.unlink(tmp_fn)
		lib_utils.msgout('notice',msg); self.logger.info(msg)
	
	def cleanup(self):
		for fn in self.to_delete_fns:
			if os.path.exists(fn):
				os.unlink(fn)

def main():
	parser = argparse.ArgumentParser(description="Divine (v%s) [author:%s]"%(VERSION,author_email))
	parser.add_argument('-q','--hpo', dest='hpo_query', required=False, default=None, help='Input patient HPO file. A file contains HPO IDs (e.g., HP:0002307), one entry per line. Refer to http://compbio.charite.de/phenomizer, https://hpo.jax.org, or https://mseqdr.org/search_phenotype.php')
	
	parser.add_argument('-v','--vcf', dest='vcf', required=False, default=None, help='input vcf file')
	parser.add_argument('-o','--out_dir', action='store', dest='out_dir', required=False, default=None, help='output directory without white space. If not exist, the directory will be created.')
	
	parser.add_argument('-c','--vcf_filter_cfg', dest='vcf_filter_cfg', required=False, default=None, help='vcf filter configuration file [None]')
	parser.add_argument('-f','--family_fn', dest='ped', required=False, default=None, help='family pedigree file [None]')
	parser.add_argument('-p','--proband_id', dest='proband_id', required=False, default=None, help='proband sample ID [None]')
	
	parser.add_argument('-d','--exp_tag', action='store', dest='exp_tag', required=False, default=None, help='specify experiment tag without white space. The tag will be contained in the output file name.[None] ')
	parser.add_argument('-i','--indel', action='store', dest='indel_fidel', required=False, type=int, default=1, help='the level of fidelity of indell call in VCF, [1]:low (e.g., samtools), 2:high (GATK haplotype caller)')
	
	parser.add_argument('-K', action='store', dest='top_k_disease', required=False, default=0, type = int, help='focus on top-K disease associated with the input phenotypes [0], set 0 to consider all')
	
	parser.add_argument('-r','--go_seed_k', action='store', dest='go_seed_k', required=False, type=float, default=3, help='the number of top-k diseases for GO enrichment [3]; set to 0 to disable')
	
	parser.add_argument('-e','--ref_exon_only', action='store', dest='ref_exon_only', required=False, type=int, default=1, help='the annotation process only runs on RefSeq coding regions 0:No, [1]:Yes')
	parser.add_argument('-C','--cadd', action='store', dest='cadd', required=False, type=int, default=1, help='use CADD prediction score, 0:No, [1]:Yes')
	parser.add_argument('-j','--cosmic', action='store', dest='cosmic', required=False, type=int, default=0, help='enable COSMIC, [0]:No, 1:Yes')
	parser.add_argument('-D','--dblink', action='store', dest='dblink', required=False, type=int, default=0, help='enable dblink, [0]:No, 1:Yes')
	parser.add_argument('-H','--hgmd', action='store', dest='hgmd', required=False, type=int, default=0, help='enable HGMD (requires a license), [0]:No, 1:Yes')
	parser.add_argument('-k','--vknown', action='store', dest='vknown', required=False, type=int, default=1, help='apply variant-level pathogenic annotation (e.g., either ClinVar or HGMD) to prioritization strategy, 0:No, [1]:Yes')

	parser.add_argument('-t', dest='capkit', required=False, default=None, help='capture kit symbol [None],SureSelect_V6,SeqCapEZ_Exome')

	parser.add_argument('--reuse', action='store_const', dest='reuse', required=False, default=False, const=True, help='Reuse previous annotation file (divine.vcf) if it is available [False]')
	
	args = parser.parse_args()
	
	lib_utils.msgout('banner','Divine (v%s) is running on [HPO:%s ,VCF:%s]'%\
									(VERSION,args.hpo_query,args.vcf))
	
	# get a Divine instance/configure program condition and inputs 
	dv = Divine(args)
	
	# analyze phenotype if avaiable
	if dv.hpo_query:
		#to store omim-hpo-gene flat file
		dv.omim = lib_disease.Omim(\
			dv.entries['ext_disease_to_gene'],\
			dv.entries['disease_desc'])
		
		dv.omim.load(dv.logger)
		
		msg = 'analyzing query phenotypes on [%s] ...'%dv.hpo_query
		lib_utils.msgout('notice',msg); dv.logger.info(msg)
	
		# to get disease score associated with the given phenotypes (syscall)
		dv.hpo_to_diseases(args.top_k_disease)

		msg = 'done. [phenotype analysis]'
		lib_utils.msgout('notice',msg); dv.logger.info(msg)

	# analyze genotype if available
	if dv.vcf:
		msg = 'analyzing variants on [%s] ...'%dv.vcf
		lib_utils.msgout('notice',msg); dv.logger.info(msg)

		# to create an instance of varant
		dv.vannotate(args.reuse)
		
		# to apply basic filter/tagging variant class
		dv.vfilter()
		
		# (TODO) to phase
		
		# to run genmod
		if dv.proband_id:
			dv.annotate_comphet_inherit()
		
	# to enrich phenogenes and also predict genetic damage score
	dv.preprocess_dmg_scores()

	# to normalize both pheno and gt dmg scores
	dv.normalize2()

	# to combine two damage scores
	_ = dv.combine_pheno_gt_dmg()

	pagerank.run_heatdiffusion(dv,dv.logger)

	if dv.vcf:
		## to generate an excel file report with a ranking score
		dv.run_vcf2xls()
		msg = 'Done.\nCheck vcf[%s]\nxls[%s]\nranked_gene_fn[%s]\nranked_disease_fn[%s].'%\
			(dv.vcf, dv.xls, dv.gene_rank_fn, dv.disease_rank_fn)
		
	elif args.hpo_query:
		# print ranked phenotype score per gene
		dv.rank_pheno_gene()
		msg = 'Done.\nCheck HPO-to-disease_similarity_fn[%s]\nranked_gene_fn[%s].'%\
			(dv.hpo2disease_fn,dv.gene_rank_fn)
	else:
		lib_utils.msgout('notice', msg, '[WARNING] Nothing to run');
		dv.logger.info(msg)

	lib_utils.msgout('notice',msg,'Divine'); dv.logger.info(msg)
	
	#to cleanup
	dv.cleanup()
	
	lib_utils.msgout('banner','Divine (v%s) is finished for [HPO:%s, VCF:%s].\nContact to %s for any question/error report/feedback'%\
		(VERSION,args.hpo_query,args.vcf,author_email))

if __name__ == '__main__':
	main()
