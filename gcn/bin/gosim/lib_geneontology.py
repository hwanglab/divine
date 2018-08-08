#!/usr/bin/env python
'''
author: changjin.hong@gmail.com
'''
import os,random
import lib_utils
import itertools

class Fastsemsim:
	def __init__(self,resource_dir,work_dir='/tmp',go_bin=os.path.join(os.environ['DIVINE'],'python_libs','bin','fastsemsim'),debug=False):
		self.resource_dir = resource_dir
		
		self.obo_fn = '%s/go.obo'%self.resource_dir
		if not lib_utils.check_if_file_valid(self.obo_fn):
			lib_utils.gen_msg_time('error','check if the file[%s] exists'%self.obo_fn,'lib_geneontology.py')
		
		self.goa_fn = '%s/gene_association.goa_human'%self.resource_dir
		if not lib_utils.check_if_file_valid(self.goa_fn):
			lib_utils.gen_msg_time('error','check if the file[%s] exists'%self.goa_fn,'lib_geneontology.py')
			
		self.work_dir = work_dir
		
		#rfix = random.randint(1,1e7)
		#rfix = 0
		
		#to prepare input_fn
		self.query_fn = None
		
		self.gobin = go_bin
		
		self.go_categories = ['biological_process','molecular_function','cellular_component']
		self.gosim_fns = {}
		self.uniprot2RefGenes = {}
		self.debug = debug
	
	def map_uniprot_to_refgene(self):
		tmp_fn = '%s_tmp%d'%(self.query_fn,random.randint(1,1e7))
		cmd = 'cut -f1 %s | sort | uniq > %s'%(self.query_fn,tmp_fn)
		lib_utils.gen_msg_time('cmd',cmd,'map_uniprot_to_refgene')
		Query = lib_utils.file2list(tmp_fn,False)
		os.unlink(tmp_fn)
		
		tmp_fn = '%s_tmp%d'%(self.query_fn,random.randint(1,1e7))
		cmd = 'cut -f2 %s | sort | uniq > %s'%(self.query_fn,tmp_fn)
		lib_utils.gen_msg_time('cmd',cmd,'map_uniprot_to_refgene')
		Query.extend(lib_utils.file2list(tmp_fn,False))
		os.unlink(tmp_fn)

		fp = lib_utils.open2(self.goa_fn,'r')
		for i in fp:
			if i[0]=='!':continue
			j=i.split('\t')
			refGene = j[2]
			uniprot = j[1]
			if uniprot in Query:
				if refGene not in self.uniprot2RefGenes:
					self.uniprot2RefGenes[uniprot]=refGene
		fp.close()
	
	def list_to_query_fn(self,Query,Target,to_uniprot=True):
		'''
		convert refSeqGene to uniprot and generate a pair of uniprot ids into a text file
		'''
		Genes=[Query,Target]
		if to_uniprot:
			Uniprots=[[],[]]

			fp = lib_utils.open2(self.goa_fn,'r')
			for i in fp:
				if i[0]=='!':continue
				j=i.split('\t')
				refGene = j[2]
				uniprot = j[1]
				for gi in range(2):
					if refGene in Genes[gi]:
						if uniprot not in Uniprots[gi]:
							Uniprots[gi].append(uniprot)
							if refGene not in self.uniprot2RefGenes:
								self.uniprot2RefGenes[uniprot]=refGene
			fp.close()
		else:
			Uniprots=[Query,Target]

		fp2 = lib_utils.open2(self.query_fn,'w')
		for query,target in list(itertools.product(Uniprots[0],Uniprots[1])):
			fp2.write('%s\t%s\n'%(query,target))
		fp2.close()
		return self.query_fn
	
	def run_sim(self,out_dir,sim_method='SimRel'):
		
		for root_node in self.go_categories:
			gocmd = self.gobin
			gocmd += " --ontology_file %s"%self.obo_fn
			gocmd += " --ontology_file_format obo"
			gocmd += " --ac_file %s"%self.goa_fn
			gocmd += " --ac_sep '\\t'"
			gocmd += " --task SS"
			gocmd += " -o GeneOntology"
			
			if root_node not in self.go_categories:
				lib_utils.gen_msg_time('error','%s not defined'%root_node,'lib_geneontology.py')
			gocmd += " --root %s"%root_node
			
			gocmd += " --ac fly"
			gocmd += " --tss %s"%sim_method
			gocmd += " --tmix max"
			gocmd += " --remove_nan"
			gocmd += " --ignore_EC TAS"
			gocmd += " --query_ss_type obj"
			gocmd += " --query_mode pairs"
			gocmd += " --query_input file"
			gocmd += " --query_file %s" % self.query_fn
			
			gosim_fn = "%s/%s_%s.tsv" % (out_dir,self.query_fn.split('/')[-1],root_node)
			gocmd += " --output_file %s"%gosim_fn
			
			if self.debug and lib_utils.check_if_file_valid(gosim_fn):
				lib_utils.gen_msg_time('notice','reuse %s'%gosim_fn,'lib_geneontology.run_sim')
			else:
				lib_utils.gen_msg_time('cmd',gocmd,'lib_geneontology.run_sim')
				#print 'reuse...' #debug

			self.gosim_fns[root_node]=gosim_fn
	
	def get_funsim_mat(self,out_fn,min_sim_score=1.):
		Pairs = {}
		Pairs_cnt = {}
		for go_cat, gosim_fn in self.gosim_fns.iteritems():
			fp = lib_utils.open2(gosim_fn,'r')
			fp.next() #skip head
			for i in fp:
				uQuery,uTarget,score = i.split('\t')
				score = score.strip()
				if score:
					upair = (uQuery,uTarget)
					if upair not in Pairs:
						Pairs[upair] = 0.
						Pairs_cnt[upair] = 0.
					Pairs[upair]+=float(score)
					Pairs_cnt[upair]+=1.
			fp.close()
		
		goSimScores = {}
		fp2 = lib_utils.open2(out_fn,'w')
		for upair,score in Pairs.iteritems():
			if Pairs_cnt[upair]>0:
				score /= Pairs_cnt[upair]
				fp2.write('%s\t%s\t%g\t%d\n'%(upair[0],upair[1],score,Pairs_cnt[upair]))
			
				if score>min_sim_score:
					#convert to refGene (HGNC) symbol
					pair=(self.uniprot2RefGenes[upair[0]],self.uniprot2RefGenes[upair[1]])
					goSimScores[pair]=score

		fp2.close()
		return goSimScores
	
	def gen_all_pairs_file(self,upair_file):
		uniprot_list_fn = '%s/go_%d.tmp'%(self.work_dir,random.randint(1,1e7))
		cmd = "cut -f2 %s | grep -vP '^!' | sort | uniq > %s"%(self.goa_fn,uniprot_list_fn)
		lib_utils.gen_msg_time('cmd',cmd,'lib_geneontology.run_sim')
		
		Uniprot = lib_utils.file2list(uniprot_list_fn)
		fp2 = lib_utils.open2(upair_file,'w')
		U = len(Uniprot)
		print 'total number of queries[%d]'%U
		for i,u1 in enumerate(Uniprot):
			print '[%d/%d]'%(i,U)
			for u2 in Uniprot[i+1:]:
				fp2.write('%s\t%s\n'%(u1,u2))
		fp2.close()
		print 'done.'