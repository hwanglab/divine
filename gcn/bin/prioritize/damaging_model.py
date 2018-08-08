#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''
import math
import scipy.stats as st

class DmgCoeff:
	
	def __init__(self,indel_fidel=1,go_seed_k=3,logger=None):
		msg='setting up damaging factors ...'
		print msg
		if logger: logger.info(msg)
		
		self.vncoding = 0.050
		self.vintronic = 0.055
		self.vexonic = 0.35
		self.vsplice = 0.45
		self.vsplice_syn = 0.060
		self.warning = 0.070
		
		# damage factor for known pathogenic variants 
		self.cgene = 0.04
		self.cintronic = 0.25
		self.cexonic = 0.50
		
		# logistic regression coeff for damage w.r.t MAF
		self.beta1= 0.0125
		self.beta2 = -0.093
		self.beta3 = -0.298
		self.maf_offset_max = 0.8
		self.maf_offset_min = 0.1

		# hetero/homozygous/compound het allele impact for rare SNVs
		self.het_rare = 0.10
		self.het2_rare = 0.20
		self.hom_rare = 0.75

		self.het_denovo_snp = 0.40
		self.het2_denovo_snp = 0.50
		self.hom_denovo_snp = 0.85
		self.cmp_mnp_penalty = 0.5
		
		# hetero/homozyous/compound het allele impact w.r.t. the mutation type for rare/denovo SNVs
		
		if indel_fidel==1:
			self.het_denovo_indel = 0.10
			self.het2_denovo_indel = 0.15
			self.hom_denovo_indel = 0.25

			# hetero/homozygous/compound het allele impact for indel SNVs
			self.het_rare_indel = 0.20
			self.het2_rare_indel = 0.25
			self.hom_rare_indel = 0.30

		elif indel_fidel==2:
			self.het_denovo_indel = 0.40
			self.het2_denovo_indel = 0.50
			self.hom_denovo_indel = 0.85

			# hetero/homozygous/compound het allele impact for indel SNVs
			self.het_rare_indel = 0.45
			self.het2_rare_indel = 0.55
			self.hom_rare_indel = 0.85

		# default average transcript length when its length is unknown
		# https://www.quora.com/Whats -the-average-size-of-a-human-protein-in-kDa
		self.avg_protein_len = 480
		self.min_dmg_prior = 1.25
		self.go_seed_k = go_seed_k
		self.go_penalty = 0.4
		self.gosim_min = 0.95
		self.ptwt = 0.5
		#self.prwt = 0.3
		self.prwt = 0.5
		
		self.conf_pheno_rank = 100
		self.conf_pheno_rank_dom = 100

	def get_maf_xoffset(self,maf):
		
		#return (1. - self.beta1 * math.exp(1000.*maf)) / self.beta2
		if maf > 0.:
			mafL = math.log10(maf)
			offset = self.beta1 * math.pow(mafL,2) + self.beta2 * mafL + self.beta3
			if offset < 0.:
				offset = self.maf_offset_min
			elif offset > self.maf_offset_max:
				offset = self.maf_offset_max
		else:
			offset = self.maf_offset_max
		return offset

def get_z(confidence):
	return st.norm.ppf(1 - (1 - confidence) / 2)

def ci_lower_bound(phat, n, z=None,confidence=None):

	if n == 0.:
		return 0
	if z is None:
		z = get_z(confidence)
	#Binomial proportion confidence using Wilson score interval
	return (phat + z * z / (2 * n) - z * math.sqrt((phat * (1 - phat) + z * z / (4 * n)) / n)) / (1 + z * z / n)
