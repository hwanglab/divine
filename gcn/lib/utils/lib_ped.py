#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''

from ped_parser import FamilyParser
from gcn.lib.utils import lib_utils, filter_cj
import gcn.lib.io.vcf as vcf

UNKNOWN,MALE,FEMALE = range(3)
HEALTHY,AFFECTED = range(1,3)
	
def check_consistency_ped_vcf(vcf,family_file,proband_id,family_type='ped',logger=None):
	
	family_parser = FamilyParser(family_file, family_type)

	family_ids = family_parser.families.keys()
	F = len(family_ids)
	if F==1:
		family_obj = family_parser.families[family_ids[0]]
		ind_ids = family_obj.individuals.keys()
		if proband_id not in ind_ids:
			msg = "the proband_id [%s] cannot be found in the pedigree file [%s]" %(proband_id,family_file)
			print msg
			if logger: logger.error(msg)
			raise RuntimeError(msg)
		v = vcf.VCFParser(vcf)
		
		for ind_id in ind_ids:
			if ind_id not in v.samples:
				msg = "individual name [%s] in the pedigree file [%s] does not exist in VCF [%s]" % \
					(ind_id,family_file,vcf)
				print msg
				if logger: logger.error(msg)
				raise RuntimeError(msg)

		return v.samples.index(proband_id)
	
	else:
		msg = "check the pedigree file [%s] you provide and it must have only one family ID" % \
			(family_file)
		print msg
		if logger: logger.error(msg)
		raise RuntimeError(msg)

def create_proband_ped(single_vcf,out_ped,gender=None,phenotype=2):
	job_name = "create_proband_ped"
	msg = "opening VCF [%s] and parse heads ..." % single_vcf
	print msg
	
	v = vcf.VCFParser(single_vcf)
	S = 0
	for rec in v:
		S = len(v.samples)
		break
	
	if S == 1:
		sample_id = v.samples[0]
		if not gender:
			gender = filter_cj.predict_gender_from_VCF(\
								single_vcf,sample_id)
		else:
			gender = UNKNOWN
	else:
		raise RuntimeError('specify sample id to work with in VCF [%s] or provide a VCF file containing only one sample ID'%single_vcf)
	
	v.stream.close()
	
	msg = "create a ped file [%s] for VCF [%s]" % (out_ped,single_vcf)
	fpw = open(out_ped,'w')
	head = "FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype"
	fpw.write('#%s\n' % head)
	recordStr = "1\t%s\t0\t0\t%d\t%d" % (sample_id,gender,phenotype)
	fpw.write('%s\n'%recordStr)
	fpw.close()
	
	msg = "Sample ID [%s] is identified for a proband analysis!"%sample_id
	print msg
	
	return sample_id

def parse_genmod_inherit_model(genmod_tag):
	
	tag = ''
	if 'hom' in genmod_tag:
		tag += 'H'
	elif 'comp' in genmod_tag:
		tag += 'h'
	
	if 'dn' in genmod_tag:
		tag += 'd'
	return tag
