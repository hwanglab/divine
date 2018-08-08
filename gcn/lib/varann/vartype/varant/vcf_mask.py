"""
.. module:: annotpipe
    :platform: Unix, Windows, MacOSX
    :synopsis: A wraper to call VARANT

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is a wrapper to call VARANT. The inputs are -
1. Unannotated VCF file path
2. Path to create annotated vcf file (Option)
"""
import numpy as np
from gcn.lib.databases.refgene import get_ucsc_chrom, to_ucsc_chrom
from gcn.lib.utils import lib_utils
from gcn.lib.varann.vartype.varant import annotateRegion

def vcf_get_max_sv_len(ref,alt):
	
	#looking for a deletion
	alts = alt.split(',')
	altSvLen = 1e9
	for alt1 in alts:
		if alt1=='.': alen = 0
		else: alen = len(alt1)
		
		if alen<altSvLen: altSvLen = alen
			
	refs = ref.split(',')
	refSvLen = 0
	for ref1 in refs:
		if ref1 == '.': rlen = 0
		else: rlen = len(ref1)
		
		if rlen>refSvLen: refSvLen = rlen
	
	offset = 0
	if refSvLen>altSvLen: #deletion
		offset = refSvLen - altSvLen

	return offset
	
class BedMaskingVCF():
	def __init__(self,bed_fn,logger=None):
		self.bed = bed_fn
		self.logger=logger

	def run(self,vcf_fn,masked_vcf_fn):
		
		job_name = 'BedMaskingVCF.run'
		cRef = annotateRegion.RefGeneUcscTB(logger=self.logger)
		cRef.bed_fn = self.bed
		eBnds = cRef.get_boundary()

		msg = 'masking the vcf file [%s] by the bed file [%s] @ %s'%(vcf_fn,cRef.bed_fn,job_name)
		lib_utils.msgout('notice',msg)
		if self.logger: self.logger.info(msg)
		
		fp = open(vcf_fn,'r')
		fp2 = open(masked_vcf_fn,'w')
		cmd_head = '##original_vcf=%s\n##masking_bed=%s'%(vcf_fn,self.bed)
		head_written = False

		for i in fp:
			if i[0]=='#':
				if i.startswith('##contig'):
					if not head_written:
						fp2.write('%s\n'%cmd_head)
						head_written = True
				fp2.write('%s'%i)
			else:
				j = i.split('\t')
				chrom = j[0]
				chrom = to_ucsc_chrom(chrom)
				if chrom in eBnds:
					pos1 = int(j[1])-1 #adjust VCF to BED coordinate
					idx = np.nonzero((eBnds[chrom][:,0]<=pos1) & (pos1<eBnds[chrom][:,1]))
					if idx[0].size>0:
						fp2.write('%s'%i)
					else:
						pos2 = pos1 + vcf_get_max_sv_len(j[3],j[4])
						if pos2>pos1:
							idx = np.nonzero((eBnds[chrom][:,0]<=pos2) & (pos2<eBnds[chrom][:,1]))
							if idx[0].size>0:
								fp2.write('%s'%i)
		fp.close()
		fp2.close()
		msg = 'done. @ %s'%job_name
		lib_utils.msgout('notice',msg)
		if self.logger: self.logger.info(msg)

def by_bed(vcf_fn,bed_fn,masked_vcf_fn,logger=None):
	cM = BedMaskingVCF(bed_fn,logger=logger)
	cM.run(vcf_fn,masked_vcf_fn)
	return masked_vcf_fn
