#!/usr/bin/env python
'''
author: changjin.hong@gmail.com
'''
import os
import datetime, time
from gcn.etc.fileconfig import FILECONFIG
from gcn.lib.io import anyopen
from gcn.lib.utils import lib_utils
from gcn.bin.hgConvert.hgvs_resource import Hgvs2
from collections import namedtuple
import csv

class ClinVitae:
	def __init__(self,out_dir):
		if out_dir is None:
			raise RuntimeError('specify output directory!')
		self.out_dir=out_dir
		self.tsv = os.path.join(self.out_dir,'clinvitae.tsv')
		if os.path.exists(FILECONFIG['CLINVITAE']):
			self.vcf = FILECONFIG['CLINVITAE']

	def download(self):
		url="http://s3-us-west-2.amazonaws.com/clinvitae/clinvitae_download.zip"

		cmd = 'cd %s; wget -c %s'%(self.out_dir,url)
		lib_utils.runcmd(cmd,call_from='ClinVitae.download()')
		
		fnz = '%s/clinvitae_download.zip'%self.out_dir
		cmd = 'cd %s; unzip %s'%(self.out_dir,fnz)
		print cmd
		lib_utils.runcmd(cmd,call_from='ClinVitae.download()')
		fn = '%s/variant_results.tsv'%self.out_dir
		print fn
		self.tsv = '%s/clinvitae.tsv'%self.out_dir
		print self.tsv
		os.rename(fn,self.tsv)
		os.unlink(fnz)
		
	def _write_vcf_head(self,fp2):
		import datetime
		
		fp2.write('##fileformat=VCFv4.2\n')
		
		ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d')
		fp2.write('##fileDate=%s\n'%ts)
		
		source_url='http://clinvitae.invitae.com/'
		fp2.write('##source=%s\n'%source_url)
		fp2.write('##reference=GRCh37\n')
		fp2.write('##INFO=<ID=cDNA,Number=.,Type=String,Description="nt change in HGVS format">\n')
		fp2.write('##INFO=<ID=VC,Number=.,Type=String,Description="reported clinical significance classification">\n')
		fp2.write('##INFO=<ID=SRC,Number=.,Type=String,Description="source">\n')
		fp2.write('##INFO=<ID=UPD,Number=.,Type=String,Description="date updated lastly">\n')
		fp2.write('##INFO=<ID=URL,Number=.,Type=String,Description="URL of the interpretation source">\n')
		fp2.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

	def _iterfile(self):
		fields = "gene,nt_change,pt_change,other_mapping,alias,tx,region,rep_class,inf_class,source,last_eval,last_upd,url,comment".split(",")
		F = len(fields)
		clinvt = namedtuple('clinvt',fields)
		
		stream = anyopen.openfile(self.tsv, 'rt')
		stream.next()
		for rec in csv.reader(stream,delimiter='\t'):
			if len(rec)==F:
				yield clinvt._make(rec)
		stream.close()
		
	def may_pass(self,rec):
		if 'clinvar' in rec.source.lower() or \
		'clinvar' in rec.url.lower() or \
		':' not in rec.nt_change or \
		not rec.tx.strip() or \
		rec.tx == '-':
			return True
		else:
			return False
	
	def decompose_cdna(self,nt_change_str):
		tx_cdnas = []
		if ';' in nt_change_str:
			tx,nt_change = nt_change_str.split(';')
			tx = tx.strip()
			cdnas = nt_change.split(';')
			for cdna in cdnas:
				tx_cdnas.append('%s:%s'%(tx,cdna.strip()))
		elif ',' in nt_change_str:
			tx,nt_change = nt_change_str.split(',')
			tx = tx.strip()
			cdnas = nt_change.split(',')
			for cdna in cdnas:
				tx_cdnas.append('%s:%s'%(tx,cdna.strip()))
		else:
			tx_cdnas.append(nt_change_str)
		return tx_cdnas
		
	def extract_cDNA(self):
		cache_tx = {}
		self.cDNA_fn = lib_utils.file_tag(self.tsv,'cDNA',None)
		fpw = open(self.cDNA_fn,'w')
		for cvt in self._iterfile():
			if self.may_pass(cvt): continue
			for tx_cdna in self.decompose_cdna(cvt.nt_change):
				fpw.write('%s\n'%tx_cdna)
		fpw.close()
		
		print 'upload Clinvitae variant file [%s] in cDNA HGVS format to [%s] to obtain chromosome contig bp location' % (self.cDNA_fn,'https://mutalyzer.nl/position-converter')
		
		return self.cDNA_fn

	def gdna_to_vcf(self,mutalyzer_batch_outfn):
		
		if not os.path.exists(mutalyzer_batch_outfn):
			raise RuntimeError('check if input file [%s] exists'%\
												mutalyzer_batch_outfn)

		cHgvs = Hgvs2()
		cHgvs.load_resource()
		
		fp = open(mutalyzer_batch_outfn,'r')
		fp.next()
		gdna_cache = {}
		for mutalyzer in fp:
			mut = mutalyzer.split('\t')
			if mut[1].strip(): continue
			gdna = mut[2].strip()
			variants = cHgvs.gdna_to_vcf(gdna)
			if variants:
				gdna_cache[mut[0].strip()] = variants
		fp.close()

		self.out_vcf = lib_utils.file_tag(self.tsv,None,'vcf')
		
		tmp_vcf = self.out_vcf+'.tmp'
		fpw = open(tmp_vcf,'w')
		self._write_vcf_head(fpw)
		qual = 100
		filter = 'PASS'
		rsid = '.'

		for cvt in self._iterfile():
			if self.may_pass(cvt): continue
			if cvt.nt_change not in gdna_cache: continue
			
			for chrom,pos,ref,alt in gdna_cache[cvt.nt_change]:
				if len(ref)>100 or len(alt)>100:continue
				info = 'cDNA=%s;'%cvt.nt_change
				info +='VC=%s;'%self.determine_vclass(cvt.rep_class)
				info +='SRC=%s;'%cvt.source
				info +='UPD=%s;'%cvt.last_upd
				info +='URL=%s'%cvt.url
				if chrom.startswith('chr'):
					if chrom.startswith('chrM'): chrom = 'MT'
					else: chrom = chrom[3:] 
				cols = [chrom,pos,rsid,ref,alt,qual,filter,info]
				fpw.write('%s\n' % lib_utils.joined(cols,'\t'))
		fpw.close()
		lib_utils.sort_tsv_by_col2(tmp_vcf,[1,2],\
			['V','n'],False,self.out_vcf)
		os.unlink(tmp_vcf)

	def determine_vclass(self,rep_class):
		vc_in_db = rep_class.lower()
		if 'conflicting' in vc_in_db:
			cvt_vc = 'vus'
		elif 'benign' in vc_in_db:
			cvt_vc = 'benign'
		elif 'pathogenic' in vc_in_db:
			cvt_vc = 'pathogenic'
		else:
			cvt_vc = 'vus'
		return cvt_vc
	
if __name__ == '__main__':
	print 'to add unittest'

