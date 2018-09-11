#!/usr/bin/env python
'''
author: changjin.hong@gmail.com
'''
import datetime, time
import os, argparse
import gcn.lib.io.vcf as vcf
from gcn.lib.utils import lib_utils

NONCODING='0'
CODING = '1'

class Cosmic:
	def __init__(self,version=None):
		self.release_ver = version
	
	def reformat_to_lite(self,infile,vtype,outfile,min_cnt = 1):
		jobname = "reformat_to_lite"
		msg = "working on vcf file [%s] ..." % infile
		print msg
		
		infoKeys = ['GENE','STRAND','CDS','AA','SNP']
		
		v = vcf.VCFParser(infile)
		
		v.add_meta_info('COSMIC_ID','.','String','cosmic ID')
		v.add_meta_info('REG','2','Integer','1:coding, 0:noncoding')
		
		if not os.path.exists(outfile):
			ostream = open(outfile, 'w')
			v.writeheader(ostream,to_del_info = infoKeys)
		else:
			ostream = open(outfile, 'a')

		pk0 = 'NA'
		cosmics = []
		cnts = []
		prev_rec = None
		
		for rec in v:
			v.parseinfo(rec)
			pk = lib_utils.joined([rec.chrom,rec.pos,rec.ref,rec.alt],'_')
			
			if pk!=pk0:
				pk0 = pk
				if prev_rec:
					prev_rec.id = '.'
					prev_rec.info['COSMIC_ID'] = cosmics
					prev_rec.info['REG'] = vtype
					for info_key in infoKeys:
						v.delete_info(prev_rec,info_key)
					if vtype == NONCODING:
						prev_rec.info['CNT'] = '1'
					v.write(ostream, prev_rec)
					
				cosmics = [rec.id[0]]
				prev_rec = rec
			else:
				pk0 = pk
				cosmics.append(rec.id[0])
				
		if prev_rec:
			prev_rec.id = '.'
			prev_rec.info['COSMIC_ID'] = cosmics
			prev_rec.info['REG'] = vtype
			for info_key in infoKeys:
				v.delete_info(prev_rec,info_key)
			
			if vtype == NONCODING:
				prev_rec.info['CNT'] = '1'
			v.write(ostream, prev_rec)

		ostream.close()
		v.stream.close()
		
		msg = "Done [%s]." % jobname
		print msg

def main():
	"""to reformat cosmic VCF file to transfer to database model"""

	desc = 'Reformatting cosmic VCF file'
	parser = argparse.ArgumentParser(description=desc)
	
	parser.add_argument('-C', dest='cosmic_coding_vcf', required=True,\
		default=None, help='cosmic coding VCF file')
	parser.add_argument('-c', dest='cosmic_noncoding_vcf', required=True,\
		default=None, help='cosmic coding VCF file')
	parser.add_argument('-i', dest='min_cnt_ncoding', required=False,\
		default=2, help='minimum number of CNT to include non-coding cosmic varianat sites [2]')
	parser.add_argument('-o', dest='cosmic_vcf4db', required=True,\
		default=None, help='output cosmic VCF file to upload database')

	args = parser.parse_args()
	
	cC = Cosmic()
	
	if os.path.exists(args.cosmic_vcf4db):
		os.unlink(args.cosmic_vcf4db)
		
	cC.reformat_to_lite(args.cosmic_coding_vcf,CODING,\
										args.cosmic_vcf4db)
	
	cC.reformat_to_lite(args.cosmic_noncoding_vcf,NONCODING,\
										args.cosmic_vcf4db,min_cnt = args.min_cnt_ncoding)

if __name__ == '__main__':
	main()
	
	"""
	./reformat_cosmic,py
	-C gcndata/cosmic/CosmicCodingMuts.vcf.gz
	-c gcndata/cosmic/CosmicNonCodingVariants.vcf.gz
	-i 2
	-o gcndata/cosmic/cosmic.vcf
	"""
	