#!/usr/bin/env python
'''
author: changjin.hong@gmail.com
'''
import os, argparse
from gcn.config import lib_config
import lib_clinvitae

def main():
	"""to generate Clinvitae VCF file"""

	desc = 'Generating Clinvitae VCF file'
	parser = argparse.ArgumentParser(description=desc)
	
	parser.add_argument('--preprocess', action='store_const', dest='preprocess',\
		 required=False, default=False, const=True, help='download Clinvitae')

	parser.add_argument('-g', dest='mutalyzer_fn', required=False,\
		default=None, help='a file generated from https://mutalyzer.nl/position-converter')
	
	parser.add_argument('-o', '--out_dir', dest='out_dir', required=False,\
		type=str, default = lib_config.gcn_path('GCN_DATA_DIR'),\
			help='output directory')
	
	args = parser.parse_args()
	
	if args.preprocess:
		clinvt = lib_clinvitae.ClinVitae(args.out_dir)
		clinvt.download()
		_ = clinvt.extract_cDNA()
	elif args.mutalyzer_fn:
		clinvt = lib_clinvitae.ClinVitae(args.out_dir)
		if os.path.exists(clinvt.tsv) and os.path.exists(args.mutalyzer_fn):
			clinvt.gdna_to_vcf(args.mutalyzer_fn)
		else:
			raise RuntimeError('check if both [%s] and [%s] exists'%\
												(clinvt.tsv,args.mutalyzer_fn))
	else:
		raise RuntimeError('it does not do anything!')
	
if __name__ == '__main__':
	main()