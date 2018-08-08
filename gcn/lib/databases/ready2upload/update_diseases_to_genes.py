#!/usr/bin/env python

import os, shutil
import argparse
from gcn.lib.databases import lib_hpo, lib_disgenet
VERSION = "0.1.1"
author_email = "hongc2@ccf.org"

def main():

	parser = argparse.ArgumentParser(description="Divine (v%s) [author:%s]" % (VERSION, author_email))
	parser.add_argument('-q', '--hpo', dest='hpo_query', required=False, default=None,
											help='Input patient HPO file. A file contains HPO IDs (e.g., HP:0002307), one entry per line. Refer to http://compbio.charite.de/phenomizer or https://mseqdr.org/search_phenotype.php')

	parser.add_argument('-v', '--vcf', dest='vcf', required=False, default=None, help='input vcf file')
	parser.add_argument('-o', '--out_dir', action='store', dest='out_dir', required=False, default=None,
											help='output directory without white space. If not exist, the directory will be created.')

	#locate two default disease-to-gene file locations (HPO and disgenet)

	hpo_db_dir=os.path.join(os.environ.get("DIVINE"),'gcndata','hpo')
	
	cH = lib_hpo.Hpo(db_dir=hpo_db_dir)
	diseases_to_add = cH.get_diseases_no_gene_assoc()
	
	cDgnet = lib_disgenet.Disgenet()
	umls2did = cDgnet.get_umls_ids(diseases_to_add)
	disGenes = cDgnet.get_genes_by_umls(umls2did)

	flat_ext_fn = os.path.join(hpo_db_dir,'ALL_SOURCES_disgnet_TYPICAL_FEATURES_diseases_to_genes_to_phenotypes.txt')
	
	shutil.copyfile(cH.flat_fn,flat_ext_fn)

	fp1=open(flat_ext_fn,'a')
	for dis,Genes in disGenes.iteritems():
		for gene in Genes:
			for hpo in cH.diseases[dis].hpos:
				fp1.write('%s\t%s\t%s\t%s\t%s\n'%(dis,gene,'-1',hpo,'NA'))
	fp1.close()
	
if __name__ == "__main__":
	main()