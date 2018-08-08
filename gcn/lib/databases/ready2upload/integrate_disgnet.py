#!/usr/bin/env python

import os, shutil
from gcn.lib.databases import lib_hpo, lib_disgenet

def main():
	
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