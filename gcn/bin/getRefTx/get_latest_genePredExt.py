#!/usr/bin/env python

'''

http://ftp.ensembl.org/pub/grch37/release-85/gtf/homo_sapiens/Homo_sapiens.GRCh37.85.gtf.gz

SELECT concat(transcript.stable_id,'.',cast(transcript.version as char(2))) as ensembl_tx, xref.display_label
FROM transcript, object_xref, xref, external_db
WHERE transcript.transcript_id = object_xref.ensembl_id
AND object_xref.ensembl_object_type = 'Transcript'
AND object_xref.xref_id = xref.xref_id
AND xref.external_db_id = external_db.external_db_id
AND external_db.db_name = 'RefSeq_mRNA';

http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
'''
import argparse, os, sys

def map_enst2mrna(enst_to_ref_fn):
	fp = open(enst_to_ref_fn,'r')
	fp.next() #skip header
	enst2mrna = {}
	enstv2mrna = {}
	
	for i in fp:
		enstv,mrna = i.rstrip().split('\t')
		enstv2mrna[enstv] = mrna
		
		enst,version2 = enstv.split('.')
		if enst not in enst2mrna:
			enst2mrna[enst] = []
		enst2mrna[enst].append([int(version2),mrna])
		
	fp.close()
	return enstv2mrna, enst2mrna

def replace_tx_by_mrna(gp_in,enstv2mrna,enst2mrna,gp_out):
	
	fp = open(gp_in,'r')
	fpw= open(gp_out,'w')
	for i in fp:
		itms = i.split('\t')
		enstv = itms[0]
		
		printKey = False
		if enstv in enstv2mrna:
			itms[0] = enstv2mrna[enstv]
			printKey = True
		else:
			enst,version2 = enstv.split('.')
			if enst in enst2mrna:
				max_ver = 0
				for ver2,mrna in enst2mrna[enst]:
					if ver2>max_ver:
						approx_mrna = mrna
				itms[0] = approx_mrna
				printKey = True

		if printKey:
			fpw.write('%s'%'\t'.join(itms))

	fp.close()
	fpw.close()

def main():

	parser = argparse.ArgumentParser(description="get latest RefGene transcript coordinates in genePredExt format")
	parser.add_argument('-i', action='store', dest='gp_in', required=True, help='')
	parser.add_argument('-d', action='store', dest='enst_to_ref_fn', required=True, help='')
	parser.add_argument('-o', action='store', dest='gp_out', required=True, help='')
	args=parser.parse_args()
	
	#read dictionary
	enstv2mrna,enst2mrna = map_enst2mrna(args.enst_to_ref_fn)
	
	#replace enst by mrna
	replace_tx_by_mrna(args.gp_in,enstv2mrna,enst2mrna,args.gp_out)

if __name__ == '__main__':
	main()