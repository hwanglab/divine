#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan and modified by changjin.hong@gmail.com
#
"""
.. module:: anyopen
		:platform: Unix, Windows, MacOSX
		:synopsis: Transparent opening of compressed and uncompressed files

.. moduleauthor:: ; changjin.hong@gmail.com

select
	am.acc_num as 'hgmd_accession_id',
	coalesce(hg19.chromosome,'') AS 'chromosome',
	coalesce(hg19.coordSTART,'') AS 'hg19_start',
	coalesce(hg19.coordEND,'') AS 'hg19_end',
	coalesce(hg19.strand,'') AS 'hg19_strand',
	coalesce(mtn.wildBASE,'') AS 'ref',
	coalesce(mtn.mutBASE,'') AS 'alt',
	coalesce(am.tag,'') as 'hgmd_classification',
	coalesce(am.omimid,'') as 'omim_id',
	coalesce(ag.mut_total,'') as 'count_reported_mutations_per_gene',
	coalesce(ag.entrezid,'') as 'entrez_id',
	coalesce(ag.gene,'') as 'gene_symbol',
	-- coalesce(ag.altsymbol,'') as 'alt_gene_symbol',
	coalesce(ag.genename,'') as 'gene_definition',
	coalesce(ag.refseq,'') as 'refseq_tx_id',
	coalesce(concat('c.', am.hgvs),'') AS 'hgvs_cdna',
	coalesce(concat(mtn.protCORE,'.', mtn.protver),'') as 'refseq_prot_id',
	coalesce(am.descr,'') as 'aa_change_desc',
	coalesce(am.dbsnp,'') as 'rsid',
	coalesce((select group_concat(hp.phenotype separator ', ')
		from hgmd_phenbase.hgmd_mutation as hm
		inner join hgmd_phenbase.hgmd_phenotype as hp on hm.phen_id = hp.phen_id
		where am.acc_num = hm.acc_num),'') as 'phenotype',
	coalesce(am.pmid,'') as 'primary_pmid',
	coalesce(am.new_date,'') as 'appear_date',
	coalesce((select group_concat(er.pmid separator ', ')
		from hgmd_pro.extrarefs as er
		where am.acc_num = er.acc_num),'') as 'secondary_pmid',
	coalesce(fsite.name,'') as 'func_site_name',
	coalesce(fsite.type,'') as 'func_site_type',
	coalesce(sp.location,'') as 'splice_loc'
from hgmd_pro.allmut as am
	left outer join hgmd_pro.allgenes as ag on ag.gene = am.gene
	left outer join hgmd_pro.mutnomen as mtn on am.acc_num = mtn.acc_num
		left outer join hgmd_pro.splice as sp on am.acc_num = sp.acc_num
		left outer join hgmd_pro.hg19_coords as hg19 on am.acc_num = hg19.acc_num
		left outer join hgmd_pro.func_anotat as fant on am.acc_num = fant.acc_num
		left outer join hgmd_pro.func_sites as fsite on fant.site_id = fsite.id
		order by hg19.chromosome, hg19.coordSTART, hg19.coordEND
	INTO OUTFILE '/var/lib/mysql-files/hgmd_pro_v2016.4.tsv' FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';
"""

import re

BENIGN, VUS, PATHOGENIC = range(-1, 2)
FP, UNKNOWN, SAME_GENE, PREDICTED, LIKELY_DM, DM = range(-1, 5)
NA, GENE, POS, ALT, OVR = range(5)

def determine_vclass(hgmd_vc,vc_found):
	'''
	there are 6 different hgmd class
	DFP (disease-associated polymorphisms with additional supporting functional evidence. These are reported to be in significant association with disease (p<0.05) and to have evidence for being of direct functional importance (e.g. as a consequence of altered expression, mRNA studies etc))
	DM (disease-causing mutations, pathological mutations reported to be disease causing in the original literature report. The other three tags are used for polymorphisms. )
	DM? (likely DM)
	DP (disease-associated polymorphisms. These are reported to be in significant association with disease (p<0.05) and are assumed to be functional (e.g. as a consequence of location, evolutionary conservation, replication studies etc), although there may be as yet no direct evidence (e.g. from an expression study) of function)
	FP (in vitro/laboratory or in vivo functional polymorphisms. These are reported to affect the structure, function or expression of the gene (or gene product), but with no disease association reported as yet)
	R (Retired entry)
	'''

	hgmd_pos = ['DM']
	hgmd_likely_pos = ['DM?']
	hgmd_predicted_pos = ['DFP', 'DP']
	hgmd_neg = ['FP', 'R']

	vc = hgmd_vc.split('|')[1]
	if vc in hgmd_pos:
		vc_found = DM
	elif vc in hgmd_likely_pos:
		if vc_found < LIKELY_DM:
			vc_found = LIKELY_DM
	elif vc in hgmd_predicted_pos:
		if vc_found < PREDICTED:
			vc_found = PREDICTED
	elif vc in hgmd_neg:
		if vc_found == UNKNOWN:
			vc_found = FP
	return vc_found

def concat_hgmd_desc(rec):
	NotAvail = '.'
	tmp = NotAvail
	if rec.HGMD_ACC: tmp = rec.HGMD_ACC
	hgmd_desc = '%s' % tmp

	tmp = NotAvail
	if rec.VC: tmp = rec.VC
	hgmd_desc += '|%s' % tmp

	tmp = NotAvail
	if rec.PMID: tmp = '_'.join(rec.PMID.split(','))
	hgmd_desc += '|%s' % tmp

	tmp = NotAvail
	if rec.PHENO: tmp = rec.PHENO.replace(',', '.')
	hgmd_desc += '|%s' % tmp

	return hgmd_desc

def get_hgmd_class(hgmd_descriptions):

	hgmd_match_types = [POS, ALT]
	vc_found = UNKNOWN
	for var in hgmd_descriptions:
		mObj = re.search(r'(\d+)\((.+)\)', var)
		if mObj:
			hgmd_match = int(mObj.group(1))
			if hgmd_match in hgmd_match_types:
				vc_found = determine_vclass(mObj.group(2),vc_found)
				if vc_found == DM:
					break
		else:
			print 'check HGMD_DESC[%s]' % var

	return vc_found

def pathogenic_gene(hgmd_descriptions):
	has_gene = False
	for var in hgmd_descriptions:
		# print var #debug
		mObj = re.search(r'(\d+)\((.+)\)', var)
		if mObj:
			hgmd_match = int(mObj.group(1))
			if hgmd_match in [GENE]:
				has_gene = True
				break
		else:
			print 'check HGMD_DESC[%s]' % var
	return has_gene

