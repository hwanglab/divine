#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan and modified by changjin.hong@gmail.com
#
"""
.. module:: preprocess_clinvar
		:platform: Unix, Windows, MacOSX
		:synopsis: Transparent opening of compressed and uncompressed files

.. moduleauthor:: ; changjin.hong@gmail.com
- 2018/01: ALLELEID is not available from clinvar.vcf.gz
"""

import os, dill, re, argparse
from gcn.lib.databases.snpdb import ClinvarDB
from collections import namedtuple
from gcn.lib.databases.ready2upload import lib_hgmd
from gcn.lib.databases.refgene import Refgene
from gcn.etc import fileconfig
from gcn.lib.utils.lib_utils import py_struct, month_to_num, msgout, open2, runcmd
from gcn.lib.io import anyopen
import gcn.lib.io.vcf as vcf

def store_variant_citations(variant_citation_fn):

	print 'storing variant citations ...'
	# linked_ids = py_struct(allele_id=[],
	# 											 variation_id=[],
	# 											 rs=[])
	linked_ids = {}
	tmp_fn = '%s.tmp'%variant_citation_fn
	cmd = "cut -f1,2,3 %s | sort -r -k1,1 -k2,2n -k3,3n | uniq > %s" % (variant_citation_fn,tmp_fn)
	runcmd(cmd)
	fp = anyopen.openfile(tmp_fn,'rt')
	head = fp.next()[:-1]
	if head.startswith('#'):
		head = head[1:]
	ntuple = namedtuple('ntuple', head.split('\t'))

	for i in fp:
		rec = i[:-1]
		linked_id = ntuple._make(rec.split('\t'))
		linked_ids[linked_id.AlleleID] = linked_id.VariationID

	fp.close()
	os.unlink(tmp_fn)
	print 'Done.'
	return linked_ids

def store_variant_summary(variant_summary_fn, linked_ids):
	'''
	objective: parse variant summary file from Clinvar
		and generate a dictionary (rcvaccess) of
		refSeq tx
		HGVSc, aa change in hgvs
		last evaluate date
		REV status
	:return:
	'''

	print 'storing variant summary ...'
	fp = anyopen.openfile(variant_summary_fn, 'rt')
	head_col = fp.next().split('\t')
	alleleid_i = head_col.index('#AlleleID')
	name_i = head_col.index('Name')
	rcvacc_i = head_col.index('RCVaccession')
	date_i = head_col.index('LastEvaluated')
	assembly_i = head_col.index('Assembly')
	review_i = head_col.index('ReviewStatus')
	# rsid_i = 9

	vars_to_summuary = {}
	for i in fp:
		itms = i.split('\t')
		if itms[assembly_i]=='GRCh37':

			if not itms[rcvacc_i].strip():
				continue

			allele_id = itms[alleleid_i]
			variation_id = None
			if allele_id in linked_ids:
				variation_id = linked_ids[allele_id]

			rcv_ids = itms[rcvacc_i].split(';')

			# if itms[rsid_i] == '-1':
			# 	rcv_ids = itms[rcvacc_i].split(';')
			# else:
			# 	rcv_ids = ['rs%s'%itms[rsid_i]]

			for rcv_id in rcv_ids:
				if rcv_id not in vars_to_summuary:
					vars_to_summuary[rcv_id] = py_struct(
						name=None,
						REFTX=None,
						HGVSc=None,
						HGVSp=None,
						DATE=None,
						REV=None,
						CLNMETHOD=None,
						allele_id=None,
						variation_id=None)

					vars_to_summuary[rcv_id].allele_id = allele_id
					if variation_id:
						vars_to_summuary[rcv_id].variation_id = variation_id

					mObj = re.search(r'(.+)\([\w]+\):c\.(.+)\s+\(p\.(.+)\)', itms[name_i])

					if mObj:
						vars_to_summuary[rcv_id].REFTX = mObj.group(1)
						vars_to_summuary[rcv_id].HGVSc = 'c.%s'%mObj.group(2)
						vars_to_summuary[rcv_id].HGVSp = 'p.%s'%mObj.group(3)
					else: #to handle a case where there is no aa hgvs
						mObj = re.search(r'(.+)\([\w]+\):c\.(.+)', itms[name_i])
						if mObj:
							vars_to_summuary[rcv_id].REFTX = mObj.group(1)
							vars_to_summuary[rcv_id].HGVSc = 'c.%s' % mObj.group(2)

					if itms[date_i] != '-':
						[mon_date,year] = itms[date_i].split(',')
						mon2, date2 = mon_date.split()
						mon2 = month_to_num(mon2)
						if mon2: month2 = mon2
						else: month2 = '01'
						vars_to_summuary[rcv_id].DATE = '%s-%s-%s'%(year.strip(),month2,date2.strip())

					vars_to_summuary[rcv_id].REV = itms[review_i].replace(' ','_')

	fp.close()
	print 'Done.'
	return vars_to_summuary

def store_submission_summary_fn2(submit_fn):

	print 'storing submission summary file ...'
	fp= anyopen.openfile(submit_fn, 'rt')
	submissions = {}
	read_heads = False
	for rec in fp:
		rec = rec[:-1]
		# print 'rec:',rec
		if not read_heads and rec.endswith('SubmittedGeneSymbol'):
			rec = rec[1:]
			heads = rec.split('\t')
			# print 'heads:',heads
			subm = namedtuple('subm', heads)
			read_heads = True
		elif read_heads:
			subm_rec = subm._make(rec.split('\t'))
			if subm_rec.VariationID not in submissions:
				submissions[subm_rec.VariationID] = py_struct(collection_methods=[])

			for cmethod in subm_rec.CollectionMethod.split(';'):
				submissions[subm_rec.VariationID].collection_methods.append(cmethod.replace(' ','_'))

	fp.close()

	print 'Done.'
	return submissions

def append_annotation_to_vcf(vcf_fn, vars_to_summuary, submissions, out_vcf):

	print 'appending annotation to clinvar VCF file ...'
	v = vcf.VCFParser(vcf_fn)
	ostream = open2(out_vcf, 'w')

	v.add_meta_info("REFTX", "1", "String", "RefSeq Transcript Name")
	v.add_meta_info("HGVSc", "1", "String", "HGVSc change in HGVS nomenclature")
	v.add_meta_info("HGVSp", "1", "String", "AA change in HGVS nomenclature")
	v.add_meta_info("SPLOC", "1", "Integer", "Distance from the predicted splice site")
	v.add_meta_info("DATE", "1", "String", "Last evaluated date")
	v.add_meta_info("REV", "1", "String", "Review status")
	v.add_meta_info("CLNMETHOD", "1", "String", "Collection methods")
	v.writeheader(ostream)

	for rec in v:
		v.parseinfo(rec)
		found = False
		for j, rcv_ids in enumerate([rec.id,rec.info.ALLELEID]):
			if not found:
				for rcv_id in rcv_ids:
					if j == 1: #strip off version
						rcv_id = rcv_id.split('.')[0]

					if rcv_id in vars_to_summuary:
						rec.info.REFTX = vars_to_summuary[rcv_id].REFTX
						if vars_to_summuary[rcv_id].HGVSc:
							rec.info.HGVSc = vars_to_summuary[rcv_id].HGVSc
							mObj = re.search(r'c\.(.*)([\+\-]\d+)\D+', rec.info.HGVSc)
							if mObj:
								SPLOC = mObj.group(2)
								if abs(int(SPLOC)) < 3:
									rec.info.SPLOC = SPLOC

						if vars_to_summuary[rcv_id].HGVSp:
							rec.info.HGVSp = vars_to_summuary[rcv_id].HGVSp
						if vars_to_summuary[rcv_id].DATE:
							rec.info.DATE = vars_to_summuary[rcv_id].DATE
						if vars_to_summuary[rcv_id].REV:
							rec.info.REV = vars_to_summuary[rcv_id].REV
						if vars_to_summuary[rcv_id].variation_id in submissions:
							cmethods = list(set(submissions[vars_to_summuary[rcv_id].variation_id].collection_methods))
							rec.info.CLNMETHOD = '|'.join(cmethods)

						found = True
						break

		for j,clndbn in enumerate(rec.info.CLNDN):
			rec.info.CLNDN[j] = clndbn.replace('\\x2c_', ',').replace('\\x2c', ',')

		v.write(ostream, rec)

	ostream.close()
	v.stream.close()
	print 'Done.'

def append_annotation_to_vcf2(vcf_fn, vars_to_summuary, submissions, out_vcf):

	print 'appending annotation to clinvar VCF file ...'
	v = vcf.VCFParser(vcf_fn)
	ostream = open2(out_vcf, 'w')

	v.add_meta_info("REFTX", "1", "String", "RefSeq Transcript Name")
	v.add_meta_info("HGVSc", "1", "String", "HGVSc change in HGVS nomenclature")
	v.add_meta_info("HGVSp", "1", "String", "AA change in HGVS nomenclature")
	v.add_meta_info("SPLOC", "1", "Integer", "Distance from the predicted splice site")
	v.add_meta_info("DATE", "1", "String", "Last evaluated date")
	v.add_meta_info("REV", "1", "String", "Review status")
	v.add_meta_info("CLNMETHOD", "1", "String", "Collection methods")
	v.writeheader(ostream)

	for rec in v:
		v.parseinfo(rec)

		# clnacc = re.split('[|,]', rec.info.ALLELEID)
		# rec.info.ALLELEID = '|'.join(list(set(clnacc)))

		uniq_rcv_ids = []
		for rcv_id_str in rec.info.ALLELEID:
			for rcv_id in rcv_id_str.split('|'):
				if rcv_id in uniq_rcv_ids: continue
				uniq_rcv_ids.append(rcv_id)

		# print 'rec.info.ALLELEID:',rec.info.ALLELEID #cj_debug
		for rcv_id in uniq_rcv_ids:

			rcv_id = rcv_id.split('.')[0]
			if rcv_id in vars_to_summuary:
				rec.info.REFTX = vars_to_summuary[rcv_id].REFTX
				if vars_to_summuary[rcv_id].HGVSc:
					rec.info.HGVSc = vars_to_summuary[rcv_id].HGVSc
					mObj = re.search(r'c\.(.*)([\+\-]\d+)\D+', rec.info.HGVSc)
					if mObj:
						SPLOC = mObj.group(2)
						if abs(int(SPLOC)) < 3:
							rec.info.SPLOC = SPLOC

				if vars_to_summuary[rcv_id].HGVSp:
					rec.info.HGVSp = vars_to_summuary[rcv_id].HGVSp
				if vars_to_summuary[rcv_id].DATE:
					rec.info.DATE = vars_to_summuary[rcv_id].DATE
				if vars_to_summuary[rcv_id].REV:
					rec.info.REV = vars_to_summuary[rcv_id].REV
				if vars_to_summuary[rcv_id].variation_id in submissions:
					cmethods = list(set(submissions[vars_to_summuary[rcv_id].variation_id].collection_methods))
					# print 'cmethods:',cmethods #cj_debug
					rec.info.CLNMETHOD = '|'.join(cmethods)

				found = True
				break

		rec.info.ALLELEID = uniq_rcv_ids
		for j,clndbn in enumerate(rec.info.CLNDN):
			rec.info.CLNDN[j] = clndbn.replace('\\x2c_', ',').replace('\\x2c', ',')

		v.write(ostream, rec)

	ostream.close()
	v.stream.close()
	print 'Done.'

def pathogenic_per_gene(cds_len_per_gene,hgmd_on=False):

	clnStat = ClinvarDB()

	vartypes = clnStat.count_lofs()

	# search for hgmd
	if hgmd_on:
		hgmdStat = lib_hgmd.HgmdRegion()
		vartypes = hgmdStat.count_lofs(vartypes = vartypes)

	# count vevent per gene
	known_patho_profs = {}
	for csigs in vartypes.itervalues():
		gene = csigs[-1]
		if gene not in known_patho_profs:
			cds_len = 480*1.25 #average CDS length in one transcript = 480
			if gene in cds_len_per_gene:
				cds_len = cds_len_per_gene[gene]

			# (benign, vus, pathogenic) x (lof,nsm)
			known_patho_profs[gene] = [[0, 0], [0, 0], [0, 0], cds_len]

		for i, csig in enumerate(csigs[:-1]):
			for j, vtype in enumerate(csig):
				known_patho_profs[gene][i][j] += 100.*vtype/cds_len

	return known_patho_profs

def known_pathov_stats(reuse=True, has_hgmd_license=False):
	"""
	to retrieve variant types (LOF, missense, etc) from known pathogenic mutation database (clinvar or HGMD)
	:return:
	"""
	pathog_prof_pyv = fileconfig.FILECONFIG['PATHOG_PROF']
	if reuse and os.path.exists(pathog_prof_pyv):
		msg = 'loading some statistics on known pathogenic variants (%s) ...' % pathog_prof_pyv
		msgout('notice', msg)
		fp = open(pathog_prof_pyv, 'rb')
		pathov_prof_gene = dill.load(fp)
		fp.close()
	else:
		refgene = Refgene()
		cds_len_per_gene = refgene.get_cds_len_per_gene()
		pathov_prof_gene = pathogenic_per_gene(cds_len_per_gene, hgmd_on=has_hgmd_license)
		fpw = open(pathog_prof_pyv, 'wb')
		dill.dump(pathov_prof_gene, fpw)
		fpw.close()

	#TODO: use SVM to infer optimal variables to classify benign vs. pathogenic

	return pathov_prof_gene

def main():
	
	desc = 'add RefSeq Transcript name, cDNA change, AA change in HGVS nomenclature, and review/date into Clinvar VCF file.'
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('-i', '--input', dest='vcf_fn', type=str, required=True, help='input clinvar VCF')
	parser.add_argument('-c', '--citation', dest='citation_fn', type=str, required=True,
											help='input clinvar citation')  # Clinvar db provides a separate tab delimited file
	parser.add_argument('-s', '--summary', dest='summary_fn', type=str, required=True, help='input clinvar summary') #Clinvar db provides a separate tab delimited file
	parser.add_argument('-S', '--submit', dest='submit_fn', type=str, required=True, help='input clinvar submit') #Clinvar db provides a separate tab delimited file
	parser.add_argument('-o', '--output', dest='out_fn', type=str, required=True, help='output clinvar VCF')
	args = parser.parse_args()

	linked_ids = store_variant_citations(args.citation_fn)
	dAnnotation = store_variant_summary(args.summary_fn, linked_ids)
	dSubmission = store_submission_summary_fn2(args.submit_fn)
	
	append_annotation_to_vcf2(args.vcf_fn, dAnnotation, dSubmission, args.out_fn)


if __name__ == "__main__":
	main()

# class Clinvar(snpdb.SNPDB):
#
# 	def __init__(self):
#
# 		#open db connection
# 		name = 'CLINVARDB'
# 		super(Clinvar, self).__init__(name)
# 		if name in dbconfig.DBCONFIG:
# 			self.load(name=name)
# 		elif os.path.exists(name):
# 			self.load(db=name)
# 		else:
# 			raise ValueError('No such database %s' % name)