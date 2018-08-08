"""
.. module:: annotator
	:platform: Unix, Windows, MacOSX
	:synopsis: SNP Annotation

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is the base module for annotating a given SNP
(Chr, Position, Ref allele, Obs allele). It calls other supporting
modules (annotateRegion, annotateSequence, annotateSite, annotateMutation)
to generate the annotaion.
"""

from annotateRegion import RegionAnnotation
from annotateSequence import SequenceAnnotation
from gcn.data.complement import COMPLEMENT
from annotateMutation import MutationAnnotation
from annotateLCR import LCR
from gcn.lib.databases.refgene import Refgene
from gcn.lib.databases.refgene import get_ucsc_chrom, to_ucsc_chrom
from bisect import bisect
from gcn.lib.databases.omim import Omim
from gcn.lib.databases.ready2upload.interpro_region import Interpro
from gcn.lib.databases.clnphe import Clnsnp
from annotateGerp import Gerp
from gcn.lib.databases.nsfp import Effpred
from gcn.lib.databases.utr import UTRdb
from gcn.lib.databases.regulomedb import Regulomedb
from gcn.lib.databases.mirna import MiRna
from gcn.lib.databases.splicedb import Splicedb
from gcn.lib.databases.snpdb import ClinvarDB
from gcn.lib.databases.snpdb import ClinvitaeDB
from gcn.lib.databases.snpdb import CosmicDB
from gcn.lib.databases.snpdb import HgmdDB
from gcn.lib.databases.snpdb import DBSNP
from gcn.lib.databases.snpdb import KGDB
from gcn.lib.databases.snpdb import EspDB
from gcn.lib.databases.snpdb import ExACDB
from gcn.lib.databases.annotator_to_db import DbLink
from annotateNimblegenCap import NimblegenCapture
from gcn.data.pseudoautosomal_genes import PSEUDO_AUTO_GENES
from gcn.etc.dbconfig import DBCONFIG
from gcn.lib.io.vcf import VCFParser
from gcn.lib.io.vcfutils import normalize_variant
from gcn.bin.hgConvert.hgvs_resource import Hgvs2
from gcn.lib.utils.lib_utils import joined
from collections import namedtuple
import os, sys, math
import time

CHROMOSOMES, ICHROMOSOMES = get_ucsc_chrom()
exacs = ['EXACAC', 'EXACAC_AFR', 'EXACAC_AMR', 'EXACAC_EAS', 'EXACAC_FIN', 'EXACAC_NFE', 'EXACAC_OTH', 'EXACAC_SAS', 'EXACAF', 'EXACAF_AFR', 'EXACAF_AMR', 'EXACAF_EAS', 'EXACAF_FIN', 'EXACAF_NFE', 'EXACAF_OTH', 'EXACAF_SAS']

def reformat_comma_string(invalue, adjlen, value='.'):
	if invalue is None:
		extraList = ['.'] * adjlen
		cmstring2 = joined(extraList, ',')
	else:
		if isinstance(invalue,str):
			myList = invalue.split(',')
		elif isinstance(invalue, (int, float, long, complex)):
			invalue = str(invalue)
			myList = [str(invalue)]
		if value is None:
			value = myList[0]
		L = len(myList)
		M = adjlen - L
		if M > 0:
			extraList = [value] * M
			cmstring2 = '%s,%s' % (invalue, joined(extraList, ','))
		elif M < 0:
			cmstring2 = joined(myList[:L],',')
		else:
			cmstring2 = invalue
	return cmstring2

class SNPAnnotation():

	def __init__(self, options):

		(capture_kit_name, probe_flanking_bp, cosmic_on, hgmd_on, max_exac_af, dblink) = options

		self.refgene = Refgene()
		self.mimdb = Omim()
		self.omim = self.mimdb.load_omim()
		self.pdomain = Interpro()
		self.clnsnpdb = Clnsnp()
		self.clnsnp = self.clnsnpdb.load_clnsnp()
		self.gerp = Gerp()
		self.genelist = self.load_refgene()
		self.effpred = Effpred()
		self.utrdb = UTRdb()
		self.utrdb.load_utrdb()
		self.mirnadb = MiRna()
		self.splicedb = Splicedb()
		self.regulomedb = Regulomedb()
		self.clinvardb = ClinvarDB()
		self.hgmd_genes = []
		self.hgmddb = None
		self.max_exac_af = max_exac_af
		self.dblink = dblink

		if hgmd_on:
			if os.path.exists(DBCONFIG['HGMDDB']['name']):
				self.hgmddb = HgmdDB()  # added by CJ
				self.hgmd_on = True
			else:
				print '[Warning] user attempts to use HGMD but it seems not configured properly (a license requires!)'
				self.hgmd_on = False
		else:
			self.hgmd_on = False

		self.clinvtdb = ClinvitaeDB()

		self.cosmicdb = None
		if cosmic_on:
			self.cosmicdb = CosmicDB()
		self.cosmic_on = cosmic_on

		self.NO_HGMD, self.HGMD_GENE, self.HGMD_POS, self.HGMD_ALT = range(4)
		self.kgdb = KGDB()
		self.dbsnp = DBSNP()
		self.espdb = EspDB()
		self.exacdb = ExACDB()  # added by CJ

		if capture_kit_name:
			self.nimbcap_on = True
			self.nimbcap = NimblegenCapture(capture_kit_name, probe_flanking_bp)
		else:
			self.nimbcap_on = False
			self.nimbcap = None

		self.lcr = LCR()

		# setting for HGVS nomenclature
		self.hgvs = Hgvs2()

	def _getpos(self, pos, ref, alt):
		"""Internal Method. Since in VCF file the coordinate for the indels
		are mentioned with reference to the previous neucleotide, so before
		annoation this correction needs to be made. e.g. -
		For Insertion case -
		In VCF : 1	   15903   .	   G	   GCAC
		is reframed as :  pos1=15903, pos2=15903, ref='-', alt='CAC'
		For Deletion case -
		In VCF : 1	   17961   .	   TGAGA	  T
		is reframed as : pos1=17962, pos2=17965, ref='GAGA',alt='-'
		Args:
			- pos:	Genomic Position
			- ref:	Reference Allele as mentioned in VCF
			- alt:	Alternate Allele as mentioned in VCF
		Returns:
			- spos:	Start genomic position for the variant
			- epos:	End genomic position for the variant
			- c_ref:	Corrected Reference Allele
			- c_alt:	Corrected Alternate Allele
		"""
		c_ref = None
		c_alt = None
		spos = None
		epos = None
		pos = int(pos)
		# Deletion
		if alt == '<DEL>':
			c_ref = ref
			c_alt = '-'
			spos = pos
			epos = spos + len(c_ref) - 1
		# Deletion
		elif len(ref) > len(alt):
			c_ref = ref[len(alt):]
			c_alt = '-'
			spos = pos + len(alt)
			epos = spos + len(c_ref) - 1
		# Insertion
		elif len(ref) < len(alt):
			c_ref = '-'
			c_alt = alt[len(ref):]
			spos = pos + len(ref) - 1
			epos = spos + 1
		# SNP
		else:
			c_ref = ref
			c_alt = alt
			spos = pos
			epos = spos + len(c_alt) - 1
		return spos, epos, c_ref, c_alt

	def load_refgene(self):
		genelist = [[] for c in CHROMOSOMES]
		for record in self.refgene.iterate_db():
			if record.chrom in ICHROMOSOMES:
				cid = ICHROMOSOMES[record.chrom]
			else:
				cid = len(CHROMOSOMES)
				ICHROMOSOMES[record.chrom] = cid
				CHROMOSOMES.append(record.chrom)
				genelist.append([])

			try:
				genelist[cid].append((int(record.str_transcript),
										int(record.end_transcript), record))
			except:
				print record.chrom, cid, CHROMOSOMES
				raise
		for key in genelist:
			key.sort()
		return genelist

	def retrieve_hgmdgenes(self):
		self.hgmd_genes = self.hgmddb.get_hgmd_genes()

	def retrieve_refgene(self, chrom, spos, epos):
		refgene_list = set([])
		chrom = to_ucsc_chrom(chrom)
		pos = int(spos)
		if chrom in ICHROMOSOMES:
			chrom = ICHROMOSOMES[chrom]
			sposidx = bisect(self.genelist[chrom], (spos, spos + 1))
			if sposidx != 0:
				for j in range(sposidx, -1, -1):
					i1, i2, record = self.genelist[chrom][j - 1]
					if pos >= i1 and pos <= i2:
						refgene_list.add(record)
			if spos != epos:
				eposidx = bisect(self.genelist[chrom], (epos, epos + 1))
				if eposidx != 0:
					for j in range(eposidx, -1, -1):
						i1, i2, record = self.genelist[chrom][j - 1]
						if pos >= i1 and pos <= i2:
							refgene_list.add(record)
		refgene_list = list(refgene_list)
		return refgene_list

	def extract_updwn_genes(self, chrom, spos, epos, updwnrange=5000):
		upgene = 'NONE(dist=NONE)'
		dwngene = 'NONE(dist=NONE)'
		if not chrom.startswith('chr'):
			 chrom = 'chr' + chrom

		if chrom in ICHROMOSOMES:
			chrom = ICHROMOSOMES[chrom]
			idx = bisect(self.genelist[chrom], (spos, spos + 1))
			if idx != 0:
				ui1, ui2, urecord = self.genelist[chrom][idx - 1]
				if idx == len(self.genelist[chrom]):
					udist = spos - ui2
					if udist <= updwnrange:
						upgene = urecord.gene + '(dist=' + str(spos - ui2)\
						+ ')'
				else:
					di1, di2, drecord = self.genelist[chrom][idx]
					if ui2 < spos and epos < di1:
						udist = spos - ui2
						if udist <= updwnrange:
							upgene = urecord.gene + '(dist=' + str(udist) + ')'
						ddist = di1 - epos + 1
						if ddist <= updwnrange:
							dwngene = drecord.gene + '(dist=' + \
												str(di1 - epos + 1) + ')'

			elif spos < self.genelist[chrom][0][0]:
				di1, di2, drecord = self.genelist[chrom][0]
				ddist = di1 - spos + 1
				if ddist <= updwnrange:
					dwngene = drecord.gene + '(dist=' + \
												str(di1 - spos + 1) + ')'

		updwngenes = ':'.join([upgene, dwngene])

		return updwngenes

	def get_revcmplmnt(self, allele):
		'''Returns reverse complemented allele'''
		if allele != '-':
			comp_allele = ''
			for neu in allele:
				comp_allele += COMPLEMENT[neu]
			comp_allele = comp_allele.upper()
			comp_allele = comp_allele[::-1]
		else:
			comp_allele = '-'
		return comp_allele

	def annotate(self, chrom, pos, ref, alt):
		'''Returns a list of annotation for the given Variant position and
		dbSNP Id'''

		annotations = []
		intron_entry = {}

		# back up the original input for HGVS representation input
		pos0 = pos
		ref0 = ref
		alt0 = alt

		self.chrom = chrom
		self.spos, self.epos, self.ref, self.alt = self._getpos(pos, ref, alt)

		if self.dblink:
			self.cDbLink.add_key(chrom, pos, ref, alt)

		# ExAc data
		exac_ant = ['EXACDB', False]
		exac_acs = [['EXACAC', []], \
					 ['EXACAC_AFR', []], \
					 ['EXACAC_AMR', []], \
					 ['EXACAC_EAS', []], \
					 ['EXACAC_FIN', []], \
					 ['EXACAC_NFE', []], \
					 ['EXACAC_OTH', []], \
					 ['EXACAC_SAS', []]]

		exac_afs = [['EXACAF', []], \
					 ['EXACAF_AFR', []], \
					 ['EXACAF_AMR', []], \
					 ['EXACAF_EAS', []], \
					 ['EXACAF_FIN', []], \
					 ['EXACAF_NFE', []], \
					 ['EXACAF_OTH', []], \
					 ['EXACAF_SAS', []]]

		E = len(exac_afs)

		flag = False
		if self.exacdb.has_snp(chrom, pos, ref):
			for rec in self.exacdb.snp_by_location(chrom, pos):
				alts = rec.alt.split(',')
				if rec.ref == ref and alt in alts:
					idx = alts.index(alt)
					exac_acs, exac_afs = self._get_exac_info(rec, idx, exac_acs, exac_afs)
					exac_ant[1] = True
					for j in range(E):
						exac_acs[j][1] = exac_acs[j][1][0]
						exac_afs[j][1] = exac_afs[j][1][0]
					if self.dblink:
						self.cDbLink.add_key(chrom,pos,ref,alt,db_name='EXAC',primary_idx=rec.idx)
					flag = True
					break

		if flag:
			#prefiltering
			if float(exac_afs[0][1]) > self.max_exac_af:
				return annotations
		else:
			varlist = normalize_variant(chrom, pos, ref, alt)  # update primary key for variant locations
			for var in varlist:
				if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
					for j in range(E):
						exac_acs[j][1].append('0')
						exac_afs[j][1].append('0.0')
					continue

				flag = False
				if self.exacdb.has_snp(chrom, var[1], var[2]):
					for rec in self.exacdb.snp_by_location(chrom, var[1]):
						alts = rec.alt.split(',')
						if rec.ref == var[2] and var[3] in alts:
							idx = alts.index(var[3])
							exac_acs, exac_afs = self._get_exac_info(rec, idx, exac_acs, exac_afs)
							exac_ant[1] = True
							if self.dblink:
								self.cDbLink.add_key(chrom, pos, var[2], var[3], db_name='EXAC', primary_idx=rec.idx)
							flag = True
							break

				if flag is False:
					for j in range(E):
						exac_afs[j][1].append('0.0')
						exac_acs[j][1].append('0')

			#prefiltering
			if len(exac_afs[0][1]) == 1 and float(exac_afs[0][1][0]) > self.max_exac_af:
				return annotations

			#converting the list in 2nd element into string
			for j in range(E):

				if len([e for e in exac_acs[j][1] if e == '0.0']) == len(exac_acs[j][1]):
					exac_acs[j][1] = '0'
				else:
					exac_acs[j][1] = ':'.join(exac_acs[j][1])

				if len([e for e in exac_afs[j][1] if e == '0.0']) == len(exac_afs[j][1]):
					exac_afs[j][1] = '0.'
				else:
					exac_afs[j][1] = ':'.join(exac_afs[j][1])

		exac = {'EXAC':[exac_ant, exac_acs, exac_afs]}

		# dbSNP
		dbsnp_version = self.dbsnp.get_version()[0].version[1:]
		dbsnp_ant = ['DB%s' % dbsnp_version, False]
		dbsnp_build = ['dbSNPBuildID', []]
		dbsnpids = []
		flag = False
		if self.dbsnp.has_snp(chrom, pos, ref):
			for rec in self.dbsnp.snp_by_location(chrom, pos):
				if rec.ref == ref and alt in rec.alt.split(','):
					dbsnp_ant[1] = True
					dbsnp_build[1] = str(rec.dbSNPBuildID)
					dbsnpids = rec.id.split(';')[-1]
					flag = True
					if self.dblink:
						self.cDbLink.add_key(chrom, pos, ref, alt,\
						                      db_name='SNPDB', primary_idx=rec.idx,\
						                      rsid=dbsnpids)
					break

		if flag is False:
			varlist = normalize_variant(chrom, pos, ref, alt)
			for var in varlist:
				if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
					dbsnp_build[1].append('.')
					dbsnpids.append('.')
					continue
				flag = False
				if self.dbsnp.has_snp(chrom, var[1], var[2]):
					for rec in self.dbsnp.snp_by_location(chrom, var[1]):
						if rec.ref == var[2] and var[3] in rec.alt.split(','):
							dbsnp_ant[1] = True
							dbsnp_build[1].append(str(rec.dbSNPBuildID))
							dbsnpids.append(rec.id.split(';')[-1])
							flag = True
							if self.dblink:
								self.cDbLink.add_key(chrom, pos, var[2], var[3],\
								                      db_name='SNPDB', primary_idx=rec.idx, rsid=dbsnpids)
							break
				if flag is False:
					dbsnp_build[1].append('.')
					dbsnpids.append('.')
			if len([e for e in dbsnp_build[1] if e == '.']) == \
											len(dbsnp_build[1]):
				dbsnp_build[1] = '.'
			else:
				dbsnp_build[1] = ':'.join(dbsnp_build[1])
			if len([e for e in dbsnpids if e == '.']) == len(dbsnpids):
				dbsnpids = '.'
			else:
				dbsnpids = ':'.join(dbsnpids)

		# 1000Genome
		kgdb_ant = ['KGDB', False]
		kgaf = ['KGAF', []]
		flag = False
		if self.kgdb.has_snp(chrom, pos, ref):
			for rec in self.kgdb.snp_by_location(chrom, pos):
				if rec.ref == ref and alt in rec.alt.split(','):
					kgdb_ant[1] = True
					kgaf[1] = str(rec.AF)
					flag = True
					if self.dblink:
						self.cDbLink.add_key(chrom, pos, ref, alt, db_name='KGDB', primary_idx=rec.idx)
					break
		if flag is False:
			varlist = normalize_variant(chrom, pos, ref, alt)
			for var in varlist:
				if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
					kgaf[1].append('0.0')
					continue
				flag = False
				if self.kgdb.has_snp(chrom, var[1], var[2]):
					for rec in self.kgdb.snp_by_location(chrom, var[1]):
						if rec.ref == var[2] and var[3] in rec.alt.split(','):
							kgdb_ant[1] = True
							kgaf[1].append(str(rec.AF))
							flag = True
							if self.dblink:
								self.cDbLink.add_key(chrom, pos, var[2], var[3], db_name='KGDB', primary_idx=rec.idx)
							break
				if flag is False:
					kgaf[1].append('0.0')
			if len([e for e in kgaf[1] if e == '0.0']) == len(kgaf[1]):
				kgaf[1] = '0.0'
			else:
				kgaf[1] = ':'.join(kgaf[1])

		# Exome Sequencing Project (ESP)
		esp_ant = ['ESPDB', False]
		espaf = ['ESPAF', []]
		flag = False
		if self.espdb.has_snp(chrom, pos, ref):
			for rec in self.espdb.snp_by_location(chrom, pos):
				if rec.ref == ref and alt in rec.alt.split(','):
					esp_ant[1] = True
					espaf[1] = str(float(rec.MAF.split(',')[-1]) / 100)
					flag = True
					if self.dblink:
						self.cDbLink.add_key(chrom, pos, ref, alt, db_name='ESP', primary_idx=rec.idx)
					break
		if flag is False:
			varlist = normalize_variant(chrom, pos, ref, alt)
			for var in varlist:  # 0:chrom, 1:pos, 2:ref, 3:alt, 4.variant_type
				if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
					espaf[1].append('0.0')
					continue
				flag = False
				if self.espdb.has_snp(chrom, var[1], var[2]):
					for rec in self.espdb.snp_by_location(chrom, var[1]):
						if rec.ref == var[2] and var[3] in rec.alt.split(','):
							esp_ant[1] = True
							espaf[1].append(str(float(
											rec.MAF.split(',')[-1]) / 100))
							flag = True
							if self.dblink:
								self.cDbLink.add_key(chrom, pos, var[2], var[3], db_name='ESP', primary_idx=rec.idx)
							break
				if flag is False:
					espaf[1].append('0.0')
			if len([e for e in espaf[1] if e == '0.0']) == len(espaf[1]):
				espaf[1] = '0.0'
			else:
				espaf[1] = ':'.join(espaf[1])

		# ClinvarDB
		cols = ['CLN' + ele for ele in ['ACC','DBN','METHOD','ORIGIN','SIG']]
		cols.extend(['REFTX','HGVSc','SPLOC', 'HGVSp','VLD', 'DATE'])

		clnacc,clndbn,clnmethod,clnorigin,clnsig \
		,clnreftx,clnhgvsc,clnsploc, clnhgvsp,clnvld,clndate \
			= [[ele, []] for ele in cols]

		clnant = [clnacc,clndbn,clnmethod,clnorigin,clnsig \
			,clnreftx,clnhgvsc,clnsploc, clnhgvsp,clnvld,clndate]

		flag = False
		if self.clinvardb.has_snp(chrom, pos, ref):
			for rec in self.clinvardb.snp_by_location(chrom, pos):
				numOfAlts = len(rec.alt.split(','))
				if rec.ref == ref and alt in rec.alt:
					if self.dblink:
						self.cDbLink.add_key(chrom, pos, ref, alt, db_name='CLINVARDB', primary_idx=rec.idx)

					cvars = [[rec.CLNACC,'.'],
									 [rec.CLNDBN,None],
									 [rec.CLNMETHOD,None],
									 [rec.CLNORIGIN,None],
									 [rec.CLNSIG,None],
									 [rec.REFTX,None],
									 [rec.HGVSc,'.'],
									 [rec.SPLOC,None],
									 [rec.HGVSp,'.'],
									 [rec.VLD,None],
									 [rec.DATE,None]]

					cvars2 = [rec.alt]
					for cvar,defVal in cvars:
						#print cvar
						cvars2.append(reformat_comma_string(cvar,numOfAlts,value=defVal))

					for cln in zip(*[e.split(',') for e in cvars2]):

						if alt == cln[0]:
							for f, v in zip(clnant, cln[1:]):
								f[1] = v
							flag = True
							break

		if flag is False:
			varlist = normalize_variant(chrom, pos, ref, alt)
			for var in varlist:
				if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
					for f in clnant:
						f[1].append('.')
					continue
				flag = False
				if self.clinvardb.has_snp(chrom, var[1], var[2]):
					for rec in self.clinvardb.snp_by_location(chrom, var[1]):
						numOfAlts = len(rec.alt.split(','))
						if rec.ref == var[2] and var[3] in rec.alt:
							if self.dblink:
								self.cDbLink.add_key(chrom, var[1], var[2], var[3], db_name='CLINVARDB', primary_idx=rec.idx)

							cvars = [rec.CLNACC,
											 rec.CLNDBN,
											 rec.CLNMETHOD,
											 rec.CLNORIGIN,
											 rec.CLNSIG,
											 rec.REFTX,
											 rec.HGVSc,
											 rec.SPLOC,
											 rec.HGVSp,
											 rec.VLD,
											 rec.DATE]

							cvars2 = [rec.alt]
							for cvar in cvars:
								#print cvar
								cvars2.append(reformat_comma_string(cvar, numOfAlts, value=None))

							for cln in zip(*[e.split(',') for e in cvars2]):
								if var[3] == cln[0]:
									for f, v in zip(clnant, cln[1:]):
										f[1].append(v)
									flag = True
									break
				if flag is False:
					for f in clnant:
						f[1].append('.')
			for f in clnant:
				if len([e for e in f[1] if e == '.']) == len(f[1]):
					f[1] = '.'
				else:
					f[1] = '__'.join(f[1])

		# Check for Nimblegen Capture Array
		nimb_capture = None
		if self.nimbcap_on:
			nimb_capture = self.nimbcap.get_info(self.chrom, self.spos, self.epos)

		# Gerp Conservation
		gerp_consv = self.gerp.get_conserve(self.chrom, self.spos, self.epos)

		# LCR Annotation
		lcr_ant = self.lcr.get_info(self.chrom, self.spos, self.epos)

		# Interpro Domain
		# TODO: this process may take so long and we'd better move to acmg.pm1()? and search for only highly suspect pathov
		# self.pdomain.get_domain_desc(self.chrom, self.spos, self.epos)
		pgenes,pdenoms,pbdns,pvdns,ppdns = self.pdomain.get_patho_density_by_pos(self.chrom, self.spos, self.epos)

		if self.dblink and self.pdomain.annot.indices:
			self.cDbLink.add_primary_idx('INTERPRO', self.pdomain.annot.indices)

		# domain annotation
		pdomain = None
		if pgenes:
			if pdenoms[0] > 0.:
				pdomain = ['PATHO_DOMAIN','%d,%g,%g,%g' % (pdenoms[0], pbdns[0], pvdns[0], ppdns[0])]

		# RegulomeDB
		regulome_score = []
		for rsid in dbsnpids.split(':'):
			rec = self.regulomedb.retrieve(rsid)
			if rec:
				key = rec[0]
				regulome_score.append(rec[2])
				if self.dblink:
					self.cDbLink.add_snpid(rsid,'REGULOME',rec[1])

		if len(regulome_score) > 0:
			regulome_score = [key, ':'.join(regulome_score)]
		else:
			regulome_score = None

		# NHGRI-GWAS
		gwas_phen = []
		for rsid in dbsnpids.split(':'):
			if rsid in self.clnsnp.db['NHGRI']:
				gwas_phen.append('__'.join(self.clnsnp.db['NHGRI'][rsid]))

				if self.dblink and rsid in self.clnsnp.rsid2idx:
					self.cDbLink.add_snpid(rsid,'CLNPHESNP',self.clnsnp.rsid2idx[rsid])

		if gwas_phen:
			gwas_phen = ':'.join(gwas_phen)

		# ------------------------------------
		# HGMD
		hgmd = ['HGMD', [], []]  # Header, string_to_printout, hgmd_acc

		NotAvail = '.'

		if self.hgmd_on:
			hgmd_alt_match_cnt = 0
			if self.hgmddb.has_snp(chrom, pos, ref):
				for rec in self.hgmddb.snp_by_location(chrom, pos):
					# to store hgmd hits per chrom_pos
					match_type = self.HGMD_POS
					hgmd_alt_match, hgmd_desc = self._hgmd_has_alt(ref, alt, rec)
					if hgmd_alt_match:
						match_type = self.HGMD_ALT
						hgmd_alt_match_cnt += 1
						if self.dblink:
							self.cDbLink.add_key(chrom, pos, ref, alt, db_name='HGMDDB', primary_idx=rec.idx)
					else:
						hgmd_desc = self._concat_hgmd_desc(rec)

					hgmd_desc = '%d(%s)' % (match_type, hgmd_desc)
					hgmd[1].append(hgmd_desc)
					hgmd[2].append(rec.HGMD_ACC)

			if hgmd_alt_match_cnt == 0:
				varlist = normalize_variant(chrom, pos, ref, alt)

				for var in varlist:  # 0:chrom, 1:pos, 2:ref, 3:alt, 4.variant_type
					if self.hgmddb.has_snp(chrom, var[1], var[2]):
						for rec in self.hgmddb.snp_by_location(chrom, var[1]):
							if not rec.HGMD_ACC in hgmd[2]:  # was it printed prevly?
								match_type = self.HGMD_POS
								hgmd_alt_match, hgmd_desc = self._hgmd_has_alt(var[2], var[3], rec)
								if hgmd_alt_match:
									match_type = self.HGMD_ALT
									hgmd_alt_match_cnt += 1
									if self.dblink:
										self.cDbLink.add_key(chrom, var[1], var[2], var[3], db_name='HGMDDB', primary_idx=rec.idx)
								else:
									hgmd_desc = self._concat_hgmd_desc(rec)

								hgmd_desc = '%d(%s)' % (match_type, hgmd_desc)
								hgmd[1].append(hgmd_desc)
								hgmd[2].append(rec.HGMD_ACC)

			hgmd_pos_match = False
			if len(hgmd[2]) > 0:
				hgmd_pos_match = True
				hgmd[1] = ','.join(hgmd[1])

		clinvt = ['CLINVT', []]
		flag = False
		if self.clinvtdb.has_snp(chrom, pos, ref):
			for rec in self.clinvtdb.snp_by_location(chrom, pos):
				clinvt_desc = self._clinvt_has_alt(ref, alt, rec)
				if clinvt_desc:
					clinvt[1].append(clinvt_desc)

					if self.dblink:
						self.cDbLink.add_key(chrom, pos, ref, alt, db_name='CLINVITAE', primary_idx=rec.idx)

					if 'pathogenic' in clinvt_desc:
						flag = True
						break

		if not flag:
			varlist = normalize_variant(chrom, pos, ref, alt)
			for var in varlist:
				if flag: break
				if self.clinvtdb.has_snp(chrom, var[1], var[2]):
					for rec in self.clinvtdb.snp_by_location(chrom, var[1]):
						clinvt_desc = self._clinvt_has_alt(ref, alt, rec)
						if clinvt_desc:
							clinvt[1].append(clinvt_desc)

							if self.dblink:
								self.cDbLink.add_key(chrom, var[1], var[2], var[3], db_name='CLINVITAE', primary_idx=rec.idx)

							if 'pathogenic' in clinvt_desc:
								flag = True
								break

		cosmic = ['COSMIC', []]
		if self.cosmic_on:
			flag = False
			if self.cosmicdb.has_snp(chrom, pos, ref):
				for rec in self.cosmicdb.snp_by_location(chrom, pos):
					alts = rec.alt.split(',')
					if rec.ref == ref and alt in alts:
						cosmic[1].append(rec.COSMIC_ID)
						if self.dblink:
							self.cDbLink.add_key(chrom, pos, ref, alt, db_name='COSMIC', primary_idx=rec.idx)
						flag = True
						break

			if flag is False:
				varlist = normalize_variant(chrom, pos, ref, alt)  # update primary key for variant locations

				for var in varlist:
					if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
						continue

					flag = False
					if self.cosmicdb.has_snp(chrom, var[1], var[2]):
						for rec in self.cosmicdb.snp_by_location(chrom, var[1]):
							alts = rec.alt.split(',')
							if rec.ref == var[2] and var[3] in alts:
								cosmic[1].append(rec.COSMIC_ID)
								if self.dblink:
									self.cDbLink.add_key(chrom, var[1], var[2], var[3], db_name='COSMIC', primary_idx=rec.idx)
								flag = True
								break
			cosmic[1] = list(set(cosmic[1]))

		# CADD raw and phred score
		cadd_raw, cadd_phred, cadd_aa, nsfp_idx = (None, None, None, None)
		if self.ref != '-' and self.alt != '-':
			rl = []
			pl = []
			ral = []
			for idx, vpos in enumerate(range(self.spos, self.epos + 1)):
				nsfp_idx, raw, phred, refaa, altaa = self.effpred.get_cadd(\
					self.chrom, vpos, self.ref[idx], self.alt[idx])
				if raw is not None:
					rl.append(str(raw))
					pl.append(str(phred))
					ral.append('%s/%s' % (refaa, altaa))

					if self.dblink:
						self.cDbLink.add_key(chrom, vpos, self.ref[idx], self.alt[idx], \
						                     db_name='NSFP', primary_idx=nsfp_idx)

			if rl:
				cadd_raw = ','.join(rl)
				cadd_phred = ','.join(pl)
				cadd_aa = ','.join(ral)

		#-------------------------------------------------------
		#gene level annotation or region-specific annotation
		refgene_entries = self.retrieve_refgene(self.chrom, self.spos, self.epos)

		# Pseudoautosomal region annotation variable
		par = ['PAR', False]

		if refgene_entries:
			hgmd_gene_printed = {}  # to avoid duplicated output
			for entry in refgene_entries:
				if self.dblink:
					self.cDbLink.add_genes(entry.gene)

				warning = ''

				# Allele changes based on strands
				if entry.strand == '+':
					alt = self.alt
					ref = self.ref
				elif entry.strand == '-':
					alt = self.get_revcmplmnt(self.alt)
					ref = self.get_revcmplmnt(self.ref)

				# Transcript and Gene
				refseq_acc = entry.refseq_acc
				gene = entry.gene

				error_state = int(entry.error)

				# Checks if gene is a pseudoautosomal gene
				if not par[1]:
					if gene in PSEUDO_AUTO_GENES:
						par[1] = True

				# Gene -Level Disease Association
				if gene in self.omim:
					omim_phen, omim_ids, omim_idx = self.omim[gene]
					if self.dblink:
						self.cDbLink.add_primary_idx('OMIM',omim_idx)
				else:
					omim_phen, omim_ids = [None] * 2

				if gene in self.clnsnp.db['GAD']:
					site_gad = self.clnsnp.db['GAD'][gene]
					site_gad = '__'.join(site_gad)
					if self.dblink:
						self.cDbLink.add_primary_idx('CLNPHESNP',self.clnsnp.gene2idx[gene])
				else:
					site_gad = None

				# Interpro Domain
				idx_domains, domains = self.effpred.get_domains(self.chrom, self.spos,
																												self.epos)

				if self.dblink and not nsfp_idx:
					self.cDbLink.add_primary_idx('NSFP',idx_domains)

				# Region
				region_annot = RegionAnnotation(self.spos, self.epos,
												self.ref, self.alt, entry)
				region_annot.annotate()
				region = region_annot.region
				splice_site = region_annot.splice_site
				exons = region_annot.exons
				if region == 'Intergenic':
					if entry.strand == '+':
						gene = self.extract_updwn_genes(self.chrom,
										self.spos - 1, self.spos - 1)
					else:
						gene = self.extract_updwn_genes(self.chrom,
											self.spos, self.spos)

				if self.hgmd_on and hgmd_pos_match is False:
					if gene in self.hgmd_genes:
						if gene not in hgmd_gene_printed:
							hgmd[1] = '%s(%s)' % (self.HGMD_GENE, gene)
							hgmd_gene_printed[gene] = True

				# if region_annot.cdna_poslist:
				if 'Intergenic' not in region:
					# cdna_pos = '_'.join([str(e) for e in region_annot.cdna_poslist])
					if alt0 == '<DEL>':
						cdna_pos = None
					else:
						cdna_pos = self.hgvs.to_cDNA(self.chrom, pos0, ref0, alt0, refseq_acc)
				else:
					cdna_pos = None

				protein_len = None
				# For Intronic Variants retaining only one transcript
				if region[-8:] == 'Intronic':
					if '__' in splice_site:
						sddist = int(splice_site.split('__')[0].split('_')[-1])
						sadist = int(splice_site.split('__')[1].split('_')[-1])
					else:
						sddist, sadist = [0, 0]

					if refseq_acc.startswith('NM'):
						rrank = 1
					else:
						rrank = 2

					annt = [gene, refseq_acc, region, exons, cdna_pos, splice_site] + \
							[None] * 8 + [None] + [omim_phen, omim_ids, site_gad,
														gerp_consv, lcr_ant, regulome_score, gwas_phen,
														[dbsnpids, dbsnp_ant, dbsnp_build,
														 kgdb_ant, kgaf, esp_ant, espaf, exac,
														 hgmd[:2], clinvt, cosmic,
														 clnacc, clndbn, clnmethod, clnorigin, clnsig \
															, clnreftx, clnhgvsc, clnsploc, clnhgvsp, clnvld, clndate],
													nimb_capture, pdomain, par, cadd_raw, cadd_phred, cadd_aa, domains]
					if gene not in intron_entry:
						intron_entry[gene] = [annt, rrank, sddist, sadist]
					else:
						if rrank < intron_entry[gene][1]:
							intron_entry[gene] = [annt, rrank, sddist, sadist]
						elif sddist < intron_entry[gene][2] or \
								sadist < intron_entry[gene][3]:
							intron_entry[gene] = [annt, rrank, sddist, sadist]

				# For the refgene entries whose CDS is not multiple of 3
				elif error_state == 1:
					annt = [gene, refseq_acc, region, exons, cdna_pos,
									splice_site] + [None] * 8 + ['CDS_NOT_MULTIPLE_OF_3'] + \
								 [omim_phen, omim_ids, site_gad, gerp_consv,
									lcr_ant, regulome_score, gwas_phen, [dbsnpids, dbsnp_ant,
																											 dbsnp_build, kgdb_ant, kgaf, esp_ant, espaf, exac,
																											 hgmd[:2], clinvt, cosmic,clnacc, clndbn, clnmethod,
																																clnorigin, clnsig, clnreftx, clnhgvsc,
																																clnsploc, clnhgvsp, clnvld, clndate],
									nimb_capture, pdomain, par, cadd_raw, cadd_phred, cadd_aa,domains]
					annotations.append(annt)

				else:
					# For Exonic and UTR variants
					# Sequence and Mutation
					seq_annot = SequenceAnnotation(region, entry,
									self.spos, self.epos, ref, alt, warning)
					status, warning = seq_annot.annotate()
					if status is True:
						if seq_annot.cds:
							protein_len = str((len(seq_annot.cds) / 3) - 1)
						mut_annot = MutationAnnotation(region, seq_annot)
						mut_annot.annotate()
						mut_type = mut_annot.mut_type
						codon_usage = mut_annot.codon_usage
						if type(seq_annot.wt_codon) == list:
							aa_change, codon_change = ([], [])
							for waa, maa, wcdn, mcdn, aapos in zip(
										 seq_annot.wt_aa, seq_annot.mut_aa,
										 seq_annot.wt_codon,
										 seq_annot.mut_codon,
										 seq_annot.wt_aa_position):
								codon_change.append(wcdn + '/' + mcdn)
								aa_change.append(waa + str(aapos) + maa)
								if mut_type in ['SynStop', 'StopLoss']:
									if int(aapos) <= int(protein_len):
										warning = 'NOT_ACTUAL_STOP_CODON__' \
										+ 'TRANSCRIPT_WITH_MULTIPLE_STOP_CODON'
							aa_change = '__'.join(aa_change)
							codon_change = '__'.join(codon_change)
						else:
							if seq_annot.wt_codon:
								codon_change = seq_annot.wt_codon + \
											'/' + seq_annot.mut_codon
								aa_change = seq_annot.wt_aa + \
											str(seq_annot.wt_aa_position) \
											+ seq_annot.mut_aa
							else:
								codon_change, aa_change = (None, None)
					else:
						mut_type, codon_usage, codon_change, aa_change, \
										wt_aa, mut_aa, aa_pos = [None] * 7

					# UTRdb
					utrannt = set([])
					if 'UTR5' in region or 'UTR3' in region:
						for pos in range(self.spos, self.epos + 1):
							utra = self.utrdb.retrieve(self.chrom,
											self.spos, refseq_acc)
							for e in list(utra):
								utrannt.add(e)

					# miRNA Binding Site Annotation
					mirnas = []
					mirnas_idx = None
					if region == 'UTR3':
						mirnas,mirnas_idx = self.mirnadb.get_annot(
							self.chrom,self.spos, refseq_acc)

					# UTR Functional Site : Combine annotation from utrdb
					# and mirna database
					utr_funsites = None
					func_sites = list(utrannt) + mirnas
					if func_sites:
						utr_funsites = '__'.join(func_sites)
						if self.dblink and mirnas_idx:
							self.cDbLink.add_primary_idx('MIRNA', mirnas_idx)

					# Splice Enhancer/Silencer Site Annotation
					if region in ['CodingExonic', 'NonCodingExonic']:
						splice_reg_site,splice_reg_site_idx = \
							self.splicedb.get_annot(self.chrom, self.spos, refseq_acc)

						if self.dblink:
							self.cDbLink.add_primary_idx('SPLICE', splice_reg_site_idx)

						if splice_reg_site:
							splice_site = splice_reg_site

					# Effect Prediction
					sift_pred, pphen_pred = (None, None)
					if type(seq_annot.wt_aa) == list:
						sp = []
						pp = []
						for wt_aa, mut_aa, aa_pos in zip(seq_annot.wt_aa,
												seq_annot.mut_aa,
												seq_annot.wt_aa_position):
							if wt_aa != mut_aa:
								if self.spos == self.epos:
									self.effpred(self.chrom, self.spos, wt_aa,
												 mut_aa.replace('*', 'X'),
												 aa_pos, self.ref, self.alt)
								else:
									self.effpred(None, None, wt_aa,
												 mut_aa.replace('*', 'X'),
												 aa_pos, None, None, gene)

								# Sift
								if self.effpred.sift_pred:
									sp.append(self.effpred.sift_pred + '_' +
											str(self.effpred.sift_score))

								# Polypen2
								if self.effpred.pp2_pred:
									pp.append(self.effpred.pp2_pred + '_' +
												str(self.effpred.pp2_score))
						sift_pred = '__'.join(sp)
						pphen_pred = '__'.join(pp)

					annt = [gene, refseq_acc, region, exons, cdna_pos,
							splice_site,
							utr_funsites, mut_type, codon_change,
							aa_change, protein_len, codon_usage, sift_pred,
							pphen_pred, warning, omim_phen,
							omim_ids, site_gad, gerp_consv, lcr_ant,
							regulome_score, gwas_phen,
									[dbsnpids, dbsnp_ant, dbsnp_build, kgdb_ant,
									 kgaf, esp_ant, espaf, exac,
									 hgmd[:2], clinvt, cosmic,
									 clnacc, clndbn, clnmethod, clnorigin, clnsig,
									 clnreftx, clnhgvsc, clnsploc, clnhgvsp, clnvld, clndate],
									nimb_capture, pdomain, par, cadd_raw,
							cadd_phred, cadd_aa, domains]

					annotations.append(annt)

			if intron_entry:
				for ele in intron_entry.values():
					annotations.append(ele[0])
		else:
			updwngenes = self.extract_updwn_genes(chrom, self.spos, self.epos)
			annt = [updwngenes, None, 'Intergenic'] + 15 * [None]\
								+ [gerp_consv, lcr_ant, regulome_score, gwas_phen,
									 [dbsnpids, dbsnp_ant,
										dbsnp_build, kgdb_ant, kgaf, esp_ant,
										espaf, exac, hgmd[:2], clinvt, cosmic,
										clnacc, clndbn, clnmethod, clnorigin, clnsig,
										clnreftx, clnhgvsc, clnsploc, clnhgvsp, clnvld, clndate],
									 nimb_capture, pdomain, par, cadd_raw, cadd_phred, cadd_aa, None]
			annotations.append(annt)

		return annotations

	def _get_exac_info(self, rec, idx, exac_acs, exac_afs):

		ac = rec.AC_Adj.split(',')[idx]
		an = float(rec.AN_Adj)
		exac_acs[0][1].append(ac)
		if an > 0.:
			exac_afs[0][1].append(str(float(ac) / an))
		else:
			exac_afs[0][1] = rec.AF.split(',')[idx]

		# to get AF for each ethnic group
		ac = rec.AC_AFR.split(',')[idx]
		an = float(rec.AN_AFR)
		exac_acs[1][1].append(ac)
		if an > 0.:
			exac_afs[1][1].append(str(float(ac) / an))
		else:
			exac_afs[1][1].append('0.')

		ac = rec.AC_AMR.split(',')[idx]
		exac_acs[2][1].append(ac)
		an = float(rec.AN_AMR)
		if an > 0.:
			exac_afs[2][1].append(str(float(ac) / an))
		else:
			exac_afs[2][1].append('0.')

		ac = rec.AC_EAS.split(',')[idx]
		exac_acs[3][1].append(ac)
		an = float(rec.AN_EAS)
		if an > 0.:
			exac_afs[3][1].append(str(float(ac) / an))
		else:
			exac_afs[3][1].append('0.')

		ac = rec.AC_FIN.split(',')[idx]
		exac_acs[4][1].append(ac)
		an = float(rec.AN_FIN)
		if an > 0.:
			exac_afs[4][1].append(str(float(ac) / an))
		else:
			exac_afs[4][1].append('0.')

		ac = rec.AC_NFE.split(',')[idx]
		exac_acs[5][1].append(ac)
		an = float(rec.AN_NFE)
		if an > 0.:
			exac_afs[5][1].append(str(float(ac) / an))
		else:
			exac_afs[5][1].append('0.')

		ac = rec.AC_OTH.split(',')[idx]
		exac_acs[6][1].append(ac)
		an = float(rec.AN_OTH)
		if an > 0.:
			exac_afs[6][1].append(str(float(ac) / an))
		else:
			exac_afs[6][1].append('0.')

		ac = rec.AC_SAS.split(',')[idx]
		exac_acs[7][1].append(ac)
		an = float(rec.AN_SAS)
		if an > 0.:
			exac_afs[7][1].append(str(float(ac) / an))
		else:
			exac_afs[7][1].append('0.')

		return (exac_acs, exac_afs)

	def _concat_hgmd_desc(self, rec):
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

	def _hgmd_has_alt(self, ref, alt, rec):
		"""
		input1: ref (nt base in reference sequence)
		input2: alt (nt base in query sample)
		input3: rec (a row in the VCF record)
		input4: hgmd (hgmd_desc ['HGMD_DESC', hgmd annotation, matching type]

		output1: flag (whether both ref and alt are matched exactly?)
		output2: hgmd_desc ()
		"""

		flag = False
		hgmd_desc = None

		alts = rec.alt.split(',')
		if rec.ref == ref and alt in alts:
			hgmd_desc = self._concat_hgmd_desc(rec)
			flag = True

		return flag, hgmd_desc

	def _concat_clinvt_desc(self, rec):
		vc_in_db = rec.VC.lower()
		cvt_vc = 'vus'
		if 'conflicting' in vc_in_db:
			cvt_vc = 'vus'
		elif 'benign' in vc_in_db:
			cvt_vc = 'benign'
		elif 'pathogenic' in vc_in_db:
			cvt_vc = 'pathogenic'
		return "%s(%s|%s)" % (cvt_vc, rec.cDNA, rec.URL)

	def _clinvt_has_alt(self, ref, alt, rec):

		clinvt_desc = None
		alts = rec.alt.split(',')
		if rec.ref == ref and alt in alts:
			clinvt_desc = self._concat_clinvt_desc(rec)
		return clinvt_desc


	def get_header(self):
		""" Intializing vcf formatted header information for populating VARANT
		annotation data"""

		Varant_h = namedtuple('Varant_h', 'id count type desc')

		header = [Varant_h('VARANT_INTERGENIC', 1, 'String', "Variant in "
							 "intergenic region and so reporting upstream and "
							 "downstream genes which are 5000bp away from "
							 "variant position and their exact (TSS/TES) "
							 "distance from variant position. Format - "
							 "VARANT_INTERGENIC"
							 "=UpstreamGene(dist=XYZ):DownstreamGene(dist=XYZ)"),
						Varant_h('VARANT_GENIC', '.', 'String', "Variant in genic "
						"region and so reporting the effects on genes. "
						"Format - VARANT_GENIC="
						"Gene1(Transcript1_id|Region|Exon_Number|AltId|cDNAPos|"
						"SpliceSite|UTRSignal|Mutation|Codon_Change|"
						"Amino_Acid_Change|Ref_Protein_Length|Codon_Usage|"
						"SIFT(pred_score)|Polyphen2(pred_score)|Warning: "
						"OMIM_Disease:OMIM_Ids:GAD_Disease). If there are "
						"more than one transcripts for the gene, the effect "
						"on them are appended after the annotation of first "
						"transcript by \':\'"),
						Varant_h('CADD_raw', '.', 'String', 'CADD raw score for '
							 'funtional prediction of a SNP. Please '
							 'refer to Kircher et al.(2014) Nature Genetics '
							 '46(3):310-5 for details. The larger the score '
							 'the more likely the SNP has damaging effect. '
							 'Scores are reported in the order in which ALT '
							 'alleles are reported.'),
						Varant_h('CADD_phred', '.', 'String', 'CADD phred-like score.'
							 ' This is phred-like rank score based on whole '
							 'genome CADD raw scores. Please refer to Kircher '
							 'et al. (2014) Nature Genetics 46(3):310-5 for '
							 'details. The larger the score the more likely '
							 'the SNP has damaging effect. Scores are reported'
							 ' in the order in which ALT alleles are reported'),
						Varant_h('CADD_aa', '.', 'String', 'Amino acid change in CADD prediction.'
							 ' This is change of amino acid based on whole '
							 'genome CADD raw scores. Please refer to Kircher '
							 'et al. (2014) Nature Genetics 46(3):310-5 for '
							 'details. The larger the score the more likely '
							 'the SNP has damaging effect. Scores are reported'
							 ' in the order in which ALT alleles are reported'),
						Varant_h('Interpro_domains', '.', 'String', 'domain or '
							 'conserved site on which the variant locates. '
							 'Domain annotations come from Interpro database.'),
						Varant_h('Pathogenicity_in_Interpro', '.', 'String', 'Pathogenicity in the Interpro domain'
							'pathogenicity density divided by sqrt of benign plus pathogenicity.'),
						]
		if self.dblink:
			header.append(Varant_h('VARIANT_IDX', 1, 'Integer', 'Variant index in raw VCF file (Temporary).'))

		return header


def add_meta(vcf, snpa):
	'''Add Annotation related meta data to VCF'''
	dbobjects = [snpa.clinvardb, snpa.clnsnpdb, snpa.dbsnp, snpa.effpred,
				 snpa.espdb, snpa.exacdb, snpa.gerp, snpa.lcr, snpa.kgdb, snpa.mirnadb,
				 snpa.mimdb, snpa.refgene,
				 snpa.regulomedb, snpa.utrdb, snpa.splicedb]

	if snpa.nimbcap_on:
		dbobjects.append(snpa.nimbcap)

	if snpa.hgmd_on:
		dbobjects.append(snpa.hgmddb)

	dbversion = None
	for obj in dbobjects:
		for ele in obj.get_version():
			if dbversion:
				dbversion += ',' + ele.dbname + '(' + ele.version + ')'
			else:
				dbversion = ele.dbname + '(' + ele.version + ')'
	vcf.add_meta('VARANTDB_VERSIONS', '"' + dbversion + '"')

	if snpa.nimbcap_on:
		for header in NimblegenCapture(snpa.nimbcap.kit_symbol, snpa.nimbcap.ext_bp).get_header():
			vcf.add_meta_info(header.id, header.count, header.type, header.desc)

	for header in Gerp().get_header():
		vcf.add_meta_info(header.id, header.count, header.type, header.desc)
	for header in snpa.lcr.get_header():
		vcf.add_meta_info(header.id, header.count, header.type, header.desc)
	for header in Clnsnp().get_header():
		vcf.add_meta_info(header.id, header.count, header.type, header.desc)
	for header in Regulomedb().get_header():
		vcf.add_meta_info(header.id, header.count, header.type, header.desc)
	for header in snpa.get_header():
		vcf.add_meta_info(header.id, header.count, header.type, header.desc)
	dbsnp_version = snpa.dbsnp.get_version()[0].version[1:]
	vcf.add_meta_info('PAR', 0, 'Flag', "Variant is in Pseudoautosomal" +
						" Region")
	vcf.add_meta_info('DB%s' % dbsnp_version, 0, 'Flag', "Variant" +
						" is present in dbSNP%s" % dbsnp_version)

	vcf.add_meta_info('dbSNPBuildID', 'A', 'String', "First dbSNP " +
						"Build for RS ordered by ALT allele")

	vcf.add_meta_info('KGDB', 0, 'Flag', "Variant" +
						" is present in 1000Genome Database")

	vcf.add_meta_info('KGAF', 'A', 'String', "Minor Allele Frequency" +
						" reported in 1000Genome Database" +
						" ordered by ALT allele")

	vcf.add_meta_info('ESPDB', 0, 'Flag', "Variant" +
							" is present in Exome Sequencing Project" +
							"(ESP) Database")
	vcf.add_meta_info('ESPAF', 'A', 'String', "The minor-allele frequency "
						"computed based on 6503 exomes samples in Exome "
						"Sequencing Project Database ordered by ALT allele")
	vcf.add_meta_info('EXACDB', 0, 'Flag', "Variant" +
							" is present in Exome Aggregation Consortium" +
							"(EXAC) Database")
	for exac in exacs:
		vcf.add_meta_info(exac, 'A', 'String', "The minor-allele frequency of %s " % exac +
							"computed based on 60,706 unrelated individuals in Exome " +
							"Aggregation Consortium Database ordered by ALT allele")

	if snpa.hgmd_on:
		vcf.add_meta_info('HGMD', '.', 'String', "[0:NO_HGMD,1:GENE_MATCH,2:POS_MATCH,3:POS_ALT_MATCH](ACC_ID_or_GeneSymbol|Variant_Class|PubMed|Phenotype)")

	if snpa.cosmic_on:
		vcf.add_meta_info('COSMIC', '.', 'String', "cosmic IDs")

	vcf.add_meta_info('CLINVT', '.', 'String', "class(cDNA|source_url)")

	clnvir = ClinvarDB()
	for field in ['CLNACC','CLNDBN','CLNMETHOD','CLNORIGIN','CLNSIG','REFTX','HGVSc','SPLOC','HGVSp','VLD','DATE']:

		meta_field = clnvir.get_meta_info(field)
		vcf.add_meta_info(meta_field.id, 'A',
							meta_field.type, meta_field.description)

	vcf.add_meta_info('PATHO_DOMAIN', '1', 'String', "total_number_reported_variants_in_the_domain,benign_density,vus_density,pathogenic_density")

	if snpa.dblink:
		vcf.add_meta_info('VARIANT_IDX', 1, 'Integer', "Variant index in raw VCF file")

def add_annotation(vcf, rec, annotations):
	'''Add annotation to the info field of the vcf record'''
	if len(annotations) > 0:
		for i in range(0, len(annotations), 2):
			vcf.add_info(rec, annotations[i], annotations[i + 1])

def get_min_maf(mafStr):
	"""Choose minimum MAF from MAF scores concatenated by :"""
	maf_min = 1.0
	for maf in mafStr.split(':'):
		if float(maf)<maf_min:
			maf_min =float(maf)
	return maf_min

def main(vcffile, outfile, logger, options):

	msg = "Annotation on [%s] in progress. Be patient (30 min+) ..." % vcffile
	print msg

	#tick
	timecnt = 0
	t1 = time.time()

	vcf = VCFParser(vcffile)
	outstream = open(outfile, 'w')

	snpa = SNPAnnotation(options)
	logger.info('Instantiated SNP annotation class..')

	add_meta(vcf, snpa)
	snpa.hgvs.load_resource() #hgvs

	# print refgene_entries
	if snpa.hgmd_on:
		snpa.retrieve_hgmdgenes()

	vcf.writeheader(outstream)

	cDbLink = None
	dblink_fpw = None
	outfile_dblink = None

	if snpa.dblink:
		outfile_dblink = outfile+'.tsv'
		d = 'creating dblink tsv file [%s].' % outfile_dblink
		logger.info(d)
		snpa.cDbLink = DbLink(outfile_dblink)

	snpcnt = 0
	for rec in vcf:
		# print 'chrm[%s],pos[%d],ref[%s]'%(rec['chrom'],rec['pos'],rec['ref']) #debug
		vcf.parseinfo(rec)

		annotations = []
		chrom, pos, ref, alts = (rec['chrom'], rec['pos'], rec['ref'], rec['alt'])

		snpcnt += 1 #assign a new variant index

		annot_info = {'Intergenic': None, 'Genic': {}}
		snpdb_ant = {}
		snpids = []
		regulome_score = []
		gwas_phens = []
		cadd_raw_list = []
		cadd_phred_list = []
		cadd_aa_list = []
		domainlist = []
		par_ant = None
		alt_cnt = 0
		valid_rec = False

		if snpa.dblink:
			snpa.cDbLink.vKeys = {}

		for alt in alts:
			valid_alt = False
			alt_cnt += 1

			for annt in snpa.annotate(chrom, pos, ref, alt):
				if not annt:
					continue
				elif not valid_alt:
					valid_alt = True

				gene = annt.pop(0)
				region = annt[1]
				domains = annt.pop(28)
				if domains:
					domainlist += domains
				cadd_aa = annt.pop(27)
				cadd_phred = annt.pop(26)
				cadd_raw = annt.pop(25)
				par = annt.pop(24)

				if par[1]:
					par_ant = par
				pdomain = annt.pop(23)
				nimb_capture = annt.pop(22)
				snpdbinfo = annt.pop(21)
				dbsnpids = snpdbinfo.pop(0)

				gwas = annt.pop(20)
				regu_score = annt.pop(19)
				lcr_ant = annt.pop(18)
				conserve = annt.pop(17)
				gad = annt.pop(16)
				omim_ids = annt.pop(15)
				omim_phen = annt.pop(14)
				data = ["" if ele is None else str(ele) for ele in annt]
				data.insert(3, str(alt_cnt))
				if region == 'Intergenic':
					annot_info['Intergenic'] = gene
				else:
					if gene in annot_info['Genic']: #create a new holder
						ele = annot_info['Genic'][gene]
						ele.append('|'.join(data))
					else: #keep updating
						annot_info['Genic'][gene] = []
						annot_info['Genic'][gene] = ["" if ele is None else
													 str(ele) for ele in
													 [omim_phen, omim_ids,
													 gad, '|'.join(data)]]
			#end for alt(annt)

			if not valid_alt:#is there at least one valid alt=?
				continue
			elif not valid_rec:
				valid_rec = True

			snpids.append(dbsnpids)
			if cadd_raw is not None:
				cadd_raw_list.append(cadd_raw)
			else:
				cadd_raw_list.append('.')
			if cadd_phred is not None:
				cadd_phred_list.append(cadd_phred)
			else:
				cadd_phred_list.append('.')
			if cadd_aa is not None:
				cadd_aa_list.append(cadd_aa)
			else:
				cadd_aa_list.append('.')
			if not regu_score:
				regulome_score.append('.')
			else:
				regulome_score.append(regu_score[1])
			if not gwas:
				gwas_phens.append('.')
			else:
				gwas_phens.append(gwas)

			for ele in snpdbinfo:
				if isinstance(ele, dict):
					if 'EXAC' in ele:
						exac_ant, exac_acs, exac_afs = ele['EXAC']
						if exac_ant[0] not in snpdb_ant:
							snpdb_ant[exac_ant[0]] = [exac_ant[1]]
						else:
							snpdb_ant[exac_ant[0]].append(exac_ant[1])

						for exac_cf in [exac_acs, exac_afs]:
							for exac in exac_cf: #for EXAC_AC's and EXAC_AF's
								if exac[0] not in snpdb_ant: #for each ethnic group
									snpdb_ant[exac[0]] = [exac[1]]
								else:
									snpdb_ant[exac[0]].append(exac[1])
				else:
					isList = isinstance(ele[1], list)
					if ele[0] not in snpdb_ant:
						if isList:
							snpdb_ant[ele[0]] = ele[1]
						else:
							snpdb_ant[ele[0]] = [ele[1]]
					else:
						if isList:
							snpdb_ant[ele[0]].extend(ele[1])
						else:
							snpdb_ant[ele[0]].append(ele[1])

		#end for alts

		if not valid_rec: continue

		# to write cDbLink to file
		if snpa.dblink:
			snpa.cDbLink.print_vkey(snpcnt)

		if annot_info['Intergenic']:
			annotations += ['VARANT_INTERGENIC', annot_info['Intergenic']]
		if annot_info['Genic']:
			eff = []
			for gene, trans in annot_info['Genic'].iteritems():
				eff.append(gene + '(' + ':'.join(trans[3:] + trans[:3]) + ')')
			annotations += ['VARANT_GENIC', eff]
		if conserve:
			annotations += conserve
		if lcr_ant:
			annotations += lcr_ant
		if len([e for e in regulome_score if e == '.']) != len(regulome_score):
			annotations += ['RegulomeScore', regulome_score]
		if len([e for e in gwas_phens if e == '.']) != len(gwas_phens):
			annotations += ['GWASPhenotype', gwas_phens]
		if nimb_capture:
			annotations += nimb_capture
		if pdomain:
			annotations += pdomain

		if par_ant:
			annotations += par_ant
		dbsnp_version = snpa.dbsnp.get_version()[0].version[1:]
		if True in snpdb_ant['DB%s' % dbsnp_version]:
			annotations += ['DB%s' % dbsnp_version, True]
			annotations += ['dbSNPBuildID', snpdb_ant['dbSNPBuildID']]
			rec['id'] = snpids
		if True in snpdb_ant['KGDB']:
			annotations += ['KGDB', True]
			annotations += ['KGAF', snpdb_ant['KGAF']]
		if True in snpdb_ant['ESPDB']:
			annotations += ['ESPDB', True]
			annotations += ['ESPAF', snpdb_ant['ESPAF']]
		if True in snpdb_ant['EXACDB']:
			annotations += ['EXACDB', True]
			for exac in exacs:
				if exac in snpdb_ant:
					annotations += [exac, snpdb_ant[exac]]

		if snpa.hgmd_on:
			if snpdb_ant['HGMD']:
				annotations += ['HGMD', snpdb_ant['HGMD']]

		if snpdb_ant['CLINVT']:
			annotations += ['CLINVT', snpdb_ant['CLINVT']]

		if snpa.cosmic_on:
			if snpdb_ant['COSMIC']:
				annotations += ['COSMIC', snpdb_ant['COSMIC']]

		if len([e for e in snpdb_ant['CLNACC'] if e == '.']) != \
											len(snpdb_ant['CLNACC']):
			cols = ['CLN' + ele for ele in ['ACC', 'DBN', 'METHOD', 'ORIGIN', 'SIG']]
			cols.extend(['REFTX', 'HGVSc', 'SPLOC', 'HGVSp', 'VLD', 'DATE'])
			for key in cols:
				if len([e for e in snpdb_ant[key] if e == '.']) != len(snpdb_ant[key]):
					annotations += [key, snpdb_ant[key]]

		if len([e for e in cadd_raw_list if e == '.']) != len(cadd_raw_list):
			annotations += ['CADD_raw', cadd_raw_list]
		if len([e for e in cadd_phred_list if e == '.']) != len(cadd_phred_list):
			annotations += ['CADD_phred', cadd_phred_list]
		if len([e for e in cadd_aa_list if e == '.']) != len(cadd_aa_list):
			annotations += ['CADD_aa', cadd_aa_list]
		if domainlist:
			domainlist = list(set(domainlist))
			annotations += ['Interpro_domains', domainlist]

		if snpa.dblink:
			annotations += ['VARIANT_IDX',snpcnt]

		# print annotations
		add_annotation(vcf, rec, annotations)
		vcf.write(outstream, rec, False)
		del annotations

		if time.time() - t1 > 60:
			timecnt += time.time() - t1
			d = 'Annotated - ' + str(snpcnt) + ' variants in ' + \
					'%.2f' % ((timecnt) / 60) + 'min,' + \
					' Currently processing chrom-%s' % chrom.strip('chr')
			logger.info(d)
			t1 = time.time()
	#end for rec
	snpa.hgvs.close_resource() #hgvs
	vcf.stream.close()

	if snpa.dblink:
		snpa.cDbLink.fpw.close()
		d = 'dblink tsv file [%s] is created.' % outfile_dblink
		logger.info(d)

	if time.time() - t1 < 60:
		timecnt += time.time() - t1
		d = 'Annotated - ' + str(snpcnt) + ' variants in ' + \
					'%.2f' % ((timecnt) / 60) + 'min,' + \
					' Currently processing chrom-%s' % chrom.strip('chr')
		logger.info(d)
	if outstream != sys.stdout:
		outstream.close()


if __name__ == '__main__':

	# variant = ['15', 25425643, 'CGG', 'TGG']
	variant = ['17', 2541615, 'G', 'A']
	#variant = ['1', 878230, 'G', 'C']
	# variant = ['1', 246703860, 'ATCT', 'A']
	# variant = ['7', 117307123, 'AGAG', 'A']
	# variant  = ['1',145103947,'TTTTATTTA','TTTTA,TTTTATTTATTTA,T']
	#variant = ['1', 980824, 'G', 'A,C']
	#variant = ['2', 152417133, 'GGC', 'GTT']
	info = {'Genic': {}, 'Intergenic': None}
	p = time.time()
	options = [None, 200, True, True, 0.05, True]
	snpa = SNPAnnotation(options)
	snpa.hgvs.load_resource() #hgvs

	q = time.time()
	print 'Time Taken to instantiate SNPAnnotation=', q - p, 'sec'

	dbobjects = [snpa.clinvardb, snpa.clnsnpdb, snpa.dbsnp, snpa.effpred,
				 snpa.espdb, snpa.exacdb, snpa.gerp, snpa.kgdb, snpa.mirnadb,
				 snpa.nimbcap, snpa.mirnadb, snpa.mimdb, snpa.refgene,
				 snpa.regulomedb, snpa.utrdb, snpa.splicedb]
	dbversion = None
	for obj in dbobjects:
		for ele in obj.get_version():
			if dbversion:
				dbversion += ', ' + ele.dbname.upper() + '(' \
							 + ele.version + ')'
			else:
				dbversion = ele.dbname.upper() + '(' + \
							ele.version + ')'
	print dbversion

	alt_cnt = 1
	par_ant = False
	for annt in snpa.annotate(*variant):
		# print annt
		gene = annt.pop(0)
		region = annt[1]
		domains = annt.pop(27)
		cadd_aa = annt.pop(26)
		cadd_phred = annt.pop(25)
		cadd_raw = annt.pop(24)
		par = annt.pop(23)
		if par[1]:
			par_ant = par
		nimb_capture = annt.pop(22)
		kwnsnp_dbinfo = annt.pop(21)
		gwas = annt.pop(20)
		regulome_score = annt.pop(19)
		lcr_ant = annt.pop(18)
		conserve = annt.pop(17)
		gad = annt.pop(16)
		omim_ids = annt.pop(15)
		omim_phen = annt.pop(14)
		data = ["" if ele is None else str(ele) for ele in annt]
		data.insert(3, str(alt_cnt))

		if region == 'Intergenic':
			info['Intergenic'] = gene
		else:
			if gene in info['Genic']:
				ele = info['Genic'][gene]
				ele.append('|'.join(data))
			else:
				info['Genic'][gene] = [], True
				info['Genic'][gene] = ["" if ele is None else str(ele) for
															ele in [omim_phen, omim_ids, gad,
																		'|'.join(data)]]

	for key, val in info.iteritems():
		if key == 'Genic':
			temp = []
			for gene, trans in info['Genic'].iteritems():
				temp.append(gene + '(' + ':'.join(trans[3:] + trans[:3]) + ')')
			info['Genic'] = ','.join(temp)

	snpa.hgvs.close_resource() #hgvs

	print info
	print cadd_raw
	print cadd_phred
	print cadd_aa
	print domains
	print kwnsnp_dbinfo
	print nimb_capture
	print gwas
	print regulome_score
	print lcr_ant
	print conserve
	print par_ant
