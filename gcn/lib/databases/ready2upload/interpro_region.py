#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: mkomim
	:platform: Unix, Windows, MacOSX
	:synopsis: Script to build a omim phenotype sqlite database.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a omim phenotype sqlite database. To use this script

usage: mkomim.py [-h] [-i INNAME1,INNAME2] [-o OUTNAME] [--force]

Create a OMIM sqlite database

optional arguments:
  -h, --help			show this help message and exit
  -i INNAME1,INNAME2, --input INNAME1,INNAME2
						sprot,trembl files in Uniprot format
  -o OUTNAME, --output OUTNAME
						Name of sqlite database to create
  --force			   Create new database even if the existing one is newer
						than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has two tables - *mimdis*, *release_info*
The *mimdis* table contains the following 4 columns

	- unipacc:	   Uniprot Accession
	- gene:	  Gene
	- mimphe_id: OMIM Phenotype Id
	- phenotype:	OMIM Phenotype

The *release_info* table contains the following 3 columns

	- file_name:	File Name
	- file_date:	The date when the file was downloaded from website
	- version:	Version Number
	- entry_count:	Number of entries in the mimdis table of sqlite DB
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import anyopen, db
from collections import namedtuple
import csv
from gcn.lib.databases import lib_mysql
from gcn.lib.databases.snpdb import ClinvarDB, HgmdDB
from gcn.lib.databases.refgene import to_grach
from gcn.lib.utils.lib_utils import py_struct
import argparse
import sys
import os
import time
from datetime import datetime as dt

db_fields = ('idx','chrom','cdsStart','cdsEnd','gene','interpro_id','domain_desc','benign_lof','benign_nsm','vus_lof','vus_nsm','pathogenic_lof','pathogenic_nsm')

db_fields_norm = ('idx','domain_idx','benign_lof','benign_nsm','vus_lof','vus_nsm','pathogenic_lof','pathogenic_nsm')

SCHEMA = ["""create table domain
				 		(%s integer PRIMARY KEY AUTOINCREMENT,
							%s text not null,
							%s integer not null,
							%s integer not null,
							%s text not null,
				 			%s text,
				 			%s text,
				 			%s integer,
				 			%s integer,
				 			%s integer,
				 			%s integer,
				 			%s integer,
				 			%s integer)
            """%db_fields,
            """create index index1 on domain (gene,interpro_id,domain_desc)""",
            """create index index2 on domain (chrom,cdsStart,cdsEnd)""",
						"""create table density
				 		(%s integer PRIMARY KEY AUTOINCREMENT,
				 			%s integer NOT NULL,
				 			%s float,
				 			%s float,
				 			%s float,
				 			%s float,
				 			%s float,
				 			%s float,
				 			FOREIGN KEY (domain_idx) REFERENCES domain (idx))
            """%db_fields_norm,
            """create index index3 on density (domain_idx)"""
          ]

class InterproReport(py_struct):
	def __init__(self):
		super(py_struct, self).__init__()
		self.name = 'interpro'

		self.found = False
		self.chrom = None
		self.spos = None
		self.epos = None
		self.genes = []
		self.denom = []
		self.benign_dens = []
		self.vus_dens = []
		self.patho_dens = []
		self.benign_lof_rates = []
		self.benign_nsm_rates = []
		self.vus_lof_rates = []
		self.vus_nsm_rates = []
		self.patho_lof_rates = []
		self.patho_nsm_rates = []
		self.desc = []
		self.indices = []

	def to_vcf_info(self):

		if self.desc:
			return [self.name, list(set(self.desc))]
		else:
			return None

class Interpro(db.DB):
	def __init__(self):
		"""Class initialization

		Argument
				name (string): Name of the database or filename for the database
		"""
		self.name = 'INTERPRO'
		super(Interpro, self).__init__()
		if self.name in dbconfig.DBCONFIG:
			self.load(name=self.name)
		elif os.path.exists(self.name):
			self.load(db=self.name)
		else:
			raise ValueError('No such database %s' % self.name)

		self.report = InterproReport()

	def get_patho_density_by_pos(self,chrom,spos,epos):
		# collect func domain density

		#self.report.reset_attr()
		# self.report.set_attr(chrom=chrom,
		# 											spos=spos,
		# 											epos=epos)

		stmt = """
			select
				dom.gene as gene,
				(sum(dom.benign_lof) + sum(dom.benign_nsm) + sum(dom.vus_lof) + sum(dom.vus_nsm) + sum(dom.pathogenic_lof) + sum(dom.pathogenic_nsm)) as denom,
				(sum(den.benign_lof) + sum(den.benign_nsm)) as benign,
				(sum(den.vus_lof) + sum(den.vus_nsm)) as vus,
				(sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as patho
			from domain as dom
				inner join density as den on den.domain_idx = dom.idx
			where dom.chrom = (?) and ((dom.cdsStart<=(?) and (?)<dom.cdsEnd) or (dom.cdsStart<=(?) and (?)<dom.cdsEnd) or ((?)<=dom.cdsStart and dom.cdsEnd<(?)))
			group by gene """

		fields = "gene denom benign vus patho".split()

		denstup = namedtuple('denstup', fields)

		results = self.execute(stmt,(chrom.strip('chr'),spos,spos,epos,epos,spos,epos))
		genes = []
		denoms = []
		benign_dens = []
		vus_dens = []
		patho_dens = []

		for row in results:
			rec = denstup._make(row)

			if int(rec.denom) > 0:
				genes.append(rec.gene)
				denoms.append(rec.denom)
				benign_dens.append(rec.benign)
				vus_dens.append(rec.vus)
				patho_dens.append(rec.patho)

		return(genes,denoms,benign_dens,vus_dens,patho_dens)

	def get_patho_density_per_gene(self):
		# collect func domain density

		self.report.reset_attr()

		stmt = """
			select
				dom.gene as gene,
					(sum(dom.benign_lof) + sum(dom.benign_nsm) + sum(dom.vus_lof) + sum(dom.vus_nsm) + sum(dom.pathogenic_lof) + sum(dom.pathogenic_nsm)) as denom,
					sum(den.benign_lof)/(sum(den.benign_lof) + sum(den.benign_nsm) + sum(den.vus_lof) + sum(den.vus_nsm) + sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as benign_lof_r,
					sum(den.benign_nsm)/(sum(den.benign_lof) + sum(den.benign_nsm) + sum(den.vus_lof) + sum(den.vus_nsm) + sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as benign_nsm_r,
					sum(den.vus_lof)/(sum(den.benign_lof) + sum(den.benign_nsm) + sum(den.vus_lof) + sum(den.vus_nsm) + sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as vus_lof_r,
					sum(den.vus_nsm)/(sum(den.benign_lof) + sum(den.benign_nsm) + sum(den.vus_lof) + sum(den.vus_nsm) + sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as vus_nsm_r,
					sum(den.pathogenic_lof)/(sum(den.benign_lof) + sum(den.benign_nsm) + sum(den.vus_lof) + sum(den.vus_nsm) + sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as patho_lof_r,
					sum(den.pathogenic_nsm)/(sum(den.benign_lof) + sum(den.benign_nsm) + sum(den.vus_lof) + sum(den.vus_nsm) + sum(den.pathogenic_lof) + sum(den.pathogenic_nsm)) as patho_nsm_r
				from domain as dom
					inner join density as den on den.domain_idx = dom.idx
				group by gene
		"""

		fields = "gene denom benign_lof_r benign_nsm_r vus_lof_r vus_nsm_r patho_lof_r patho_nsm_r".split()

		denstup = namedtuple('denstup', fields)

		results = self.execute(stmt)

		for row in results:
			rec_tup = denstup._make(row)

			if int(rec_tup.denom) > 0:
				self.report.genes.append(rec_tup.gene)
				self.report.benign_lof_rates.append(rec_tup.benign_lof_r)
				self.report.benign_nsm_rates.append(rec_tup.benign_nsm_r)
				self.report.vus_lof_rates.append(rec_tup.vus_lof_r)
				self.report.vus_nsm_rates.append(rec_tup.vus_nsm_r)
				self.report.patho_lof_rates.append(rec_tup.patho_lof_r)
				self.report.patho_nsm_rates.append(rec_tup.patho_nsm_r)

	def _iter_pathogenicty_per_gene(self):

		stmt = """
			select
				dom.gene as gene,
				sum(den.benign_lof) as benign_lof_dens,
				sum(den.benign_nsm) as benign_nsm_dens,
				sum(den.vus_lof) as vus_lof_dens,
				sum(den.vus_nsm) as vus_nsm_dens,
				sum(den.pathogenic_lof) as patho_lof_dens,
				sum(den.pathogenic_nsm) as patho_nsm_dens
			from domain as dom
				inner join density as den on den.domain_idx = dom.idx
			group by gene;
		"""

		fields = "gene total_region_len benign_lof_pct benign_nsm_pct vus_lof_pct vus_nsm_pct patho_lof_pct patho_nsm_pct".split()

		domtup = namedtuple('domtup', fields)

		results = self.execute(stmt)
		for row in results:
				yield domtup._make(row)

	def get_pathogenicity_per_gene(self):

		if self.report.prof:
			return

		self.report.set_attr(prof={})

		for domtup in self._iter_pathogenicty_per_gene():
			self.report.prof[domtup.gene] = domtup


	def get_domain_desc(self,chrom,spos,epos):
		"""Retrieves Interpro domains that overlaps with the given variant
		positions"""
		self.report.reset_attr()
		self.report.set_attr(chrom=chrom,
												 spos=spos,
												 epos=epos)

		stmt = "SELECT idx, domain_desc FROM domain where chrom = '%s' and ((cdsStart<=%d and %d<cdsEnd) or (cdsStart<=%d and %d<cdsEnd) or (%d<=cdsStart and cdsEnd<%d))" % (chrom.strip('chr'),spos,spos,epos,epos,spos,epos)

		results = self.execute(stmt)
		domains = []

		for row in results:
			e = row[1]
			if e:
				self.report.indices.append(row[0])
				domains += e.split(',')

		dom_desc = [e.replace('/', '_').replace(';', '_').replace('-', '_') for e in domains]

		self.report.desc.extend(dom_desc)

class InterproDB(db.DB):
	"""Class to insert OMIM data into a sqlite3 database"""

	def createdb(self, inname, outname, force, hgmd_on = 0):
		"""Create the database

		Args:
			inname (str): interpro file name from ucsc SQL results
			outname (str): Name of sqlite3 database
			force (bool): If True overwrite existing database even if it
						  is newer than any of  `inname`
		"""
		self.logger = fileconfig.getlogger()
		if os.path.exists(outname):
			newer = fileutils.file_newer(inname, outname)
		else:
			newer = True

		if not newer and not force:
			self.logger.info('Not Updating. INTERPRO database already uptodate')
		else:
			self.inname = inname
			self.outname = outname
			self.hgmd_on = hgmd_on
			t1 = time.time()
			self._makedb()
			time_taken = (time.time() - t1) / 60
			self.logger.info("Time Taken for creating %s is %f min"
							 % (self.outname, time_taken))

		return

	def _iterfile(self):

		"""Internal method. Do not use"""

		fields = "ucsc_gene chrom start_bp end_bp gene domain_id domain_desc".split()
		domtup = namedtuple('domtup', fields)

		clnStat = ClinvarDB()
		hgmdStat = HgmdDB()

		stream = anyopen.openfile(self.inname)
		stream.next()  # skip header

		processed_regions = {}

		for rec in csv.reader(stream, delimiter='\t'):

			rec[1] = to_grach(rec[1])
			if rec[1] is None: continue

			domain = domtup._make(rec)

			region_key = '%s_%s_%s' % (domain.chrom, domain.start_bp, domain.end_bp)

			if region_key in processed_regions: continue
			processed_regions[region_key] = True

			# query to clinvar
			vartypes = clnStat.count_lofs(goi_tup=domain)

			# search for hgmd
			if self.hgmd_on:
				vartypes = hgmdStat.count_lofs(domain, vartypes)

			# count
			vcounts = [[0, 0], [0, 0], [0, 0]]
			for csigs in vartypes.itervalues():
				for i, csig in enumerate(csigs[:-1]):
					for j, vtype in enumerate(csig):
						vcounts[i][j] += vtype

			yield domain, vcounts

		stream.close()

		clnStat.conn.close()
		if self.hgmd_on:
			hgmdStat.conn.close()

	def _makedb_domain_norm(self):

		curs = self.conn.cursor()

		denom_stmt = """
			select
				sum(1.*benign_lof/(cdsEnd-cdsStart)) as benign_lof_denom,
				sum(1.*benign_nsm/(cdsEnd-cdsStart)) as benign_nsm_denom,
				sum(1.*vus_lof/(cdsEnd-cdsStart)) as vus_lof_denom,
				sum(1.*vus_nsm/(cdsEnd-cdsStart)) as vus_nsm_denom,
				sum(1.*pathogenic_lof/(cdsEnd-cdsStart)) as patho_lof_denom,
				sum1.*(pathogenic_nsm/(cdsEnd-cdsStart)) as patho_nsm_denom
			from domain
		"""

		curex = curs.execute(denom_stmt)
		psum = [row for row in lib_mysql.iter_tupfetchall(curex)][0]

		dens_stmt = """
			select
				idx,
				1.*benign_lof/(cdsEnd-cdsStart)/(?) as benign_lof,
				1.*benign_nsm/(cdsEnd-cdsStart)/(?) as benign_nsm,
				1.*vus_lof/(cdsEnd-cdsStart)/(?) as vus_lof,
				1.*vus_nsm/(cdsEnd-cdsStart)/(?) as vus_nsm,
				1.*pathogenic_lof/(cdsEnd-cdsStart)/(?) as patho_lof,
				1.*pathogenic_nsm/(cdsEnd-cdsStart)/(?) as patho_nsm
			from domain
		"""

		curex = curs.execute(dens_stmt, (
			psum.benign_lof_denom,
			psum.benign_nsm_denom,
			psum.vus_lof_denom,
			psum.vus_nsm_denom,
			psum.patho_lof_denom,
			psum.patho_nsm_denom))

		query = "insert into density (%s) values (%s)" % (','.join(db_fields_norm[1:]),','.join('?'*7))

		n = 0
		entry_cnt = 0
		entries = []

		for row in lib_mysql.iter_tupfetchall(curex):

			entries.append((row.idx,\
											row.benign_lof,row.benign_nsm,\
											row.vus_lof,row.vus_nsm,\
											row.patho_lof,row.patho_nsm))

			n += 1
			if n == 100:
				# print entries #debug
				curs.executemany(query, entries)
				n = 0
				entry_cnt += 100
				entries = []

		if entries:
			curs.executemany(query, entries)
			entry_cnt += len(entries)

		curs.close()

	def _makedb(self):
		"""Internal method. Do not use"""

		self.logger.info('Creating INTERPRO database ...')
		self.logger.info('Input files: %s' % (self.inname))

		if not os.path.exists(self.inname):
			self.logger.error('%s: No such file' % self.inname)
			self.logger.error('Database not created')
			sys.exit(1)

		self.load(db=self.outname)

		for s in SCHEMA:
			print s #debug
			self.createtable(s, True)

		curs = self.conn.cursor()

		query = "insert into domain (%s) values (%s)" % (','.join(db_fields[1:]), ','.join('?' * 12))

		n = 0
		entry_cnt = 0
		entries = []
		for dom,vc in self._iterfile():
			#[[benign_lof,benign_nsm],[vus_lof,vus_nsm],[pathogenic_lof,pathogenic_nsm]]
			entries.append((dom.chrom, dom.start_bp, dom.end_bp, dom.gene, dom.domain_id, dom.domain_desc, vc[0][0],vc[0][1],vc[1][0],vc[1][1],vc[2][0],vc[2][1]))

			n += 1
			if n == 100:
				#print entries #debug
				curs.executemany(query, entries)
				n = 0
				entry_cnt += 100
				entries = []

		if entries:
			curs.executemany(query, entries)
			entry_cnt += len(entries)

		curs.close()

		self._makedb_domain_norm()

		# Add version details
		fd = dt.fromtimestamp(os.path.getmtime(self.inname)).strftime('%Y-%m-%d')
		version = "v%s_%s" % tuple(fd.split('-')[:2])
		self.set_version('Interpro', fd, version, entry_cnt)

		self.conn.commit()
		self.conn.close()
		self.logger.info('... INTERPRO database created')

def main():
	"""Main script to create the OMIM database"""

	infile = fileconfig.FILECONFIG['INTERPRO']
	outfil = dbconfig.DBCONFIG['INTERPRO']['name']

	desc = 'Create INTERPRO sqlite database'
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('-i', '--input', dest='inname', type=str,
						default=infile,
						help='Interpro tab delimited file')
	parser.add_argument('-o', '--output', dest='outname', type=str,
						default=outfil,
						help='Name of sqlite database to create')
	parser.add_argument('-H', '--hgmd', dest='hgmd',
											required=False, type=int, default=0,
											help='enable HGMD (requires a license), [0]:No, 1:Yes')
	parser.add_argument('--force', action='store_true',
						help='Create new database even if the existing one ' +
							 'is newer than *inname*')
	args = parser.parse_args()
	import time
	t1 = time.time()
	InterproDB().createdb(args.inname, args.outname, args.force, hgmd_on=args.hgmd)
	print 'TIME TAKEN : ', (time.time() - t1) / 60, 'min'
	sys.exit(0)


if __name__ == "__main__":
	main()
