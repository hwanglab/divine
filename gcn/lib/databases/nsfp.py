#
# COPYRIGHT (C) 2002-2012 Kunal Kundu
#
"""
.. module:: NSFP
		:platform: Unix, Windows, MacOSX
		:synopsis: Class for accessing  information from the NSFP sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for accessing POLYPHEN2/SIFT predictions for given hg19 genomic
coordinate from the NSFP sqlite database
"""
from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
from collections import namedtuple
from gcn.etc import dbconfig
import os


class Effpred(db.DB):
		"""Class to retrieve polyphen/sift effect predictions for given
		hg19 genomic coordinate from NSFP database"""

		def __init__(self):
				"""Class initialization

				Argument

						name (string): Name of the database or filename for the database
				"""
				name = 'NSFP'
				super(Effpred, self).__init__()
				if name in dbconfig.DBCONFIG:
						self.load(name=name)
				elif os.path.exists(name):
						self.load(db=name)
				else:
						raise ValueError('No such database %s' % name)

		def __call__(self, chrom, pos, wtaa, altaa, aapos, ref=None,
								 alt=None, gene=None):
				if chrom:
						self.chrom = chrom.strip('chr')
				else:
						self.chrom = chrom
				self.pos = pos
				self.ref = ref
				self.alt = alt
				self.wtaa = wtaa
				self.altaa = altaa
				self.aapos = aapos
				self.gene = gene
				self.pp2_names = {'B': 'PP2B', 'P': 'PP2PD', 'D': 'PP2D'}
				self.sift_names = {'T': 'T', 'D': 'D'}
				self.pp2_pred = None
				self.pp2_score = None
				self.sift_score = None
				self.sift_pred = None
				self.rec_idx = None
				self._predict()

		def _predict(self):
				"""Retrieve the polyphen/sift predicted effect and scores
				and assigns to the the class variable for the given
				isoform. (Internal Method)"""

				if self.ref and self.alt:
						stmt = 'SELECT idx,sift_score,sift_pred,pp2_hvar_score,pp2_hvar_pred \
						from nsfp where chr = (?) and pos = (?) and ref = (?) and alt = (?)\
						 and refaa = (?) and altaa = (?) and aapos = (?) limit 1'

						results = self.execute(stmt, (self.chrom, self.pos, self.ref,
											self.alt, self.wtaa, self.altaa, self.aapos,
											)).fetchall()
				else:
						stmt = 'SELECT idx,sift_score,sift_pred,pp2_hvar_score,pp2_hvar_pred \
						from nsfp where gene = (?) \
						 and refaa = (?) and altaa = (?) and aapos = (?) order by \
						 pp2_hvar_score desc limit 1'
						results = self.execute(stmt, (self.gene,
														self.wtaa, self.altaa, self.aapos, )).fetchall()
				if len(results) > 0:
						self.sift_score = results[0][1]
						if self.sift_score is not None:
								self.sift_pred = self.sift_names[str(results[0][2])]
						self.pp2_score = results[0][3]
						if self.pp2_score is not None:
								self.pp2_pred = self.pp2_names[results[0][4].strip()]
						self.rec_idx = results[0][0]

		def get_domains(self, chrom, start_pos, end_pos):
				"""Retrieves Interpro domains that overlaps with the given variant
				positions"""
				domains = []
				indices = []
				stmt = "SELECT idx, domains FROM nsfp WHERE chr='%s' AND pos "\
							 "BETWEEN %d and %d" % (chrom.strip('chr'), start_pos, end_pos)
				results = self.execute(stmt)
				for row in results:
						e = row[1]
						if e:
								indices.append(row[0])
								domains += e.split(',')
				domains = list(set(domains))
				domains = [e.replace('/', '_').replace(';', '_').replace('-', '_')
									 for e in domains]
				return list(set(indices)),domains

		def get_cadd(self, chrom, pos, ref, alt):
				"""Retrieves CADD raw and phred scores for given variant position,
				reference and alternate allele"""
				idx, cadd_raw, cadd_phred, refaa,  altaa = (None, None, None, None, None)
				stmt = "SELECT idx, cadd_raw, cadd_phred, refaa, altaa FROM nsfp WHERE chr='%s' \
								AND pos=%d AND ref='%s' AND alt='%s' limit 1" \
								% (chrom.strip(), pos, ref, alt)
				results = self.execute(stmt).fetchall()
				if len(results) > 0:
						idx, cadd_raw, cadd_phred, refaa, altaa = results[0]
				return idx, cadd_raw, cadd_phred, refaa, altaa

		def get_amino_acid(self, chrom, pos, ref, alt):
				"""Retrieves refaa, altaa for given variant position,
				reference and alternate allele"""

				refaa, altaa = (None, None)
				stmt = "SELECT refaa  altaa from nsfp WHERE chr='%s' \
									AND pos=%d AND ref='%s' AND alt='%s' limit 1" \
									% (chrom.strip(), pos, ref, alt)
				results = self.execute(stmt).fetchall()
				if len(results) > 0:
						refaa, altaa = results[0]
				return refaa, altaa

if __name__ == '__main__':
		eff = Effpred()
		print eff.get_version()
		eff('chr1', 69115, 'G', 'S', 9, 'G', 'A')
		print eff.pp2_pred, eff.pp2_score, eff.sift_pred, eff.sift_score
		eff('chr1', 69119, 'L', 'H', 10, 'T', 'A')
		print eff.pp2_pred, eff.pp2_score, eff.sift_pred, eff.sift_score
		eff('22', 16266932, 'H', 'P', 506, None, None, 'POTEH')
		print eff.pp2_pred, eff.pp2_score, eff.sift_pred, eff.sift_score
		eff(None, None, 'P', 'S', 245, None, None, 'ACTA2')
		print eff.pp2_pred, eff.pp2_score, eff.sift_pred, eff.sift_score
		domains = eff.get_domains('11', 200000, 270000)
		print domains
		domains = eff.get_domains('1', 879297, 879297)
		print domains
		raw, phred, refaa, altaa = eff.get_cadd('11', 212428, 'A', 'G')
		print raw, phred, refaa, altaa
