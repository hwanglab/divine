"""
.. module:: splicedb
		:platform: Unix, Windows, MacOSX
		:synopsis: Class for accessing  information from an sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for accessing splicedb information from an sqlite database
"""
from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
from gcn.etc import dbconfig
import os


class Splicedb(db.DB):
		"""Class to splice annotation from SPLICE database"""

		def __init__(self):
				"""Class initialization

				Argument

						name (string): Name of the database or filename for the database
				"""
				name = 'SPLICE'
				super(Splicedb, self).__init__()
				if name in dbconfig.DBCONFIG:
						self.load(name=name)
				elif os.path.exists(name):
						self.load(db=name)
				else:
						raise ValueError('No such database %s' % name)

		def get_header(self):
				info = "##INFO=<ID=SpliceESE,Number=0,Type=Flag,Description=\"If \
				variant is present in Exonic Splice Enhancer Site \">\n"
				info += "##INFO=<ID=SpliceESS,Number=0,Type=Flag,Description=\"If \
				variant is present in Exonic Splice Silencer Site \">\n"
				return info

		def get_annot(self, chrom, position, refmrna_acc):
				"""Retrieve lod name and score for the given genomic coordinate.
				Args:
						chr (str):   Chromosome Number (starting with 'chr' Eg. chr1,
												 chr2 etc)
						position (integer):    Chromosome Position
						refmrna_acc(str):    RefSeq mRNA Accession (Starts with NM_)

				Returns:
						- annotation:    If the variant is in splice enhancer/silencer
														 site then returns either SpliceESS or SpliceESE
														 else None
				"""
				if not chrom.startswith('chr'):
						chrom = 'chr' + chrom
				stmt = 'select idx, annotation from splice where chr=(?) and pos=(?)\
				 and refseq_acc=(?)'
				results = self.execute(stmt, (chrom, position, \
																			refmrna_acc, )).fetchall()
				idx = None
				if results:
						if results[0][1] == 'ESE':
								annotation = 'ESE'
								idx = results[0][0]

						elif results[0][1] == 'ESS':
								annotation = 'ESS'
								idx = results[0][0]
				else:
						annotation = None
				return annotation,idx

if __name__ == '__main__':
		print Splicedb().get_version()
		print Splicedb().get_annot('10', 23605935, 'NM_153714')
