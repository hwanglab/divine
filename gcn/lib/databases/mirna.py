#
# COPYRIGHT (C) 2002-2011 Kunal Kundu
#
"""
.. module:: mirna
		:platform: Unix, Windows, MacOSX
		:synopsis: Class for accessing mirna information from an sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for accessing mirna information from an sqlite database
""" 
from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
from collections import namedtuple
from gcn.etc import dbconfig
import os

mirna = namedtuple('mirna', "idx,mirbase_acc,mirna,\
refseq_acc,gene_align_seq,coord_info,energy,mirsvr_score")


class MiRna(db.DB):
		"""Class to retrieve gene information from REFGENE database"""

		def __init__(self):
				"""Class initialization

				Argument

						name (string): Name of the database or filename for the database
				"""
				name = 'MIRNA'
				super(MiRna, self).__init__()
				if name in dbconfig.DBCONFIG:
						self.load(name=name)
				elif os.path.exists(name):
						self.load(db=name)
				else:
						raise ValueError('No such database %s' % name)

		def retrieve(self, refseq_acc):
				"""Retrieve details of all genes associated with a given
				chromosomal location.

				Args:
						refseq_acc (str):   Refseq Accession which starts with 'NM_'

				Returns:
						List of mirna associated with the refseq accession. Each mirna is
						 a named tuple with the following fields
								- idx               index
								- mirbase_acc       mirbase accession
								- mirna:            miRNA Name
								- refseq_acc:       Refseq accession
								- gene_align_seq:   Aligned sequence
								- coord_info:       Genome coordinate
								- energy:           Binding energy
								- mirsvr_score:     MirSvr Algo score
				"""

				stmt = 'SELECT * FROM mirna WHERE refseq_acc=(?)'
				results = self.execute(stmt, (refseq_acc,)).fetchall()
				mirnadata_list = []
				for r in results:
						mirnadata_list.append(mirna(r[0], r[1], r[2], r[3], r[4], \
																				r[5], r[6], r[7]))
				return mirnadata_list

		def get_annot(self, chrom, pos, refseq_acc):
				mirna = []
				indices = []
				if chrom[:3] == 'chr':
						chrom = chrom[3:]
				mirnadata_list = self.retrieve(refseq_acc)
				for data in mirnadata_list:
						info = data.coord_info.split(':')
						coordlist = info[2].split(',')
						poslist = []
						for coord in coordlist:
								x, y = coord.split('-')
								poslist += range(int(x), int(y) + 1)
						if pos in poslist:
								fltr_poslist = []
								for neu, pos in zip(data.gene_align_seq, poslist):
										if neu.isupper():
												fltr_poslist.append(pos)
								if pos in fltr_poslist:
										mirna.append(str(data.mirna))
										indices.append(data.idx)

				return mirna,indices


if __name__ == '__main__':
		mr = MiRna()
		print mr.get_version()
		mirna = mr.get_annot('8', 39775404, 'CR591703')
		print mirna
		#for e in mr.retrieve('NM_147191'):
		#    print e
