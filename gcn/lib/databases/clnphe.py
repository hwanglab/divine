"""
.. module:: clnphe
		:platform: Unix, Windows, MacOSX
		:synopsis: Class for accessing  information from an sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for retrieving associations information between gene and
phenotypes from GAD and NHGRI-GWAS Catalog database.
"""
from gcn.lib.io import db
from gcn.etc import dbconfig
from collections import namedtuple
import os

class Clnsnp3:
	def __init__(self):
		self.db = {'GAD': {}, 'NHGRI': {}}
		self.rsid2idx = {}
		self.gene2idx = {}

class Clnsnp(db.DB):
		"""Class to retrieve gene and clinical phenotype associations from
		CLNPHESNP database"""

		def __init__(self):
				"""Class initialization

				Argument

						name (string): Name of the database or filename for the database
				"""
				name = 'CLNPHESNP'
				super(Clnsnp, self).__init__()
				if name in dbconfig.DBCONFIG:
						self.load(name=name)
				elif os.path.exists(name):
						self.load(db=name)
				else:
						raise ValueError('No such database %s' % name)

		def load_clnsnp(self):

				clnsnp = Clnsnp3()

				stmt = "SELECT idx, snp, gene, phenotype, refdb from clnsnp where "\
																										"asnstatus = 'Y'"
				results = self.execute(stmt).fetchall()
				for row in results:
						idx, snp, gene, phen, refdb = row[0], row[1], row[2], row[3], str(row[4])

						if snp not in clnsnp.rsid2idx:
							clnsnp.rsid2idx[snp] = []
						clnsnp.rsid2idx[snp].append(idx)

						if gene not in clnsnp.gene2idx:
							clnsnp.gene2idx[gene] = []
						clnsnp.gene2idx[gene].append(idx)

						phen = phen.strip().replace(',', '').replace(' ', '_')\
										.replace('-', '').replace('.', '').replace(')', '')\
										.replace('(', '').replace(';', '')\
										.replace('*', '').replace('|', '').replace('\n','')
						gwas = clnsnp.db[refdb]

						if refdb == 'GAD':
								if snp:
										try:
												data = gwas[(gene, snp)]
												data.add(phen)
										except:
												gwas[(gene, snp)] = set()
												gwas[(gene, snp)].add(phen)
								else:
										try:
												data = gwas[gene]
												data.add(phen)
										except:
												gwas[gene] = set()
												gwas[gene].add(phen)
						elif refdb == 'NHGRI':
								if snp:
										try:
												data = gwas[snp]
												data.add(phen)
										except:
												gwas[snp] = set()
												gwas[snp].add(phen)

				for refid, val in clnsnp.db.iteritems():
						for id, info in val.iteritems():
								val[id] = list(info)

				for snp,idx in clnsnp.rsid2idx.iteritems():
					clnsnp.rsid2idx[snp] = list(set(idx))

				for gene,idx in clnsnp.gene2idx.iteritems():
					clnsnp.gene2idx[gene] = list(set(idx))

				return clnsnp

		def get_gadinfo(self, gene, rsid=None):
				"""Retrieve the clinical phenotype name, GAD accession
				and PubmedID that reported the association for a given gene.
				Args:
						- gene:    Gene Name (HUGO Gene Nomenclature)
						- rsid:    dbSNP accession (This is optional. If snpid not
												given then it retrieves association at gene level)
				Returns:
						- gadinfo:    '|' delimited string. If rsId given output is like -
												(GadSnp or NotGadSnp|phenotype1,phenotype2,|pubmedid1,
												pubmedid2,|gadid1,gadid2,)
												and if rsId not given then output is like -
												(phenotype1,phenotype2,|pubmedid1,pubmedid2,
												|gadid1,gadid2,)
				"""
				stmt = "SELECT snp, phenotype, pubmedid, refdbid from clnsnp where "\
												"refdb = 'GAD' and asnstatus = 'Y' and gene = (?)"
				results = self.execute(stmt, (gene, )).fetchall()
				rsids = set()
				gadids = set()
				phen = set()
				pubmed = set()
				for row in results:
						if row[0]:
								rsids.add(row[0])
						phenotype = row[1].replace(',', '').replace(' ', '_')\
										.replace('-', '').replace('.', '').replace(')', '')\
										.replace('(', '').replace(';', '').replace('|', '__')\
										.replace('*', '').strip().replace('\n','')
						phen.add(phenotype)
						pubmed.add(str(row[2]))
						gadids.add(str(row[3]))
				if len(gadids) == 0:
						gadinfo = None
				else:
						phen = list(phen)
						pubmed = list(pubmed)
						gadids = list(gadids)
						if rsid:
								if rsid in rsids:
										gadinfo = list(phen)
								else:
										gadinfo = list(phen)
						else:
								gadinfo = list(phen)
				return gadinfo

		def get_gwasinfo(self, rsid):
				"""Retrieve the clinical phenotype name
				and PubmedID that reported the association for a given SNP.
				Args:
						- rsid:    dbSNP accession
				Returns:
						- gadinfo:    '|' delimited string.
												(phenotype1,phenotype2,|pubmedid1,pubmedid2,)
				"""
				stmt = "SELECT phenotype, pubmedid from clnsnp where refdb = 'NHGRI' "\
																														"and snp = (?)"
				results = self.execute(stmt, (rsid, )).fetchall()

				phen = set()
				pubmed = []
				for row in results:
						phenotype = row[0].replace(',', '').replace(' ', '_')\
										.replace('-', '').replace('.', '').replace(')', '')\
										.replace('(', '').replace(';', '').strip().replace('\n','')
						phen.add(phenotype)
						pubmed.append(str(row[1]))

				if len(phen) == 0:
						gwasinfo = None
				else:
						gwasinfo = list(phen)
				return gwasinfo

		def get_header(self):
				""" Intializing vcf formatted header information for GWAS data"""

				Gwas_h = namedtuple('Gwas_h', 'id count type desc')

				header = [Gwas_h('GWASPhenotype', '.', 'String', "NHGRI-GWAS"\
												 " phenotypes associated with the Variant.")]
				return header

if __name__ == '__main__':
		clnsnp = Clnsnp()
		print clnsnp.get_version()
		print clnsnp.get_gadinfo('AKT1')
		print clnsnp.get_gadinfo('APP')
		print clnsnp.get_gadinfo('CETP', 'rs708272')
		print clnsnp.get_gwasinfo('rs1926203')
		print clnsnp.get_gwasinfo('rs13069000')
		lcsnp,rec = clnsnp.load_clnsnp()
		print lcsnp['GAD']['AKT1']
		print lcsnp['GAD']['CETP', 'rs708272']
		print lcsnp['NHGRI']['rs13069000']
		print lcsnp['GAD']['EPHA2']
		print clnsnp.get_header()
