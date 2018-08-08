"""
from gcn.etc.dbconfig import DBCONFIG
clinvardb  clinvitae  clnphesnpdb  cosmic  dbsnp  esp  exac  geneontology  hgmd  kgdb  mimdb  mirna  nsfpdb  refgene  refmrna  regulomedb  splicedb  utrdb
"""

from collections import namedtuple
from gcn.lib.utils import lib_utils

def get_vkey(chrom,pos,ref,alt):
	#assembly = 'b37'
	return '%s_%s_%s_%s' % (chrom, pos, ref, alt)

class VariantKey:
	def __init__(self, chrom, pos, ref, alt):
		#self.assembly = 'b37'
		self.chrom = chrom
		self.pos = pos
		self.ref = ref
		self.alt = alt

		self.rsids = []
		self.genes = []
		self.dbs = {}
		self.initialize_db()

	def initialize_db(self):
		DBINDEX = [['EXAC', 'snps', 'idx', []],
		           ['SNPDB', 'snps', 'idx', []],
		           ['KGDB', 'snps', 'idx', []],
		           ['ESP', 'snps', 'idx', []],
		           ['CLINVARDB', 'snps', 'idx', []],
		           ['CLINVITAE', 'snps', 'idx', []],
		           ['REGULOME', 'regulome', 'idx', []],
		           ['CLNPHESNP', 'clnsnp', 'idx', []],
		           ['HGMDDB', 'snps', 'idx', []],
		           ['COSMIC', 'snps', 'idx', []],
		           ['NSFP', 'nsfp', 'idx', []],
		           ['SPLICE', 'splice', 'idx', []],
		           ['MIRNA', 'mirna', 'idx', []],
		           ['OMIM', 'omim', 'idx', []]]
		for tbFmt in DBINDEX:
			self.dbs[tbFmt[0]] = tbFmt[-1]

class DbLink:
	def __init__(self,out_fn):
		self.vKeys = {}
		self.tsv = out_fn
		self.fpw = self.print_header(out_fn)

	def add_key(self,chrom,pos,ref,alt,\
	                  db_name=None,primary_idx=None,\
	                  rsid = None):

		vkey = get_vkey(chrom,pos,ref,alt)

		if vkey not in self.vKeys:
			self.vKeys[vkey] = VariantKey(chrom,pos,ref,alt)

		if db_name and primary_idx:
			self.vKeys[vkey].dbs[db_name]=primary_idx

		if rsid:
			self.vKeys[vkey].rsids.append(rsid)

	def add_snpid(self,rsid,db_name,primary_index):
		for vkey in self.vKeys.iterkeys():
			if rsid in self.vKeys[vkey].rsids:
				self.vKeys[vkey].dbs[db_name]=primary_index
				break

	def add_genes(self,gene):
		for vkey in self.vKeys.iterkeys():
			self.vKeys[vkey].genes.append(gene)

	def add_primary_idx(self,db_name,primary_index):
		for vkey in self.vKeys.iterkeys():
			if isinstance(self.vKeys[vkey].dbs[db_name],list):
				if isinstance(primary_index,list):
					self.vKeys[vkey].dbs[db_name].extend(list(set(primary_index)))
				else:
					self.vKeys[vkey].dbs[db_name].append(primary_index)
			else:
				self.vKeys[vkey].dbs[db_name]=primary_index

	def print_header(self,tsv):
		DBINDEX = [['EXAC', 'snps', 'idx', []],
		           ['SNPDB', 'snps', 'idx', []],
		           ['KGDB', 'snps', 'idx', []],
		           ['ESP', 'snps', 'idx', []],
		           ['CLINVARDB', 'snps', 'idx', []],
		           ['CLINVITAE', 'snps', 'idx', []],
		           ['REGULOME', 'regulome', 'idx', []],
		           ['CLNPHESNP', 'clnsnp', 'idx', []],
		           ['HGMDDB', 'snps', 'idx', []],
		           ['COSMIC', 'snps', 'idx', []],
		           ['NSFP', 'nsfp', 'idx', []],
		           ['SPLICE', 'splice', 'idx', []],
		           ['MIRNA', 'mirna', 'idx', []],
		           ['OMIM', 'omim', 'idx', []]]

		heads = ['variant_index','chrom','pos','ref','alt','rsid']

		fpw = open(tsv,'w')
		for dbindex in DBINDEX:
			heads.append('%s'%('.'.join(dbindex[:-1])))
		fpw.write('#%s\n'%(lib_utils.joined(heads,'\t')))

		return fpw

	def print_vkey(self,snpcnt):
		DBINDEX = ['EXAC','SNPDB','KGDB','ESP','CLINVARDB','CLINVITAE',\
		           'REGULOME','CLNPHESNP','HGMDDB',\
		           'COSMIC','NSFP','SPLICE','MIRNA','OMIM']

		if self.vKeys:
			for cVkey in self.vKeys.itervalues():

				cols = [snpcnt,cVkey.chrom,cVkey.pos,cVkey.ref,cVkey.alt,lib_utils.joined(cVkey.rsids,',')]

				for tb in DBINDEX:
					tbIdx = cVkey.dbs[tb]
					if isinstance(tbIdx,list):
						tbIdx = lib_utils.joined(list(set(tbIdx)), ',')
						if not tbIdx:
							tbIdx = 'NULL'
						cols.append(tbIdx)
					else:
						if not tbIdx:
							tbIdx = 'NULL'
						cols.append(tbIdx)

				self.fpw.write('%s\n'%(lib_utils.joined(cols,'\t')))
				del cols

		self.cleanup()

	def cleanup(self):
		DBINDEX = ['EXAC','SNPDB','KGDB','ESP','CLINVARDB','CLINVITAE',\
		           'REGULOME','CLNPHESNP','HGMDDB',\
		           'COSMIC','NSFP','SPLICE','MIRNA','OMIM']

		if self.vKeys:
			for cVkey in self.vKeys.itervalues():
				for tb in DBINDEX:
					del cVkey.dbs[tb]
				cVkey.dbs.clear()
				del cVkey.rsids
				del cVkey.genes
				cVkey.rsids = []
				cVkey.genes = []

			self.vKeys.clear()
			self.vKeys = {}


if __name__ == '__main__':

	#initialize
	outfile = '/tmp/test'
	outfile_dblink = outfile + '.tsv'
	d = 'creating dblink tsv file [%s].' % outfile_dblink
	print d

	cDbLink = DbLink(outfile_dblink)

	chrom,pos,ref,alt = ['chr1', '1234', 'A', 'C']
	cDbLink.add_key(chrom, pos, ref, alt)

	cDbLink.add_key(chrom, pos, ref, alt, 'EXAC', '4')

	chrom, pos, ref, alt = ['chr1', '1234', 'A', 'G']
	cDbLink.add_key(chrom, pos, ref, alt, 'SNPDB', '4',rsids = ['rs1234','rs3456'])

	cDbLink.add_snpid('rs1234', 'REGULOME', '6')
	cDbLink.add_snpid('rs3456', 'CLNPHESNP', '9')

	cDbLink.add_primary_idx('SPLICE', '3')

	cDbLink.add_genes('CA3AR1')

	cDbLink.print_vkey(1)

	cDbLink.fpw.close()
	d = 'Done [%s].' % outfile_dblink
	print d

