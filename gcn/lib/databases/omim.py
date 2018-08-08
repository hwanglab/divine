#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: omim
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing  information from an sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for accessing association among gene, omim phenotype id and omim
phenotype name from a sqlite database
"""
from gcn.lib.io import db
from gcn.etc import dbconfig
from collections import namedtuple
import os


class Omim(db.DB):
    """Class to retrieve gene and omim phenotype associations from
    OMIM database"""

    def __init__(self):
        """Class initialization

        Argument

            name (string): Name of the database or filename for the database
        """
        name = 'OMIM'
        super(Omim, self).__init__()
        if name in dbconfig.DBCONFIG:
            self.load(name=name)
        elif os.path.exists(name):
            self.load(db=name)
        else:
            raise ValueError('No such database %s' % name)

    def load_omim(self):
        omim = {}
        gene2idx = {}
        stmt = 'SELECT idx,gene,mimphe_id,phenotype from mimdis'
        results = self.execute(stmt).fetchall()
        for row in results:

            gene = row[1]
            if gene not in gene2idx:
                gene2idx[gene] = []
            gene2idx[gene].append(row[0])

            phen_id = str(row[2])
            phen = row[3].split(';')[0].replace(',', '').replace(' ', '_')\
                    .replace('-', '').replace('.', '').replace(')', '')\
                    .replace('(', '').replace(':', '_').strip().replace('\n','')
            try:
                phen_list = omim[gene]
                phen_list[0] += '__' + phen
                phen_list[1] += '__' + phen_id
            except:
                omim[gene] = [phen, phen_id]

        #for key, val in omim.iteritems():
        #    omim[key] = ":".join(val)

        for gene in omim.iterkeys():
            omim[gene].append(list(set(gene2idx[gene])))

        return omim

    def get_miminfo(self, gene):
        """Retrieve the omim phenotype id and omim phenotype name for a
        given gene.
        Args:
            - gene:    Gene Name (HUGO Gene Nomenclature)
            -format:    Specify the output format. Default is list of tuple.
                        But when calling this method for input to vcf INFO
                        field one may specify the format = 'vcfinfo'
        Returns:
            - miminfo:    List whose elements are in form -
                          [_ delimited phenotype, _ delimited mimids]
        """
        stmt = 'SELECT mimphe_id,phenotype from mimdis where gene = (?)'
        results = self.execute(stmt, (gene, )).fetchall()
        ids = set()
        phe = set()
        for row in results:
            pheno = row[1].split(';')[0].replace(',', '').replace(' ', '_')\
                    .replace('-', '').replace('.', '').replace(')', '')\
                    .replace('(', '').replace(':', '_').strip().replace('\n','')
            ids.add(str(row[0]))
            phe.add(pheno)
        if len(ids) == 0 and len(phe) == 0:
            miminfo = None
        elif len(ids) > 0:
            miminfo = ['__'.join(list(phe)), '__'.join(list(ids))]
        return miminfo

    def get_mimgene(self, mimpheid):
        """Retrieve the gene names for the given omim phenotype id.

        Args:
            - mimpheid:    OMIM Phenotype Id

        Returns:
            - genelist(list):    A Gene list i.e all the genes clinically
                                associated with the given omim phenotype.
        """
        genes = set()
        stmt = 'SELECT gene from mimdis where mimphe_id = (?)'
        results = self.execute(stmt, (mimpheid, )).fetchall()
        for row in results:
            genes.add(str(row[0]))
        if len(genes) == 0:
            genes = None
        else:
            genes = list(genes)
        return genes


if __name__ == '__main__':
    mim = Omim()
    print mim.tablenames()
    print mim.get_version()
    print mim.get_miminfo('AKT1')
    print mim.get_miminfo('APP')
    print mim.get_miminfo('ACTA2')
    print mim.get_mimgene(608967)
    lmim = mim.load_omim()
    print lmim['ACTA2']
    print lmim['APP']
    print lmim['OXCT1']
