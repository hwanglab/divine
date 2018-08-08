#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: refmrna
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing  information from an sqlite database

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Class for accessing refmrna information from an sqlite database
"""
from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
from collections import namedtuple
import os

mrna = namedtuple('mrna', 'idx, refseq_acc, sequence')


class Refmrna(db.DB):
    """Class to retrieve seq information from REFMRNA database"""

    def __init__(self, candgene=[]):
        """Class initialization

        Argument

            name (string): Name of the database or filename for the database
        """
        name = 'REFMRNA'
        self.candgene = candgene
        super(Refmrna, self).__init__()
        if name in DBCONFIG:
            self.load(name=name)
        elif os.path.exists(name):
            self.load(db=name)
        else:
            raise ValueError('No such database %s' % name)

    def retrieve(self, refseq_acc, gid):
        """Retrieve mrna seq for the given refseq accession.

        Args:
            refseq_acc (str):   RefSeq Accession (starts with 'NM_')
            gid (integer):    Unique identifier of each transcript
                              This is added to the database beacuse there
                              are cases where different gene definition
                              exists with same refseq accession. For eg.
                              for NM_001005277 there exits 3 transcript
                              definition in refGene file.

        Returns:
            Named tuple with the following fields

                - refseq_acc:       refseq accession
                - sequence:        mrna sequence


        """
        stmt = 'SELECT idx,name,sequence FROM refmrna WHERE gid=(?) and name=(?)'
        results = self.execute(stmt, (gid, refseq_acc,)).fetchall()
        if len(results) > 0:
            mrna_info = mrna(results[0][0],results[0][1],results[0][2])
        else:
            mrna_info = None
        return mrna_info


if __name__ == '__main__':
    print Refmrna().get_version()
    print Refmrna().retrieve('NM_001005277', 28117)
    print Refmrna().retrieve('NM_002714', 16486)
