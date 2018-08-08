"""
.. module:: regulomedb
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing regulome database information

.. moduleauthor:: Kunal Kundu (kunal.kundu@gmail.com)

Class for accessing regulome database information from its distributed
tsv files.
The tsv files are downloaded from http://regulome.stanford.edu/downloads
"""

from collections import namedtuple
from gcn.lib.io import db
from gcn.etc import dbconfig
import os


class Regulomedb(db.DB):
    """Class to retrieve Regulomedb data"""

    def __init__(self):
        """Class initialization

        Argument

            name (string): Name of the database or filename for the database
        """
        name = 'REGULOME'
        super(Regulomedb, self).__init__()
        if name in dbconfig.DBCONFIG:
            self.load(name=name)
        elif os.path.exists(name):
            self.load(db=name)
        else:
            raise ValueError('No such database %s' % name)

    def retrieve(self, snpid):
        rec = None
        stmt = 'SELECT idx, score from regulome where snpid = (?) limit 1'
        results = self.execute(stmt, (snpid, )).fetchall()
        if results:
            rec = ['RegulomeScore', str(results[0][0]), str(results[0][1])]
        return rec

    def get_header(self):
        """ Intializing vcf formatted header information for Regulomedb
        data"""

        Regulome_h = namedtuple('Regulome_h', 'id count type desc')

        header = [Regulome_h('RegulomeScore', '.', 'String', "Regulome Score "\
                    "as described here - http://regulome.stanford.edu/help")]
        return header


if __name__ == '__main__':
    regdb = Regulomedb()
    print regdb.get_version()
    print regdb.retrieve('rs2887286')
    print regdb.retrieve('rs116904365')
    print regdb.retrieve('rs74479209')
