"""
.. module:: utrdb
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing utrdb information from an sqlite database

.. moduleauthor:: Kunal Kundu (kunal.kundu@gmail.com)

Class for accessing utrdb information from an sqlite database
"""
from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
from gcn.lib.databases.refgene import get_ucsc_chrom
from bisect import bisect
from gcn.etc import dbconfig
import os

CHROMOSOMES, ICHROMOSOMES = get_ucsc_chrom()

class UTRdb(db.DB):
    """Class to retrieve UTR information from UTRDB database"""

    def __init__(self):
        """Class initialization

        Argument

            name (string): Name of the database or filename for the database
        """
        name = 'UTRDB'
        super(UTRdb, self).__init__()
        if name in dbconfig.DBCONFIG:
            self.load(name=name)
        elif os.path.exists(name):
            self.load(db=name)
        else:
            raise ValueError('No such database %s' % name)

    def load_utrdb(self):
        self.anntlist = [[] for c in CHROMOSOMES]
        stmt = 'SELECT start_pos, end_pos, refseq_acc, std_name,repeat_type,\
        repeat_family FROM utrdb WHERE chrom=(?)'
        import time
        t1 = time.time()
        for chrno in CHROMOSOMES:
            results = self.execute(stmt, (chrno, )).fetchall()
            if results:
                for row in results:
                    #print row
                    std_name = row[3]
                    rep_type = row[4]
                    rep_fam = row[5]
                    if std_name or rep_type or rep_fam:
                        cid = ICHROMOSOMES[chrno]
                        if row not in self.anntlist[cid]:
                            self.anntlist[cid].append(row)
        for key in self.anntlist:
            key.sort()

    def retrieve(self, chrno, position, refseq_acc):
        utrannt = set([])
        temp = set()
        entry_flag = False
        if not chrno.startswith('chr'):
            chrno = 'chr' + chrno
        pos = int(position)
        if chrno in ICHROMOSOMES:
            chrom = ICHROMOSOMES[chrno]
            idx = bisect(self.anntlist[chrom], (position + 1, position + 2))
            if idx != 0:
                temp = set()
                for j in range(idx, -1, -1):
                    i1, i2, acc, std_name, rep_type, rep_fam = \
                                            self.anntlist[chrom][j - 1]
                    if i1 <= pos <= i2 and acc == refseq_acc:
                        if std_name:
                            temp.add(std_name.split('(')[0].strip().\
                                                replace(' ', ''))
                        if rep_type:
                            temp.add('Repeat-' + rep_type.replace('(', '').\
                                     replace(')', '').strip().\
                                     replace(' ', '').replace('_', '') + '-' \
                                     + rep_fam.\
                                     replace(' ', '').strip().replace('_', ''))
                        entry_flag = True
                    else:
                        if entry_flag == True:
                            break
        if len(temp) > 0:
            utrannt = temp
        return utrannt

#    def retrieve(self, chromosome, position, refseq_acc):
#        """Retrieve details of all utr signals associated with a given
#        chromosomal location.
#        Args:
#            chromosome (str):   Chromosome name in the form chrN where
#                                N=1-23,X,Y
#            position (integer): Location on the chromosome, using 0-based
#                                indexing
#        Returns:
#            List of utrsignals associated with the position. Each gene is a
#            named tuple with the following fields

#                - site:         utr signal motif name
#                - desc:         description of the signal, if known
#                - std_name:     standard name of the signal
#                - rep_type:     repeat type, if signal is a repeat region
#                - rep_fam:      repeat family, if signal is a repeat region
#       """
#        info = None
#        self.load(DBCONFIG['UTRDB']['name'])
#        c = self.conn.cursor()
#        if 'chr' not in chromosome:
#            chromosome = 'chr'+chromosome
#        stmt = "SELECT std_name,repeat_type,repeat_family FROM utrdb WHERE \
#                chrom=(?) AND (?)>=start_pos AND (?)<=end_pos \
#                AND refseq_acc=(?)"
#        results = c.execute(stmt, (chromosome, position, position, \
#                                   refseq_acc)).fetchall()
#        utr_signals = set()
#        if results:
#            for row in results:
#                if row[0]:
#                    utr_signals.add(row[0].split('(')[0].strip()\
#                                    .replace(' ', '_'))
#                elif row[1]:
#                    utr_signals.add('Repeat-' + row[1].replace('(', '').\
#                                replace(')', '').strip().replace(' ', '_') +\
#                                     '-' + row[2].replace(' ', '_').strip())
#            if len(utr_signals) > 0:
#                info = "|".join(list(utr_signals))
#
#        return info


if __name__ == '__main__':
    udb = UTRdb()
    print udb.get_version()
    udb.load_utrdb()
    print udb.retrieve('chr12', 125627229, 'NM_023928')
    print udb.retrieve('chr2', 189877298, 'NM_000090')
    print udb.retrieve('chrX', 107397617, 'NM_052936')
    print udb.retrieve('1', 55505447, 'NM_174936')
    print udb.retrieve('10', 70968818, 'NM_003171')
