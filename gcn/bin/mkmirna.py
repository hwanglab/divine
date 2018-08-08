#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: mkmirna
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a mirna sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a mirna sqlite database. To use this script

usage: mkmirna.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a MIRNA sqlite database

optional arguments:
  -h, --help            show this help message and exit
  -i INNAME, --input INNAME
                        Fasta file of sequences
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
  --force               Create new database even if the existing one is newer
                        than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has two tables - *mirna*, *release_info*.
The *mirna* table contains the following columns which are all of type text
unless otherwise indicated.

    - mirbase_acc:             mirbase accession number (Eg. MIMAT0000062)
    - mirna_name:              mirna name (Eg. hsa-let-7a)
    - ext_transcript_id :      RefSeq Accession (Eg. NM_005504)
    - gene_alignment :         Alignment of the gene sequence (\
                                Eg. guagGUAAAGGAAACUACCUCa)
                               where the caps nucleotide are point of \
                                contacts for binding.
    - genome_coordinates :     Alignment genome coordinates (\
                               Eg. hg19:3:119536478-119536497:+)
    - energy :                 (FLoat)Binding energy
    - mirsvr_score :           (Float)Mirsvr (Mirand Algorithm) score

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the mirna table of sqlite DB

The input file (Good mirSVR score,Conserved miRNA) for the program is \
available from the mirna.org [1] site

.. [1] http://www.microrna.org/microrna/getDownloads.do
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import anyopen, db
from collections import namedtuple
import csv
import argparse
import sys
import os
import time
from datetime import datetime as dt

db_fields = ['idx','mirbase_acc','mirna','refseq_acc','gene_align_seq','coord_info','energy','mirsvr_score']
SCHEMA = ["""create table mirna
            (%s integer PRIMARY KEY AUTOINCREMENT,
             %s text not null,
             %s text not null,
             %s text not null,
             %s text not null,
             %s text not null,
             %s FLOAT not null,
             %s FLOAT not null)
            """%db_fields,
            """create index mrindex1 on mirna (refseq_acc)""",
            ]


class MiRnaDB(db.DB):
    """Class to insert REFGENE data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Name of the csv file containing mirna data
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database even if it
                          is newer than `inname`
        """
        if os.path.exists(outname):
            newer = fileutils.file_newer(inname, outname)
        else:
            newer = True

        self.logger = fileconfig.getlogger()

        if not newer and not force:
            self.logger.info('Not Updating. MIRNA database already uptodate')
        else:
            self.inname = inname
            self.outname = outname
            t1 = time.time()
            self._makedb()
            time_taken = (time.time() - t1) / 60
            self.logger.info("Time Taken for creating %s is %f min" \
                             % (self.outname, time_taken))
        return

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating MIRNA database ...')
        self.logger.info('Input file: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such file' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)

        for s in SCHEMA:
            self.createtable(s, True)

        self.curs = self.conn.cursor()
        entry_cnt = 0
        for entry in self._iterfile():
            self._insertmirna(entry)
            entry_cnt += 1

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(self.inname)\
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('miRNAdb', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... MIRNA database created')

    def _iterfile(self):
        """Internal method. Do not use"""

        fields = """mir_acc mirna gene_id gene_sym trans_id ext_trans_id \
                mirna_alig alig gene_alig mirna_start mirna_end gene_start \
                gene_end genome_coord conserv align_score
                     seed_cat energy mirsvr_score """.split()

        mirna = namedtuple('mirna', fields)

        with anyopen.openfile(self.inname) as stream:
            for rec in csv.reader(stream, delimiter='\t'):
                if '#' not in rec[0]:
                    yield mirna._make(rec)

    def _insertmirna(self, entry):
        """Internal method. Do not use"""
        rginsert = "insert into mirna (%s) values (%s)" % (','.join(db_fields[1:]), ','.join('?' * 7))

        self.curs.execute(rginsert, (entry.mir_acc, entry.mirna, \
                                     entry.ext_trans_id, entry.gene_alig,
                                     entry.genome_coord.strip('[]'), \
                                     float(entry.energy),\
                                     float(entry.mirsvr_score)))


def main():
    """Main script to create the MIRNA database"""

    infile = fileconfig.FILECONFIG['MIRNA']
    outfil = dbconfig.DBCONFIG['MIRNA']['name']
    desc = 'Create a MIRNA sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='miRNA DB raw text file')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true', default=True,
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    MiRnaDB().createdb(args.inname, args.outname, args.force)
    sys.exit(0)


if __name__ == "__main__":
    main()
