#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
# !/usr/bin/env python
"""
.. module:: mkregulomedb
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a RegulomeDB sqlite database.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a regulomedb sqlite database. To use this script

usage: mkclnsnpdb.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a RegulomeDB sqlite database

optional arguments:
  -h, --help    show this help message and exit
  -i INNAME, --input INNAME
                        Directory name that contains the Regulomedb files.
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
  --force               Create new database even if the existing one is newer
                        than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has two tables - *regulome*, *release_info*.
The *regulome* table contains the following columns which are all of type text
unless otherwise indicated.

    - chrom:    Chromosome
    - pos:      Genomic position
    - snpid:    dbSNP Accession
    - exp_info: Experiment Information
    - score:    Score based on experiment information.

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the regulome table of sqlite DB

The input directory for the program contains five .txt file -
downloaded from  - http://www.regulomedb.org/downloads
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import anyopen, db
import csv
import argparse
import sys
import os
import time
from datetime import datetime as dt

db_fields = ['idx','chrom','pos','snpid','exp_info','score']

SCHEMA = ["""create table regulome
               ( idx integer PRIMARY KEY AUTOINCREMENT,
                 chrom text not null,
                 pos integer not null,
                 snpid text,
                 exp_info text not null,
                 score text not null)
            """,
            """create index reg1 on regulome (snpid)""",
            """create index reg2 on regulome (chrom, pos, ref, alt)""",
            ]


class RegulomeDB(db.DB):
    """Class to insert RegulomeDB data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Name of the directory containing
                            regulomedb data
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database
                            even if it is newer than `inname`
        """
        self.inname = inname
        if os.path.exists(outname):
            for fn in os.listdir(inname):
                filepath = os.path.join(inname, fn)
                newer = fileutils.file_newer(filepath, outname)
                if newer is True:
                    break
        else:
            newer = True

        self.logger = fileconfig.getlogger()

        if not newer and not force:
            self.logger.info('Not Updating. RegulomeDB' +
                            'database already uptodate')
        else:
            self.inname = inname
            self.outname = outname
            t1 = time.time()
            self._makedb()
            time_taken = (time.time() - t1) / 60
            self.logger.info("Time Taken for creating %s is %f min"
                             % (self.outname, time_taken))
        return

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating Regulome database ...')
        self.logger.info('Input directory: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such directory' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)

        for s in SCHEMA:
            self.createtable(s, True)

        self.curs = self.conn.cursor()

        temp = []
        entry_cnt = 0
        for entry in self._iterfile():
            temp.append(entry)
            if len(temp) == 10000:
                self._insert(temp)
                entry_cnt += 10000
                temp = []
        if temp:
            self._insert(temp)
            entry_cnt += len(temp)
            temp = []

        for fn in os.listdir(self.inname):
            infile = os.path.join(self.inname, fn)

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(infile)
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('Regulome', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... Regulome database created')

    def _iterfile(self):
        """Internal method. Do not use"""
        for filename in os.listdir(self.inname):
            self.logger.info('Processing %s..' % filename)
            infile = os.path.join(self.inname, filename)
            with anyopen.openfile(infile) as stream:
                for rec in csv.reader(stream, delimiter='\t'):
                    yield tuple(rec)

    def _insert(self, entrylist, sid=None):
        """Internal method. Do not use"""
        rginsert = "insert into regulome (%s) values (%s)" % (','.join(db_fields[1:]), ','.join('?'*5))
        self.curs.executemany(rginsert, entrylist)


def main():
    """Main script to create the Regulomedb sqlite database"""

    infile = fileconfig.FILECONFIG['REGULOME']
    outfil = dbconfig.DBCONFIG['REGULOME']['name']
    desc = 'Create a Regulomedb sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Directory path containing regulomedb data \
                        files')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    RegulomeDB().createdb(infile, outfil, args.force)
    sys.exit(0)


if __name__ == "__main__":
    main()
