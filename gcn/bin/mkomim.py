#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: mkomim
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a omim phenotype sqlite database.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a omim phenotype sqlite database. To use this script

usage: mkomim.py [-h] [-i INNAME1,INNAME2] [-o OUTNAME] [--force]

Create a OMIM sqlite database

optional arguments:
  -h, --help            show this help message and exit
  -i INNAME1,INNAME2, --input INNAME1,INNAME2
                        sprot,trembl files in Uniprot format
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
  --force               Create new database even if the existing one is newer
                        than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has two tables - *mimdis*, *release_info*
The *mimdis* table contains the following 4 columns

    - unipacc:       Uniprot Accession
    - gene:      Gene
    - mimphe_id: OMIM Phenotype Id
    - phenotype:    OMIM Phenotype

The *release_info* table contains the following 3 columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the mimdis table of sqlite DB
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import anyopen, db
import argparse
import sys
import os
import time
from datetime import datetime as dt

db_fields = ['idx','unipacc','gene','mimphe_id','phenotype']

SCHEMA = ["""create table mimdis
               ( %s integer PRIMARY KEY AUTOINCREMENT,
                 %s text not null,
                 %s text not null,
                 %s integer not null,
                 %s text not null)
            """%db_fields,
            """create index index1 on mimdis (gene)""",
            """create index index2 on mimdis (mimphe_id)"""
          ]


class OmimDB(db.DB):
    """Class to insert OMIM data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Swissprot file. Note - The trEMBL file does not
                            contain OMIM annotation and even if it exists
                            the information is partial like only mimid is
                            specified but does not say anything is the id
                            is gene omim id or phenotype omim id.
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database even if it
                          is newer than any of  `inname`
        """
        self.logger = fileconfig.getlogger()
        if os.path.exists(outname):
            newer = fileutils.file_newer(inname, outname)
        else:
            newer = True

        if not newer and not force:
            self.logger.info('Not Updating. OMIM database already uptodate')
        else:
            self.inname = inname
            self.outname = outname
            t1 = time.time()
            self._makedb()
            time_taken = (time.time() - t1) / 60
            self.logger.info("Time Taken for creating %s is %f min"
                             % (self.outname, time_taken))

        return

    def getphenotype(self, mimid):
        """For the given OMIM Phenotype ID this methods retrieves
        the OMIM Phenotype using the EUtils service provides by NCBI"""
        phenotype = None
        URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi" + \
                    "?db=omim&id=%s" % (mimid)
        stream = anyopen.openfile(URL)
        for line in stream:
            line = line.strip()
            if line:
                if 'Item Name="Title"' in line:
                    phenotype = line.split('>')[1].split('<')[0].strip()
                    if phenotype:
                        break
                elif 'Item Name="AltTitles"' in line:
                    phenotype = line.split('>')[1].split('<')[0].strip()
                    break
        return phenotype

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating OMIM database ...')
        self.logger.info('Input files: %s' % (self.inname))

        if not os.path.exists(self.inname):
            self.logger.error('%s: No such file' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)

        for s in SCHEMA:
            self.createtable(s, True)

        self.curs = self.conn.cursor()
        query = "insert into mimdis (%s) values (%s)" % (','.join(db_fields[1:]), ','.join('?' * 4))

        n = 0
        entry_cnt = 0
        entries = []
        for acnum, gene, mimids in self._iterfile():
            if len(mimids) > 0:
                for mimid in mimids:
                    phenotype = self.getphenotype(mimid)
                    if phenotype:
                        entries.append((acnum, gene, mimid, phenotype))
                        n += 1
                        if n == 100:
                            self.curs.executemany(query, entries)
                            n = 0
                            entry_cnt += 100
                            entries = []
        if entries:
            self.curs.executemany(query, entries)
            entry_cnt += len(entries)
            entries = []

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(self.inname)
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('Omim', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... OMIM database created')

    def _iterfile(self):
        """Internal method. Do not use"""
        mimids = set()
        stream = anyopen.openfile(self.inname)
        acnum = None
        gene = None
        for line in stream:
            if line[:2] == '//':
                yield acnum, gene, list(mimids)
                mimids = set()
            elif line[:2] == 'ID':
                acnum = None
            elif line[:2] == 'AC':
                if acnum is not None:
                    continue
                acnum = line[2:].split(';')[0].strip()
            elif line[:10] == 'GN   Name=':
                gene = line.split(';')[0].split('=')[1].strip()
            elif line[:9] == 'DR   MIM;':
                if 'phenotype' in line:
                    line = line.strip()
                    id = line.split(';')[1].strip()
                    mimids.add(int(id))


def main():
    """Main script to create the OMIM database"""

    infile = fileconfig.FILECONFIG['OMIM']
    outfil = dbconfig.DBCONFIG['OMIM']['name']
    desc = 'Create a OMIM sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Uniprot files')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    import time
    t1 = time.time()
    OmimDB().createdb(args.inname, args.outname, args.force)
    print 'TIME TAKEN : ', (time.time() - t1) / 60, 'min'
    sys.exit(0)


if __name__ == "__main__":
    main()
