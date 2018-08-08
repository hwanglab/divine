#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: mkrefgene
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a refgene sqlite database.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a refgene sqlite database. To use this script

usage: mkrefgene.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a REFGENE sqlite database

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


The created database has four tables, called *refgene*, *exons*, *defect*,
*release_info*.
The *refgene* table contains the following columns which are all of type text
unless otherwise indicated.

    - name :        Gene name. Usually starts with NM_ or NR_
    - chrom:        Chromosome on which the gene is present (e.g chr1, chr2)
    - strand:       '+' or '-' depending on which strand the gene is
    - txStart:      Starting position of the transcript (integer)
    - txEnd:        Ending position of the transcript (integer)
    - cdsStart:     Coding region start position (integer)
    - cdsEnd:       Coding region end position (integer)
    - exonCount:    Number of exons in the transcript (integer)
    - score:        An integer of unknown meaning (integer)
    - name2:        Alternative name for the gene
    - cdsStartStat: An enum from one of ('none','unk','incmpl','cmpl')
    - cdsEndStat:   An enum from one of ('none','unk','incmpl','cmpl')
    - error:    Error flag. Values are 0 or 1. This flag denotes 0 if the CDS
                derived from gene definition is multiple of 3 else 1

The *exons* table contains the following columns

    - name:         Gene name. Usually starts with NM_ or NR_
    - number:       The numerical id of the exon, first exon is 1, second is
                    2 and so on (integer)
    - chrom:        Chromosome on which the exon is present (e.g chr1, chr2)
    - start:        Starting position of the exon (integer)
    - end:          Ending position of the exon (integer)
    - frame:        The frame of this exon, 0, 1, 2 or -1 (integer)

The *defect* table contains the following columns

    - gid:    Unique identifier for the entry in refgene table
    - name:   Gene name. Usually starts with NM_ or NR_
    - remark:    Type of defect found in the transcript of the gene.
                 For Eg. CDS not multiple of 3 etc

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the refgene table of sqlite DB

The input file for the program is available from the UCSC [1] site

.. [1] ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import anyopen, db
from gcn.lib.varann.vartype.varant import annotateRegion
from collections import namedtuple
import csv
import argparse
import sys
import os
from datetime import datetime as dt
import time

fields_gene = ['gid','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd',\
               'exonCount','score','name2','cdsStartStat','cdsEndStat','error']

fields_exon = ['idx','gid','name','number','start','end','frame']

fields_defect=['idx','gid','name','remark']

SCHEMA = ["""create table refgene
               ( gid integer PRIMARY KEY,
                 name text not null,
                 chrom text not null,
                 strand text not null,
                 txStart integer not null,
                 txEnd integer not null,
                 cdsStart integer not null,
                 cdsEnd integer not null,
                 exonCount integer not null,
                 score integer not null,
                 name2 text not null,
                 cdsStartStat text not null,
                 cdsEndStat text not null,
                 error integer not null)
            """,
            """ create table exons
                (idx integer PRIMARY KEY AUTOINCREMENT,
                 gid integer not null,
                 FOREIGN KEY (gid) REFERENCES refgene(gid),
                 name text not null,
                 number integer not null,
                 start integer not null,
                 end integer not null,
                 frame integer not null)
            """,
            """ create table defect
                (idx integer PRIMARY KEY AUTOINCREMENT,
                 FOREIGN KEY (gid) REFERENCES refgene(gid),
                 gid integer not null,
                 name text not null,
                 remark text not null)
            """,
            """create index rgindex1 on refgene (name)""",
            """create index rgindex2 on refgene (chrom)""",
            """create index rgindex3 on refgene (name2)""",
            """create index rgindex4 on refgene (txStart)""",
            """create index rgindex5 on refgene (txEnd)""",
            """create index rgindex6 on refgene (gid)""",
            """create index exonindex1 on exons (name)""",
            """create index exonindex2 on exons (gid)""",
            """create index defindex1 on defect (gid) """,
            """create index defindex2 on defect (name)""",
            ]


class RefgeneDB(db.DB):
    """Class to insert REFGENE data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Name of the csv file containing refgene data
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
            self.logger.info('Not Updating. REFGENE database already \
                            uptodate')
        else:
            self.inname = inname
            self.outname = outname
            t1 = time.time()
            self._makedb()
            cRef = annotateRegion.RefGeneUcscTB()
            cRef.create_bed()
            time_taken = (time.time() - t1) / 60
            self.logger.info("Time Taken for creating %s is %f min"
                             % (self.outname, time_taken))
        return

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating REFGENE database ...')
        self.logger.info('Input file: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such file' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)

        for s in SCHEMA:
            self.createtable(s, True)

        self.curs = self.conn.cursor()

        self._gid = 0
        for entry in self._iterfile():
            err_flag = self._insertdefect(entry)
            self._insertgene(entry, errorflag=err_flag)
            self._insertexons(entry)
            self._gid += 1

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(self.inname)
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('refgene', fd, version, self._gid)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... REFGENE database created')

    def _iterfile(self):
        """Internal method. Do not use"""

        fields = """bin name chrom strand txStart txEnd cdsStart cdsEnd
                    exonCount exonStarts exonEnds score name2 cdsStartStat
                    cdsEndStat exonFrames""".split()

        refgene = namedtuple('Refgene', fields)

        with anyopen.openfile(self.inname) as stream:
            for rec in csv.reader(stream, delimiter='\t'):
                if rec[0].isdigit():
                    yield refgene._make(rec)

    def get_cds_exonnum(self, frames):
        """Internal method. Do not use"""
        str_cds_exon_number = 0
        end_cds_exon_number = 0
        for e in frames:
            if e != -1:
                break
            else:
                str_cds_exon_number += 1
        for e in frames:
            if e == -1:
                continue
            else:
                end_cds_exon_number += 1
        end_cds_exon_number = str_cds_exon_number + end_cds_exon_number - 1
        return str_cds_exon_number, end_cds_exon_number

    def _insertdefect(self, entry, gid=None):
        """Internal method. Do not use"""
        error_flag = 0
        if gid is None:
            gid = self._gid
        if not entry.name.startswith('NR_'):
            str_cds = int(entry.cdsStart)
            end_cds = int(entry.cdsEnd)
            exon_str = [int(e) for e in entry.exonStarts.split(',')[:-1]]
            exon_end = [int(e) for e in entry.exonEnds.split(',')[:-1]]
            exon_lengths = [x - y for (x, y) in zip(exon_end, exon_str)]
            frames = [int(e) for e in entry.exonFrames.split(',')[:-1]]
            cds_len = 0
            if entry.strand == '+':
                str_cds_exon_number, end_cds_exon_number = \
                                                self.get_cds_exonnum(frames)
                if str_cds_exon_number == end_cds_exon_number:
                    cds_len = end_cds - str_cds
                else:
                    for idx in range(str_cds_exon_number,
                                     end_cds_exon_number + 1):
                        if idx == str_cds_exon_number:
                            cds_len += exon_end[idx] - str_cds
                        elif idx == end_cds_exon_number:
                            cds_len += end_cds - exon_str[idx]
                        else:
                            cds_len += exon_end[idx] - exon_str[idx]
            elif entry.strand == '-':
                exon_str.reverse()
                exon_end.reverse()
                exon_lengths.reverse()
                frames.reverse()
                str_cds_exon_number, end_cds_exon_number = \
                                                self.get_cds_exonnum(frames)
                if str_cds_exon_number == end_cds_exon_number:
                    cds_len = end_cds - str_cds
                else:
                    for idx in range(str_cds_exon_number,
                                                end_cds_exon_number + 1):
                        if idx == str_cds_exon_number:
                            cds_len += end_cds - exon_str[idx]
                        elif idx == end_cds_exon_number:
                            cds_len += exon_end[idx] - str_cds
                        else:
                            cds_len += exon_end[idx] - exon_str[idx]
            if divmod(cds_len, 3)[1] != 0:
                error_flag = 1

                definsert = 'insert into defect (%s) values (?,?,?)' % (','.join(fields_defect[1:]))

                self.curs.execute(definsert, (gid, entry.name,
                                        'CDS not multiple of 3'))
        return error_flag

    def _insertgene(self, entry, errorflag, gid=None):
        """Internal method. Do not use"""
        if gid is None:
            gid = self._gid

        rginsert = "insert into refgene (%s) values (%s)" % (','.join(fields_gene), ','.join('?' * 14))

        self.curs.execute(rginsert, (gid, entry.name, entry.chrom,
                                     entry.strand,
                                    int(entry.txStart), int(entry.txEnd),
                                    int(entry.cdsStart), int(entry.cdsEnd),
                                    int(entry.exonCount), int(entry.score),
                                    entry.name2, entry.cdsStartStat,
                                    entry.cdsEndStat, errorflag))

    def _insertexons(self, entry, gid=None):
        """Internal method. Do not use"""

        if gid is None:
            gid = self._gid

        exinsert = 'insert into exons (%s) values (?,?,?,?,?,?)'%(','.join(fields_exon[1:]))
        exonstarts = [int(e) for e in entry.exonStarts.split(',')[:-1]]
        exonends = [int(e) for e in entry.exonEnds.split(',')[:-1]]
        exonframes = [int(e) for e in entry.exonFrames.split(',')[:-1]]
        n = len(exonstarts)
        if len(exonends) != n or len(exonframes) != n:
            self.logger.error('Incorrect exon information in %s' % entry.name)
            sys.exit(1)
        idx = 1
        for es, ee, ef in zip(exonstarts, exonends, exonframes):
            self.curs.execute(exinsert, (gid, entry.name, idx,
                                         es, ee, ef))
            idx += 1


def main():
    """Main script to create the REFGENE database"""

    infile = fileconfig.FILECONFIG['REFGENE']
    outfil = dbconfig.DBCONFIG['REFGENE']['name']
    desc = 'Create a REFGENE sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Refgene File')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    RefgeneDB().createdb(args.inname, args.outname, args.force)
    sys.exit(0)


if __name__ == "__main__":
    main()
