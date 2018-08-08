#!/usr/bin/env python
"""
Script to build a utrdb sqlite database. To use this script

usage: mkutrdb.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a UTRDB sqlite database

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


The created database has two tables - *utrdb*, *release_info* .
The *utrdb* table contains the following columns which are all of type text
unless otherwise indicated.

    - chrom:          Chromosome on which the gene is present (e.g chr1, chr2)
    - strand:         '+' or '-' depending on which strand the gene is
    - start_pos:      Starting position of the site (integer)
    - end_pos:        Ending position of the site (integer)
    - refseq_acc:     RefSeq Accession
    - utr_typr:       3'UTR or 5'UTR
    - site:           Site name
    - evid:           Evidence for the site detection
    - desc:           Description of the site
    - std_name:       Standard name of the site
    - repeat_type:    Repeat type
    - repeat_family:  Repeat family

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the utrdb table of sqlite DB

The input file for the program is available from the Institute for Biomedical \
Technologies - Univ of Bari [1] site

.. [1] ftp://ftp.ba.itb.cnr.it/pub/Databases/UTR/data/3UTRaspic.Hum.dat.gz
.. [2] ftp://ftp.ba.itb.cnr.it/pub/Databases/UTR/data/5UTRaspic.Hum.dat.gz
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import db
from parseutr import UTR_Parser
from hgConvert.converter import hgConvert
import argparse
import sys
import os
from datetime import datetime as dt
import time

db_fields = ['idx','chrom','start_pos','end_pos','strand',\
          'refseq_acc','utr_type','site','evid','desc',\
          'std_name','repeat_type','repeat_family']

SCHEMA = ["""create table utrdb
            (idx integer PRIMARY KEY AUTOINCREMENT,
             chrom text not null,
             start_pos integer not null,
             end_pos integer not null,
             strand text not null,
             refseq_acc text not null,
             utr_type text not null,
             site text not null,
             evid text not null,
             desc text not null,
             std_name text not null,
             repeat_type text not null,
             repeat_family text not null)
            """,
            """create index utrindex1 on utrdb (refseq_acc)""",
            """create index utrindex2 on utrdb (chrom)""",
            """create index utrindex3 on utrdb (start_pos)""",
            """create index utrindex4 on utrdb (end_pos)"""
            ]


class UtrDB(db.DB):
    """Class to insert REFGENE data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Name of the csv file containing utr data
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database even if it
                          is newer than `inname`
        """
        newfiles = []
        self.inname = inname
        if os.path.exists(outname):
            for fn in os.listdir(inname):
                infile = os.path.join(inname, fn)
                newer = fileutils.file_newer(infile, outname)
                if newer is True:
                    break
        else:
            newer = True

        self.logger = fileconfig.getlogger()

        if not newer and not force:
            self.logger.info('Not Updating. UTRdb' +
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
        self.logger.info('Creating UTRDB database ...')
        self.logger.info('Input file: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such file' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)
        self.load(db=self.outname)
        for s in SCHEMA:
            self.createtable(s, True)
        self.curs = self.conn.cursor()
        tmplist = []
        entry_cnt = 0
        for entry in self._iterfile():
            for parsed_dict in entry:
                for row_data in self._build_dbrows(parsed_dict):
                    tmplist.append(row_data)
                    if len(tmplist) > 10000:
                        self._insertutr(tmplist)
                        entry_cnt += 10000
                        tmplist = []
        if tmplist:
            self._insertutr(tmplist)
            entry_cnt += len(tmplist)
            tmplist = []

        for fn in os.listdir(self.inname):
            infile = os.path.join(self.inname, fn)
            break

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(infile)
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('UTRdb', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... UTRDB database created')

    def _iterfile(self):
        """Internal method. Do not use"""
        sections = ['FT']
        for filename in os.listdir(self.inname):
            fp = os.path.join(self.inname, filename)
            self.hg19map = self._hg_convert(fp)
            p = UTR_Parser(fp, sections)
            parsed_dict = p.parse()
            yield parsed_dict

    def _hg_convert(self, filepath):
        """ Internal Method. For the given UTRdb file path
        it scans for all the lines with genome coordinates and
        creates an input coordinate list for the hg convert module.
        Then it calls the hg convert module and builds hg19 coord
        look up dictionary. This hg19map dictionary is returned
        for further use.
        """
        inlist = set([])
        cmd = "grep '/genome=' %s | grep 'FT'" % (filepath)
        stream = os.popen(cmd)
        for line in stream:
            line = line.strip()
            chrom, pos, strand = line.split('=')[1].strip('"').split(':')
            start_pos, end_pos = [int(e) for e in pos.split('-')]
            if start_pos > end_pos:
                tmp = start_pos
                start_pos = end_pos
                end_pos = tmp
            inlist.add((chrom, start_pos, end_pos))

        hgc = hgConvert()
        hg19map = hgc.hg18_to_hg19(list(inlist))
        return hg19map

    def get_convrt_coordinate(self, genome):
        """For the given genome (Format: chrom:startpos-endpos:starnd)
        this method looks up in the hg19map dictionary and returns the
        converted hg19 coordinate."""
        chrom, pos, strand = genome.split(':')
        start_pos, end_pos = [int(e) for e in pos.split('-')]
        if start_pos > end_pos:
            tmp = start_pos
            start_pos = end_pos
            end_pos = tmp
        key = chrom + ':' + str(start_pos) + '-' + str(end_pos)
        hgcoord = list(self.hg19map[key])
        if hgcoord:
            hgcoord.append(strand)
        return hgcoord

    def _build_dbrows(self, parsed_dict):
        """Internal method. Do not use"""
        dict = parsed_dict['FT']

        for m in dict['source']:
            if m.keys()[0] == 'db_xref':
                if m.values()[0].split(':')[0] == 'RefSeq':
                    refseq_acc = m.values()[0].split(':')[1]

        if "3'UTR" in dict.keys():
            utr_type = "3'UTR"
        elif "5'UTR" in dict.keys():
            utr_type = "5'UTR"

        for key, val in dict.iteritems():
            desc = ""
            evidence = ""
            std_name = ""
            repeat_type = ""
            repeat_fam = ""
            genome = ""
            hg19_c = []

            if key == 'ConservedRegion':
                for m in val:
                    if m.keys()[0] == 'genome':
                        genome = m.values()[0]
                        hg19_c = self.get_convrt_coordinate(genome)
                    elif m.keys()[0] == 'description':
                        desc = m.values()[0]
            elif key == 'repeat_region':
                for m in val:
                    if m.keys()[0] == 'evidence':
                        evidence = m.values()[0]
                    elif m.keys()[0] == 'repeat_type':
                        repeat_type = m.values()[0]
                    elif m.keys()[0] == 'repeat_family':
                        repeat_fam = m.values()[0]
                    elif m.keys()[0] == 'genome':
                        genome = m.values()[0]
                        hg19_c = self.get_convrt_coordinate(genome)
            else:
                for m in val:
                    if m.keys()[0] == 'evidence':
                        evidence = m.values()[0]
                    elif m.keys()[0] == 'standard_name':
                        std_name = m.values()[0]
                    elif m.keys()[0] == 'genome':
                        genome = m.values()[0]
                        hg19_c = self.get_convrt_coordinate(genome)

            if key not in ['source', "5'UTR", "3'UTR"] and hg19_c:
                out = (hg19_c[0], int(hg19_c[1]), int(hg19_c[2]), hg19_c[3],
                       refseq_acc, utr_type, key, evidence, desc,
                       std_name, repeat_type, repeat_fam)
                yield out
            elif key not in ['source', "5'UTR", "3'UTR"] and \
                            not hg19_c and len(genome) > 0:
                out = genome + 'could not be converted..\n'
                # prints the hg18 genome coordinates which could not be \
                # converted to hg19 coordinate.. need to log this info
                self.logger.warning(out)

    def _insertutr(self, row_data):
        """Internal method. Do not use"""

        insert = "insert into utrdb (%s) values (%s)" % (','.join(db_fields[1:]), ','.join('?' * 12))

        self.curs.executemany(insert, row_data)


def main():
    """Main script to create the UTRDB database"""

    infile = fileconfig.FILECONFIG['UTRDB']
    outfil = dbconfig.DBCONFIG['UTRDB']['name']
    desc = 'Create a UTRDB sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Directory containing the UTRdb raw files')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    UtrDB().createdb(args.inname, args.outname, args.force)


if __name__ == "__main__":
    t1 = time.time()
    main()
    print 'Time Taken to create UTRdb:', float(time.time() - t1) / 60, 'min'
