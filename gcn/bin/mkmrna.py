#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: mkmrna
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a refmrna sqlite database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a refmrna sqlite database. To use this script

usage: refmkmrna.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a REFMRNA sqlite database

optional arguments:
  -h, --help            show this help message and exit
  -i INNAME, --input INNAME
                        Directory path for human genome
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
  --force               Create new database even if the existing one is newer
                        than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has two tables - *refmrna*, *release_info*.
The *refmrna* table contains the following columns which are all of type text
unless otherwise indicated.

    - gid :    Unique identifier for each transcript
    - name :    Gene name. Usually starts with NM_ or NR_
    - seq_len:    Sequence length
    - sequence :    Neucleotide sequence

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the refmrna table of sqlite DB

The input directory of the latest human genome (hg19) for the program
is available from the UCSC [1] site
.. [1]
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import db, anyopen
from gcn.lib.databases.refgene import Refgene
from gcn.lib.io.entrystream import FastaStream
from gcn.lib.io import fasta
import argparse
import sys
import os
import time
from datetime import datetime as dt

db_fields = ['idx','gid','name','seq_len','sequence']

SCHEMA = ["""create table refmrna
               ( %s integer PRIMARY KEY AUTOINCREMENT,
                 %s integer not null,
                 %s text not null,
                 %s integer not null,
                 %s text not null)
            """%db_fields,
            """create index rmindex1 on refmrna (gid)""",
            """create index rmindex2 on refmrna (name)""",
            ]


class RefmrnaDB(db.DB):
    """Class to insert REFMRNA sequence data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Path to human genome reference file
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database even if it
                          is newer than `inname`
        """
        if os.path.exists(outname):
            if fileutils.file_newer(inname, outname) or \
                    fileutils.file_newer(dbconfig.DBCONFIG['REFGENE']['name'],\
                                          outname):
                newer = True
            else:
                newer = False
        else:
            newer = True

        self.logger = fileconfig.getlogger()

        if not newer and not force:
            self.logger.info('Not Updating. REFMRNA database already uptodate')
        else:
            self.inname = inname
            self.outname = outname
            t1 = time.time()
            self._makedb()
            time_taken = (time.time() - t1) / 60
            self.logger.info("Time Taken for creating %s is %f min" \
                             % (self.outname, time_taken))
        return

    def set_input_dir(self, inname):
        """This method sets the input directory path. This method is added
        explicitly to access the other methods like load_hg(), getseq() which
        requires the inname variable to be in an assigned state"""
        self.inname = inname

    def load_hg(self, chrom, inname=None):
        """Internal method. This method loads the human genome in a dictionary
        where the key is chromosome name and value is sequence"""
        self.chrom_seq = None
        if not inname:
            inname = self.inname
        if chrom[:3] == 'chr':
            chrom = chrom[3:]
        stream = anyopen.openfile(inname)
        es = FastaStream(stream)
        for seqi in es:
            header = seqi[0].rstrip().split()[0]
            if 'chr' in header:
              header = header[4:]
            else:
              header = header[1:]
            if header == chrom:
                self.chrom_seq = ''.join([s[:-1] for s in seqi[1:]]).lower()
                break
        return self.chrom_seq

    def getseq(self, chrom, start_positions, end_positions, strand):
        """Extracts the sequence for the given gene coordinate
        Args:
            chrom (str):    Chromosome name in the form chrN where
                                N=1-23,X,Y
            start_positions (list):    List of start exon positions. 0-based
                                indexing. For Eg. [34610, 35276, 35720]
            end_positions (list):    List of end exon positions. 0-based
                                indexing. For Eg. [35174,35481,36081]
            strand (str):    Strand is either + or -

        Returns:
            gene_seq:    mRNA Sequence of the gene
        """
        gene_seq = None
        try:
            chrom_seq = self.chrom_seq
            for sidx, eidx in zip(start_positions, end_positions):
                try:
                    gene_seq += chrom_seq[sidx:eidx]
                except:
                    gene_seq = chrom_seq[sidx:eidx]
            if strand == '-':
                gene_seq = fasta.revcomplement(gene_seq)
        except:
            self.logger.error('Chromosome sequence not found for %s' % (chrom))
        return gene_seq

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating REFMRNA database ...')
        self.logger.info('Input file: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such file' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)

        for s in SCHEMA:
            self.createtable(s, True)

        self.curs = self.conn.cursor()

        refg = Refgene()

        human_chrom = ['chr' + str(e) for e in range(1, 23)] + \
                                ['chrX', 'chrY', 'chrM']
        records = []
        entry_cnt = 0
        for chrom in human_chrom:
            self.logger.info('Processing chrom-%s' % (chrom))
            self.load_hg(chrom)
            for entry in refg.iterate_db(chrom):
                chrom = entry.chrom
                if '_' not in chrom:
                    strand = entry.strand
                    start_positions = [int(pos) for pos in \
                                       entry.str_exon_pos.split(',')]
                    end_positions = [int(pos) for pos in \
                                     entry.end_exon_pos.split(',')]
                    seq = self.getseq(chrom, start_positions, \
                                       end_positions, strand)
                    if seq:
                        records.append((entry.gid, entry.refseq_acc, \
                                        len(seq), seq.lower()))
                        if len(records) == 1000:
                            #print 'Inserted %d records [seqLen:%d]' % (1000,len(seq))
                            self._insert(records)
                            entry_cnt += 1000
                            records = []
        if records:
            self._insert(records)
            entry_cnt += len(records)
            records = []

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(\
                                dbconfig.DBCONFIG['REFGENE']['name'])\
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('refmrna', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... REFMRNA database created')

    def _insert(self, records):
        """Internal method. Do not use"""
        rginsert = "insert into refmrna (%s) values (%s)" % (','.join(db_fields[1:]), ','.join('?'* 4))

        self.curs.executemany(rginsert, records)


def main():
    """Main script to create the REFMRNA database"""
    infile = fileconfig.FILECONFIG['REFGENOME']
    outfil = dbconfig.DBCONFIG['REFMRNA']['name']
    force = True
    desc = 'Create a REFMRNA sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Directory path for human genome')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true', default=force,
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    t1 = time.time()
    RefmrnaDB().createdb(args.inname, args.outname, args.force)
    print 'Time Taken to create the DB:', float(time.time() - t1) / 60, 'min'
    sys.exit(0)


if __name__ == "__main__":
    main()
