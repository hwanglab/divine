#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: mknsfpdb
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a sqite db for the NSFP Database

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a sqite db for the NSFP Database
for hg19 exome coordinates. To use this script

usage: mknsfpdb.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a NSFP sqlite database

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


The created database has two tables - *nsfp*, *release_info*.
The *nsfp* table contains the following columns which are all of type text
unless otherwise indicated.

    - chr:    chromosome
    - pos:    position (integer)
    - ref:    Reference Allele
    - alt:    Observed Allele
    - refaa:    Wild type Amino acid
    - altaa:    Mutant type Amino acid
    - gene:    Gene Name
    - embl_transcriptid:    Ensembl transcript ids (separated by ";")
    - aapos:    Amino acid position as to the protein (Integer)
    - uniprot_acc_polyphen2:    Uniprot accession
    - sift_score:    SIFT Score (Float)
    - sift_pred:    If a score is smaller than 0.05 the corresponding NS is
                    predicted as "D(amaging)"; otherwise it is
                    predicted as "T(olerated)"
    - pp2_hvar_score:    Polyphen2 score based on HumVar, i.e. hvar_prob.
                        The score ranges from 0 to 1, and the corresponding
                        prediction is "probably damaging" if it is in
                        [0.909,1]; "possibly damaging" if it is in
                        [0.447,0.908]; "benign" if it is in [0,0.446].(FLoat)
    - pp2_hvar_pred:    Polyphen2 prediction based on HumVar, "D"
                        ("porobably damaging"),"P" ("possibly damaging")
                        and "B" ("benign").

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the regulome table of sqlite DB

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

fields = """chr pos ref alt refaa altaa gene transid aapos unpacc """\
         """siftscr siftpred hvarscr hvarpred caddraw caddphred """\
         """domains""".split()

db_fields = """chr pos ref alt refaa altaa gene embl_transcriptid aapos uniprot_acc_polyphen2 """\
         """sift_score sift_pred pp2_hvar_score pp2_hvar_pred cadd_raw cadd_phred """\
         """domains""".split()

SCHEMA = ["""create table nsfp
               ( idx integer PRIMARY KEY AUTOINCREMENT,
                 chr text not null,
                 pos integer not null,
                 ref text not null,
                 alt text not null,
                 refaa text not null,
                 altaa text not null,
                 gene text not null,
                 embl_transcriptid text not null,
                 aapos integer not null,
                 uniprot_acc_polyphen2 text,
                 sift_score float,
                 sift_pred text,
                 pp2_hvar_score float,
                 pp2_hvar_pred text,
                 cadd_raw float,
                 cadd_phred float,
                 domains text)
            """,
            """create index index1 on nsfp (chr,pos,ref,alt,refaa,altaa,aapos)""",
            """create index index2 on nsfp (gene, refaa, altaa, aapos)""",
        ]


class NsfpDB(db.DB):
    """Class to insert NSFP data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Name of the directory containing NSFP Database files
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database even if it
                          is newer than `inname`
        """

        self.indir = os.path.dirname(inname)
        if os.path.exists(outname):
            newer = False
            for filename in os.listdir(self.indir):
                if 'chr' in os.path.splitext(filename)[1]:
                    infile = os.path.join(self.indir, filename)
                    newer = fileutils.file_newer(infile, outname)
                    if newer:
                        break
        else:
            newer = True

        self.logger = fileconfig.getlogger()

        if not newer and not force:
            self.logger.info('Not Updating. NSFP database already uptodate')
        else:
            self.outname = outname
            t1 = time.time()
            self._makedb()
            time_taken = (time.time() - t1) / 60
            self.logger.info("Time Taken for creating %s is %f min"
                             % (self.outname, time_taken))
        return

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating NSFP database ...')
        self.logger.info('Input directory: %s' % self.indir)
        if not os.path.exists(self.indir):
            self.logger.error('%s: No such directory' % self.indir)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)

        for s in SCHEMA:
            self.createtable(s, True)

        self.curs = self.conn.cursor()

        query = "insert into nsfp (%s) values (%s)" % (','.join(db_fields), ','.join('?' * 17))

        entry_cnt = 0
        entries = []
        for entry in self._iterfile():
            # print entry
            entries.append((entry.chr, int(entry.pos), entry.ref, entry.alt,
                            entry.refaa, entry.altaa, entry.gene,
                            entry.transid, int(entry.aapos), entry.unpacc,
                            entry.siftscr, entry.siftpred, entry.hvarscr,
                            entry.hvarpred, entry.caddraw, entry.caddphred,
                            entry.domains))
            
            entry_cnt += 1
            if not entry_cnt % 4000000:
                self.curs.executemany(query, entries)
                entries = []
            
        if entries:
            self.curs.executemany(query, entries)
            entries = []

        for fn in os.listdir(self.indir):
            if 'chr' in os.path.splitext(fn)[1]:
                infile = os.path.join(self.indir, fn)
                break

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(infile)
                              ).strftime('%Y-%m-%d')
        version = 'v' + os.path.split(infile)[1].split('_')[0].strip('dbNSFP')
        self.set_version('dbNSFP', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... NSFP database created')

    def _iterfile(self):
        """Internal method. Do not use"""

        nsfp_data = namedtuple('Nsfp', fields)
        for filename in os.listdir(self.indir):
            if 'chr' in os.path.splitext(filename)[1]:
            #if '.chr1' == os.path.splitext(filename)[1]: #debug
                self.logger.info('Processing File : %s' % filename)
                stream = anyopen.openfile(os.path.join(self.indir, filename))
                for rec in csv.reader(stream, delimiter='\t'):
                    if rec[0] == '#chr':
                        h = []
                        for idx, e in enumerate(rec):
                            if idx == 0:
                                h.append('chr')
                            elif idx == 8:
                                h.append('pos')
                            else:
                                e = e.strip()
                                e = e.replace('(', '_')
                                e = e.replace(')', '_')
                                e = e.replace('-', '_')
                                e = e.replace('+', '')
                                e = e.replace('1000', 'K1')
                                h.append(e)
                        rt = namedtuple('rt', h)
                    else:
                        rec_tup = rt._make(tuple(rec))
#                         print 'pos1:%s|hg19-pos1:%s'%(rec_tup.pos_1_based_,rec_tup.pos)
#                         if rec_tup.pos == '949523':
#                           debug = 1
                        # CADD raw and phred scores
                        if rec_tup.CADD_raw != '.':
                            cadd_raw = float(rec_tup.CADD_raw)
                            cadd_phred = float(rec_tup.CADD_phred)
                        else:
                            cadd_raw, cadd_phred = (None, None)

                        # Interpro domain
                        if rec_tup.Interpro_domain != '.':
                            domains = [e.split('(')[0].strip()
                                       for e in
                                       rec_tup.Interpro_domain.split(';') if e]
                            domains = [e.replace(' ', '').replace(',', '_')
                                       for e in domains]
                            domains = ','.join(domains)
                        else:
                            domains = None

                        # ENSEMBL TRANSCRIPT IDS
                        trans_ids = [trans.strip() for trans in
                                    rec_tup.Ensembl_transcriptid.split(';')]

                        # AMINO ACID POSITIONS IN TRANSCRIPTS
                        trans_aapos = [int(pos.strip()) for pos in
                                            rec_tup.aapos.split(';')]
                        if len(trans_ids) > 1 and len(trans_aapos) == 1:
                            trans_aapos = trans_aapos * len(trans_ids)

                        # GENE NAME
                        if rec_tup.genename != '.':
                            gene = rec_tup.genename.strip()
                        else:
                            gene = None

                        # UNIPROT ACC AND AMINO ACID POSITION IN UNIPROT
                        if rec_tup.Uniprot_acc_Polyphen2.strip() != '.':
                            unip_acs = [acc.strip() for acc in
                                            rec_tup.Uniprot_acc_Polyphen2.split(';')]
                            unp_aapos = [int(pos.strip()) for pos in
                                         rec_tup.Uniprot_aapos_Polyphen2.split(';')]
                        else:
                            unip_acs = None

                        # SIFT SCORE / PREDICTION
                        sifteff = {}
                        if rec_tup.SIFT_score.strip() != '.':
                            sift_score = [float(score.strip()) for score in
                                            rec_tup.SIFT_score.split(';') if score!='.']
                            
#                             sift_pos = [int(e.split(':')[1][1:-1])
#                                         for e in rec_tup.aapos_SIFT.split(';')]

                            for ss, key in zip(sift_score, trans_aapos):
                                flag = False
                                if key not in sifteff:
                                    sifteff[key] = [ss]
                                    flag = True
                                else:
                                    if ss < sifteff[key][0]:
                                        sifteff[key] = [ss]
                                        flag = True
                                if flag is True:
                                    if ss < 0.05:
                                        # Damaging
                                        sifteff[key].append('D')
                                    else:
                                        # Tolerated
                                        sifteff[key].append('T')

                        # POLYPHEN2(HUMVAR) SCORE / PREDICTION
                        if rec_tup.Polyphen2_HVAR_score.strip() != '.':
                            pp2_hvar_score = [score.strip() for score in
                                             rec_tup.Polyphen2_HVAR_score.split(';')]
                            pp2_hvar_pred = [pred.strip() for pred in
                                            rec_tup.Polyphen2_HVAR_pred.split(';')]
                        else:
                            pp2_hvar_score = None
                            pp2_hvar_pred = None

                        # UNIPROT ACC EFFECT PREDICTION MAP
                        hvareff = {}
                        if pp2_hvar_pred:
                            for pos, acc, hvscore, hvpred in \
                                                    zip(unp_aapos,
                                                        unip_acs,
                                                        pp2_hvar_score,
                                                        pp2_hvar_pred):
                                if hvscore != '.':
                                    hvscore = float(hvscore)
                                    if pos not in hvareff:
                                        hvareff[pos] = [hvscore, hvpred,
                                                            acc]
                                    else:
                                        if hvscore > hvareff[pos][0]:
                                            hvareff[pos][0] = hvscore
                                            hvareff[pos][1] = hvpred

                        # YIELD THE ROW FOR DB INSERTION
                        for tid, pos in zip(trans_ids, trans_aapos):
                            if pos in sifteff:
                                ss, sp = sifteff[pos]
                            else:
                                ss, sp = (None, None)
                                
                            if pos in hvareff:
                                eff = hvareff[pos]
                                yield nsfp_data._make((rec_tup.chr,
                                                       int(rec_tup.pos),
                                                       rec_tup.ref,
                                                       rec_tup.alt,
                                                       rec_tup.aaref,
                                                       rec_tup.aaalt, gene,
                                                       tid, pos, eff[2],
                                                       ss, sp, eff[0], eff[1],
                                                       cadd_raw, cadd_phred,
                                                       domains))
                            else:
                                if ss is not None or cadd_raw is not None:
                                    yield nsfp_data._make((rec_tup.chr,
                                                           int(rec_tup.pos),
                                                           rec_tup.ref,
                                                           rec_tup.alt,
                                                           rec_tup.aaref,
                                                           rec_tup.aaalt, gene,
                                                           tid, pos, None,
                                                           ss, sp, None, None,
                                                           cadd_raw,cadd_phred,
                                                           domains)
                                                        )


def main():
    """Main script to create the NSFP database"""

    infile = fileconfig.FILECONFIG['NSFP']
    outfil = dbconfig.DBCONFIG['NSFP']['name']
    desc = 'Create a NSFP sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Directory containing dbNSFP text files')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    t1 = time.time()
    NsfpDB().createdb(args.inname, args.outname, args.force)
    print 'The Time taken for creating the db =', \
            (time.time() - t1) / 60, 'min'
    sys.exit(0)


if __name__ == "__main__":
    main()
