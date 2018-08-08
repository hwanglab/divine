"""
Script to build a splicedb sqlite database. To use this script

usage: mksplicedb.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a SPLICEDB sqlite database

optional arguments:
  -h, --help            show this help message and exit
  -i INNAME, --input INNAME
                        RefGene File
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
  --force               Create new database even if the existing one is newer
                        than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has four tables - *splice*, *ese_stat*, *ess_stat* and
*release_info*. The *splice* table contains the following columns which are all
of type text unless otherwise indicated.

    - chr:    Chromosome on which the gene is present (e.g chr1, chr2)
    - pos:    Genomic position-hg19 (integer)
    - refseq_acc:    Refseq mRMA accession
    - strand:       '+' or '-' depending on which strand the gene is
    - gene:    Gene name
    - annotation:    Values can be [ESE/ESS] where in ESE is Exonic Splice
                     Enhance Site and ESS is Exonic Splice Silencer Site.

The *ese_stat* table contains the following columns

    - chr:    Chromosome on which the gene is present (e.g chr1, chr2)
    - refseq_acc:    Refseq mRMA accession
    - gene:    Gene name
    - mrna_len:    mRNA Sequence Length
    - exon_cnt:    Exon count for the gene
    - ese_cnt:    Exonic Splice Enhancer Site count for the gene

The *ess_stat* table contains the following columns

    - chr:    Chromosome on which the gene is present (e.g chr1, chr2)
    - refseq_acc:    Refseq mRMA accession
    - gene:    Gene name
    - mrna_len:    mRNA Sequence Length
    - exon_cnt:    Exon count for the gene
    - ese_cnt:    Exonic Splice Silencer Site count for the gene

The *release_info* table contains the following columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the utrdb table of sqlite DB
"""

from gcn.etc import fileconfig, dbconfig
from gcn.lib.io import db
from gcn.lib.utils import fileutils
from gcn.lib.databases.refgene import Refgene
from gcn.lib.databases.refmrna import Refmrna
from gcn.data.ese_human import ESE_MOTIFS
from gcn.data.ess_human import ESS_MOTIFS
from multiprocessing import Pool
import argparse
import re
import os
import sys
import time
from datetime import datetime as dt

ESE_MER = 6
ESS_MER = 6

fields_spl = ['idx','chr','pos','refseq_acc','strand','gene','annotation']
fields_ese = ['idx','chrom','refseq_acc','gene','mrna_len','exon_cnt','ese_cnt']
fields_ess = ['idx','chrom','refseq_acc','gene','mrna_len','exon_cnt','ess_cnt']

SCHEMA = ["""create table splice
               ( idx integer PRIMARY KEY AUTOINCREMENT,
                 chr text not null,
                 pos integer not null,
                 refseq_acc text not null,
                 strand text not null,
                 gene text not null,
                 annotation text not null)
            """,
            """create table ese_stat
                ( idx integer PRIMARY KEY AUTOINCREMENT,
                  chrom text not null,
                  refseq_acc text not null,
                  gene text not null,
                  mrna_len integer not null,
                  exon_cnt integer not null,
                  ese_cnt integer not null)
            """,
            """create table ess_stat
                ( idx integer PRIMARY KEY AUTOINCREMENT,
                  chrom text not null,
                  refseq_acc text not null,
                  gene text not null,
                  mrna_len integer not null,
                  exon_cnt integer not null,
                  ess_cnt integer not null)
            """,
            """create index index1 on splice (chr,pos)""",
            """create index index2 on splice (chr,pos,refseq_acc)""",
            """create index index3 on splice (gene)""",
            """create index index4 on ese_stat (chrom,gene)""",
            """create index index5 on ese_stat (chrom,refseq_acc)""",
            """create index index6 on ess_stat (chrom,gene)""",
            """create index index7 on ess_stat (chrom,refseq_acc)""",
            ]


class SpliceRegSite(db.DB):
    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): The UCSC's refgene text file. Because the splice
                            regulating sites are computed for the
                            transcripts defined in the refgene file.
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database even if it
                          is newer than any of  `inname`
        """
        self.logger = fileconfig.getlogger()
        if os.path.exists(outname):
            newer = False
            newer = fileutils.file_newer(inname, outname)
        else:
            newer = True

        if not newer and not force:
            self.logger.info('Not Updating. SPLICE database already uptodate')
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

        self.logger.info('Creating SPLICE database ...')
        self.logger.info('Input directory: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such directory' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)
        self.load(db=self.outname)
        for s in SCHEMA:
            self.createtable(s, True)
        self.curs = self.conn.cursor()
        query1 = 'insert into splice (%s) values (?,?,?,?,?,?)'%(','.join(fields_spl[1:]))
        query2 = 'insert into ese_stat (%s) values (?,?,?,?,?,?)'%(','.join(fields_ese[1:]))
        query3 = 'insert into ess_stat (%s) values (?,?,?,?,?,?)'%(','.join(fields_ess[1:]))
        ese_prev_acc = None
        ess_prev_acc = None
        entries1 = []
        entries2 = []
        entries3 = []
        entry_cnt = 0
        for entry in self._get_splice_data():
            seq_len = entry.pop()
            site_cnt = entry.pop()
            exon_cnt = entry.pop()
            entry_cnt += 1
            acc = entry[2]
            gene = entry[4]
            chrom = entry[0]
            annt = entry[5]
            if annt == 'ESE':
                if ese_prev_acc is None:
                    entries2.append((chrom, acc, gene, seq_len, exon_cnt,
                                     site_cnt))
                    if len(entries2) > 10000:
                        self.curs.executemany(query2, entries2)
                        entries2 = []
                    ese_prev_acc = acc
                elif ese_prev_acc != acc:
                    entries2.append((chrom, acc, gene, seq_len, exon_cnt,
                                     site_cnt))
                    if len(entries2) > 10000:
                        self.curs.executemany(query2, entries2)
                        entries2 = []
                    ese_prev_acc = acc
            elif annt == 'ESS':
                if ess_prev_acc is None:
                    entries3.append((chrom, acc, gene, seq_len, exon_cnt,
                                     site_cnt))
                    if len(entries3) > 10000:
                        self.curs.executemany(query3, entries3)
                        entries3 = []
                    ess_prev_acc = acc
                elif ess_prev_acc != acc:
                    entries3.append((chrom, acc, gene, seq_len, exon_cnt,
                                     site_cnt))
                    if len(entries3) > 10000:
                        self.curs.executemany(query3, entries3)
                        entries3 = []
                    ess_prev_acc = acc
            entries1.append(tuple(entry))
            if len(entries1) > 50000:
                self.curs.executemany(query1, entries1)
                entries1 = []
        if entries1:
            self.curs.executemany(query1, entries1)
        if entries2:
            self.curs.executemany(query2, entries2)
        if entries3:
            self.curs.executemany(query3, entries3)

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(self.inname)
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('Splicedb', fd, version, entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... SPLICE database created')

    def _get_exons(self, rgene, seq):
        """Internal Method. Builds exon boundary and sequence dictionary.
        Args:
            - rgene:    Named Tuple (chrom,refseq_acc,strand,str_transcript,
                                    end_transcript,str_cds,end_cds,exon_number,
                                    str_exon_pos,end_exon_pos,gene,frames)
            - seq:    mRNA Sequence
        Returns:
            - exons_map:    A dictionary where key is exon boundary
                             (start position, end position) and value is
                            the exon sequence. For eg. {
                            (100010,100220):'GATCAGACAGTA...',
                            (100900,101400):'ATAGACAGT..'}
        """
        exons = {}
        prev = 0
        str_exon = rgene.str_exon_pos.split(',')
        end_exon = rgene.end_exon_pos.split(',')
        if rgene.strand == '+':
            for sp, ep in zip(str_exon, end_exon):
                exon_len = int(ep) - int(sp)
                exons[(int(sp), int(ep))] = seq[prev:prev + exon_len].upper()
                prev += exon_len
        elif rgene.strand == '-':
            str_exon.reverse()
            end_exon.reverse()
            for sp, ep in zip(end_exon, str_exon):
                exon_len = int(sp) - int(ep)
                exons[(int(sp), int(ep))] = seq[prev:prev + exon_len].upper()
                prev += exon_len
        return exons

    def _get_ese(self, (exons_map, rgene_info)):
        """Internal Method. Computes the ESE sites present in the given exon.
        Args:
            - (exons_map, rgene_info): 'exons_map' is a dictionary where key is
                                        exon boundary (start position,
                                        end position) and value is the exon
                                        sequence. For eg. {
                                        (100010,100220):'GATCAGACAGTA...',
                                         (100900,101400):'ATAGACAGT..'}
                                        'rgene_info' is a named tuple for the
                                        refgene table.
        Returns:
            - ese_sites:    List of genomic position that are part of ESE site
            - ese_cnt:    Number of ESE hexamers present in the given exon
        """
        ese_cnt = 0
        ese_sites = set()
        for pos, seq in exons_map.iteritems():
            str_pos, end_pos = pos
            for ese in ESE_MOTIFS:
                idxlist = [s.start() for s in re.finditer(ese, seq)]
                if rgene_info[2] == '+':
                    for idx in idxlist:
                        for coord in range(str_pos + idx, str_pos + idx +
                                            ESE_MER):
                            ese_sites.add(coord)
                elif rgene_info[2] == '-':
                    for idx in idxlist:
                        for coord in range(str_pos - idx - ESE_MER + 1,
                                            str_pos - idx + 1):
                            ese_sites.add(coord)
                if len(idxlist) > 0:
                    ese_cnt += 1
        return (list(ese_sites), ese_cnt, rgene_info)

    def _get_ess(self, (exons_map, rgene_info)):
        """Internal Method. Computes the ESS sites present in the given exon.
        Args:
            - (exons_map, rgene_info): 'exons_map' is a dictionary where key is
                                        exon boundary (start position,
                                        end position) and value is the exon
                                        sequence. For eg. {
                                        (100010,100220):'GATCAGACAGTA...',
                                        (100900,101400):'ATAGACAGT..'}
                                        'rgene_info' is a named tuple for the
                                        refgene table.
        Returns:
            - ess_sites:    List of genomic position that are part of ESS site
            - ess_cnt:    Number of ESS hexamers present in the given exon
        """
        ess_cnt = 0
        ess_sites = set()
        for pos, seq in exons_map.iteritems():
            str_pos, end_pos = pos
            for ess in ESS_MOTIFS:
                idxlist = [s.start() for s in re.finditer(ess, seq)]
                if rgene_info[2] == '+':
                    for idx in idxlist:
                        for coord in range(str_pos + idx, str_pos + idx +
                                            ESS_MER):
                            ess_sites.add(coord)
                elif rgene_info[2] == '-':
                    for idx in idxlist:
                        for coord in range(str_pos - idx - ESS_MER + 1,
                                           str_pos - idx + 1):
                            ess_sites.add(coord)
                if len(idxlist) > 0:
                    ess_cnt += 1
        return (list(ess_sites), ess_cnt, rgene_info)

    def _get_splice_data(self):
        """Iterates over the RefGene DB and returns a tuple containing
            the predicted ESE site details
        Returns:
            - data:    A List containing (chrom, pos, refseq_acc,
                        strand, gene, ESE/ESS, Exon_count, Number of
                        ESE/ESS_sites, mRNA_length)
        """
        CHROM = ['chr' + str(e) for e in range(1, 23) + ['X', 'Y']]
        core = 5
        p = Pool(processes=core)
        inlist = []
        self.refgene = Refgene()
        self.refmrna = Refmrna()
        for rgene in self.refgene.iterate_db():
            if rgene.chrom not in CHROM:
                continue
            refacc = rgene.refseq_acc
            gid = rgene.gid
            seq = self.refmrna.retrieve(refacc, gid)
            if seq:
                exons_map = self._get_exons(rgene, seq.sequence)
                inlist.append((exons_map, [rgene.chrom, refacc, rgene.strand,
                                           rgene.gene, rgene.exon_number,
                                           len(seq.sequence)]))
                if len(inlist) == 10000:
                    for ese_sites, ese_cnt, rgene_info in \
                                            p.map_async(unwrap_self_get_ese,
                                                        inlist).get():
                        for site in ese_sites:
                            data = [rgene_info[0]] + [int(site)] +\
                            rgene_info[1:4] + ['ESE'] + [rgene_info[4]]\
                            + [ese_cnt] + [rgene_info[-1]]
                            yield data

                    for ess_sites, ess_cnt, rgene_info in\
                                            p.map_async(unwrap_self_get_ess,
                                                        inlist).get():
                        for site in ess_sites:
                            data = [rgene_info[0]] + [int(site)] +\
                            rgene_info[1:4] + ['ESS'] + [rgene_info[4]] +\
                            [ess_cnt] + [rgene_info[-1]]
                            yield data
                    inlist = []
            else:
                self.logger.info('SEQ NOT FOUND for ' + refacc + '\t'
                                 + str(gid))
        if inlist:
            for ese_sites, ese_cnt, rgene_info in\
                                            p.map_async(unwrap_self_get_ese,
                                                        inlist).get():
                for site in ese_sites:
                    data = [rgene_info[0]] + [int(site)] + rgene_info[1:4] +\
                    ['ESE'] + [rgene_info[4]] + [ese_cnt] +\
                            [rgene_info[-1]]
                    yield data

            for ess_sites, ess_cnt, rgene_info in\
                                            p.map_async(unwrap_self_get_ess,
                                                        inlist).get():
                for site in ess_sites:
                    data = [rgene_info[0]] + [int(site)] + rgene_info[1:4]\
                            + ['ESS'] + [rgene_info[4]] + [ess_cnt] +\
                            [rgene_info[-1]]
                    yield data
            inlist = []


def unwrap_self_get_ese(kwarg):
    """Method to call class method _get_ese() during multiprocessing"""
    srs = SpliceRegSite()
    return srs._get_ese(kwarg)


def unwrap_self_get_ess(kwarg):
    """Method to call class method _get_ess() during multiprocessing"""
    srs = SpliceRegSite()
    return srs._get_ess(kwarg)


def main():
    """Main script to create the SPLICE database"""

    infile = fileconfig.FILECONFIG['REFGENE']
    outfil = dbconfig.DBCONFIG['SPLICE']['name']
    desc = 'Create a SPLICE sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Refgene files')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    SpliceRegSite().createdb(args.inname, args.outname,
                             args.force)
    sys.exit(0)


if __name__ == "__main__":
    t1 = time.time()
    main()
    print 'Time Taken to create Splicedb:', (time.time() - t1) / 60, 'min'
