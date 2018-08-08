#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
#!/usr/bin/env python
"""
.. module:: mkclnsnpdb
    :platform: Unix, Windows, MacOSX
    :synopsis: Script to build a heterogenous clinical assocated Phenotype SNP
               sqlite database.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to build a heterogenous clinical assocated Phenotype SNP
sqlite database. To use this script

usage: mkclnsnpdb.py [-h] [-i INNAME] [-o OUTNAME] [--force]

Create a CLNPHESNP sqlite database

optional arguments:
  -h, --help    show this help message and exit
  -i INNAME, --input INNAME
                        Directory name that contains the GAD and NHGRI-GWAS
                        files.
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
  --force               Create new database even if the existing one is newer
                        than the *inname*

**NOTE**: Any existing database with the specified outname will be overwritten.
To be safe one should write the new database to a temporary file which, upon
successful completion, can be renamed to the desired database name.


The created database has two tables - *clnsnp*, *release_info*.
The *clnsnp* table contains the following columns which are all of type text
unless otherwise indicated.

    - snp :        dbSNP rs Ids
    - gene:        Gene name. (Note:- Entry with 'Gene1-Gene2' types denotes
                    that the variant was in an intergenic region and the Gene1
                    & Gene2 are the upstream and downstream genes.)
    - phenotype:   Phenotype/Disease/Trait name associated with the SNP
    - asnstatus: 'Y' if the gene associated with a phenotype and
                   'U' if the author did not specify about it.
    - pubmedid:    Pubmed ID that defines the association between the snp
                   and phenotype(integer)
    - refdb:       Name of the reference database from which the association
                   was extracted
    - refdbid:     Id of the record in the ref db.

The *release_info* table contains the following 3 columns

    - file_name:    File Name
    - file_date:    The date when the file was downloaded from website
    - version:    Version Number
    - entry_count:    Number of entries in the cln table of sqlite DB

The input directory for the program contains two .txt file -
.. [1] all.txt - downloaded from
                                http://geneticassociationdb.nih.gov/all.zip
.. [2] gwascatalog.txt - downloaded from
                                http://www.genome.gov/admin/gwascatalog.txt
"""

from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import fileutils
from gcn.lib.io import anyopen, db
from collections import namedtuple
import csv
import argparse
import sys
import os
import re
import time
from datetime import datetime as dt

fields = ['idx','clnid','snp','gene','phenotype','asnstatus','pubmedid','refdb','reddbid']

SCHEMA = ["""create table clnsnp
               ( %s integer PRIMARY KEY AUTOINCREMENT,
                 %s integer not null,
                 %s text,
                 %s text not null,
                 %s text not null,
                 %s text not null,
                 %s integer not null,
                 %s text not null,
                 %s text not null)
            """%fields,
            """create index clnsnp1 on clnsnp (snp)""",
            """create index clnsnp2 on clnsnp (gene)""",
            """create index clnsnp3 on clnsnp (phenotype)""",
            """create index clnsnp4 on clnsnp (asnstatus)"""
            ]


class ClnSNPDB(db.DB):
    """Class to insert CLNPHESNP data into a sqlite3 database"""

    def createdb(self, inname, outname, force):
        """Create the database

        Args:
            inname (str): Name of the directory containing
                            clinical phenotype data
            outname (str): Name of sqlite3 database
            force (bool): If True overwrite existing database
                            even if it is newer than `inname`
        """
        
        if os.path.exists(outname):
            for fn in os.listdir(inname):
                infile = os.path.join(inname, fn)
                newer = fileutils.file_newer(infile, outname)
                if newer == True:
                    break
        else:
            newer = True

        self.logger = fileconfig.getlogger()

        if not newer and not force:
            self.logger.info('Not Updating. CLNPHESNP' + \
                            'database already uptodate')
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

        self.logger.info('Creating CLNPHESNP database ...')
        self.logger.info('Input directory: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such directory' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)
        
        update = True
        if update:
            self.curs = self.conn.cursor()
            last_clnid = int(self.execute("SELECT max(clnid) FROM clnsnp").fetchone()[0])
        else:
            for s in SCHEMA:
                self.createtable(s, True)

            self.curs = self.conn.cursor()
            last_clnid = 0

        self._clnid = last_clnid
        for entry in self._iterfile():
            self._insertclnsnp(entry)
            self._clnid += 1
        self._deletersids(list(set(self.gadrs.keys()) \
        & set(self.nhgrirs.keys())))  # delete redundent rsId

        self.conn.commit()

        # Add version details
        gad_entry_cnt = self.execute("SELECT COUNT(clnid) FROM \
                        clnsnp WHERE refdb='GAD'").fetchall()[0][0]
        gwas_entry_cnt = self.execute("SELECT COUNT(clnid) FROM \
                        clnsnp WHERE refdb='NHGRI'").fetchall()[0][0]
        if not update:
            fd = dt.fromtimestamp(os.path.getmtime(os.path.join(self.inname, \
                                                                'all.txt'))\
                                  ).strftime('%Y-%m-%d')
            version = "v%s_%s" % tuple(fd.split('-')[:2])
            self.set_version('NCBI_GAD', fd, version, gad_entry_cnt)
            
        fd = dt.fromtimestamp(os.path.getmtime(os.path.join(self.inname, \
                                                        'gwascatalog.txt'))\
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        self.set_version('NHGRI_GWAS', fd, version, gwas_entry_cnt)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... CLNPHESNP database created')

    def _parsersId(self, data):
        """Internal method. Do not use"""
        snp = set()
        if data:
            if data[2:].isdigit() == False:
                for ele in re.split('\/|\s|\,|\;|\(|\position [0-9]+', data):
                    ele = ele.lower()
                    if ele.startswith('rs'):
                        sid = re.split('\*|\-|\)|[a,t,g,c]', ele)[0]
                        if sid.strip('rs'):
                            snp.add(sid)
                    elif ele.isdigit():
                        if 'haplotype' not in data:
                            snp.add('rs' + ele)
            else:
                snp.add(data)
        return list(snp)

    def _iterfile(self):
        """Internal method. Do not use"""

        self.gadrs = {}
        self.nhgrirs = {}
        fields = """snp gene phen asnsta pubid refdb refdbid""".split()
        clnsnp = namedtuple('Clnsnp', fields)
        for filename in os.listdir(self.inname):
            #print filename
            infile = os.path.join(self.inname, filename)
            if filename == 'all.txt':
                with anyopen.openfile(infile) as stream:
                    for rec in csv.reader(stream, delimiter='\t'):
                        if len(rec) > 1:
                            if 'Association(Y/N)' in rec:
                                continue
                            if rec[1] != 'N' and rec[8] and rec[2] and \
                                                    rec[13].isdigit():
                                refdbid = rec[0]
                                genes = rec[8].split(',')
                                phen = rec[2]
                                snps = self._parsersId(rec[28])
                                if rec[1]:
                                    asnsta = rec[1]
                                else:
                                    asnsta = 'U'
                                pubid = rec[13]

                                for gene in genes:
                                    gene = gene.strip()
                                    if gene:
                                        if snps:
                                            for snp in snps:
                                                self.gadrs[snp] = [\
                                        gene.upper(), pubid, \
                                        phen.lower().strip()]
                                                yield clnsnp._make((snp,\
                                 gene, phen, asnsta, pubid, 'GAD', refdbid))
                                        else:
                                            yield clnsnp._make((None, gene,\
                                    phen, asnsta, pubid, 'GAD', refdbid))
            elif filename == 'gwascatalog.txt':
                with anyopen.openfile(infile) as stream:
                    for rec in csv.reader(stream, delimiter='\t'):
                        if len(rec) > 1:
                            if rec[14] != ' - ':
                                if 'Date Added to Catalog' in rec or 'DATE ADDED TO CATALOG' in rec:
                                    continue
                                for snpid in rec[21].split(','):
                                    self.nhgrirs[snpid] = [rec[14].upper(),\
                                                    rec[1], rec[7].lower().strip()]
                                    yield clnsnp._make((snpid, rec[14], rec[7],\
                                                'Y', rec[1], 'NHGRI', rec[0]))

    def _insertclnsnp(self, entry, sid=None):
        """Internal method. Do not use"""
        
        cntsql = "select count(distinct snp) from clnsnp where snp=? and refdb = ?"
        prev  = self.curs.execute(cntsql,(entry.snp, entry.refdb)).fetchone()[0]
        
        if not prev:
            print 'adding SNP_ID [%s]' % entry.snp
            if sid is None:
                sid = self._clnid + 1

            rginsert = "insert into clnsnp (%s) values (%s)"%(','.join(fields[1:]),','.join('?'*8))

            self.curs.execute(rginsert, (sid, entry.snp, entry.gene, entry.phen, \
                                     entry.asnsta, int(entry.pubid), \
                                     entry.refdb, entry.refdbid))

    def _deletersids(self, rsidlist, refdb='GAD'):
        """Internal method. Do not use"""
        for rs in rsidlist:
            delstm = "delete from clnsnp where snp=? and refdb = ?"
            self.curs.execute(delstm, (rs, refdb))


def main():
    """Main script to create the Clinical Phenotype of SNPs database"""

    infile = fileconfig.FILECONFIG['CLNPHESNP']
    outfil = dbconfig.DBCONFIG['CLNPHESNP']['name']
    desc = 'Create a CLNPHESNP sqlite database'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        default=infile,
                        help='Directory that contains the GAD and GWAS files')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default=outfil,
                        help='Name of sqlite database to create')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')
    args = parser.parse_args()
    ClnSNPDB().createdb(args.inname, args.outname, args.force)
    sys.exit(0)


if __name__ == "__main__":
    main()
