#!/usr/bin/env python
"""
Script to build an sqlite database from a dbsnp or 1000 genomes vcf file.
The script ignores all sample information and only processes those INFO fields
that are described in the meta INFO section of the vcf file.  Additionally
QUAL and FILTER fields are also ignored since these tend to be superfluous
in the files that are intended to be processed by this script

To use this script

usage: mksnpdb.py [-h] [-i INNAME [-n DBNAME]] [-o OUTNAME]

Create an sqlite database with the input vcf file

optional arguments:
  -h, --help            show this help message and exit
  -i INNAME, --input INNAME
                        Input VCF file
  -n DBNAME  --name DBNAME
                        Name of the database (either -i or -n should be
                                              specified)
  -o OUTNAME, --output OUTNAME
                        Name of sqlite database to create
 --force               Create new database even if the existing one is newer
                        than the *inname*

"""
from __future__ import absolute_import
from gcn.lib.io import db
from gcn.lib.io.vcf import VCFParser, toint, tofloat
import argparse
import sys
import os
from datetime import datetime as dt

TYPEMAP = {
    'Integer': 'int',
    'Float': 'float',
    'String': 'text',
    'Flag': 'int',
    'text': 'text',
}

INFOMAP = {'Integer': toint,
           'Float': tofloat,
           'Flag': bool,
           'Character': lambda e: e,
           'String': lambda e: e
           }

RTYPEMAP = dict((value, key) for key, value in INFOMAP.items())

INFOTABLE = """create table info
                (id text PRIMARY KEY not null,
                 type text not null,
                 count text not null,
                 description text)"""

class SNPDB(db.DB):

    def _getschema(self):
        schema = ["""create table snps
                   (idx integer PRIMARY KEY not null,
                    chrom text not null,
                    pos int not null,
                    id text not null,
                    ref text not null,
                    alt text not null,
                  """,
                  """create index chromposindex on snps (chrom, pos)""",
                  """create index snpid on snps (id)"""
                  ]
        return schema

    def _makedb(self):
        """Internal method. Do not use"""

        self.logger.info('Creating SNP database from ...')
        self.logger.info('Input file: %s' % self.inname)
        if not os.path.exists(self.inname):
            self.logger.error('%s: No such file' % self.inname)
            self.logger.error('Database not created')
            sys.exit(1)

        self.load(db=self.outname)
        self.conn.text_factory = str
        self._vcf = VCFParser(self.inname)
        self._infokeys, self._schema = self._createschema()
        for s in self._schema:
            self.createtable(s, True)

        stmt = "insert into snps values (%s)" % (','.join('?' * (6 + len(self._infonum))))

        self.curs = self.conn.cursor()

        cache = []
        n = 0
        next = self._vcf.__next__
        while 1:

            try:
                rows = self._insertentry(next())
                rows.insert(0,n)
                cache.append(rows)
            except StopIteration:
                break

            n += 1
            if not n % 5000000:
                self.curs.executemany(stmt, cache)
                self.logger.info('Processed %d entries' % n)
                cache = []

        if cache:
            self.curs.executemany(stmt, cache)
        self.logger.info('Processed %s entries' % n)

        # Add version details
        fd = dt.fromtimestamp(os.path.getmtime(self.inname)
                              ).strftime('%Y-%m-%d')
        if 'dbsnp' in self.inname:
            version = 'v' + str(self._vcf.meta['dbSNP_BUILD_ID'])
        elif 'ESP6500' in self.inname:
            fn = os.path.splitext(os.path.split(self.inname)[1])[0]
            version = fn.split('.')[0]
        elif 'ExAC' in self.inname:
            fn = os.path.splitext(os.path.split(self.inname)[1])[0]
            version = fn.split('.')[1]
        else:
            version = "v%s_%s" % tuple(fd.split('-')[:2])

        self.set_version(os.path.split(self.outname)[1], fd, version, n)

        self.conn.commit()
        self.curs.close()
        self.conn.close()
        self.logger.info('... SNP database created')

    def _createschema(self):
        """Internal method to create schema based on the INFO fields in the VCF
        file"""
        schema = self._getschema()
        info = self._vcf.meta['INFO']
        keys = list(info.keys())
        keys.sort()
        infocols = []
        cursor = self.conn.cursor()
        cursor.execute(INFOTABLE)
        stmt = 'insert into info values (?,?,?,?)'
        self._infonum = {}
        for key in keys:
            knumb, ktype, kdesc = info[key]
            if knumb not in '01':
                ktype = 'text'
                self._infonum[key] = True
            else:
                ktype = TYPEMAP[RTYPEMAP.get(ktype, 'text')]
                self._infonum[key] = False
            cursor.execute(stmt, ('info_' + key, ktype, knumb, kdesc))
            infocols.append('info_%s %s' % (key, ktype))
        self.conn.commit()
        cursor.close()
        schema[0] += ','.join(infocols)
        schema[0] += ')'
        self.logger.info(schema)
        return keys, schema

    def _insertentry(self, rec):
        self._vcf.parseinfo(rec)
        args = [rec.chrom, rec.pos, ';'.join(rec.id),
                rec.ref, ','.join(rec.alt)]
        infn = self._infonum
        info = rec.info.get
        for key in self._infokeys:
            v = info(key, None)
            if v is not None:
                if infn[key]:
                    v = ','.join(str(el) for el in v)
            args.append(v)
        return args


def main():
    """Main script to create the database from a vcf file"""
    from gcn.etc import fileconfig, dbconfig

    desc = 'Create a sqlite database from a vcf file'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', dest='inname', type=str,
                        help='Input vcf file')
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        help='Name of sqlite database to create')
    parser.add_argument('-n', '--name', dest='dbname', type=str,
                        help='Name of the db to be created - ' +
                        'CLINVARDB, HGMDDB, DBSNP, KGDB, ESP, EXAC, etc')
    parser.add_argument('--force', action='store_true',
                        help='Create new database even if the existing one ' +
                             'is newer than *inname*')

    args = parser.parse_args()
    dbname = None
    if args.dbname:
        dbname = args.dbname.upper()

    if not args.inname:
        if dbname:
            args.inname = fileconfig.FILECONFIG[dbname]
        else:
            raise RuntimeError('Either input file or database name should ' +
                                   'be specified')
    if not args.outname:
        if dbname:
            args.outname = dbconfig.DBCONFIG[dbname]['name']
        else:
            raise RuntimeError('Either output file or database name should ' +
                                   'be specified')

    print args.inname, args.outname
    SNPDB().createdb(args.inname, args.outname, args.force)
    sys.exit(0)

if __name__ == "__main__":
    main()
