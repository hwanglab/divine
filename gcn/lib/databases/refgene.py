#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: refgene
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing refgene information from an sqlite database

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Class for accessing refgene information from an sqlite database. A refgene
is an instance of the `Gene` namedtuple
"""

from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
from collections import namedtuple
import os

Gene = namedtuple('Gene', ['gid', 'name', 'chrom', 'strand',
                           'txStart', 'txEnd', 'cdsStart',
                           'cdsEnd', 'exonCount', 'score',
                           'name2', 'cdsStartStat',
                           'cdsStartEnd', 'exons'])
refgene = namedtuple('refgene', 'gid,chrom,refseq_acc,strand,\
                                str_transcript,end_transcript,\
                                str_cds,end_cds,exon_number,\
                                str_exon_pos,end_exon_pos,gene,\
                                frames,error')


#Gene.__doc__ = \
"""
Gene('name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
    'exonCount', 'score', 'name2', 'cdsStartStat', 'cdsStartEnd', 'exons')

The fields are as follows:
    - gid:          Unique transcript identifier
    - name:         gene name
    - chrom:        chromosome number
    - strand:       '+' or '-'
    - txStart:      transcript start position
    - txEnd:        transcript end position
    - cdsStart:     Coding region start position (integer)
    - cdsEnd:       Coding region end position (integer)
    - exonCount:    Number of exons in the transcript (integer)
    - score:        An integer of unknown meaning (integer)
    - name2:        Alternative name for the gene
    - cdsStartStat: An enum from one of
                    ('none','unk','incmpl','cmpl')
    - cdsEndStat:   An enum from one of
                    ('none','unk','incmpl','cmpl')
    - exons:        List of all exons associated with the gene.
                    Each item in the list is a 3-tuple of
                    start position, end position and frame
"""

def get_ucsc_chrom():
  CHROMOSOMES = ['chr%d' % s for s in range(1, 23)] + ['chrX', 'chrY', 'chrM']
  ICHROMOSOMES = dict((s, i) for (i, s) in enumerate(CHROMOSOMES))
  return CHROMOSOMES, ICHROMOSOMES

def to_ucsc_chrom(chrom):
    if not chrom.startswith('chr'):
       if chrom=='MT':
           chrom = 'chrM'
       else:
           chrom = 'chr' + chrom
    return chrom

def to_grach(chrom):
    chromosomes, _ = get_ucsc_chrom()
    if chrom in chromosomes:
        if chrom == 'chrM':
            chrom2 = 'MT'
        else:
            chrom2 = chrom[3:]
    else:
        chrom2 = None
    return chrom2

class Refgene(db.DB):
    """Class to retrieve gene information from REFGENE database"""

    def __init__(self, candgene=[]):
        """Class initialization

        Argument

            name (string): Name of the database or filename for the database
        """
        name = 'REFGENE'
        self.candgene = candgene
        super(Refgene, self).__init__()
        if name in DBCONFIG:
            self.load(name=name)
        elif os.path.exists(name):
            self.load(db=name)
        else:
            raise ValueError('No such database %s' % name)

    def gene_by_location(self, chromosome, position):
        """Retrieve details of all genes associated with a given
        chromosomal location.

        Args:
            chromosome (str):   Chromosme name in the form chrN where
                                N=1-23,X,Y
            position (integer): Location on the chromosome, using 0-based
                                indexing

        Returns:
            List of genes associated with the position. Each gene is an
            instance of a `Gene` namedtuple
        """
        if not chromosome.startswith('chr'):
            chromosome = 'chr' + chromosome
        stmt = 'SELECT * FROM refgene WHERE chrom=(?) AND ' + \
               '(?)>=txStart AND (?)<=txEnd AND error = 0'
        results = self.execute(stmt, (chromosome, position, position))
        genes = []
        for r in results:
            exons = self._exons_from_gid(r[0])
            genes.append(Gene._make(r[0:13] + (exons,)))
        return genes

    def _check_defect(self, name):
        """Internal method. Do not use """
        stmt = 'SELECT count(name) FROM defect WHERE name=(?)'
        results = self.execute(stmt, (name,)).fetchall()
        if results[0][0] > 0:
            defect = 'Mismatch in mRNA seq length'
            status = True
        else:
            defect = None
            status = False
        status = False
        return status

    def gene_by_name(self, name):
        """Retrieve details of genes with a specified name
        Args:
            name (str): Gene name

        Returns:
            List of genes with the requested name. Each gene is an
            instance of a `Gene` namedtuple
        """
        stmt = 'SELECT * FROM refgene WHERE name2=(?)'
        genes = []
        for r in self.execute(stmt, (name,)):
            exons = self._exons_from_gid(r[0])
            genes.append(Gene._make(r[1:] + (exons,)))

        return genes

    def get_cds_length(self, trans_acc):
        """Returns length of CDS of the transcript
        Args:
            trans_acc (str):    Transcript Accession (NM_)

        Retuns:
            size (int):    CDS Length
        """
        size = 0
        stmt = "SELECT r.cdsStart, r.cdsEnd, e.end, e.start, e.frame \
                FROM refgene r, exons e WHERE r.gid = e.gid AND r.name = '%s' \
                AND e.frame >= 0 ORDER BY e.start" % trans_acc
        results = self.execute(stmt).fetchall()
        for idx, row in enumerate(results):
            if idx == 0:
                size += row[2] - row[0]
            elif idx + 1 == len(results):
                size += row[1] - row[3]
            else:
                size += row[2] - row[3]
        return size

    def get_max_cds_length(self):
        stmt = "select gene, tx, MAX(cds_len) \
                from \
                (select \
                    refgene.name2 as gene, \
                    exons.name as tx, \
                    sum(exons.end-exons.start) as cds_len\
                from exons \
                    inner join refgene on refgene.gid = exons.gid \
                    group by exons.name)\
                group by gene"
        results = self.execute(stmt).fetchall()
        genes_by_cds_len = {}
        for r in results:
            genes_by_cds_len[r[0]] = int(r[2])
        return genes_by_cds_len

    def iterate_db(self, chrom=None, gene=None):
        """Iterates over the databases and retrieves the columns"""
        if chrom:
            stmt = "SELECT * FROM refgene where chrom = '%s'" % (chrom)
        elif gene:
            stmt = "SELECT * FROM refgene where name2 = '%s'" % (gene)
        else:
            stmt = "SELECT * FROM refgene"
        results = self.execute(stmt)
        for row in results:
            gid = row[0]
            exons = self._exons_from_gid(gid)
            exon_starts = None
            exon_ends = None
            frame = None
            for ele in exons:
                if exon_starts and exon_ends:
                    exon_starts += ',' + str(ele[0])
                    exon_ends += ',' + str(ele[1])
                    frame += ',' + str(ele[2])
                else:
                    exon_starts = str(ele[0])
                    exon_ends = str(ele[1])
                    frame = str(ele[2])
            rgene = refgene(gid, str(row[2]), str(row[1]), str(row[3]), \
                            str(row[4]), str(row[5]), str(row[6]), \
                            str(row[7]), str(row[8]), exon_starts, \
                            exon_ends, str(row[10]), frame, str(row[13]))

            yield rgene

    def get_gene_exons(self, chrom, gene):
        """Returns the exon coodinates along with its number for the gene
        Args:
            chrom(str):    Chromosome
            gene(str):    Gene Name

        Returns:
            exoninfo(dict):    Exon coordinates along with the exon number
                               where key=exon coordinates and value is
                               exon number.
        """
        stmt = 'SELECT DISTINCT(e.start), e.end FROM refgene r, exons e \
        WHERE r.gid = e.gid AND r.name2 = (?) AND r.chrom = (?) \
        ORDER BY e.start'
        exoninfo = {}
        cnt = 0
        for ele in self.execute(stmt, (gene, chrom, )).fetchall():
            cnt += 1
            key = (ele[0], ele[1])
            exoninfo[key] = cnt
        return exoninfo

    def get_genes(self, chrom, spos, epos):
        """Returns list of genes that resides in the given coordinate
        interval along with exon details like exon numbers, exon coordinates.
        Args:
            chrom(str):    Chromosome
            spos(int):    Start position
            epos(int):    End position

        Returns:
            geneinfo(dictionary):    Genes in the interval along with
                                     the exon details.
        """
        stmt = 'SELECT DISTINCT(e.start), e.end, r.name2 FROM refgene r, \
        exons e WHERE r.gid = e.gid AND ((start BETWEEN (?) AND (?)) OR \
        (end BETWEEN (?) AND (?))) AND r.chrom = (?)'
        if 'chr' not in chrom:
            chrom = 'chr' + str(chrom)
        geneinfo = {}
        for ele in self.execute(stmt, (spos, epos, spos, epos, \
                                    chrom, )).fetchall():
            if ele[2] not in geneinfo:
                geneinfo[ele[2]] = [[ele[0], ele[1]]]
            else:
                geneinfo[ele[2]].append([ele[0], ele[1]])
        for gene in geneinfo:
            exoninfo = self.get_gene_exons(chrom, gene)
            exons = geneinfo[gene]
            exnum = [exoninfo[tuple(ele)] for ele in exons]
            val = zip(exons, exnum)
            val.sort()
            geneinfo[gene] = val
        return geneinfo

    def _exons_from_gid(self, gid):
        """Internal method"""
        stmt = 'SELECT * FROM exons WHERE gid=(?) ORDER BY start ASC'
        exons = []
        for ex in self.execute(stmt, (gid,)):
            exons.append((ex[4], ex[5], ex[6])) #start,end,frame
        return exons

    def get_strand(self, transcipt_id):
        """For the given transcript it retrieves the strand"""
        stmt = "SELECT strand from refgene where name = '%s' limit 1"\
                                                 % (transcipt_id)
        result = self.execute(stmt).fetchall()
        if result:
            return result[0][0]
        else:
            None

    def _compute_coding_exons(self, rec):
        """Internal Method. This returns the coding exonic boundary
        for the given record"""
        cds_exons = []
        for str_pos, end_pos, frame in zip(rec.str_exon_pos.split(','), \
                                           rec.end_exon_pos.split(','), \
                                           rec.frames.split(',')):
            str_pos = int(str_pos)
            end_pos = int(end_pos)
            frame = int(frame)
            str_cds = int(rec.str_cds)
            end_cds = int(rec.end_cds)
            if frame == -1:
                continue
            if str_pos <= str_cds < end_pos:
                cds_exons.append([rec.chrom, str_cds, end_pos - 1, \
                                  rec.strand, rec.gene, rec.refseq_acc])
            elif str_pos <= end_cds < end_pos:
                cds_exons.append([rec.chrom, str_pos, end_cds - 1, \
                                  rec.strand, rec.gene, rec.refseq_acc])
            else:
                cds_exons.append([rec.chrom, str_pos, end_pos - 1, \
                                  rec.strand, rec.gene, rec.refseq_acc])
        return cds_exons

    def get_exons(self, outbedfile, genelist=None):
        """For the given gene list extract the exact coding exonic boundary
        and if genelist is not provided compute for all human genes and writes
        in bed format to the output file"""
        bedfile = open(outbedfile, 'w')
        bedfile.write('\t'.join(['#CHROM', 'CDS_START', 'CDS_END', 'STRAND', \
                                 'GENE', 'TRANSCRIPT']) + '\n')
        chromlist = ['chr' + str(ele) for ele in range(1, 23)] + \
                            ['chrX', 'chrY']
        if genelist:
            for gname in genelist:
                for rec in self.iterate_db(gene=gname):
                    if rec.chrom in chromlist and rec.refseq_acc[:2] == 'NM':
                        cds_exons = self._compute_coding_exons(rec)
                        for cds in cds_exons:
                            cds = [str(ele) for ele in cds]
                            bedfile.write('\t'.join(cds) + '\n')
        else:
            for chrno in chromlist:
                for rec in self.iterate_db(chrom=chrno):
                    if rec.refseq_acc[:2] == 'NM':
                        cds_exons = self._compute_coding_exons(rec)
                        for cds in cds_exons:
                            cds = [str(ele) for ele in cds]
                            bedfile.write('\t'.join(cds) + '\n')

    def get_gid(self, chrom, pos, refseq_acc):
        """Returns the gid"""
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        stmt = "SELECT gid FROM refgene WHERE chrom = '%s' and %d >= txStart" \
               % (chrom, pos) + " and %d <= txEnd and "\
               "name = '%s'" % (pos, refseq_acc)
        gid = self.execute(stmt).fetchall()[0][0]
        return gid

    def get_record(self, chromosome, position, name):
        """Retrieve details of all genes associated with a given
        chromosomal location.

        Args:
            chromosome (str):   Chromosme name in the form chrN where
                                N=1-23,X,Y
            position (integer): Location on the chromosome, using 0-based
                                indexing

        Returns:
            List of genes associated with the position. Each gene is an
            instance of a `Gene` namedtuple
        """
        if not chromosome.startswith('chr'):
            chromosome = 'chr' + chromosome
        stmt = 'SELECT * FROM refgene WHERE chrom=(?) AND ' + \
               '(?)>=txStart AND (?)<=txEnd AND name=(?) AND error = 0'
        results = self.execute(stmt, (chromosome, position, position, name))
        genes = []
        for r in results:
            exons = self._exons_from_gid(r[0])
            genes.append(Gene._make(r[0:13] + (exons,)))
        return genes


if __name__ == '__main__':
    import sys
    rg = Refgene()
    print rg.get_version()
    print rg.gene_by_location('chr14', 47120219)
    print rg.gene_by_name('APP')
    print rg.get_strand('NM_201414')
    print rg.get_cds_length('NM_001141945')

    region = ('X', 119504443, 119658996)
    print 'Following are the genes and their exons that '\
    'overlaps with %s -' % (str(region))
    for gene, exinfo in rg.get_genes(*region).items():
        print gene, exinfo
    #rg.get_exons(sys.argv[1], )
