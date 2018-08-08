"""
.. module:: genemap
    :platform: Unix, Windows, MacOSX
    :synopsis: Maps Gene coordinate to its translated protein position

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules maps the gene coordinate to its translates protein position
and a returns a dictionary of format -
map_dict[Refseq_acc]={aa_pos: [aa, codon, cds_numbers, chrom, coordinates, \
                                strand]}
"""

from gcn.lib.databases.refmrna import Refmrna
from gcn.lib.databases.refgene import Refgene
from gcn.data.codontable import CODONTABLE
from gcn.lib.databases.splicedb import Splicedb


class GeneMap():

    def __init__(self, gene):
        self.gene = gene
        self.refgene = Refgene()
        self.refmrna = Refmrna()
        self.genelist = self.load_refgene(self.gene)
        self.splicedb = Splicedb()

    def load_refgene(self, gene):
        '''Loads the refgene data for given Gene'''
        genelist = []
        for record in self.refgene.iterate_db(gene=self.gene):
            genelist.append(record)
        return genelist

    def get_cds_exonnum(self, frames):
        '''Returns the Start coding exon number and End Coding exon number'''
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

    def translate(self, mrna):
        protein = ""
        codon = ""
        for neu in mrna[:-3]:
            codon += neu
            if len(codon) == 3:
                aa = CODONTABLE['codons'][codon.upper()]
                protein += aa
                codon = ""
        return protein

    def _expand_exon(self, entry):
        '''Internal Method'''
        chrom = entry.chrom
        strand = entry.strand
        refseq_acc = entry.refseq_acc
        seq = self.refmrna.retrieve(entry.refseq_acc, entry.gid).sequence
        str_cds = int(entry.str_cds)
        end_cds = int(entry.end_cds)
        exon_str = [int(e) for e in entry.str_exon_pos.split(',')]
        exon_end = [int(e) for e in entry.end_exon_pos.split(',')]
        frames = [int(e) for e in entry.frames.split(',')]
        cds_flag = False
        cds_cnt = 0
        neu_cnt = 0
        cds_seq = ""
        gene_info = []
        if strand == '+':
            exon_cnt = 0
            for start, end, frame in zip(exon_str, exon_end, frames):
                exon_cnt += 1
                if frame != -1:
                    for i in range(start, end):
                        if i == str_cds:
                            cds_flag = True
                        elif i == end_cds:
                            cds_flag = False
                        if cds_flag == True:
                            splice_reg_site = self.splicedb.get_annot(chrom, \
                                                                i + 1, \
                                                                refseq_acc)
                            if not splice_reg_site:
                                splice_reg_site = ''
                            cds_cnt += 1
                            gene_info.append([cds_cnt, chrom, i + 1, strand,
                                              str(seq[neu_cnt]),
                                              'E-' + str(exon_cnt),
                                              splice_reg_site])
                            cds_seq += seq[neu_cnt]
                        neu_cnt += 1
                else:
                    neu_cnt += end - start
        elif strand == '-':
            exon_cnt = len(exon_str) + 1
            exon_str.reverse()
            exon_end.reverse()
            frames.reverse()
            for start, end, frame in zip(exon_str, exon_end, frames):
                exon_cnt -= 1
                if frame != -1:
                    for i in reversed(range(start + 1, end + 1)):
                        if i == end_cds:
                            cds_flag = True
                        elif i == str_cds:
                            cds_flag = False
                        if cds_flag == True:
                            splice_reg_site = self.splicedb.get_annot(chrom, \
                                                                      i + 1, \
                                                                refseq_acc)
                            if not splice_reg_site:
                                splice_reg_site = ''
                            cds_cnt += 1
                            gene_info.append([cds_cnt, chrom, i, strand, \
                                              str(seq[neu_cnt]), \
                                              'E' + str(exon_cnt),
                                              splice_reg_site])
                            cds_seq += seq[neu_cnt]
                        neu_cnt += 1
                else:
                    neu_cnt += end - start
        protein_seq = self.translate(cds_seq)
        return gene_info, protein_seq

    def mapcoord(self):
        '''Maps the gene coordinates to protein position and
        return a dictionary of format - map_dict[Refseq_acc]={aa_pos: [aa, \
        codon, cds_numbers, Exon_numers, chrom, coordinates, strand]}'''
        refgene_entries = self.genelist
        map_dict = {}
        for entry in refgene_entries:
            acc = entry.refseq_acc
            map_dict[acc] = {}
            gene_info, protein_seq = self._expand_exon(entry)
            aa_pos = 0
            for idx in range(0, len(gene_info), 3):
                if idx == len(gene_info) - 3:
                    aa = 'STOP_CODON'
                else:
                    aa = protein_seq[aa_pos]
                aa_pos += 1
                d = zip(gene_info[idx], gene_info[idx + 1], gene_info[idx + 2])
                map_dict[acc][aa_pos] = [aa, ''.join(list(d[4])), list(d[0]), \
                                         list(d[5]), d[1][0], list(d[2]), \
                                         d[3][0], list(d[6])]
        return map_dict


def main(gene):
    header = ['Refseq_acc', 'AA_pos', 'AA', 'Codon', 'CSD_numbers', \
              'Exon_number', 'Chrom', 'Coordinates', 'Strand', \
              'Exon_Splice_Regulation_Site']
    print '\t'.join(header)
    genec = GeneMap(gene)
    map_dict = genec.mapcoord()
    for acc, val in map_dict.iteritems():
        for aa_pos, map_data in val.iteritems():
            d = [acc, str(aa_pos), map_data[0], map_data[1], \
                 ','.join([str(e) for e in map_data[2]]), \
                 ','.join([str(e) for e in map_data[3]]), \
                 map_data[4], ','.join([str(e) for e in \
                                        map_data[5]]), \
                 map_data[6], ','.join(map_data[7])]
            print '\t'.join(d)


if __name__ == '__main__':
    import sys
    gene = sys.argv[1] # Gene Name like CFTR etc
    main(gene)
