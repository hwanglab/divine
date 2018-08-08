"""
.. module:: annotateSequence
    :platform: Unix, Windows, MacOSX
    :synopsis: Sequence Annotation

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

For a given SNP, this modlues derives annotations like -
if position is in 'Exon', it computes for gene name, cds,
exon sequences, wild codon, wild amino acid, mutant codon,
mutant amino acid, amino acid position and protein sequence.

If position is in '5'UTR', it computes for 5'utr sequence
"""

from gcn.data.codontable import CODONTABLE
from gcn.lib.databases.refmrna import Refmrna
from gcn.data.pseudoautosomal_genes import PSEUDO_AUTO_GENES


class SequenceAnnotation():

    def __init__(self, region, refgene_entry, spos, epos, ref, alt, warning):
        self.region = region
        self.warning = warning
        self.refgene_entry = refgene_entry
        self.ref = ref
        self.alt = alt
        self.strand = refgene_entry.strand
        self.poslist = range(spos, epos + 1)
        if self.strand == '+':
            self.position = self.poslist[0] - 1
        else:
            self.position = self.poslist[-1]
        self.refseq_acc = refgene_entry.refseq_acc
        self.seq = None
        self.initiate_variables()
        self.cds = None
        self.utr5_seq = None
        self.codon_idx = None
        self.wt_codon = None
        self.wt_aa_position = None
        self.wt_aa = None
        self.mut_codon = None
        self.mut_aa_position = None
        self.mut_aa = None
        self.posidx = None
        self.utr5_codons = set()
        self.exon_number = 0
        self.idx = 0

    def load_file(self):
        mrna_info = Refmrna().retrieve(self.refseq_acc, self.refgene_entry.gid)
        return mrna_info

    def initiate_variables(self):
        mrna_info = self.load_file()
        if mrna_info != None:
            self.idx = mrna_info.idx
            self.seq = mrna_info.sequence
            self.str_cds = int(self.refgene_entry.str_cds)
            self.end_cds = int(self.refgene_entry.end_cds)
            self.exon_str = [int(e) for e in \
                             self.refgene_entry.str_exon_pos.split(',')]
            self.exon_end = [int(e) for e in \
                             self.refgene_entry.end_exon_pos.split(',')]
            self.exon_lengths = [x - y for (x, y) in \
                                  zip(self.exon_end, self.exon_str)]
            self.frames = [int(e) for e in \
                           self.refgene_entry.frames.split(',')]

    def get_cds_exonnum(self, frames):
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

    def get_cds_posidx(self, position):
        str_cds_exon_number, end_cds_exon_number = \
                                self.get_cds_exonnum(self.frames)
        exon_starts = self.exon_str[str_cds_exon_number:\
                                        end_cds_exon_number + 1]
        exon_ends = self.exon_end[str_cds_exon_number:\
                                        end_cds_exon_number + 1]
        l = len(exon_ends)
        exon_ends[l - 1] = self.end_cds
        exon_starts[0] = self.str_cds
        idx = 0
        if self.strand == '+':
            for x, y in zip(exon_starts, exon_ends):
                if x <= position < y:
                    idx += position - x
                    break
                else:
                    idx += y - x
        elif self.strand == '-':
            exon_starts.reverse()
            exon_ends.reverse()
            for x, y in zip(exon_starts, exon_ends):
                if y >= position > x:
                    idx += y - position
                    break
                else:
                    idx += y - x
        return idx

    def compute_cdsseq(self):
        if self.strand == '+':
            str_cds_exon_number, end_cds_exon_number = \
                                        self.get_cds_exonnum(self.frames)
            str_mrna_idx = sum(self.exon_lengths[:str_cds_exon_number])
            end_mrna_idx = sum(self.exon_lengths[:end_cds_exon_number + 1])
            cds_exon_seq = self.seq[str_mrna_idx:end_mrna_idx]
            start_offset = self.str_cds - self.exon_str[str_cds_exon_number]
            end_offset = len(cds_exon_seq) - \
                            (self.exon_end[end_cds_exon_number] - self.end_cds)
        elif self.strand == '-':
            self.exon_str.reverse()
            self.exon_end.reverse()
            self.exon_lengths.reverse()
            self.frames.reverse()
            str_cds_exon_number, end_cds_exon_number = \
                                        self.get_cds_exonnum(self.frames)
            str_mrna_idx = sum(self.exon_lengths[:str_cds_exon_number])
            end_mrna_idx = sum(self.exon_lengths[:end_cds_exon_number + 1])
            cds_exon_seq = self.seq[str_mrna_idx:end_mrna_idx]
            start_offset = self.exon_end[str_cds_exon_number] - self.end_cds
            end_offset = len(cds_exon_seq) - (self.str_cds - \
                                        self.exon_str[end_cds_exon_number])
            self.exon_str.reverse()
            self.exon_end.reverse()
            self.exon_lengths.reverse()
            self.frames.reverse()
        self.cds = cds_exon_seq[start_offset:end_offset]

    def compute_utr5seq(self):
        if self.refgene_entry.strand == '+':
            str_cds_exon_number, end_cds_exon_number = \
                                self.get_cds_exonnum(self.frames)
            idx = 0
            for exon_num in range(len(self.exon_lengths)):
                if exon_num < str_cds_exon_number:
                    idx += self.exon_lengths[exon_num]
                elif exon_num == str_cds_exon_number:
                    idx += self.str_cds - self.exon_str[exon_num]
                    break
        elif self.refgene_entry.strand == '-':
            self.exon_str.reverse()
            self.exon_end.reverse()
            self.exon_lengths.reverse()
            self.frames.reverse()
            str_cds_exon_number, end_cds_exon_number = \
                                self.get_cds_exonnum(self.frames)
            idx = 0
            for exon_num in range(len(self.exon_lengths)):
                if exon_num < str_cds_exon_number:
                    idx += self.exon_lengths[exon_num]
                elif exon_num == str_cds_exon_number:
                    idx += self.exon_end[exon_num] - self.end_cds
                    break
            self.exon_str.reverse()
            self.exon_end.reverse()
            self.exon_lengths.reverse()
            self.frames.reverse()
        self.utr5_seq = self.seq[:idx]

    def mutate_cds(self, poslist, alleles):
        self.mut_cds = None
        seq = list(self.cds)
        for pos, allele in zip(poslist, alleles):
            seq[pos] = allele
        self.mut_cds = "".join(seq).lower()

    def mutate_5utr(self, poslist, alleles):
        self.mut_5utr = None
        seq = list(self.utr5_seq)
        for pos, allele in zip(poslist, alleles):
            seq[pos] = allele
        self.mut_5utr = "".join(seq).lower()

    def compute_codon_info(self, position_idx, mut_allele=None):
        if mut_allele != None:
            seq = self.mut_cds
        elif mut_allele == None:
            seq = self.cds
        codon_idx = None
        codon = ''
        aa_position = ''
        aa = ''
        position_idx += 1
        if position_idx != 1:
            # check for last base of the codon
            if divmod(position_idx, 3)[1] == 0:
                codon = seq[position_idx - 3:position_idx].upper()
                aa_position = position_idx / 3
                aa = CODONTABLE['codons'][codon]
                codon_idx = 3
            # check for middle base of the codon
            elif divmod(position_idx + 1, 3)[1] == 0:
                codon = seq[position_idx - 2:position_idx + 1].upper()
                aa_position = (position_idx + 1) / 3
                aa = CODONTABLE['codons'][codon]
                codon_idx = 2
            # check for first base of the codon
            elif divmod(position_idx - 1, 3)[1] == 0:
                codon = seq[position_idx - 1:position_idx + 2].upper()
                aa_position = (position_idx + 2) / 3
                aa = CODONTABLE['codons'][codon]
                codon_idx = 1
        else:
            codon = (seq[:3]).upper()
            aa_position = 1
            aa = CODONTABLE['codons'][codon]
            codon_idx = 1
        return codon, aa_position, aa, codon_idx

    def get_5utr_posidx(self, position):
        pos_idx = 0
        if self.strand == '+':
            for x, y in zip(self.exon_end, self.exon_str):
                if y <= position < x:
                    pos_idx += position - y
                    break
                else:
                    pos_idx += x - y
        elif self.refgene_entry.strand == '-':
            self.exon_str.reverse()
            self.exon_end.reverse()
            for x, y in zip(self.exon_end, self.exon_str):
                if x >= position > y:
                    pos_idx += x - position
                    break
                else:
                    pos_idx += x - y
            self.exon_str.reverse()
            self.exon_end.reverse()
        return pos_idx

    def compute_utr_codon(self, poslist):
        for pos in poslist:
            for codon in [self.mut_5utr[pos:pos + 3], \
                      self.mut_5utr[pos - 1:pos + 2], \
                      self.mut_5utr[pos - 2:pos + 1]]:
                if codon and len(codon) == 3:
                    self.utr5_codons.add(codon)

    def annotate(self):
        if self.seq == None:
            self.warning = 'mRNA_SEQUENCE_NOT_FOUND'
            return False, self.warning
        elif list(set(self.seq))[0].upper() in ['N', ' ']:
            if self.refgene_entry.gene in PSEUDO_AUTO_GENES:
                self.warning = 'GENE_IS_IN_PSEUDOAUTOSOMAL_REGION'
            else:
                self.warning = \
                'mRNA_SEQUENCE_DOES_NOT_COMPRISE_OF_NEUCLEOTIDE_BASES'
            return False, self.warning
        elif self.alt == 'N':
            self.warning = 'ALT_ALLELE_IS_NOT_BASE'
            return False, self.warning
        else:
            if self.region == 'CodingExonic':
                self.compute_cdsseq()
                if '-' != self.ref and '-' != self.alt:
                    mutposlist = [(idx, ele[1]) \
                                       for idx, ele in enumerate([(x == y, y) \
                                        for x, y in zip(self.ref, self.alt)]) \
                                       if not ele[0]]
                    posidxlist = []
                    if self.strand == '+':
                        for ele in mutposlist:
                            posidxlist.append((self.get_cds_posidx(\
                                                self.position + ele[0]), \
                                                ele[1]))
                    elif self.strand == '-':
                        for ele in mutposlist:
                            posidxlist.append((self.get_cds_posidx(\
                                                self.position - ele[0]), \
                                                ele[1]))

                    if self.cds[posidxlist[0][0]:posidxlist[-1][0] \
                                + 1].upper() != self.ref[\
                            mutposlist[0][0]:mutposlist[-1][0] + 1]:
                            self.warning = 'MISMATCH_IN_REF_ALLELE'

                    altlist = [ele[1] for ele in mutposlist]

                    self.mutate_cds([e[0] for e in posidxlist], altlist)

                    for posidx, aale in posidxlist:
                        wt_codon, wt_aa_pos, wt_aa, codon_idx = \
                                            self.compute_codon_info(posidx)
                        mut_codon, mut_aa_pos, mut_aa, codon_idx = \
                                            self.compute_codon_info(\
                                                        posidx, aale)
                        try:
                            if wt_aa_pos not in self.wt_aa_position:
                                self.wt_codon.append(wt_codon)
                                self.wt_aa_position.append(wt_aa_pos)
                                self.wt_aa.append(wt_aa)
                                self.codon_idx.append(codon_idx)
                                self.mut_codon.append(mut_codon)
                                self.mut_aa_position.append(mut_aa_pos)
                                self.mut_aa.append(mut_aa)
                        except:
                            self.wt_codon = [wt_codon]
                            self.wt_aa_position = [wt_aa_pos]
                            self.wt_aa = [wt_aa]
                            self.codon_idx = [codon_idx]
                            self.mut_codon = [mut_codon]
                            self.mut_aa_position = [mut_aa_pos]
                            self.mut_aa = [mut_aa]
                elif '-' == self.ref:
                    self.posidx = self.get_cds_posidx(self.position)
                    self.wt_codon, self.wt_aa_position, self.wt_aa, \
                            self.codon_idx = \
                            self.compute_codon_info(self.posidx)
                    self.mut_codon = 'Insert'
                    self.mut_aa = 'Insert'
                elif '-' == self.alt:
                    self.posidx = self.get_cds_posidx(self.position)
                    computed_ref = self.cds[self.posidx:self.posidx + \
                                                len(self.ref)].upper()
                    if computed_ref != self.ref:
                        self.warning = 'MISMATCH_IN_REF_ALLELE'
                    self.wt_codon, self.wt_aa_position, self.wt_aa, \
                            self.codon_idx = \
                            self.compute_codon_info(self.posidx)
                    self.mut_codon = 'Delete'
                    self.mut_aa = 'Delete'
            elif self.region == 'UTR5':
                if '-' != self.ref and '-' != self.alt:
                    self.compute_utr5seq()
                    mutposlist = [(idx, ele[1]) \
                                       for idx, ele in enumerate([(x == y, y) \
                                        for x, y in zip(self.ref, self.alt)]) \
                                       if not ele[0]]
                    posidxlist = []
                    if self.strand == '+':
                        for ele in mutposlist:
                            posidxlist.append((self.get_5utr_posidx(\
                                                self.position + ele[0]), \
                                                ele[1]))
                    elif self.strand == '-':
                        for ele in mutposlist:
                            posidxlist.append((self.get_5utr_posidx(\
                                                self.position - ele[0]), \
                                                ele[1]))

                    if self.utr5_seq[posidxlist[0][0]: posidxlist[-1][0] \
                                + 1].upper() != self.ref[\
                            mutposlist[0][0]:mutposlist[-1][0] + 1]:
                            self.warning = 'MISMATCH_IN_REF_ALLELE'

                    altlist = [ele[1] for ele in mutposlist]
                    self.mutate_5utr([e[0] for e in posidxlist], altlist)
                    self.compute_utr_codon([e[0] for e in posidxlist])
                elif self.ref == '-':
                    self.compute_utr5seq()
                    self.posidx = self.get_5utr_posidx(self.position)
            elif 'boundary' in self.region and 'UTR5' in self.region:
                    if self.ref == '-':
                        self.compute_utr5seq()
                        self.posidx = self.get_5utr_posidx(self.position)
            return True, self.warning
