"""
#
# COPYRIGHT (C) 2011 Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com
#
.. module:: Mutation Type Annotator
    :platform: Unix, Windows, MacOSX

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

    For a given SNP's wild type codon, mutant type codon and amino acid \
    position this module computes the mutation type which could be \
    one of the following -
    1.Synonymous Mutation
    2.Synonymous Stop Mutation
    3.Nonsynonymous Mutation
    4.Nonsynonymous Mutation (STOP LOSS)
    5.Nonsynonymous Mutation (START LOSS)
    6.Nonsense Mutation (STOP GAIN)
"""
from gcn.data.cut import CUT
from gcn.data.codontable import CODONTABLE
import re


class MutationAnnotation():

    def __init__(self, region, seq_annot):
        self.codon_idx = None
        if region == 'CodingExonic':
            self.wild_codon = seq_annot.wt_codon
            self.wild_aa = seq_annot.wt_aa
            self.mut_codon = seq_annot.mut_codon
            self.mut_aa = seq_annot.mut_aa
            self.aa_position = seq_annot.wt_aa_position
            self.codon_idx = seq_annot.codon_idx
        self.region = region
        self.ref = seq_annot.ref
        self.alt = seq_annot.alt
        self.posidx = seq_annot.posidx
        self.cds = seq_annot.cds
        self.codon_usage = None
        self.mut_type = None
        self.seq_annot = seq_annot

    def _get_insert_codons(self, seq):
        codons = []
        seq = seq.upper()
        if self.codon_idx == 1:
            codons = re.findall('...', seq[self.posidx] + self.alt + \
                                    seq[self.posidx + 1] + \
                                    seq[self.posidx + 2])
        elif self.codon_idx == 2:
            codons = re.findall('...', seq[self.posidx - 1] + \
                                    seq[self.posidx] + self.alt + \
                                    seq[self.posidx + 1])
        elif self.codon_idx == 3:
            codons = re.findall('...', self.alt)
        elif not self.codon_idx:
            if self.posidx == 0:
                f1 = seq[0] + self.alt + seq[1:3]
                f2 = self.alt + seq[1:3]
                f3 = self.alt[1:] + seq[1:3]
            elif self.posidx == len(seq) - 1:
                f1 = seq[self.posidx - 1] + seq[self.posidx] + self.alt
                f2 = seq[self.posidx] + self.alt
                f3 = self.alt
            elif self.posidx == len(seq) - 2:
                f1 = seq[self.posidx - 1] + seq[self.posidx] + self.alt + \
                                                        seq[self.posidx + 1]
                f2 = seq[self.posidx] + self.alt + seq[self.posidx + 1]
                f3 = self.alt + seq[self.posidx + 1]
            else:
                f1 = seq[self.posidx - 1] + seq[self.posidx] + self.alt + \
                                seq[self.posidx + 1] + seq[self.posidx + 2]
                f2 = seq[self.posidx] + self.alt + seq[self.posidx + 1] + \
                                                        seq[self.posidx + 2]
                f3 = self.alt + seq[self.posidx + 1] + seq[self.posidx + 2]
            for f in [f1, f2, f3]:
                codons += re.findall('...', f)
        return codons

    def mutation_type(self):
        if self.region[-8:] == 'boundary':
            reg1, reg2, _ = self.region.split('_', 2)
            if (reg1 == 'CodingExonic' and reg2 == 'UTR3')\
                or (reg1 == 'UTR3' and reg2 == 'CodingExonic'):
                if self.ref != '-':
                    self.mut_type = 'StopLoss'
            elif (reg1 == 'UTR5' and reg2 == 'CodingExonic')\
                or (reg1 == 'CodingExonic' and reg2 == 'UTR5'):
                if self.ref != '-':
                    self.mut_type = 'StartLoss'
                elif self.ref == '-':
                    utr5_seq = self.seq_annot.utr5_seq
                    codons = self._get_insert_codons(utr5_seq)
                    if 'ATG' in codons:
                        self.mut_type = 'StartGain'
        elif self.ref != '-' and self.alt != '-':
            if self.region == 'CodingExonic':
                mut_type = []
                codon_usage = []
                for wildaa, mutaa, aapos, mutcodon, wildcodon in zip(\
                            self.wild_aa, self.mut_aa, self.aa_position, \
                            self.mut_codon, self.wild_codon):
                    if wildaa == '*' and mutaa == '*':
                        mut_type.append('SynStop')
                        codon_usage.append(self.compute_codon_usage(\
                                                    wildcodon, mutcodon))
                    elif wildaa == '*' and mutaa != '*':
                        mut_type.append('StopLoss')
                    elif wildaa != '*' and mutaa == '*':
                        mut_type.append('StopGain')
                    elif wildaa != '*' and mutaa != '*':
                        codon_usage.append(self.compute_codon_usage(\
                                                        wildcodon, mutcodon))
                        if aapos == 1:
                            if mutcodon != 'ATG':
                                mut_type.append('StartLoss')
                            elif wildcodon != 'ATG' and \
                                        mutcodon == 'ATG':
                                mut_type.append('NonSynStart')
                        else:
                            if wildaa == mutaa:
                                mut_type.append('Syn')
                                codon_usage.append(self.compute_codon_usage(\
                                                        wildcodon, mutcodon))
                            elif wildaa != mutaa:
                                mut_type.append('NonSyn')
                temp = []
                for ele in mut_type:
                    if ele not in temp:
                        temp.append(ele)
                self.mut_type = '__'.join(temp)
                temp = []
                if codon_usage:
                    for ele in codon_usage:
                        if ele not in temp:
                            temp.append(ele)
                    self.codon_usage = '__'.join(temp)
            elif self.region == 'UTR5':
                self.utr5_codons = list(self.seq_annot.utr5_codons)
                if 'atg' in self.utr5_codons:
                    self.mut_type = 'StartGain'
        elif self.ref == '-':
            if self.region == 'CodingExonic':
                alt_length = len(self.alt)
                cds_len = len(self.cds)
                if self.posidx in [cds_len - 1, cds_len - 2, cds_len - 3]:
                    codons = self._get_insert_codons(self.cds)
                    if self.cds[-3:].upper() == codons[0]:
                        self.mut_type = 'NoCDSChange'
                    elif codons[0] in CODONTABLE['stop_codons']:
                        self.mut_type = 'SynStop'
                    elif codons[0] not in CODONTABLE['stop_codons']:
                        self.mut_type = 'StopLoss'
                else:
                    if divmod(alt_length, 3)[1] == 0:
                        self.mut_type = 'NonFrameShiftInsert'
                        codons = self._get_insert_codons(self.cds)
                        for stop_codons in CODONTABLE['stop_codons']:
                            if stop_codons in codons:
                                self.mut_type = 'StopGain'
                                break
                    else:
                        self.mut_type = 'FrameShiftInsert'
            elif self.region == 'UTR5':
                utr5_seq = self.seq_annot.utr5_seq
                codons = self._get_insert_codons(utr5_seq)
                if 'ATG' in codons:
                    self.mut_type = 'StartGain'
        elif self.alt == '-':
            if self.region == 'CodingExonic':
                ref_length = len(self.ref)
                if ref_length > 1:
                    endposidx = self.posidx + len(self.ref) - 1
                    end_codon, end_aa_position, end_aa, end_codon_idx = \
                        self.seq_annot.compute_codon_info(endposidx)
                else:
                    end_aa_position = self.aa_position
                    end_aa = self.wild_aa

                if divmod(ref_length, 3)[1] == 0:
                    if self.aa_position == 1 or end_aa_position == 1:
                        self.mut_type = 'StartLoss'
                    else:
                        if self.wild_aa == '*' or end_aa == '*':
                            self.mut_type = 'StopLoss'
                        else:
                            self.mut_type = 'NonFrameShiftDelete'
                else:
                    if self.aa_position == 1 or end_aa_position == 1:
                        self.mut_type = 'StartLoss'
                    else:
                        if self.wild_aa == '*' or end_aa == '*':
                            self.mut_type = 'StopLoss'
                        else:
                            self.mut_type = 'FrameShiftDelete'
        return self.mut_type

    def compute_codon_usage(self, wt_codon, mut_codon):
        wt_freq = CUT[wt_codon]
        mut_freq = CUT[mut_codon]
        if wt_freq < mut_freq:
            codon_usage = 'CodonUsageUp'
        elif wt_freq > mut_freq:
            codon_usage = 'CodonUsageDown'
        else:
            codon_usage = 'CodonUsageNoChange'
        return codon_usage

    def annotate(self):
        self.mutation_type()
