"""
.. module:: gene_pmrph
    :platform: Unix, Windows, MacOSX
    :synopsis: Computes the polymorphic score for all coding genes in human.

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules creates a tsv file comprising of polymorphic scores for all
coding genes in human.

The header of the output tsv file is -
['CHROM', 'GENE', 'EXONIC_VARIANTS_IN_GENE', 'GENE_LENGTH',
'POLYMORPHIC_SCORE']
"""
import argparse
from gcn.lib.databases.refgene import Refgene
from gcn.lib.io.vcf import VCFParser
from gcn.lib.varann.vartype.varant import varant_parser as vp
from collections import OrderedDict
from gcn.etc.fileconfig import FILECONFIG


def get_genelens():
    '''Computes gene's CDS length based on RefGene definition

    Returns:
        gene_lens(dictionary):    Returns the gene's CDS length
    '''
    gene_lens = {}
    gene_exons = {}
    gene_cds_se = {}
    rfg = Refgene()
    for rec in rfg.iterate_db():
        chrom = rec.chrom
        if 'chr' in chrom:
            chrom = chrom[3:]
        gene = rec.gene
        for exs, exe in zip(rec.str_exon_pos.split(','), \
                                            rec.end_exon_pos.split(',')):
            try:
                gene_exons[(chrom, gene)].add((int(exs), int(exe)))
                if int(rec.str_cds) < gene_cds_se[(chrom, gene)][0]:
                    gene_cds_se[(chrom, gene)][0] = int(rec.str_cds)
                if int(rec.end_cds) > gene_cds_se[(chrom, gene)][1]:
                    gene_cds_se[(chrom, gene)][1] = int(rec.end_cds)
            except:
                gene_exons[(chrom, gene)] = set([])
                gene_exons[(chrom, gene)].add((int(exs), int(exe)))
                gene_cds_se[(chrom, gene)] = [int(rec.str_cds),
                                              int(rec.end_cds)]

    for chrom_gene, exons in gene_exons.items():
        str_cds, end_cds = gene_cds_se[chrom_gene]
        cds_len = 0
        flag = False
        exons = list(exons)
        exons.sort()
        for spos, epos in exons:
            if spos <= str_cds <= epos and spos <= end_cds <= epos:
                cds_len += end_cds - str_cds
                break
            elif spos <= str_cds <= epos:
                flag = True
                cds_len += epos - str_cds
            elif spos <= end_cds <= epos:
                cds_len += end_cds - spos
                flag = False
                break
            elif flag == True:
                cds_len += epos - spos
        gene_lens[chrom_gene] = cds_len
    return gene_lens


def get_varcnt(invcf):
    '''Computes number of exonic variants per gene

    Args:
        invcf(str):    VARANT annotated VCF

    Returns:
        varcnt(dictionary):    Returns exonic variant count per gene
    '''
    varcnt = {}
    vcf = VCFParser(invcf)
    for rec in vcf:
        vcf.parseinfo(rec)
        ant = vp.parse(rec.info)
        prant = vp.prio_trans(ant)
        cache = []
        for altid, antinfo in prant.items():
            if altid != 'intergenic':
                genelist = antinfo.keys()
                for gene in genelist:
                    txant = antinfo[gene]['TRANSCRIPT']
                    key = (rec.chrom, gene)
                    if 'CodingExonic' in txant.region.split('_')\
                             and txant.mutation != 'Syn' \
                             and rec.info['ESPAF'] < 5.0 and key not in cache: #TODO (to be replaced by ExAC?)
                        cache.append(key)
                        if key not in varcnt:
                            varcnt[key] = 1
                        else:
                            varcnt[key] += 1
    return varcnt


def compute(invcf, outfile):
    '''Main method to compute the polymorphic score for human coding genes

    Args:
        invcf(str):    VARANT annotated VCF file
        outfile(str):    Output tsv file location.
                         Note - The genes in the file are listed in
                         descending order of their polymorphic score.
    '''
    varcnt = get_varcnt(invcf)
    genelens = get_genelens()
    gene_pmrp = {}
    out = open(outfile, 'w')
    out.write('#VCF File used to compute the polymorphic score is %s\n' \
                                                            % invcf)
    out.write('#Polymorphic score = Number of exonic variants in gene '\
                                                    '/ Gene\'s CDS length\n')
    out.write('\t'.join(['CHROM', 'GENE', 'EXONIC_VARIANTS_IN_GENE',
                         'GENE_LENGTH', 'POLYMORPHIC_SCORE']) + '\n')
    for key, vcnt in varcnt.items():
        if key in genelens:
            gene_pmrp[(key[0], key[1], str(vcnt), str(genelens[key]))] = \
                                    100 * (float(vcnt) / float(genelens[key]))
    sorted_data = OrderedDict(sorted(gene_pmrp.items(), key=lambda x: x[1], \
                                                            reverse=True))
    for key, val in sorted_data.items():
        od = list(key) + [str(val)]
        out.write('\t'.join(od) + '\n')

if __name__ == '__main__':
    outfile = FILECONFIG['POLYMORPH_GENES']
    desc = "Script creates a tsv file that comprise of polymorphic "\
    "scores for human coding genes. Note - The genes in the file are "\
    "listed in descending order of their polymorphic score."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--inputVCF', dest='invcf', type=str,
                        help='Path to input VARANT annotated VCF file')
    parser.add_argument('-o', '--outputVCF', dest='outtsv', type=str,
                        default=outfile, help='Path to output tsv file')
    args = parser.parse_args()
    compute(args.invcf, args.outtsv)
