"""
.. module:: varant_parser
    :platform: Unix, Windows, MacOSX
    :synopsis: Module to parse the output of Varant annotation.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

This module parses VARANT annotation present in an annotated VCF file.
"""
from gcn.lib.io.vcf import VCFParser
from collections import namedtuple
from collections import OrderedDict

Trans_ant = namedtuple('Trans_ant', ['trans_id', 'region', 'exon', 'cdna',
                                     'splice', 'utr_signal', 'mutation',
                                     'codon', 'aa', 'protein_len', 'usage',
                                     'sift', 'pp2', 'warning'])

MUT_RANKS = {'StopGain': 1, 'StartLoss': 1, 'StopLoss': 1, 'NonSyn': 3,
             'NonSynStart': 3, 'SynStop': 3, 'FrameShiftDelete': 2,
             'FrameShiftInsert': 2, 'NonFrameShiftDelete': 4,
             'NonFrameShiftInsert': 4, 'Syn': 5, 'StartGain': 5,
             'NoCDSChange': 6, '': 6}


def get_rank(mutations):
    """Given mutations type annotations returns effect rank with lower
    rank being more deleterious and higher rank less deleterious.

    Args:
        mutations(str):    Mutation Type.

    Returns:
        rank(int):    Mutation Rank
    """
    rank = []
    for mut in mutations.split('__'):
        rank.append(MUT_RANKS[mut])
    rank.sort()
    return rank[0]


def parse(info):
    """This method parses the Varant annotation and return dictionary

    Args:
        info(dictionary): Parsed INFO dictionary created by VCF parser
                           i.e. rec['info']

    Returns:
        parsed_dict(dictionary): Parsed VARANT annotation dictionary.
                                 Format of this is -
                                 {1: {Gene-1: {TRANSCRIPTS: [nt1, nt2, ],
                                            MIM_PHENS: [phen1, phen2, ],
                                            MIM_IDS: [mimid1, mimid2, ],
                                            GAD_PHENS: [phen1, phen2, ]
                                             },
                                  Gene-2: {TRANSCRIPTS: [nt1, nt2, ],
                                            MIM_PHENS: [phen1, phen2, ],
                                            MIM_IDS: [mimid1, mimid2, ],
                                            GAD_PHENS: [phen1, phen2, ]
                                             },
                                      }
                                  2: {Gene-1: {TRANSCRIPTS: [nt1, ]
                                               MIM_PHENS: [phen1, ]
                                               MIM_IDS: [mimid1, ]
                                               GAD_PHENS: [phen1, ]
                                              }
                                      }
                                 }
                                 where 1=first alt allele, 2=second alt
                                 allele and so on..
                                 nt1, nt2 etc are named tuple holding the
                                 transcipt specific annotations for each
                                 transcript and is of following format -
                                 ['trans_id', 'region', 'exon', 'splice',
                                   'utr_signal', 'mutation', 'codon',
                                   'aa', 'usage', 'sift',
                                   'pp2', 'warning']
    """
    parsed_dict = OrderedDict()
    if 'VARANT_GENIC' in info:
        for ant in info['VARANT_GENIC']:
            gene = ant.split('(')[0]
            v = ant.split('(')[1].split(')')[0].split(':')
            mim_phen, mim_ids, gad_phen = v[-3:]
            al = []
            for e in v[:-3]:
                d = e.split('|')
                alt_cnt = d.pop(3)
                if int(alt_cnt) not in parsed_dict:
                    parsed_dict[int(alt_cnt)] = OrderedDict()
                    parsed_dict[int(alt_cnt)][gene] = \
                                        {'TRANSCRIPTS': [Trans_ant._make(d)],
                                         'MIM_PHENS': mim_phen.split('__'),
                                         'MIM_IDS': mim_ids.split('__'),
                                         'GAD_PHENS': gad_phen.split('__')
                                        }
                else:
                    if gene in parsed_dict[int(alt_cnt)]:
                        parsed_dict[int(alt_cnt)][gene]['TRANSCRIPTS']\
                                                    .append(Trans_ant._make(d))
                    else:
                        parsed_dict[int(alt_cnt)][gene] = \
                                        {'TRANSCRIPTS': [Trans_ant._make(d)],
                                         'MIM_PHENS': mim_phen.split('__'),
                                         'MIM_IDS': mim_ids.split('__'),
                                         'GAD_PHENS': gad_phen.split('__')
                                        }
    if 'VARANT_INTERGENIC' in info:
        y = info['VARANT_INTERGENIC']
        parsed_dict['intergenic'] = {y: ''}
    return parsed_dict


def prio_trans(parsed_dict):
    """This method prioritizes transcripts for a gene based on its
    annotations.

    Args:
        parsed_dict(dictionary): Parsed VARANT annotation dictionary.
                                 Format is -
                                 {1: {Gene-1: {TRANSCRIPTS: [nt1, nt2, ],
                                            MIM_PHENS: [phen1, phen2, ],
                                            MIM_IDS: [mimid1, mimid2, ],
                                            GAD_PHENS: [phen1, phen2, ]
                                             },
                                  Gene-2: {TRANSCRIPTS: [nt1, nt2, ],
                                            MIM_PHENS: [phen1, phen2, ],
                                            MIM_IDS: [mimid1, mimid2, ],
                                            GAD_PHENS: [phen1, phen2, ]
                                             },
                                      }
                                  2: {Gene-1: {TRANSCRIPTS: [nt1, ],
                                               MIM_PHENS: [phen1, ],
                                               MIM_IDS: [mimid1, ],
                                               GAD_PHENS: [phen1, ]
                                              }
                                      }
                                 }
                                 where 1=first alt allele, 2=second alt
                                 allele and so on..
                                 where nt1, nt2 etc are named tuple holding the
                                 transcipt specific annotations for
                                 each transcript and is of following format -
                                 ['trans_id', 'region', 'exon', 'splice',
                                   'utr_signal', 'mutation', 'codon',
                                   'aa', 'usage', 'sift',
                                   'pp2', 'warning']

    Returns:
        prior_dict(dictionary): Transcript prioritized dictionary.
                                Format is -
                                {1: {Gene-1: {TRANSCRIPT: nt,
                                            MIM_PHENS: [phen1, phen2, ],
                                            MIM_IDS: [mimid1, mimid2, ],
                                            GAD_PHENS: [phen1, phen2, ]
                                             },
                                  Gene-2: {TRANSCRIPT: nt,
                                            MIM_PHENS: [phen1, phen2, ],
                                            MIM_IDS: [mimid1, mimid2, ],
                                            GAD_PHENS: [phen1, phen2, ]
                                             },
                                      }
                                  2: {Gene-1: {TRANSCRIPT: nt,
                                               MIM_PHENS: [phen1, ],
                                               MIM_IDS: [mimid1, ],
                                               GAD_PHENS: [phen1, ]
                                              }
                                      }
                                 }
                                 where 1=first alt allele, 2=second alt
                                 allele and so on..
                                 and nt is the prioritized transcript
                                 based on its annotation.
    """
    prior_dict = OrderedDict()
    for alt_cnt, annot in parsed_dict.items():
        prior_dict[alt_cnt] = OrderedDict()
        for gene, val in annot.items():
            if val:
                t = {}
                for e in val['TRANSCRIPTS']:
                    if e.region not in t:
                        t[e.region] = [e]
                    else:
                        t[e.region].append(e)

                if len([r for r in t.keys() if 'boundary' in r]) > 0:
                    for region in [r for r in t.keys() if 'boundary' in r]:
                        for s in t[region]:
                            if gene not in prior_dict[alt_cnt]:
                                prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                                mut_rank = get_rank(s.mutation)
                            else:
                                cur_mut_rank = get_rank(s.mutation)
                                if cur_mut_rank < mut_rank or s.utr_signal \
                                                                or s.splice:
                                    prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                                    mut_rank = cur_mut_rank
                elif 'CodingExonic' in t.keys():
                    for s in t['CodingExonic']:
                        if gene not in prior_dict[alt_cnt]:
                            prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                            mut_rank = get_rank(s.mutation)
                        else:
                            cur_mut_rank = get_rank(s.mutation)
                            if cur_mut_rank < mut_rank:
                                prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                                mut_rank = cur_mut_rank
                elif 'UTR5' in t.keys():
                    for s in t['UTR5']:
                        funct_site = s.utr_signal
                        if gene not in prior_dict[alt_cnt] or funct_site:
                            prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                elif 'UTR3' in t.keys():
                    for s in t['UTR3']:
                        funct_site = s.utr_signal
                        if gene not in prior_dict[alt_cnt] or funct_site:
                            prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                elif 'CodingIntronic' in t.keys():
                    for s in t['CodingIntronic']:
                        if gene not in prior_dict[alt_cnt]:
                            prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                            if '_' in s.splice:
                                sd = s.splice.split('_')[2]
                                sa = s.splice.split('_')[-1]
                            else:
                                sd, sa = [0, 0]
                        else:
                            if '_' in s.splice:
                                csd = int(s.splice.split('_')[2])
                                csa = int(s.splice.split('_')[-1])
                            else:
                                csd, csa = [0, 0]
                            if csd < sd or csa < sa:
                                prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                                sd, sa = [csd, csa]
                elif 'NonCodingExonic' in t.keys():
                    for s in t['NonCodingExonic']:
                        if gene not in prior_dict[alt_cnt] or s.splice:
                            prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                elif 'NonCodingIntronic' in t.keys():
                    for s in t['NonCodingIntronic']:
                        if gene not in prior_dict[alt_cnt]:
                            prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                            if '_' in s.splice:
                                sd = s.splice.split('_')[2]
                                sa = s.splice.split('_')[-1]
                            else:
                                sd, sa = [0, 0]
                        else:
                            if '_' in s.splice:
                                csd = int(s.splice.split('_')[2])
                                csa = int(s.splice.split('_')[-1])
                            else:
                                csd, csa = [0, 0]
                            if csd < sd or csa < sa:
                                prior_dict[alt_cnt][gene] = {'TRANSCRIPT': s}
                                sd, sa = [csd, csa]
                prior_dict[alt_cnt][gene]['MIM_PHENS'] = \
                                        parsed_dict[alt_cnt][gene]['MIM_PHENS']
                prior_dict[alt_cnt][gene]['MIM_IDS'] = \
                                        parsed_dict[alt_cnt][gene]['MIM_IDS']
                prior_dict[alt_cnt][gene]['GAD_PHENS'] = \
                                        parsed_dict[alt_cnt][gene]['GAD_PHENS']
            else:
                prior_dict[alt_cnt] = {gene: ''}
    return prior_dict


def get_prior_geneannot(info, alltrans=False):
    """Returns a prioritized gene with a prioritized transcript
    Args:
        info(dictionary):    Parsed INFO dictionary created by VCF parser
                             i.e. rec['info']
        alltrans(Boolean):    Boolean value where True would return all
                              the transcripts of prioritized gene and
                              False would return prioritized transcript
                              of the prioritized gene.
    Returns:
        geneannot(dictionary):    Parsed Gene annotation.
                                  Format of the dictionary is as -
                                  geneannot = {'GENE': GeneName
                                                'TRANSCRIPT': [nt],
                                                'MIM_PHENS': [phen1, ],
                                                'MIM_IDS': [id1,],
                                                'GAD_PHENS': [phen1, ]
                                                }
                                  nt is a named tuple holding the
                                  transcipt specific annotations for each
                                  transcript and is of following format -
                                  ['trans_id', 'region', 'exon', 'splice',
                                   'utr_signal', 'mutation', 'codon',
                                   'aa', 'usage', 'sift',
                                   'pp2', 'warning']
    """
    geneannot = {'GENE': ''}
    parsed_dict = parse(info)
    prior_trans_dict = prio_trans(parsed_dict)
    temp = {}
    for altid, geneinfo in prior_trans_dict.items():
        if altid != 'intergenic':
            for gene, ants in geneinfo.items():
                region = ants['TRANSCRIPT'].region
                if region not in temp:
                    temp[region] = [(gene, ants)]
                else:
                    temp[region].append((gene, ants))
        else:
            temp['Intergenic'] = (geneinfo.keys()[0], None)

    if 'CodingExonic' in temp:
        prev_mut_rank = None
        for gene, ants in temp['CodingExonic']:
            cur_mut_rank = get_rank(ants['TRANSCRIPT'].mutation)
            if cur_mut_rank < prev_mut_rank or not prev_mut_rank:
                prev_mut_rank = cur_mut_rank
                geneannot = {'GENE': gene}
                geneannot.update(ants)
    elif len([r for r in temp.keys() if 'boundary' in r]) > 0:
        for region in [r for r in temp.keys() if 'boundary' in r]:
            prev_mut_rank = None
            for gene, ants in temp[region]:
                cur_mut_rank = get_rank(ants['TRANSCRIPT'].mutation)
                if cur_mut_rank < prev_mut_rank or not prev_mut_rank or \
                                            ants['TRANSCRIPT'].utr_signal or \
                                            ants['TRANSCRIPT'].splice:
                    prev_mut_rank = cur_mut_rank
                    geneannot = {'GENE': gene}
                    geneannot.update(ants)
    elif 'CodingIntronic' in temp.keys():
        for gene, ants in temp['CodingIntronic']:
            if '_' not in ants['TRANSCRIPT'].splice:
                geneannot = {'GENE': gene}
                geneannot.update(ants)
                break
            else:
                if not geneannot['GENE']:
                    geneannot = {'GENE': gene}
                    geneannot.update(ants)
                    if '_' in ants['TRANSCRIPT'].splice:
                        sd = int(ants['TRANSCRIPT'].splice.split('_')[2])
                        sa = int(ants['TRANSCRIPT'].splice.split('_')[-1])
                    else:
                        sd, sa = [0, 0]
                else:
                    if '_' in ants['TRANSCRIPT'].splice:
                        csd = int(ants['TRANSCRIPT'].splice.split('_')[2])
                        csa = int(ants['TRANSCRIPT'].splice.split('_')[-1])
                    else:
                        csd, csa = [0, 0]
                    if csd < sd or csa < sa:
                        geneannot = {'GENE': gene}
                        geneannot.update(ants)
                        sd, sa = [csd, csa]
    elif 'UTR5' in temp:
        prev_mut_rank = None
        for gene, ants in temp['UTR5']:
            cur_mut_rank = get_rank(ants['TRANSCRIPT'].mutation)
            if cur_mut_rank < prev_mut_rank or not prev_mut_rank or \
                                ants['TRANSCRIPT'].utr_signal:
                prev_mut_rank = cur_mut_rank
                geneannot = {'GENE': gene}
                geneannot.update(ants)
    elif 'UTR3' in temp:
        for gene, ants in temp['UTR3']:
            if not geneannot['GENE'] or ants['TRANSCRIPT'].utr_signal:
                geneannot = {'GENE': gene}
                geneannot.update(ants)
    elif 'NonCodingExonic' in temp:
        for gene, ants in temp['NonCodingExonic']:
            if not geneannot['GENE'] or ants['TRANSCRIPT'].splice:
                geneannot = {'GENE': gene}
                geneannot.update(ants)
    elif 'NonCodingIntronic' in temp.keys():
        for gene, ants in temp['NonCodingIntronic']:
            if '_' not in ants['TRANSCRIPT'].splice:
                geneannot = {'GENE': gene}
                geneannot.update(ants)
                break
            else:
                if not geneannot['GENE']:
                    geneannot = {'GENE': gene}
                    geneannot.update(ants)
                    if '_' in ants['TRANSCRIPT'].splice:
                        sd = int(ants['TRANSCRIPT'].splice.split('_')[2])
                        sa = int(ants['TRANSCRIPT'].splice.split('_')[-1])
                    else:
                        sd, sa = [0, 0]
                else:
                    if '_' in ants['TRANSCRIPT'].splice:
                        csd = int(ants['TRANSCRIPT'].splice.split('_')[2])
                        csa = int(ants['TRANSCRIPT'].splice.split('_')[-1])
                    else:
                        csd, csa = [0, 0]
                    if csd < sd or csa < sa:
                        geneannot = {'GENE': gene}
                        geneannot.update(ants)
                        sd, sa = [csd, csa]
    elif 'Intergenic' in temp.keys():
        ants = Trans_ant._make(['', 'Intergenic', '', '', '', '', '',
                                '', '', '', '', '', '', ''])
        geneannot = {'GENE': temp['Intergenic'][0]}
        geneannot.update({'TRANSCRIPT': ants,
                                            'MIM_PHENS': [], 'MIM_IDS': [],
                                            'GAD_PHENS': []})
    prior_gene = geneannot['GENE']
    ta = geneannot['TRANSCRIPT']
    geneannot['TRANSCRIPT'] = [ta]
    if alltrans and \
            geneannot['TRANSCRIPT'][0].region != 'Intergenic':
        for altid, geneinfo in parsed_dict.items():
            if prior_gene in geneinfo:
                for nt in geneinfo[prior_gene]['TRANSCRIPTS']:
                    if nt not in geneannot['TRANSCRIPT']:
                        geneannot['TRANSCRIPT'].append(nt)

    return geneannot

if __name__ == '__main__':
    import sys
    vcffile = sys.argv[1]
    vcf = VCFParser(vcffile)
    for rec in vcf:
        print rec.chrom, rec.pos, rec.ref, rec.alt, rec.id
        vcf.parseinfo(rec)
        print 'Parsed Varant Annotation - '
        par_ant = parse(rec.info)
        for ac, val in par_ant.items():
            print ac, '\t', val
        print 'Prioritized transcript per gene dictionary - '
        pt_par_ant = prio_trans(par_ant)
        for ac, val in pt_par_ant.items():
            print ac, '\t', val
        print 'Prioritized gene per vcf record - '
        ga = get_prior_geneannot(rec.info, alltrans=False)
        print ga
        print '\n'
