"""
.. module:: annotateGrep
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for retrieving GERP Conservation scores for given SNP

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for retrieving GERP Conservation scores for given SNP. For each conserved
 variant the RS-Score and p-value is also extracted

The GERP Conserved Element file for hg19 is downloaded from
http://mendel.stanford.edu/SidowLab/downloads/gerp/
"""

from bisect import bisect
from gcn.lib.io import anyopen
from gcn.etc import fileconfig
import os
from collections import namedtuple
from datetime import datetime as dt


class Gerp:
    """ This Module annotates the variants
        based on Gerp dataset """

    def __init__(self):
        self.gerpfiles = fileconfig.FILECONFIG['GERP']
        self.ld_gerp, self.gerp_score = self.load_interval()

    def get_version(self):
        """Returns the version of the Gerp Files"""
        Version = namedtuple('Version', ['dbname', 'file_date',
                                         'version', 'entry_count'])
        gerp_dir = os.path.dirname(self.gerpfiles)
        for fn in os.listdir(gerp_dir):
            infile = os.path.join(gerp_dir, fn)
            break
        fd = dt.fromtimestamp(os.path.getmtime(infile)\
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        return [Version._make(('Gerp', fd, version, self.entry_cnt))]

    def get_header(self):
        """ Intializing vcf formatte dheader information for Gerp data"""

        Gerp_h = namedtuple('Gerp_h', 'id count type desc')

        header = [Gerp_h('GerpConserve', 0, 'Flag', "If variant is " \
                         "present in Gerp Element file it is tagged "\
                         "with this flag"),
                  Gerp_h('GerpRSScore', 1, 'Float', "Gerp Rejected "\
                         "Substitutions Score. Rejected substitutions "\
                         "are a natural measure of constraint that "\
                         "reflects the strength of past purifying "\
                         "selection on the element."),
                  Gerp_h('GerpPValue', 1, 'String', "Gerp P-Value")]
        return header

    def load_interval(self):
        """ Loading GERP Element interval data"""
        d = {}
        d_scores = {}
        chrom = None
        self.entry_cnt = 0
        gerp_dir = os.path.dirname(self.gerpfiles)
        for filename in os.listdir(gerp_dir):
            if '_elems.txt' in filename:
                chrom = filename.split('_')[1].strip('chr')
                gerp_file = os.path.join(gerp_dir, filename)
                for line in anyopen.openfile(gerp_file, 'r'):
                    if line[:1] == '#':
                        continue
                    self.entry_cnt += 1
                    p1, p2, len, score, pvalue = line.split()
                    d_scores[(chrom, int(p1), int(p2))] = (score, pvalue)
                    try:
                        pos = d[chrom]
                    except KeyError:
                        pos = []
                        d[chrom] = pos
                    p1 = int(p1)
                    p2 = int(p2)
                    pos.append((p1, p2))
        for key in d:
            v = d[key]
            v.sort()
            d[key] = zip(*v)
        return d, d_scores

    def check(self, chrno, pos):
        """ Checking the variants for GERP Conservation and
        return RS-score and p-value """
        rs_score, pvalue = (None, None)
        pos = int(pos)
        if chrno in self.ld_gerp:
            s, e = self.ld_gerp[chrno]
            idx = bisect(s, pos)
            if idx != 0:
                if s[idx - 1] <= pos <= e[idx - 1]:
                    rs_score, pvalue = self.gerp_score[(chrno, s[idx - 1], \
                                                        e[idx - 1])]
        return (rs_score, pvalue)

    def get_conserve(self, chrno, spos, epos=None):
        """  extracting the annotation information of GERP conservation
        Args:
            - snp:    A tuple that contains (chrom, position)

        Returns:
            - info(String):    It contains the Conserve Flag,
                                RS-Score and P-Value
        """
        if chrno[:3] == 'chr':
            chrno = chrno[3:]
        if epos == None:
            epos = spos
        consv_info = None
        rs_score = 0.0
        pvalue = None
        for pos in range(spos, epos + 1):
            tup = self.check(chrno, pos)
            if tup[0] is not None:
                if float(tup[0]) > rs_score:
                    rs_score = float(tup[0])
                    pvalue = tup[1]
        if pvalue:
            consv_info = ['GerpConserve', 'True', 'GerpRSScore',
                          rs_score, 'GerpPValue', tup[1]]
        return consv_info

if __name__ in "__main__":
    fg = Gerp()
    print fg.get_version()
    print fg.get_conserve('chr1', 2543587, 2543589)
    print fg.get_conserve('1', 1334052)
