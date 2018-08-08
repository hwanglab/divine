"""
.. module:: annotateNimblegenCap
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for retrieving Low complexity
               region annotation for given variant.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for retrieving Low complexity region annotation for given variant if it
is part of the LCR.
The LCR data used was distributed as supplement material for the paper -
'Towards Better Understanding of Artifacts in Variant Calling from
High-Coverage Samples' by Heng Li
Raw file can be downloaded from -
https://raw.githubusercontent.com/lh3/varcmp/master/scripts/LCR-hs37d5.bed.gz
"""

from bisect import bisect
from gcn.etc import fileconfig
from collections import namedtuple
import os
from datetime import datetime as dt


class LCR:
    """ This Module annotates the variants
        based on LCR data """

    def __init__(self):
        self.bedfile = fileconfig.FILECONFIG['LCR']
        self.bed_file = self.load_bed()

    def get_version(self):
        """Returns the version of the LCR File"""
        Version = namedtuple('Version', ['dbname', 'file_date',
                                         'version', 'entry_count'])
        fd = dt.fromtimestamp(os.path.getmtime(self.bedfile)
                              ).strftime('%Y-%m-%d')
        version = "v%s_%s" % tuple(fd.split('-')[:2])
        return [Version._make(('LCR', fd, version,
                               self.entry_cnt))]

    def get_header(self):
        """ Intializing vcf formatted header information for LCR
        annotation"""
        Lcr_h = namedtuple('Lcr_h', 'id count type desc')
        header = [Lcr_h('LCR', 0, 'Flag', "Variant position"
                        "is part of low-complexity region")]
        return header

    def load_bed(self):
        """ Loading LCR data"""
        d = {}
        self.entry_cnt = 0
        for line in open(self.bedfile):

            if line[:5] == 'track':
                continue
            self.entry_cnt += 1
            chrom, p1, p2 = line.split(None, 3)[:3]

            try:
                pos = d[chrom]
            except KeyError:
                pos = []
                d[chrom] = pos
            p1 = int(p1)
            p2 = int(p2)
            pos.append(((p1 - 50), (p2 + 50)))

        for key in d:
            v = d[key]
            v.sort()
            d[key] = zip(*v)
        return d

    def check(self, chrno, pos):
        """ Checking if the position is part of LCR"""
        lcr_ant = (0, None, None, None)
        if 'chr' in chrno:
            chrno = chrno[3:]
        qpos = int(pos) - 1
        if chrno in self.bed_file:
            s, e = self.bed_file[chrno]
            idx = bisect(s, qpos)
            if idx != 0:
                value = 0
                for j in range(idx, 0, -1):
                    if s[j - 1] <= qpos <= e[j - 1]:
                        if (s[j - 1] + 50) <= qpos < (e[j - 1] - 50):
                            data = 'LCR'
                            tvalue = True
                        elif qpos < (s[j - 1] + 50):
                            tvalue = (s[j - 1] + 50) - qpos
                            data = 'LCR5p'
                        elif qpos >= (e[j - 1] - 50):
                            tvalue = qpos - (e[j - 1] - 50) + 1
                            data = 'LCR3p'
                        if data == 'LCR':
                            lcr_ant = (1, data, tvalue)
                            break
                        elif tvalue < value or not value:
                            value = tvalue
                            lcr_ant = (1, data, value)
                    elif e[idx - 1] < qpos:
                        break
        return lcr_ant

    def get_info(self, chrno, spos, epos=None):
        """  Comparing the given variants within NimbleGen Capture data
             and extracting the annotation information """
        if epos is None:
            epos = spos
        info = None
        for pos in range(spos, epos + 1):
            tup = self.check(chrno, pos)
            if tup is not None and tup[0] == 1 and tup[1] == 'LCR':
                info = [tup[1], tup[2]]
                break
        return info

if __name__ in "__main__":
    lcr = LCR()
    print lcr.get_version()
    lcr_ant = lcr.get_info('1', 29055360)
    print lcr_ant
