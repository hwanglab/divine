"""
.. module:: annotateNimblegenCap
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for retrieving Nimblegen Capture
               annotation for given SNP.

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Class for retrieving Nimblegen capture annotation for given SNP if it is part
of capture array.
"""

from bisect import bisect
from gcn.etc import fileconfig
from collections import namedtuple
import os
from datetime import datetime as dt


class NimblegenCapture:
    """ This Module annotates the variants
        based on NimbleGen Capture Array dataset """

    def __init__(self, capture_kit_name='SeqCapEZ_Exome', ext_bp=50):
        self.ext_bp = ext_bp #bp length to extend from both ends
        
        if capture_kit_name == 'SeqCapEZ_Exome':
          self.kit_company = 'NimbleGen'
          self.kit_alias = 'Nimblegen_Capture_Array'
        elif capture_kit_name == 'SureSelect_V6':
          self.kit_company = 'Agilent'
          self.kit_alias = 'V6_WES'

        self.kit_symbol = capture_kit_name
        self.bedfile = fileconfig.FILECONFIG[self.kit_symbol]
        self.bed_file = self.load_bed()

    def get_version(self):
        """Returns the version of the Capture Array Files"""
        Version = namedtuple('Version', ['dbname', 'file_date',
                                         'version', 'entry_count'])
        fd = dt.fromtimestamp(os.path.getmtime(self.bedfile)\
                              ).strftime('%Y-%m-%d')
        version = os.path.basename(self.bedfile)
        return [Version._make((self.kit_alias, fd, version, \
                               self.entry_cnt))]

    def get_header(self):
        """ Intializing vcf formatted header information for NimbleGen
        Capture Array data"""
        Nimblegen_h = namedtuple('%s_h'%self.kit_company, 'id count type desc')

        header = [Nimblegen_h('CaptureCore', 0, 'Flag', 'Variant is present in %s Capture bed file'%self.kit_symbol\
                              ),
                  Nimblegen_h('Capture5p', 1, 'Integer', "Number of "\
                              "bases that the position is 5p upstream"\
                              " of capture region"),
                  Nimblegen_h('Capture3p', 1, 'Integer', "Number of "\
                              "bases that the position is 3p downstream"\
                              " of capture region"),
                  Nimblegen_h('CaptureIntervalLen', 1, 'Integer',
                              "Length of the capture interval")]
        return header

    def load_bed(self):
        """ Loading Nimblegen Capture array data"""
        d = {}
        self.entry_cnt = 0
        for line in open(self.bedfile):

            if line.startswith('browser') or line.startswith('track'):
                continue
            self.entry_cnt += 1
            chrom, p1, p2 = line.split(None, 3)[:3]
            
            chrom0 = chrom[3:]
            if chrom0 in d:
              pos = d[chrom0]
            else:
              pos = []
              d[chrom0] = pos

            p1 = int(p1)
            p2 = int(p2)
            pos.append(((p1 - self.ext_bp), (p2 + self.ext_bp)))

        for key in d:
            v = d[key]
            v.sort()
            d[key] = zip(*v)
        return d

    def check(self, chrno, pos):
        """ Checking if the position is part of NimbleGen Capture array """
        cap_ant = (0, None, None, None)
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
                        if (s[j - 1] + self.ext_bp) <= qpos < (e[j - 1] - self.ext_bp):
                            data = 'CaptureCore'
                            tvalue = True
                        elif qpos < (s[j - 1] + self.ext_bp):
                            tvalue = (s[j - 1] + self.ext_bp) - qpos
                            data = 'Capture5p'
                        elif qpos >= (e[j - 1] - self.ext_bp):
                            tvalue = qpos - (e[j - 1] - self.ext_bp)
                            data = 'Capture3p'
                        intervalLen = e[j - 1] - s[j - 1] - self.ext_bp*2
                        if data == 'CaptureCore':
                            cap_ant = (1, intervalLen, data, tvalue)
                            break
                        elif tvalue < value or not value:
                            value = tvalue
                            cap_ant = (1, intervalLen, data, value)
                    elif e[idx - 1] < qpos:
                        break
        return cap_ant

    def get_info(self, chrno, spos, epos=None):
        """  Comparing the given variants within NimbleGen Capture data
             and extracting the annotation information """
        if epos == None:
            epos = spos
        info = None
        for pos in range(spos, epos + 1):
            tup = self.check(chrno, pos)
            if tup != None and tup[0] == 1:
                info = ['CaptureIntervalLen', tup[1], tup[2], tup[3]]
                break
        return info

if __name__ in "__main__":
    nimb_cap = NimblegenCapture()
    print nimb_cap.get_version()
    nimbc_ant = nimb_cap.get_info('5', 10250430)
    print nimbc_ant
