"""
.. module:: utr_parser
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for building a parsed dictionary for the utrdb data

.. moduleauthor:: Kunal Kundu (kunal.kundu@gmail.com)

Class for building a parsed dictionary for the utrdb data
(currently only parses the 'FT' lines)
"""

from gcn.etc import fileconfig
from gcn.lib.io import anyopen


class UTR_Parser():

    def __init__(self, filepath, sections, delimiter='//'):
        self.filename = filepath
        self.delimit = delimiter
        self.sec = sections
        self.utr3_features = []
        self.utr5_features = []

    def records(self):
        """ iterates over the records of utrdb"""
        stream = anyopen.openfile(self.filename, 'r')
        temp = ""
        for line in stream:
            if self.delimit not in line[:3]:
                temp += line
            else:
                yield temp
                temp = ""

    def extract_sec(self, entry):
        """ For a given entry/record retrieves the sections"""
        sec_dict = {}
        for s in self.sec:
            sec_dict[s] = ""
            for line in entry.split('\n'):
                if line.startswith(s):
                    sec_dict[s] += line + '\n'
        sec_dict.values()[0].strip()
        return sec_dict

    def parse_FT(self, entry):
        '''Builds a parsed dictionary for the feature section
        of the utrdb record'''
        ft_dict = {}
        for line in entry.split('\n'):
            feature = line.strip('FT')[:19].strip()
            if feature:
                key = feature
                ft_dict[key] = []
            else:
                if line:
                    if '=' in line[20:].strip():
                        t = line[20:].strip().strip('/').split('=')
                        m = {}
                        m[t[0]] = t[1].strip("\"")
                        ft_dict[key].append(m)
                    else:
                        m = ft_dict[key].pop()
                        m[t[0]] += line[20:].strip("\"").strip()
                        ft_dict[key].append(m)
        return ft_dict

    def parse(self):
        'iterates the parsed dictionary for each utrdb record'
        parsed_dict = {}
        for entry in self.records():
            sec_dict = self.extract_sec(entry)
            for sec, data in sec_dict.iteritems():
                if sec == 'FT':
                    parsed_dict[sec] = self.parse_FT(data)
            yield parsed_dict


if __name__ == '__main__':
    import os
    dirpath = fileconfig.FILECONFIG['UTRDB']
    for filename in os.listdir(dirpath):
            f = os.path.join(dirpath, filename)
            sections = ['FT']
            p = UTR_Parser(f, sections)
            for pd in p.parse():
                print pd['FT']
