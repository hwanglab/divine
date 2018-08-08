"""
.. module:: converter
    :platform: Unix, Windows, MacOSX
    :synopsis: hg coordinate version convertor

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This converts the human genome coordinate from hg18 to hg19 and visa versa.
"""

import os
from gcn.config import lib_config

class hgConvert:

    def __init__(self):
        self.location = lib_config.gcn_path('GCN')
        self.liftover = os.path.join(self.location, \
                                      'gcn/bin/hgConvert/liftOver_exe', \
                                      'liftOver')
        self.chaindata = os.path.join(self.location, 'gcn/bin/hgConvert', \
                                      'chainData')
        self.tmpdir = os.path.join(self.location, 'gcn/bin/hgConvert', 'temp')
        if not os.path.exists(self.tmpdir):
            os.system('mkdir %s' % (self.tmpdir))

    def hg19_to_hg18(self, coordlist):
        """ Converts a hg19 genome coordinate to hg18 genome coordinate
        Args:
            coordlist (str):    A list of hg19 coordinates. Format of the
                                element in the list is a tuple -
                                [(chrom, start_pos, end_pos),]
                                For Eg. [('chr1', 16382877, 16382878),
                                         ('chr1', 16382890, 16382910)]
                                For the cases where only one coordinate
                                is needed to be converted the start_pos equals
                                end_pos.
                                For Eg. [('chr1', 16382877, 16382877), ]
        """
        mapdict = {}
        infile = os.path.join(self.tmpdir, 'hg19coord')
        outfile = os.path.join(self.tmpdir, 'hg18coord')
        outfile_unmapped = os.path.join(self.tmpdir, 'hg18coord_unmapped')
        fw = open(infile, 'w')
        for coord in coordlist:
            if not coord[0].startswith('chr'):
                chrom = 'chr' + coord[0]
            else:
                chrom = coord[0]
            data = str(chrom) + ':' + str(coord[1]) + '-' + \
                                str(coord[2]) + '\n'
            fw.write(data)
        fw.close()

        cmd = ' '.join([self.liftover, '-positions', infile, os.path.join(\
                        self.chaindata, 'hg19ToHg18.over.chain.gz'), \
                        outfile, outfile_unmapped])
        os.system(cmd)

        map_stream = open(outfile, 'r')
        map_data = [e.strip() for e in map_stream.readlines()]

        unmap_stream = open(outfile_unmapped, 'r')
        unmap_data = [e.strip() for e in unmap_stream.readlines() \
                      if e[0] != '#']

        stream = open(infile, 'r')
        indata = [e.strip() for e in stream.readlines()]

        for coord in indata:
            if coord in unmap_data:
                mapdict[coord] = ()
            else:
                chrom, pos = map_data.pop(0).split(':')
                start_pos, end_pos = pos.split('-')
                mapdict[coord] = (chrom, int(start_pos), int(end_pos))

        #Removing the generated intermediate files
        os.system("rm %s" % (self.tmpdir + '/*'))
        os.system("rm *.bed*")

        return mapdict

    def hg18_to_hg19(self, coordlist):
        """ Converts a hg18 genome coordinate to hg19 genome coordinate
        Args:
            coordlist (str):    A list of hg18 coordinates. Format of the
                                element in the list is a tuple -
                                [(chrom, start_pos, end_pos),]
                                For Eg. [('chr1', 16255464, 16255466),
                                         ('chr1', 16255470, 16255490)]
                                For the cases where only one coordinate
                                is needed to be converted the start_pos equals
                                end_pos.
                                For Eg. [('chr1', 16255464, 16255464), ]
        """
        mapdict = {}
        infile = os.path.join(self.tmpdir, 'hg18coord')
        outfile = os.path.join(self.tmpdir, 'hg19coord')
        outfile_unmapped = os.path.join(self.tmpdir, 'hg19coord_unmapped')
        fw = open(infile, 'w')

        for coord in coordlist:
            if not coord[0].startswith('chr'):
                chrom = 'chr' + coord[0]
            else:
                chrom = coord[0]
            data = str(chrom) + ':' + str(coord[1]) + '-' + \
                            str(coord[2]) + '\n'
            fw.write(data)
        fw.close()

        cmd = ' '.join([self.liftover, '-positions', infile, os.path.join(\
                        self.chaindata, 'hg18ToHg19.over.chain.gz'), \
                        outfile, outfile_unmapped])
        os.system(cmd)

        map_stream = open(outfile, 'r')
        map_data = [e.strip() for e in map_stream.readlines()]

        unmap_stream = open(outfile_unmapped, 'r')
        unmap_data = [e.strip() for e in unmap_stream.readlines() \
                      if e[0] != '#']

        stream = open(infile, 'r')
        indata = [e.strip() for e in stream.readlines()]

        #for coord in list(set(indata) - set(unmap_data)):
        for coord in indata:
            if coord in unmap_data:
                mapdict[coord] = ()
            else:
                chrom, pos = map_data.pop(0).split(':')
                start_pos, end_pos = pos.split('-')
                mapdict[coord] = (chrom, int(start_pos), int(end_pos))

        #Removing the generated intermediate files
        os.system("rm %s" % (self.tmpdir + '/*'))
        os.system("rm *.bed*")

        return mapdict

if __name__ == '__main__':
    chrom = '1'
    hg = hgConvert()
    pos = 16255464
    print 'hg18 to hg19 : ', hg.hg18_to_hg19(
                            [(chrom, pos, pos + 5),
                            (chrom, pos + 1000000000, pos + 1000000000),
                            (chrom, pos + 2000, pos + 20100),
                            (chrom, pos + 2000000000, pos + 2000000000)])
    pos = 16382877
    print 'hg19 to hg18 : ', hg.hg19_to_hg18(
                            [(chrom, pos, pos), \
                            (chrom, pos + 1000000000, pos + 1000000000),
                            (chrom, pos + 2000, pos + 20100),
                            (chrom, pos + 2000000000, pos + 2000000000)])
    print 'hg18 to hg19 : ', hg.hg18_to_hg19([('chr17', 59875702, 59875755)])
