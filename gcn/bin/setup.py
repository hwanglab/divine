#
# COPYRIGHT (C) 2013-2014 TCS Ltd
#
"""
.. module:: setup
    :platform: Unix
    :synopsis: Script to databases used by Varant

.. moduleauthor:: Kunal Kundu (kunal@atc.tcs.com); modified by changjin.hong@gmail.com

Script to set up all the databases used by Varant. To use this script

usage: setup.py [-h] [-p PROC]

Set up Databases for Varant

optional arguments:
  -h, --help            show this help message and exit
  -p PROC, --proc PROC  Number of processors to be used during installation

This script primarily performs 3 steps for all the databases used by Varant-
1. Creates folder where the data source file will be downloaded
2. Downloads the data
3. Calls method that uncompresses the downloaded data
4. Creates SQLite database for the downloaded data

The script uses gcn/data/data_links.txt file to look up for the download
links for different data sources. Thus this text file must be up to date with
the latest download URLs for the respective data sources.
"""

import os
from gcn.etc.fileconfig import FILECONFIG, getlogger
from gcn.etc.dbconfig import DBCONFIG
import urllib
from gcn.bin.mksnpdb import SNPDB
from gcn.bin.mkrefgene import RefgeneDB
from gcn.bin.mkmrna import RefmrnaDB
from gcn.bin.mkclnsnpdb import ClnSNPDB
from gcn.bin.mkomim import OmimDB
from gcn.bin.mkmirna import MiRnaDB
from gcn.bin.mkutrdb import UtrDB
from gcn.bin.mknsfpdb import NsfpDB
from gcn.bin.mkregulomedb import RegulomeDB
from gcn.bin.mksplicedb import SpliceRegSite
from multiprocessing import Pool
import time
import datetime
import zipfile
import tarfile
import shutil
import argparse

gcnpath = os.getenv('GCN', None)
DATA_LINKS = os.path.join(gcnpath, 'gcn', 'data', 'data_links.txt')
CREATE_OBJ = {'DBSNP': SNPDB(), 'KGDB': SNPDB(), 'ESP': SNPDB(), 'EXAC': SNPDB(),
              'CLINVARDB': SNPDB(), 'REFGENE': RefgeneDB(),
              'REFMRNA': RefmrnaDB(), 'CLNPHESNP': ClnSNPDB(),
              'OMIM': OmimDB(), 'NSFP': NsfpDB(), 'SPLICE': SpliceRegSite(),
              'REGULOME': RegulomeDB(), 'UTRDB': UtrDB(), 'MIRNA': MiRnaDB()}
CHROMLIST = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
             '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
             'X', 'Y']


def _get_download_links():
    """Private Method. Loads the data/data_links.txt to memory
    The data_links.txt comprises of download links for all the data sources
    used by Varant"""
    dlinks = {}
    s = open(DATA_LINKS)
    for line in s:
        line = line.strip().strip('\n')
        if not line or line.startswith('#'):
            continue
        d = line.split('\t')
        if d[0] not in dlinks:
            dlinks[d[0]] = [d[1]]
        else:
            dlinks[d[0]].append(d[1])
    return dlinks


def _uncompress(downloaded_files):
    """Private Method. Uncompresses files"""
    for filepath in downloaded_files:
        dirpath = os.path.dirname(filepath)
        if zipfile.is_zipfile(filepath):  # Uncompress .zip files
            zfile = zipfile.ZipFile(filepath)
            for name in zfile.namelist():
                zfile.extract(name, dirpath)
                exdir = os.path.dirname(name)
                if exdir:
                    shutil.copy(os.path.join(dirpath, name), dirpath)
        elif tarfile.is_tarfile(filepath):  # Uncompress .tar.gz files
            tar = tarfile.open(filepath, 'r:gz')
            for item in tar:
                tar.extract(item, dirpath)
        else:
            fn, ext = os.path.splitext(os.path.basename(filepath))
            if ext in ['.txt', 'tsv', '.csv', '.dat']:
                continue
            if ext == '.gz':
                os.system("gunzip %s" % filepath)


def build((dtype, durls, data_dest_dir)):
    """The key method that performs the following steps -
    1. Creates folder where the data source file will be downloaded
    2. Downloads the data
    3. Calls method that uncompresses the downloaded data
    4. Creates SQLite database for the downloaded data
    """
    t1 = time.time()
    logger.info('Processing %s' % dtype)
    print 'Processing %s' % dtype

    # Create the directory where the data-source will be downloaded
    if not os.path.exists(data_dest_dir):
        try:
            os.makedirs(data_dest_dir)
            logger.info('%s path created..' % data_dest_dir)
        except:
            logger.error('%s path cannot be created..' % data_dest_dir)
            return dtype
    else:
        logger.info('%s path already exists..' % data_dest_dir)

    # Download data-source
    testfile = urllib.URLopener()
    downloaded_files = []
    for link in durls:
        file_name = os.path.basename(link)
        local_filename = os.path.join(data_dest_dir, file_name)
        try:
            testfile.retrieve(link, local_filename)
            downloaded_files.append(local_filename)
        except:
            logger.error('Cannot download file from %s' % link)
            return dtype

    # If compressed, uncompress the data
    try:
        _uncompress(downloaded_files)
    except:
        logger.error('Encountered an error while uncompressing %s'
                     % ','.join(downloaded_files))
        return dtype

    # Create SQLite for required data sources
    if dtype in CREATE_OBJ:
        if dtype == 'ESP':
            infile = FILECONFIG[dtype]
            outfile = DBCONFIG[dtype]['name']
            for chrom in CHROMLIST:
                if chrom == '1':
                    continue
                cmd1 = "cat %s/ESP6500SI-V2-SSA137.updatedRsIds.chr%s.snps_indels.vcf "\
                "| grep -v '^#' >> %s/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf"\
                        % (data_dest_dir, chrom, data_dest_dir)
                os.system(cmd1)
            os.rename(data_dest_dir + '/ESP6500SI-V2-SSA137.updatedRsIds.chr1.snps_indels.vcf', infile)
            try:
                CREATE_OBJ['ESP'].createdb(infile, outfile, True)
            except:
                logger.error('%s creation failed..' % outfile)
                return dtype
        elif dtype == 'REFGENE':
            for t in ['REFGENE', 'REFMRNA']:
                if t == 'REFGENE':
                    infile = FILECONFIG['REFGENE']
                    outfile = DBCONFIG['REFGENE']['name']
                elif t == 'REFMRNA':
                    print t
                    infile = FILECONFIG['REFGENOME']
                    outfile = DBCONFIG['REFMRNA']['name']
                if os.path.exists(infile):
                    try:
                        CREATE_OBJ[t].createdb(infile, outfile, True)
                    except:
                        logger.error('%s creation failed..' % outfile)
                        return t
                else:
                    logger.error('%s not found' % infile)
                    return t
        else:
            infile = FILECONFIG[dtype]
            outfile = DBCONFIG[dtype]['name']
            if os.path.exists(infile):
                try:
                    CREATE_OBJ[dtype].createdb(infile, outfile, True)
                except:
                    logger.error('%s creation failed..' % outfile)
                    return dtype
            else:
                logger.error('%s not found' % infile)
                return dtype
    t2 = time.time()
    tt = float(t2 - t1) / float(60)
    logger.info('Time taken for processing %s is %f min' % (dtype, tt))
    print 'Time taken for processing %s is %f min' % (dtype, tt)


def build_splicedb():
    """Method to call the module that creates Splice DB"""
    dtype = 'SPLICE'
    t1 = time.time()
    logger.info('Processing %s' % dtype)
    print 'Processing %s' % dtype

    infile = FILECONFIG['REFGENE']
    outfile = DBCONFIG[dtype]['name']
    if os.path.exists(infile):
        try:
            CREATE_OBJ[dtype].createdb(infile, outfile, True)
        except:
            logger.error('%s creation failed..' % outfile)
            return dtype
    else:
        logger.error('%s not found' % infile)
        return dtype
    t2 = time.time()
    tt = float(t2 - t1) / float(60)
    logger.info('Time taken for processing %s is %f min' % (dtype, tt))
    print 'Time taken for processing %s is %f min' % (dtype, tt)
    return


def main(proc):
    global logger
    code_loc = os.environ.get('GCN')
    ts = datetime.datetime.fromtimestamp(
                            time.time()).strftime('%Y%m%d_%H%M')
    logfile = os.path.join(code_loc, 'varant_install_%s.log' % ts)
    FILECONFIG['LOGFILE'] = logfile
    logger = getlogger()
    logger.info('Loading Data download links from %s' % DATA_LINKS)
    dlinks = _get_download_links()
    logger.info('Data download links loaded to memory')

    p = Pool(proc)

    data_inputs = []
    not_build = []
    refgenome = []

    for ds, links in dlinks.items():
        if ds in FILECONFIG.keys():
            data_dest = os.path.dirname(FILECONFIG[ds])
            if ds == 'REFGENOME':
                refgenome.append((ds, links, data_dest))
            else:
                data_inputs.append((ds, links, data_dest))
        else:
            not_build.append(ds)
            logger.error('Download data destination folder is not set '
                         'for %s in etc/fileconfig.py and thus cannot '
                         'be build..' % ds)
    t1 = time.time()
    for e in refgenome:
        build(e)
    r = p.map(build, data_inputs)
    [not_build.append(e) for e in r if e]
    if 'REFGENE' in dlinks.keys() and 'REFGENE' not in not_build and \
                                        'REFMRNA' not in not_build:
        r = build_splicedb()
        if r is not None:
            not_build.append(r)
    t2 = time.time()
    tt = float(t2 - t1) / float(60)
    if not_build:
        logger.warning('Failed to build the following data sources - ')
        print 'Failed to build the following data sources - '
        for idx, d in enumerate(not_build):
            logger.warning(str(idx + 1) + ':' + d)
            print str(idx + 1) + ':' + d
    logger.info('(%d Processor) Total time taken to set up Varant '
                'databases is %f min' % (proc, tt))
    print '(%d Processor) Total time taken to set up is %f min' % (proc, tt)

if __name__ == '__main__':
    proc = 1
    desc = 'Set up Databases for Varant'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p', '--proc', dest='proc', type=int,
                        default=proc, help='Number of processors to be '
                        'used during installation')
    args = parser.parse_args()
    main(args.proc)
