#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: fileconfig
    :platform: Unix, Windows, MacOSX
    :synopsis: Location of various files needed for annotations etc

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com

Configuration information for the various files that are need for
annotations in the Genome Commons Navigator.

This module exports one dictionary *FILECONFIG* that contains as keys
the names of the various resources needed for GCN and as
value the full paths to where the files are stored on the local system.

"""

import os
import logging
import logging.handlers
from gcn.config import lib_config

_LOCALDIR = '/opt/data'

def set_localdir():
    """Set the local database directory from the enviornmental variable
    GCN_DATA_DIR"""
    local = lib_config.gcn_path('GCN_DATA_DIR')
    if local is not None:
        global _LOCALDIR
        _LOCALDIR = local

set_localdir()

FILECONFIG = {
    'REFGENE': os.path.join(_LOCALDIR, 'refgene', 'refseq.txt'),
    'GENCODE': os.path.join(_LOCALDIR, 'gencode',
                            'wgEncodeGencodeBasicV19.txt'),
    'FUNSIM': os.path.join(_LOCALDIR, 'geneontology',
                            'uprot20160505_funsim.tsv'),
    'GOA': os.path.join(_LOCALDIR, 'geneontology',
                            'gene_association.goa_human'),
    'GOBO': os.path.join(_LOCALDIR, 'geneontology',
                            'go.obo'),
    'MIRNA': os.path.join(_LOCALDIR, 'mirna',
                          'human_predictions_S_C_aug2010.txt'),
    'UTRDB': os.path.join(_LOCALDIR, 'utrdb', ''),
    'CLNPHESNP': os.path.join(_LOCALDIR, 'clnphesnp', ''),
    'LOGFILE': lib_config.gcn_path('GCN_LOGFILE'),
    'UNIPROT': os.path.join(_LOCALDIR, 'uniprot',
                           'uniprot_sprot_human.dat'),
    'OMIM': os.path.join(_LOCALDIR, 'omim',
                           'uniprot_sprot_human.dat'),
    'NSFP': os.path.join(_LOCALDIR, 'dbnsfp', ''),
    'GERP': os.path.join(_LOCALDIR, 'gerp', '*.txt'),
    'REFGENOME': os.path.join(_LOCALDIR, 'refgenome', 'human_g1k_v37.fasta'),
    'REFGENOME_UCSC': os.path.join(_LOCALDIR, 'refgenome', 'hg19.fasta'),
    'REGULOME': os.path.join(_LOCALDIR, 'regulome', ''),
    'SeqCapEZ_Exome': os.path.join(_LOCALDIR, 'capture_bed', 'nimblegen',
                                       'SeqCap_EZ_Exome_v3_capture.bed'),
    'SureSelect_V6': os.path.join(_LOCALDIR, 'capture_bed', 'sureselect',
                                       'S07604514_Covered.bed'),
    'CLINVARDB': os.path.join(_LOCALDIR, 'snps', 'dbsnp',
                               'clinvar-latest.vcf'),
    'DBSNP': os.path.join(_LOCALDIR, 'snps', 'dbsnp', '00-All.vcf'),
    'ESP': os.path.join(_LOCALDIR, 'snps', 'esp',
                        'ESP6500SI-V2-SSA137.updatedRsIds.snps_indels.vcf'),
    'EXAC': os.path.join(_LOCALDIR, 'snps', 'exac',
                        'ExAC.r0.3.sites.vcf'),    
    'KGDB': os.path.join(_LOCALDIR, 'snps',
                         '1000genomes',
         'ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf'),
    'LCR': os.path.join(_LOCALDIR, 'lcr', 'LCR-hs37d5.bed'),
    'CLINVITAE': os.path.join(_LOCALDIR, 'snps','clinvitae','clinvitae.vcf'),
    'COSMIC': os.path.join(_LOCALDIR, 'snps','cosmic','cosmic.vcf'),
    'INTERPRO': os.path.join(_LOCALDIR, 'interpro','interpro_hg19_coord.tsv'),
    'PATHOG_PROF': os.path.join(_LOCALDIR, 'snps','dbsnp','pathov_gene_prof.json')
    }


def set_location(resource, filename, exists=True):
    """Set the location of the file for a particular data item.

    Args:
        resource (string): Name of the resource (e.g REFGENE, REFMRNA etc.)
        filename (string): Full path to the file that contains the information
        exists (boolean):  If true then the filename should already exist

    Returns
        None. Updates the `FILECONFIG` dictionary in-place

    Raises
        IOError if the file is required to exist already but does not exist

    This method can be used as an alternative to configuring the locations
    of files on a per-user basis or when customized resources are to be used
    instead of the system defaults.  An example would be to temporarily
    change the LOGFILE to something other than the system default
    """

    if exists and not os.path.exists(filename):
        raise IOError('%s no such file' % filename)

    FILECONFIG[resource] = filename


def load_default_configuration(filename=None, logger=None):
    """Get file paths for various resources such as dbsnp, refgene etc.

    Args:
        filename (str): File containing configuration information
        logger (logging.Logger): A logger to which loading status
                                 will be written (optional)

    Returns:
        Loads default configuration into `FILECONFIG`

    If `filename` is not specified an attempt is made to load
    configuration information from pre-defined configuration files.
    Two directories are checked for the presence of configuration files.
    These are, in order, $GCNHOME/config/gcnfileconfig and
    $HOME/.gcn/config/gcnfileconfig. The configuration options from $HOME
    override the options in $GCNHOME.
    """

    if filename is not None:
        _load_file(filename)
        return

    sysdir = os.environ.get('GCNHOME', None)
    userdir = os.environ.get('HOME', None)
    if userdir is not None:
        userdir = os.path.join(userdir, '.gcn')
        if not os.path.exists(userdir):
            userdir = None

    for cfdir in (sysdir, userdir):
        if cfdir is not None:
            cfgfile = os.path.join(cfdir, 'config', 'gcnfileconfig')
            if not os.path.exists(cfgfile):
                continue
            _load_file(cfgfile)


def _load_file(cfgfile):
    with open(cfgfile, 'rU') as stream:
        for line in stream:
            if line[:1] == '#':
                continue
            resource, name = line.split('=', 1)
            set_location(resource.upper().strip(), name.strip())


def getlogger(level=logging.INFO,program_name='gcn'):
    """Retrieve a logger for the package. The logger uses the
    `LOGFILE` value in the `FILECONFIG` dictionary as the file to
    which actions are logged.

    Args:
        level (enum): sets the logging level, should be one of
                      values specified in the logging module.

    Returns:
        A logger with a rotating file handler
    """

    logger = logging.getLogger(program_name)
    if logger.handlers:
        return logger
    handler = logging.handlers.RotatingFileHandler(FILECONFIG['LOGFILE'],
                                                   maxBytes=1024000,
                                                   backupCount=5)
    f = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(f)
    handler.setFormatter(formatter)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger

load_default_configuration()
if __name__ == '__main__':
    logger = getlogger()
