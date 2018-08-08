#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: progconfig
    :platform: Unix, Windows, MacOSX
    :synopsis: Location of various programs and options for programs

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com

Configuration information for various programs used as part of GCN.

Program specific configuration information is read from files with the
same name as the program.  The config files can be system-wide or per user.
By default the program looks for files in $GCNHOME/.gcn and $HOME/.gcn
directories. An user can also specify a $GCNPROGDIR for config files.
Finally the set_configuration method in this module can be called with the name
of a file to load configuration information for a program.
"""
from __future__ import absolute_import
import os
import shlex
from gcn.etc.fileconfig import getlogger


PROGCONFIG = {}


def get_configuration(progname, logger=None, reload=False):
    """Get configuration information associated with a program

    Args:
        progname (string): Name of the program (case insenstive)
        logger (file): A logging object. If None, the system logger
                       is used
        reload (bool): If true, then the configuration will be reloaded
                       from file, even if it had been loaded previously

    Returns:
        A configuration dictionary with the configuration options
        as key
    """

    plower = progname.lower()

    config = PROGCONFIG.get(plower, None)

    if config is None or reload:
        config = PROGCONFIG[plower] = {}
        sysdir = os.environ.get('GCNHOME', None)
        userdir = os.environ.get('HOME', None)
        custdir = os.environ.get('GCNPROGDIR', None)
        if userdir is not None:
            userdir = os.path.join(userdir, '.gcn')
            if not os.path.exists(userdir):
                userdir = None

        if logger is None:
            logger = getlogger()

        for cfdir in (sysdir, userdir, custdir):
            if cfdir is not None:
                cfgfile = os.path.join(cfdir, 'config', progname)
                if os.path.exists(cfgfile):
                    if logger:
                        logger.info('Loading configuration from %s' % cfgfile)
                    config.update(load_file(cfgfile))
                    if logger:
                        logger.info('Loaded configuration successfully')

    return config


def load_file(cfgfile):
    """Load configuration for a program from a file

    Args:
        cfgfile (string): Name of file with configuration information

    Returns:
        A dictionary of options and corresponding values
    """
    cfgdict = empty_configuration()
    with open(cfgfile, 'rU') as stream:
        for line in stream:
            if line.startswith('%include'):
                fname = os.path.join(os.path.dirname(cfgfile),
                                     line.split(None, 1)[1].strip())
                dtmp = load_file(fname)
                cfgdict.update(dtmp)
                continue
            else:
                parse_line(line, cfgdict)

    return cfgdict


def set_configuration(progname, cfgfile, logger=None):
    """Load and save configuration for a program from the specified file

    Args:
        progname (string): Name of the program (case insenstive)
        cfgfile (string): Name of file with configuration information
        logger (file): A logging object. If None, the system logger
                       is used

    Returns:
        A configuration dictionary with the configuration options
        as key

    """
    if logger is None:
        logger = getlogger()
    logger.info('Loading configuration for %s from file %s' % (progname,
                                                               cfgfile))
    PROGCONFIG[progname.lower()] = cfgdict = load_file(cfgfile)
    logger.info('Loaded configuration for %s from file %s' % (progname,
                                                              cfgfile))
    return cfgdict


def empty_configuration():
    """Generate an empty configuration dictionary"""
    return {'exe': [],
            'separator': '=',
            'preopts': [],
            'shortopts': [],
            'longopts': [],
            'definitions': [],
            'preargs': [],
            'postargs': []
            }


def parse_line(line, cfgdict):
    """Parse one line from the configuration file"""
    if not line.strip():
        return
    opt, value = [v.strip() for v in  line.split('=', 1)]

    if '.' in opt:
        opttype, optname = opt.split('.', 1)
        cfgdict[opttype].append((optname, value.strip() or None))

    elif opt == 'separator':
        cfgdict[opt] = value.strip()

    else:
        if opt not in cfgdict:
            cfgdict[opt] = []
        cfgdict[opt].extend(shlex.split(value.strip()))
