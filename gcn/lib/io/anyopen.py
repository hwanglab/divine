#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: anyopen
    :platform: Unix, Windows, MacOSX
    :synopsis: Transparent opening of compressed and uncompressed files

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Helper to open files of various types.  This module defines a function
*openfile* which is passed the name of the file and returns a file object
opened for reading.  The function can open files compressed using gzip,
bzip2, unix compress in addition to plain text uncompressed files.
"""

import gzip
import bz2
import os
import tempfile
import sys
if sys.version_info[0] == 3:
    from urllib.request import urlopen
    PYVERSION = 3
else:
    from urllib2 import urlopen
    PYVERSION = 2


class _wGzip(gzip.GzipFile):

    def readline(self, size=-1, _rl=gzip.GzipFile.readline):
        return _rl(self, size).decode()


class _wBzip(bz2.BZ2File):

    def readline(self, size=-1, _rl=bz2.BZ2File.readline):
        return _rl(self, size).decode()


def openfile(filename, mode='rb'):
    """Open a file for reading

    Many of the files that are used in bioinformatics are compressed due
    of their huge size. Some are compressed using gzip and others using
    unix compress. It becomes tedious for programs that process the data in
    these files to check for the file format and then use the appropriate
    method to open the file.  Using this function obviates the need for
    such programs to check for the file format. The function handles files
    compressed using gzip, bzip2 and unix compress in addition to uncompressed
    files. As a bonus files can also be opened using http or ftp protocols

    Args:
        filename (str): Name or URL of the file to be opened

    Returns:
        A file object opened for reading

    Examples

    >>> infile = openfile('/var/data/uniprot/uniprot_sprot.dat.gz')
    >>> for line in infile:
    ...     pass

    >>> infile = openfile('/var/data/enzyme/enzyme.dat.Z')
    >>> for line in infile:
    ...     pass

    """

    if filename.startswith('http://') or filename.startswith('ftp://') or \
           filename.startswith('file://'):
        tmpname = __geturl(filename)
    else:
        tmpname = filename

    fobj = open(tmpname, 'rb')
    data = fobj.read(3)
    fobj.close()

    if data[:2] == b'\037\213':
        if PYVERSION == 2:
            return gzip.open(tmpname, 'r')
        else:
            return _wGzip(tmpname, 'rb')
    elif data == b'BZh':
        if PYVERSION == 2:
            return bz2.BZ2File(tmpname, 'r')
        else:
            return _wBzip(tmpname, 'rb')
    elif data[:2] == b'\037\235':
        return os.popen('gzip -dc %s' % tmpname)
    else:
        return open(tmpname, 'r')


def readfile(filename, splitlines=True):
    """Read data from a file.

    This function reads all data from

        -- file compressed using gzip, bzip2 or unix compress,
        -- uncompressed file

    Args:
        filename (str): Name of file or URL
        splitlines (bool): If specified and True then the data is split into
           lines

    Returns:
        A string containing all the data if *splitlines* is False (the default)
        or a list of lines if *splitlines* is True.

    """

    if filename.startswith('http://') or filename.startswith('ftp://') or \
           filename.startswith('file://'):
        tmpname = __geturl(filename)
    else:
        tmpname = filename

    fobj = openfile(tmpname)
    data = fobj.read()
    fobj.close()

    if tmpname != filename:
        os.unlink(tmpname)

    if splitlines:
        return data.splitlines(True)
    else:
        return data


def __geturl(url):

    # private method to fetch a file given an URL

    fd, tmpname = tempfile.mkstemp()
    fobj = urlopen(url)
    while 1:
        data = fobj.read(1024000)
        if not data:
            break
        os.write(fd, data)
    os.close(fd)
    fobj.close()
    return tmpname
