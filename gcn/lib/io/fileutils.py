#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: fileutils
    :platform: Unix, Windows, MacOSX
    :synopsis: Various utility functions

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com

Various utility functions that don't fit anywhere else
"""
import os
import time

def file_time(filename, format='i'):
    """Get last modified time of a file

    Args:

        filename (str): Name of the file
        format (str):   Choice of 'i' (seconds), 's' (nicely formatted date),
                        or 't' a nine tuple

    Returns:
        The last modified time of the file in the specified format

    Raises:
        IOError if the file does not exist
    """
    ftime = os.stat(filename).st_ctime
    if format == 't':
        ftime = time.localtime(ftime)
    elif format == 's':
        ftime = time.ctime(ftime)
    return ftime


def file_newer(file1, file2):
    """Check if a file is newer than another

    Args:
        file1 (str): Name of the first file
        file2 (str): Name of the second file

    Returns:
        True is file1 has a later modification time than file2 and False
        otherwise

    """
    return file_time(file1) > file_time(file2)


def normfile(filename):
    """Return the absolute path of a filename"""
    return os.path.abspath(os.path.expanduser(os.path.normpath(filename)))

