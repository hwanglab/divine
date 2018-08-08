#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: entrystream
    :platform: Unix, Windows, MacOSX
    :synopsis: Iterate over each entry in a file containing multiple entries

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com


This class can be used for those multi-entry files where each entry
is separated by a marker and the marker occurs at the *beginning* of a line.
The class can also handle those entries with distinct start of entry and
end of entry markers.  Finally, the module defines convenience functions for
iterating over several of the common biological and chemical file formats.
"""


class EntryStream(object):
    """
    Class for iterating over entries in a multi-entry filestream.

    Most biological sequence files contain multiple entries per file,
    for e.g Genbank, Swissprot, EMBL etc.  This class can be used to
    iterate over successive entries in such a file.  Each entry is
    returned as a single string (including end-of-line markers).

    **NOTE** The tags/markers that define the beginning and end of a line
    must appear at the start of a new line.

    Examples
    --------

    To iterate over each sequence in a FASTA format file::

    >>> stream = open('test.seq', 'rb')
    >>> seqstream = EntryStream(stream, '>')
    >>> n = 0
    >>> for entry in seqstream:
    ...     n = n + 1
    ...     # do something with entry
    >>> print 'Total number of entries in file =', n

    Similarly, to iterate over each sequence in a SWISS-PROT file::

    >>> stream = open('test.seq', 'rb')
    >>> seqstream = EntryStream(stream, 'ID', '//')
    >>> n = 0
    >>> for entry in seqstream:
    ...     n = n + 1
    ...     # do something with entry
    >>> print 'Total number of entries in file =', n

    Alternatively one can also iterate in the following manner

    >>> stream = open('test.seq', 'rb')
    >>> seqstream = EntryStream(stream, 'ID', '//')
    >>> n = 0
    >>> while 1:
    ...     try:
    ...         e = seqstream.next()
    ...     except StopIteration:
    ...         break
    ...     else:
    ...         n = n + 1
    ...         # do something with entry
    >>> print 'Total number of entries in file =', n

    """

    __slots__ = ['_fh', '_st', '_et', '_unget', 'shared_data', '_buffer']

    def __init__(self, filehandle, starttag, endtag=None, shareddata=False):
        """Create an instance of the EntryStream class

        Args:
            filehandle: A file object opened for reading.
            starttag(str): The marker that flags the beginning of an entry
            endtag(str): The marker that flags the end of an entry. This is
                optional and if omitted the *startag* will be used as the
                sole delimiter between entries
            shareddata(bool): If specified and True, then all all lines
                prior to the first entry will be saved as a list in the
                instance variable *shared_data*

        Returns:
            An instance of the EntryStream class
        """

        if endtag is None:
            endtag = starttag
            self._unget = 1
        else:
            self._unget = 0

        self._fh = filehandle
        self._st = starttag
        self._et = endtag
        line = self._scroll()
        if not line:
            raise StopIteration('File has no data')

        self._buffer = line

    def __iter__(self):
        # advance to the first entry
        return self

    def next(self):
        """Fetch the next entry in the file stream

        Returns:
            A list of lines corresponding to the entry

        Raises:
            StopIteration if end of file is reached.
        """

        ET = self._et
        ETLEN = len(ET)
        ST = self._st
        STLEN = len(ST)

        if self._buffer:
            entry = [self._buffer]
            self._buffer = ''
        else:
            entry = None

        for line in self._fh:
            if line[:ETLEN] == ET:
                if self._unget:
                    self._buffer = line
                    return entry
                else:
                    entry.append(line)
                    self._buffer = ''
                    return entry
            elif line[:STLEN] == ST:
                entry = [line]
            else:
                try:
                    entry.append(line)
                except AttributeError:
                    pass
        if entry:
            return entry

        raise StopIteration('End of file reached')

    __next__ = next

    def _scroll(self):
        sd = self.shared_data = []
        for line in self._fh:
            if line.startswith(self._st):
                return line
            else:
                sd.append(line)
        else:
            return ''


def FastaStream(filehandle):
    """Iterate over entries in a Fasta formatted file

    Args:
        filehandle: A file like object containing Fasta formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, '>')


def UniprotStream(filehandle):
    """Iterate over entries in a Uniprot formatted file

    Args:
        filehandle: A file like object containing Uniprot format
             entries that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

   """

    return EntryStream(filehandle, 'ID', '//')


def SwissprotStream(filehandle):
    """Iterate over entries in a SWISS_PROT formatted file

    Args:
        filehandle: A file like object containing Swissprot formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'ID', '//')


def TaxonomyStream(filehandle):
    """Iterate over entries in a TAXONOMY formatted file

    Args:
        filehandle: A file like object containing TAXONOMY formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

   """

    return EntryStream(filehandle, 'ID', '//')


def EnzymeStream(filehandle):
    """Iterate over entries in a ENZYME formatted file

    Args:
        filehandle: A file like object containing ENZYME formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

   """

    return EntryStream(filehandle, 'ID', '//')


def dbSNPStream(filehandle):

    """Iterate over entries in a dbSNP formatted file

    Args:
        filehandle: A file like object containing dbSNP formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'rs', '\n')


def GenbankStream(filehandle):
    """Iterate over entries in a Genbank flat text formatted file

    Args:
        filehandle: A file like object containing Genbank formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'LOCUS', '//')


def PrositeStream(filehandle):
    """Iterate over entries in a Prosite file

    Args:
        filehandle: A file like object containing PROSITE formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'ID', '//')


def PDOCStream(filehandle):
    """Iterate over entries in a Prosite Documentation file

    Args:
        filehandle: A file like object containing Prosite entry
            documentation formatted files
                       that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, r'{PDOC[0-9]{5,5}}', r'{END}')


def EMBLStream(filehandle):
    """Iterate over entries in a EMBL formatted file

    Args:
        filehandle: A file like object containing EMBL formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'ID', '//')


def PIRStream(filehandle):
    """Iterate over entries in a PIR CODATA formatted file

    Args:
        filehandle: A file like object containing PIR CODATA formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'ENTRY', '///')


def KEGGLigand(filehandle):
    """Iterate over entries in a KEGG LIGAND file

    Args:
        filehandle: A file like object containing KEGG formatted ligand
            files that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'ENTRY', '///')


def CIFStream(filehandle):
    """Iterate over entries in an mmcif file

    Args:
        filehandle: A file like object containing (MM)CIF formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'data_')


def OMIMTextStream(filehandle):
    """Iterate over entries in an OMIM text file

    Args:
        filehandle: A file like object containing OMIM formatted files
            that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

   """

    return EntryStream(filehandle, '\*RECORD\*')


def RebaseStream(filehandle):
    """Iterate over entries in a Rebase bairoch format file

    Args:
        filehandle: A file like object containing REBASE Bairoch format
            entries that has been opened for reading.

    Returns:
        An itertator that can be used to access each entry in the file.

    """

    return EntryStream(filehandle, 'ID', '//')


def AdvRebaseStream(filehandle):
    """Class for iteration over entries in a Rebase bairoch format file"""

    return EntryStream(filehandle, 'TY', 'ER')


def PDBHetStream(filehandle):
    """Iterate over entries in the pdb distributed het_dictionary.txt file"""

    return EntryStream(filehandle, 'RESIDUE')


def MDLStream(filehandle):
    """Iterate over entries in an MDL Mol File"""
    return EntryStream(filehandle, '\n|[A-Z,a-z,0-9]+', 'M  END')


def SDStream(filehandle):
    """Iterate over entries in an MDL SD file"""
    return EntryStream(filehandle, '\n|[A-Z,a-z,0-9]+', '\$\$\$\$')


if __name__ == "__main__":
    import anyopen
    import sys
    from time import time
    stream = anyopen.openfile(sys.argv[1])
    es = SwissprotStream(stream)
    n = 0
    t0 = time()
    elist = []
    for e in es:
        elist.append(e)
        if len(elist) == 10000:
            n += 10000
            print('Read %d entries in %s seconds' % (n, time() - t0))
            elist = []
    n += len(elist)
    print(n)

    for line in e:
        print(line),
