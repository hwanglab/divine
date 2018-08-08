"""
.. module:: bed
    :platform: Unix, Windows, MacOSX
    :synopsis: Read BED files

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Parser for BED (Browser Extensible Data)files. A detailed description of
the format can be found at the UCSC_ site.

.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
"""
import anyopen
import os

__all__ = ['BEDRecord', 'readfile', 'parseline']

CHROMOSOMES = {}


class BEDRecord:
    """Class to represent a single record in a Bed file"""

    __slots__ = ['chrom', 'start', 'end', 'name', 'score',
                 'strand', 'thickStart', 'thickEnd', 'itemRGB',
                 'blockCount', 'blockSizes', 'blockStarts']

    def __init__(self, chrom='', start=0, end=None, name=None, score=None,
                 strand=None, thickStart=None, thickEnd=None,
                 itemRGB=None, blockCount=None, blockSizes=None,
                 blockStarts=None):

        self.chrom = chrom
        if chrom not in CHROMOSOMES:
            CHROMOSOMES[chrom] = len(CHROMOSOMES) + 1
        self.start = start
        if end is None:
            end = start + 1
        if self.start < 0 or self.start >= end:
            raise ValueError('BED start should be positive and greater ' + \
                                 'than BED end')

        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRGB = itemRGB
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts

    def __str__(self):
        v = []
        for el in self.__slots__:
            w = getattr(self, el)
            if w is None:
                w = ''
            v.append(str(w))
        return '\t'.join(v)

    def __lt__(self, other):
        """Check if this record's start is less than that of the other"""

        if other.chrom == self.chrom:
            return self.start < other.start

        return CHROMOSOMES[self.chrom] < CHROMOSOMES[other.chrom]

    def __repr__(self):
        v = [self.__class__.__name__ + '(']
        for el in self.__slots__:
            w = getattr(self, el)
            if isinstance(w, str):
                w = "'" + w + "'"
            v.append('%s=%s, ' % (el, w))
        v[-1] = v[-1][:-2] + ')'
        return ''.join(v)

    def __contains__(self, other):
        """Check if another interval is wholly contained within this"""
        return (self.chrom == other.chrom) and \
               (self.start <= other.start) and \
               (self.end >= other.end)

    def __eq__(self, other):
        """Is this interval the same as another"""
        return (self.chrom == other.chrom) and \
               (self.start == other.start) and \
               (self.end == other.end)

    def overlaps(self, other):
        """Check for overlap between two intervals"""

        if self.chrom != other.chrom:
            return False
        if other.start >= self.start and other.start < self.end:
            return True
        if other.end >= self.start and other.end <= self.end:
            return True
        if other.start <= self.start and other.end >= self.end:
            return True

        return False

    def merge(self, other, check=True, intersect=0, max=max, min=min):
        """Merge 2 overlapping intervals to a single interval.
        Raises TypeError if intervals do not overlap.
        """
        if check:
            if not self.overlaps(other):
                raise TypeError('Non-overlapping intervals cannot be merged')
        if intersect:
            start = max(self.start, other.start)
            end = min(self.end, other.end)
        else:
            start= min(self.start, other.start)
            end = max(self.end, other.end)

        return self.__class__(self.chrom, start, end)

    def __unpack(self, other):
        """private method to unpack an argument or a BedRecord to
        chromosome, start and end locations"""
        try:
            chrom = other.chrom
            start = other.start
            end = other.end
        except AttributeError:
            if len(other) == 3:
                chrom, start, end = other
            elif len(other) == 2:
                chrom, start = other
                end = start
            else:
                raise ValueError('Expected a sequence of BedRecord')
        return chrom, start, end

    def write(self, stream):
        """Write this record to a stream opened for writing"""
        stream.write(self.__str__())
        stream.write(os.newline)

    def hasregion(self, region):
        """Check if a region, specified in samtools format, is within this
        interval"""
        c, s, e = _parseregion(region)
        return (c == self.chrom) and (s >= self.start) and (s < self.end)



def readfile(filename):
    """Iterate over the records in a BED file"""
    stream = anyopen.openfile(filename)
    line = stream.readline()
    if isinstance(line, bytes):
        rec = parseline(line.decode())
        if rec is not None:
            yield rec
        for line in stream:
            rec = parseline(line.decode())
            if rec is not None:
                yield rec
    else:
        rec = parseline(line)
        if rec is not None:
            yield rec
        for line in stream:
            rec = parseline(line)
            if rec is not None:
                yield rec


def parseline(line):
    """Parse a single BED record into an instance of BEDRecord"""

    args = line.strip().split('\t')
    if not args:
        return None

    # ignore lines that are comment lines or meant for annotation
    if not valid_chrom(args[0]):
        return None

    if len(args) < 2:
        raise ValueError('A valid BED record should have at least 2 fields')

    v = tuple(FIELD_MAP[idx](arg) for (idx, arg) in enumerate(args))
    return BEDRecord(*v)


def fromregion(region):
    """Create a BEDRecord from a samtools type region string, e.g.
    chr1:1212232-1221324. *Note* that the start and end positions
    are assumed to be 1 based and part of a closed interval."""

    c, s, e = _parseregion(region)
    return BEDRecord(c, s, e)


def _parseregion(region):
    """private method to parse a samtools style region"""
    pos = region.find(':')
    if pos < 0:
        raise ValueError('Improper region specified')
    chrom = region[:pos]
    interval = region[pos + 1:]
    pos = region.find('-')
    if pos < 0:
        start = int(interval) - 1
        if start < 0:
            raise ValueError('Start of a region should be 1 or higher')
        end = start + 1
    else:
        start = int(interval[:pos])
        end = int(interval[pos + 1:])

    return chrom, start, end


def no_op(arg):
    """Dummy do nothing function"""
    return arg


def tochr(arg):
    """Check if argument is a 1 character string containing + or - or !"""
    if len(arg) == 1 and arg in '+-!':
        return arg
    raise ValueError('Invalid value for strand: %s' % arg)


def rgb(arg):
    """Convert a comma separate string of 3 integers to RGB values"""
    colors = arg.split(',')
    if len(colors) == 1:
        try:
            r = int(colors[0])
        except ValueError:
            raise ValueError('Colour should be a single integer or ' + \
                             '3 comma separated integers')
        else:
            return r, r, r
    if len(colors) != 3:
        raise ValueError('Colour should be specified as 3 integers')
    try:
        colors = tuple([int(e) for e in colors])
    except ValueError:
        raise ValueError('Colour should be a single integer or ' + \
                         '3 comma separated integers')
    else:
        return colors


def int_tuple(arg):
    """Convert a string of comma separated numbers to a tuple of integers"""
    args = arg.split(',')
    try:
        return tuple([int(e) for e in args])
    except ValueError:
        raise ValueError('Expected a comma separated list of integers, ' +
                         'found %s' % arg)


def valid_chrom(arg):
    """Check if a string represents a valid chromosome"""
    arg = arg.lower()
    if arg[:3] == 'chr':
        return True
    if arg.isdigit() or arg in ('x', 'y', 'mt'):
        return True
    if arg[:8] == 'scaffold' or arg[:2] == 'gl' or arg[:3] == 'par':
        return True
    return False

FIELD_MAP = {
             0: no_op, 1: int, 2: int, 3: no_op, 4: float,
             5: tochr, 6: int, 7: int, 8: rgb, 9: int,
             10: int_tuple, 11: int_tuple,
            }

if __name__ == "__main__":
    import sys
    bedfile = sys.argv[1]
    for line in anyopen.open(bedfile):
        v = parseline(line)

    print(repr(v))
    print(v)
