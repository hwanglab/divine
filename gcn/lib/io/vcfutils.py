#
# COPYRIGHT (C) 2002-2012 Rajgopal Srinivasan
#
"""
.. module:: vcfutils
    :platform: Unix, Windows, MacOSX
    :synopsis: Parsing VCF 4.1 format files

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Utilities for reading and writing  VCF 4.1 format files. The format itself
is documented at the 1000genomes project site
"""
from __future__ import absolute_import
#from ..dtypes.bunch import Bunch
import itertools
import sys
if sys.version_info.major == 3:
    izip = itertools.zip_longest
else:
    izip = itertools.izip_longest


def toint(v):
    """Check and convert a string to an integer"""
    try:
        return int(v) if v not in '.AG' else v
    except ValueError:
        raise ValueError('Could not convert %s to int' % v)


def tofloat(v):
    """Check and convert a string to a real number"""
    try:
        return float(v) if v != '.' else v
    except ValueError:
        raise ValueError('Could not convert %s to float' % v)


def tochar(v):
    """Check that a value is a single character string"""
    if isinstance(v, str) and len(v) == 1:
        return v
    else:
        raise ValueError("Expected a character. Got %s" % v)


TYPEMAP = {
    'float': tofloat,
    'integer': toint,
    'flag': bool,
    'string': lambda e: e,
    'character': tochar,
    }


def splitwords(line, sep=','):
    """Split a line into words using `sep`. Instances of `sep` within
    quotes are ignored as separators.

    >>> splitwords('ID=GQ,Number=1,Type=Float,Description="Genotype Quality"')
    ['ID=GQ', 'Number=1', 'Type=Float', 'Description=Genotype Quality']
    >>> splitwords('ID=GQ,Number=1,Type=Float,Description="Genotype, Quality"')
    ['ID=GQ', 'Number=1', 'Type=Float', 'Description=Genotype, Quality']
    """

    words = []
    inquote = 0
    prevc = w = ''
    for c in line:
        if c == sep:
            if prevc == '\\':
                w += c
                continue
            if not inquote:
                words.append(w)
                w = ''
            else:
                w += c
        elif c == '"' and prevc != '\\':
            if inquote:
                words.append(w)
                w = ''
                inquote = False
            else:
                inquote = True
        else:
            w += c
        prevc = c
    if w:
        words.append(w)

    return words


def parsemetainfo(line):
    """Parse an INFO meta line

    >>> s = '<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
    >>> infod = parsemetainfo(s)
    >>> infod['ID']
    'AF'
    >>> infod['Number']
    'A'
    >>> infod['Type'] == tofloat
    True
    >>> infod['Description']
    'Allele Frequency'
    """
    line = line.rstrip()
    if line[:1] != '<' or line[-1:] != '>':
        raise ValueError("%s:\nNot a valid INFO line" % line)

    infod = {}
    for item in splitwords(line[1:-1]):
        key, value = item.split('=', 1)
        if key.lower() == 'type':
            value = TYPEMAP[value.lower()]
        infod[key] = value

    for key in ('ID', 'Number', 'Type', 'Description'):
        if key not in infod:
            raise ValueError("%s:Not a valid INFO/FORMAT line. Missing %s key"
                             % line)
    return infod


def parsemeta(stream):
    """Parse all meta information lines in a VCF file stream"""

    DISPATCH = {
        'FORMAT': parsemetaformat,
        'INFO': parsemetainfo,
        'FILTER': parsemetafilter,
        'CONTIG': parsemetacontig,
        'SAMPLE': parsemetasample,
        'ALT': parsemetaalt,
        }

    DMETA = {}
    for key in DISPATCH:
        DMETA[key] = {}

    headerlines = []
    next = stream.readline
    while 1:
        line = next()
        if not line:
            break
        if line[:2] == '##':
            headerlines.append(line)
            key, value = line[2:].split('=', 1)
            if key.upper() in DISPATCH:
                d = DISPATCH[key.upper()](value)
                DMETA[key.upper()][d['ID']] = d
            else:
                if key in DMETA:
                    v = DMETA[key]
                    if not isinstance(v, list):
                        v = [v]
                        DMETA[key] = v
                    v.append(value)
                else:
                    DMETA[key] = value

        elif line[:6] == '#CHROM':
            args = line[1:].split('\t')
            args[-1] = args[-1].strip()
            if args[:8] != ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                            'FILTER', 'INFO']:
                print(args[:8])
                raise ValueError("Improper VCF File")
            if len(args) > 8:
                samplelist = args[9:]
            else:
                samplelist = []
            break
    return DMETA, samplelist, headerlines


def parsemetafilter(line):
    """Parse a FILTER meta line

    >>> s = '<ID=q10,Description="Quality below 10">'
    >>> mfilt = parsemetafilter(s)
    >>> mfilt['ID']
    'q10'
    >>> mfilt['Description']
    'Quality below 10'
    """
    line = line.rstrip()
    if line[:1] != '<' or line[-1:] != '>':
        raise ValueError("%s:Not a valid FILTER line" % line)

    mfilt = {'ID': None, 'Description': None}

    for item in splitwords(line[1:-1]):
        key, value = item.split('=', 1)
        if key not in ('ID', 'Description'):
            raise ValueError('Unknown key %s in FILTER record' % key)
        mfilt[key] = value

    return mfilt


def parsemetaformat(line):
    """Parse a FORMAT meta line

    >>> s = '<ID=GT,Number=1,Type=String,Description="Genotype">'
    >>> formd = parsemetaformat(s)
    >>> formd['ID']
    'GT'
    >>> formd['Number']
    '1'
    >>> formd['Description']
    'Genotype'
    """
    return parsemetainfo(line)


def _parsemetageneric(line, lfunc, checkid=False):
    # private method
    d = {}
    line = line.rstrip()
    if line[:1] != '<' or line[-1:] != '>':
        raise ValueError('Not a valid %s line' % lfunc)

    for item in splitwords(line[1:-1]):
        key, value = item.split('=', 1)
        d[key] = value

    if checkid:
        if 'ID' not in d:
            print(line)
            print(d)
            raise ValueError("ID field not present in %s meta record" % lfunc)

    return d


def parsemetaalt(line):
    """Parse an ALT meta line"""
    return _parsemetageneric(line, 'ALT', 1)


def parsemetacontig(line):
    """Parse a contig meta line.

    According to the VCF 4.1 specification a contig line has the
    following format

    ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

    However GATK typically writes out lines like so

    ##contig=<ID=3,length=198022430,assembly=b37>

    >>> dc = parsemetacontig('<ID=3,length=198022430,assembly=b37>')
    >>> dc['ID']
    '3'
    >>> dc['length']
    '198022430'
    >>> dc['assembly']
    'b37'
    """
    return _parsemetageneric(line, 'contig', 1)


def parsemetasample(line):
    """Parse a SAMPLE meta record

    >>> s = '<ID=Blood,Genomes=Germline,Mixture=1.,'
    >>> s += 'Description="Patient germline genome">'
    >>> ds = parsemetasample(s)
    >>> ds['ID']
    'Blood'
    >>> ds['Genomes']
    'Germline'
    >>> ds['Mixture']
    '1.'
    >>> ds['Description']
    'Patient germline genome'
    """
    return _parsemetageneric(line, 'SAMPLE', 1)


def parsevcfrecord(line):
    """Parse a VCF record into a dictionary"""

    rec = line.split('\t')
    rec[-1] = rec[-1].rstrip()
    return Bunch({
                  'chrom': rec[0],
                    'pos': int(rec[1]),
                     'id': rec[2].split(';'),
                    'ref': rec[3].strip().upper(),
                    'alt': rec[4].split(','),
                   'qual': tofloat(rec[5]),
                 'filter': rec[6].split(';'),
                  '_info': rec[7],
                 'format': rec[8:9],
             '_genotypes': rec[9:],
                '_numalt': rec[4].count(',') + 1,
              '_gtparsed': False,
            '_infoparsed': False,
                  '_line': line,
            })


def parseinfo(vcf, INFO):
    """Parse a VCF info string into a dictionary"""

    if vcf._infoparsed:
        return

    dinfo = Bunch()
    for item in vcf._info.split(';'):
        if item in INFO:
            dinfo[item] = True
            continue
        key, value = item.split('=', 1)
        try:
            IKEY = INFO[key]
            n = IKEY['Number']
            Type = IKEY['Type']
        except KeyError:
            raise KeyError("%s: Unknown INFO parameter" % key)
        if n == '1':
            dinfo[key] = Type(value)
            continue
        value = [Type(v) for v in value.split(',')]
        if n != '.':
            if n == 'A':
                n = vcf['_numalt']
            elif n == 'G':
                n = 3
            else:
                n = int(n)
            if len(value) != n:
                raise ValueError("Need %d values for %s. Got %d" % (
                        (n, key, len(value))))
        dinfo[key] = value
    vcf._infoparsed = True
    return dinfo


def parsegenotypes(vcf, FORMAT, sampleids=None):
    """Parse genotypes for samples"""
    if vcf._gtparsed:
        return
    if vcf.format:
        fields = vcf.format[0].split(':')

    if sampleids is None:
        sampleids = range(len(vcf._genotypes))

    gts = {}
    for idx in sampleids:
        gt = Bunch()
        gts[idx] = gt
        for name, val in izip(fields, vcf._genotypes[idx].split(':')):
            if val is None:
                continue
            n = FORMAT[name]['Number']
            if n == '1':
                gt[name] = FORMAT[name]['Type'](val)
                continue
            val = [FORMAT[name]['Type'](v) for v in val.split(',')]
            if n != '.':
                if n == 'A':
                    n = vcf._numalt
                elif n == 'G':
                    n = 3
                else:
                    n = int(n)
                if len(val) != n:
                    raise ValueError("Need %d values for %s. Got %d" % (
                            (n, name, len(val))))
            gt[name] = val
    vcf._gtparsed = 1
    return gts


def isclinical(vcf):
    """Sample should have CLN info tag which indicates that
    SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)
    """
    return vcf.info.CLN


def isnovel(vcf):
    """SNP should not have been reported in either DBSNP or 1000Genomes"""
    return not (vcf.info.DB135 or vcf.info.KGDB or vcf.info.DB)


def isexonic(vcf):
    info = vcf.info
    return info.AnnovarExonic or info.AnnovarExonicSplicing or \
        info.AnnovarSplicing
        #info.AnnovarNCRNAExonic or info.AnnovarNCRNAExonicSplicing


def isfunctional(vcf):
    info = vcf.info
    return not (info.AnnovarSyn or info.SYN)


def issnp(sample, mingq=10.0):
    return sample.GT in ('0/1', '1/1') and sample.GQ >= mingq


def ishomoalt(sample, mingq=10.0):
    return sample.GT == '1/1' and sample.GQ >= mingq


def ishomoref(sample, mingq=10.0):
    return sample.GT == '0/0' and sample.GQ >= mingq


def ishet(sample, mingq=10.0):
    return sample.GT == '0/1' and sample.GQ >= mingq


def isnocall(sample, mingq=10.0):
    return sample.GT == './.' or sample.GQ < mingq

def expand_variant(chrom, pos, ref, alt):
    '''Returns expanded form of the variant. This method was written
    because some VCFs generated by Haplotype Callers like FreeBayes
    has variants of all 3 types (MNP, Insertion, Deletion) represented
    as one entry in VCF.
    For example - a variant called in FreeBayes caller as - 
    chr1, pos=29018316, id=., ref=TTTTATTTA, alt=TTTTA,TTTTATTTATTTA,T
    comprised of insertion as well as deletion.
    We need to expand such entries in VCF many purposes like - 
        - to compare variants in two VCF files
        - to annotate variants present in other databases like dbSNP,
          1000Genome or ESP.
    Args:
        - chrom(str):    Chromosome Number
        - pos(int):    Variant position
        - ref(str):    Reference Allele
        - alt(str):    Alternate allele
    
    Returns:
        - explist(list):    List of expanded variants
    '''
    explist = []
    if len(ref) == len(alt):
        #SNP
        if len(ref) == 1:
            explist.append((chrom, pos, ref, alt))
        #MNP/Complex
        else:
            for r, a, p in zip (ref, alt, range(pos, pos + len(ref) + 1)):
                if r != a:
                    explist.append((chrom, p, r, a))
    #Deletion
    elif len(ref) > len(alt):
        r = ref[-(len(ref) - len(alt) + 1):]
        a = alt[-1]
        p = pos + len(alt) - 1
        explist.append((chrom, p, r, a))
    #Insertion
    elif len(ref) < len(alt):
        r = ref[-1]
        a = alt[-(len(alt) - len(ref) + 1):]
        p = pos + len(ref) - 1
        explist.append((chrom, p, r, a))
    return explist

def normalize_variant(chrom, pos, ref, alt):
    '''Returns normalized form of the variant. This method was written
    because some VCFs generated by Haplotype Callers like FreeBayes
    has variants of all 3 types (MNP, Insertion, Deletion) represented
    as one entry in VCF.
    For example - a variant called in FreeBayes caller as -
    chr1, pos=29018316, id=., ref=TTTTATTTA, alt=TTTTA,TTTTATTTATTTA,T
    comprised of insertion as well as deletion.
    We need to normalize such entries in VCF for many purposes like -
        - to compare variants in two VCF files
        - to annotate variants present in other databases like dbSNP,
          1000Genome or ESP.
    Args:
        - chrom(str):    Chromosome Number
        - pos(int):    Variant position
        - ref(str):    Reference Allele
        - alt(str):    Alternate allele

    Returns:
        - explist(list):    List of normalized variants
    '''
    explist = []
    if len(ref) == len(alt):
        #SNP
        if len(ref) == 1:
            explist.append((chrom, pos, ref, alt, 'SNP'))
        #MNP/Complex
        else:
            for r, a, p in zip(ref, alt, range(pos, pos + len(ref) + 1)):
                if r != a:
                    explist.append((chrom, p, r, a, 'SNP'))
    #Deletion
    elif len(ref) > len(alt):
        r = ref[len(alt) - 1:]
        a = alt[-1]
        p = pos + len(alt) - 1
        if ref[:len(alt) - 1] == alt[:-1] and r[0] == a[0]:
            explist.append((chrom, p, r, a, 'DEL'))
        else:
            explist.append((chrom, p, r, a, 'COMPLEX_DEL'))
    #Insertion
    elif len(ref) < len(alt):
        r = ref[-1]
        a = alt[len(ref) - 1:]
        p = pos + len(ref) - 1
        if ref[:-1] == alt[:len(ref) - 1] and r[0] == a[0]:
            explist.append((chrom, p, r, a, 'INS'))
        else:
            explist.append((chrom, p, r, a, 'COMPLEX_INS'))
    return explist

def get_var_type(refNt,altNts):
  '''
  by looking ref and alt, the variant type will be determined!
  '''
  R = len(refNt)
  vtype = None
  if len(altNts) > 1:
    As = [len(alt) for alt in altNts]
    Aprev = As[0]
    for A in As:
      if A != Aprev:
        vtype = 'complex'
        break
  else:
    A = len(altNts[0])

  if not vtype:
    if R == 1:
      if A == 1:
        vtype = 'mismatch'
      elif R > A:
        vtype = 'del'
      elif R < A:
        vtype = 'ins'
    elif R == A:
      vtype = 'mnp'

  return vtype

if __name__ == "__main__":
    import doctest
    doctest.testmod()