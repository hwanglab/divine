#
# COPYRIGHT (C) 2011 Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com
#
"""
.. module:: cut
    :platform: Unix, Windows, MacOSX
    :synopsis: Codon Usage Table

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com

Codon usage table for Homo Sapiens. The table is derived from

http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606

The module exports one dictionary *CUT*, whose keys are the codons and
the values are the usage frequency per 1000 codons.
"""

CUT = {
    'AAA': 24.4, 'AAC': 19.1, 'AAG': 31.9, 'AAT': 17.0,
    'ACA': 15.1, 'ACC': 18.9, 'ACG':  6.1, 'ACT': 13.1,
    'AGA': 12.2, 'AGC': 19.5, 'AGG': 12.0, 'AGT': 12.1,
    'ATA':  7.5, 'ATC': 20.8, 'ATG': 22.0, 'ATT': 16.0,
    'CAA': 12.3, 'CAC': 15.1, 'CAG': 34.2, 'CAT': 10.9,
    'CCA': 16.9, 'CCC': 19.8, 'CCG':  6.9, 'CCT': 17.5,
    'CGA':  6.2, 'CGC': 10.4, 'CGG': 11.4, 'CGT':  4.5,
    'CTA':  7.2, 'CTC': 19.6, 'CTG': 39.6, 'CTT': 13.2,
    'GAA': 29.0, 'GAC': 25.1, 'GAG': 39.6, 'GAT': 21.8,
    'GCA': 15.8, 'GCC': 27.7, 'GCG':  7.4, 'GCT': 18.4,
    'GGA': 16.5, 'GGC': 22.2, 'GGG': 16.5, 'GGT': 10.8,
    'GTA':  7.1, 'GTC': 14.5, 'GTG': 28.1, 'GTT': 11.0,
    'TAA':  1.0, 'TAC': 15.3, 'TAG':  0.8, 'TAT': 12.2,
    'TCA': 12.2, 'TCC': 17.7, 'TCG':  4.4, 'TCT': 15.2,
    'TGA':  1.6, 'TGC': 12.6, 'TGG': 13.2, 'TGT': 10.6,
    'TTA':  7.7, 'TTC': 20.3, 'TTG': 12.9, 'TTT': 17.6,
    }
