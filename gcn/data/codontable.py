"""
#
# COPYRIGHT (C) 2011 Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com
#
.. module:: codontable
    :platform: Unix, Windows, MacOSX
    :synopsis: Codon Table from the NCBI site

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com

    Codon table from the NCBI site.

    This module exports one dictionary, *CODONTABLE*.  The keys in
    this dictionary are as follows:

    o name - name of the codon table

    o altnames - other names for this codon table

    o mnemonic - a short string mnemonic for the table

    o codons - a dictionary containing as key the codon triplet and
               as value the amino acid one letter symbol. This includes
               codons with ambiguous bases

    o unambiguous_codons - subset of the codons dictionary with the
                           restriction that every codon triplet is
                           unambiguous

    o start_codons - a list of all codons that are start sites for
                     for protein translation

    o stop_codons - a list of all codons that terminate protein
                    translation.
"""
CODONTABLE = {
    'altnames': [],
    'codons': {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAR': 'K',
               'AAT': 'N', 'AAY': 'N', 'ACA': 'T', 'ACB': 'T',
               'ACC': 'T', 'ACD': 'T', 'ACG': 'T', 'ACH': 'T',
               'ACK': 'T', 'ACM': 'T', 'ACN': 'T', 'ACR': 'T',
               'ACS': 'T', 'ACT': 'T', 'ACV': 'T', 'ACW': 'T',
               'ACY': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
               'AGR': 'R', 'AGT': 'S', 'AGY': 'S', 'ATA': 'I',
               'ATC': 'I', 'ATG': 'M', 'ATH': 'I', 'ATM': 'I',
               'ATT': 'I', 'ATW': 'I', 'ATY': 'I', 'CAA': 'Q',
               'CAC': 'H', 'CAG': 'Q', 'CAR': 'Q', 'CAT': 'H',
               'CAY': 'H', 'CCA': 'P', 'CCB': 'P', 'CCC': 'P',
               'CCD': 'P', 'CCG': 'P', 'CCH': 'P', 'CCK': 'P',
               'CCM': 'P', 'CCN': 'P', 'CCR': 'P', 'CCS': 'P',
               'CCT': 'P', 'CCV': 'P', 'CCW': 'P', 'CCY': 'P',
               'CGA': 'R', 'CGB': 'R', 'CGC': 'R', 'CGD': 'R',
               'CGG': 'R', 'CGH': 'R', 'CGK': 'R', 'CGM': 'R',
               'CGN': 'R', 'CGR': 'R', 'CGS': 'R', 'CGT': 'R',
               'CGV': 'R', 'CGW': 'R', 'CGY': 'R', 'CTA': 'L',
               'CTB': 'L', 'CTC': 'L', 'CTD': 'L', 'CTG': 'L',
               'CTH': 'L', 'CTK': 'L', 'CTM': 'L', 'CTN': 'L',
               'CTR': 'L', 'CTS': 'L', 'CTT': 'L', 'CTV': 'L',
               'CTW': 'L', 'CTY': 'L', 'GAA': 'E', 'GAC': 'D',
               'GAG': 'E', 'GAR': 'E', 'GAT': 'D', 'GAY': 'D',
               'GCA': 'A', 'GCB': 'A', 'GCC': 'A', 'GCD': 'A',
               'GCG': 'A', 'GCH': 'A', 'GCK': 'A', 'GCM': 'A',
               'GCN': 'A', 'GCR': 'A', 'GCS': 'A', 'GCT': 'A',
               'GCV': 'A', 'GCW': 'A', 'GCY': 'A', 'GGA': 'G',
               'GGB': 'G', 'GGC': 'G', 'GGD': 'G', 'GGG': 'G',
               'GGH': 'G', 'GGK': 'G', 'GGM': 'G', 'GGN': 'G',
               'GGR': 'G', 'GGS': 'G', 'GGT': 'G', 'GGV': 'G',
               'GGW': 'G', 'GGY': 'G', 'GTA': 'V', 'GTB': 'V',
               'GTC': 'V', 'GTD': 'V', 'GTG': 'V', 'GTH': 'V',
               'GTK': 'V', 'GTM': 'V', 'GTN': 'V', 'GTR': 'V',
               'GTS': 'V', 'GTT': 'V', 'GTV': 'V', 'GTW': 'V',
               'GTY': 'V', 'MGA': 'R', 'MGG': 'R', 'MGR': 'R',
               'RAC': 'B', 'RAT': 'B', 'RAY': 'B', 'SAA': 'Z',
               'SAG': 'Z', 'SAR': 'Z', 'TAA': '*', 'TAC': 'Y',
               'TAG': '*', 'TAR': '*', 'TAT': 'Y', 'TAY': 'Y',
               'TCA': 'S', 'TCB': 'S', 'TCC': 'S', 'TCD': 'S',
               'TCG': 'S', 'TCH': 'S', 'TCK': 'S', 'TCM': 'S',
               'TCN': 'S', 'TCR': 'S', 'TCS': 'S', 'TCT': 'S',
               'TCV': 'S', 'TCW': 'S', 'TCY': 'S', 'TGA': '*',
               'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 'TGY': 'C',
               'TRA': '*', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L',
               'TTR': 'L', 'TTT': 'F', 'TTY': 'F', 'YTA': 'L',
               'YTG': 'L', 'YTR': 'L'},
     'mnemonic': 'SGC0',
     'name': 'standard',
     'start_codons': ['TTG', 'CTG', 'ATG'],
     'stop_codons': ['TAA', 'TAG', 'TGA'],
     'unambiguous_codons': {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
                            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                            'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                            'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
                            'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
                            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                            'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                            'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
                            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                            'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                            'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'},
    }

for key, value in CODONTABLE['codons'].items():
    CODONTABLE['codons'][key.lower()] = value.lower()

for key, value in CODONTABLE['unambiguous_codons'].items():
    CODONTABLE['unambiguous_codons'][key.lower()] = value

