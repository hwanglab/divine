#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: complement
    :platform: Unix, Windows, MacOSX
    :synopsis: Nucleic Acid Complement

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com
The module exports one dictionary *COMPLEMENT*, whose keys are the
nucleic acid and the values are its complement.
"""


COMPLEMENT = {'A': 'T',
           'T': 'A',
           'G': 'C',
           'C': 'G',
           'a': 't',
           't': 'a',
           'c': 'g',
           'g': 'c',
           'N': 'N',
           'n': 'n    '
           }
