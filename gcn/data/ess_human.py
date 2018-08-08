"""
.. module:: ese_human
    :platform: Unix, Windows, MacOSX
    :synopsis: Exon Splicing Enhancer Sites

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

Exon Splicing Silencer Sites predicted for human genes. This data is derived
from http://genes.mit.edu/fas-ess/fas-hex3.txt

The module exports one list *ESS_MOTIFS*, which has 103 hexamers exonic
splicing silencer motifs of human genes and are likely to be highly enriched
hexamers that represent core portions of ESS motifs.
Ref Paper :
https://grail.cs.washington.edu/education/courses/cse590c/05wi/wang.pdf
"""

ESS_MOTIFS = ['ACTAGG', 'AGGGGG', 'AGGGGT', 'AGGGTA', 'AGGGTG', 'AGGGTT',
              'AGGTAA', 'AGGTAG', 'AGGTAT', 'AGGTTA', 'AGGTTT', 'AGTAGG',
              'AGTCGT', 'AGTGTA', 'AGTGTT', 'AGTTAG', 'AGTTCG', 'AGTTCT',
              'AGTTGT', 'AGTTTA', 'AGTTTG', 'ATGGGG', 'CCTGGG', 'CGTTCC',
              'CTAGGG', 'CTTAGG', 'GATGGG', 'GGATGG', 'GGGAGG', 'GGGATG',
              'GGGGAG', 'GGGGAT', 'GGGGGA', 'GGGGGG', 'GGGGGT', 'GGGGTA',
              'GGGGTG', 'GGGGTT', 'GGGTAA', 'GGGTAG', 'GGGTGG', 'GGGTTA',
              'GGGTTT', 'GGTAAG', 'GGTAGG', 'GGTAGT', 'GGTATA', 'GGTGGG',
              'GGTTAG', 'GGTTTA', 'GGTTTG', 'GGTTTT', 'GTAAGG', 'GTAAGT',
              'GTAGGG', 'GTAGGT', 'GTAGTT', 'GTATAG', 'GTCGTT', 'GTGGGG',
              'GTGGGT', 'GTGTAG', 'GTTAGG', 'GTTAGT', 'GTTCCT', 'GTTCGT',
              'GTTCTG', 'GTTGTT', 'GTTTAG', 'GTTTGG', 'GTTTGT', 'TAAGTG',
              'TAGGGG', 'TAGGGT', 'TAGGTA', 'TAGGTT', 'TAGTGG', 'TAGTGT',
              'TAGTTA', 'TAGTTT', 'TATAGG', 'TCCTGG', 'TCGTTC', 'TCTAGG',
              'TGGGGG', 'TGGGTG', 'TGGTGG', 'TGTTAG', 'TGTTCC', 'TTAGGG',
              'TTAGGT', 'TTAGTG', 'TTAGTT', 'TTCCTG', 'TTCGTT', 'TTCTAG',
              'TTGGGG', 'TTGTTC', 'TTTAGG', 'TTTAGT', 'TTTGGG', 'TTTGTT',
              'TTTTGG']
