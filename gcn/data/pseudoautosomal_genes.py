"""
.. module:: pseudoautosomal_genes
    :platform: Unix, Windows, MacOSX
    :synopsis: Genes in Pseudoautosomal Region

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

The pseudoautosomal regions, PAR1, PAR2 and PAR3 are homologous
sequences of nucleotides on the X and Y chromosomes.
The pseudoautosomal regions get their name because any genes located
within them (so far at least 29 have been found) are inherited just
like any autosomal genes.

Inheritance and function
------------------------
Normal male mammals have two copies of these genes: one in the
pseudoautosomal region of their Y chromosome, the other in the
corresponding portion of their X chromosome. Normal females also
possess two copies of pseudoautosomal genes, as each of their two
X chromosomes contains a pseudoautosomal region. Crossing over between
the X and Y chromosomes is normally restricted to the pseudoautosomal
regions; thus, pseudoautosomal genes exhibit an autosomal, rather than
sex-linked, pattern of inheritance. So, females can inherit an allele
originally present on the Y chromosome of their father and males can
inherit an allele originally present on the X chromosome of their father.

The function of these pseudoautosomal regions is that they allow the X and
Y chromosomes to pair and properly segregate during meiosis in males.

Pseudoautosomal genes are found in three different locations: PAR1, PAR2 and
PAR3.
These are believed to have evolved independently.
Ref - http://en.wikipedia.org/wiki/Pseudoautosomal_region
"""
PSEUDO_AUTO_GENES = ['ASMT', 'ASMTL', 'CD99', 'CRLF2', 'CSF2RA',
                     'SFRS17A', 'DHRSXY', 'GTPBP6', 'IL3RA', 'P2RY8',
                     'PLCXD1', 'PPP2R3B', 'SHOX', 'SLC25A6', 'XG',
                     'ZBED1', 'SPRY3', 'SYBL1', 'IL9R', 'PCDH11X',
                     'PCDH11Y', 'TGIF2LX', 'TGIF2LY']