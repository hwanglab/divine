#
# COPYRIGHT (C) 2013 Rajgopal Srinivasan
#
"""
.. module:: phase
    :platform: Unix, Windows, MacOSX
    :synopsis: Phase the genotype of a child given parents genotypes

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

This module exports 4 funtions:

        o phase(p1, p2, child): phase the genotype of a child given the
                                genotypes of both parents
        o phase1(p1, child, idx): phase the genotype of a child given the
                                  genotype of one parent
        o phaseX(mother, child): phase X chromosome genotype of a male child
                                 using mother's genotype
        o phaseY(father, child): phase Y chromosome genotype of a male child
                                 using father's genotype
"""


__all__ = ['phase', 'phase1', 'phaseX', 'phaseY']

# GTYPES holds all possible genotypes of a child given the genotypes of parent
# and child. This is a list of length 256
GTYPES = []


def set_gtypes():
    """Compute the set of all possible genotypes"""
    if GTYPES:
        return
    # all possible genotypes
    diploids = [(i, j) for i in '0123456789' for j in '0123456789']

    # all possible inheritable genotypes
    for i, j in diploids:
        for k, l in diploids:
            GTYPES.append(set(((i, k), (i, l), (j, k), (j, l))))


def get_gtypes(p1, p2, int=int):
    """Given the genotypes of two parents, get all possible genotypes
    for the child"""
    #a, _, b = p1
    a, b = p1.split('/')
    #c, _, d = p2
    c,  d = p2.split('/')
    return GTYPES[int(d) + (int(c) * 10) + (int(b) * 100) + (int(a) * 1000)]


def phase(p1, p2, child):
    """Phase a child's genotype given the genotypes of the parents.

    Arguments:
        p1 (string): Genotype of parent encoded in the form 0/1, 1/2 etc
        p2 (string): Genotype of the other parent
        child (string): Genotype of the child

    Returns:
        Phased genotype of the child.  If the genotype of one of the parents
        is unavailable (indicated by ./.) then an attempt is made to phase
        using information from only one parent. The allele corresponding to
        p1 appears to the left of the | symbol and that of p2 to the right.
        If the phased genotype is not unambiguously resolvable then the
        genotype is left unphased

    >>> phase('0/1', '0/1', '1/1')
    '1|1'
    >>> phase('0/1', '0/1', '0/1')
    '0/1'
    >>> phase('0/0', '0/1', '0/1')
    '0|1'
    >>> phase('1/2', '0/1', '0/2')
    '2|0'
    >>> phase('./.', './.', '0/2')
    '0/2'
    >>> phase('./.', '0/1', '1/2')
    '2|1'
    """
    if child == './.':
        gt = child.split('/')
    elif p1 == './.':
        gt = phase1(p2, child, 1)
    elif p2 == './.':
        gt = phase1(p1, child, 0)
    else:
        #ca, _, cb = child
        ca, cb = child.split('/')
        gtypes = get_gtypes(p1, p2)
        sep = '/'
        if (ca, cb) in gtypes:
            if (cb, ca) in gtypes:
                if cb == ca:
                    sep = '|'
            else:
                sep = '|'
        elif (cb, ca) in gtypes:
            ca, cb = cb, ca
            sep = '|'
        gt = ca + sep + cb
    return gt


def phase1(p1, child, idx):
    """Phase a genotype given only one parents genotype.

    Arguments:
        p1 (string): genotype of parent  expressed in form 0/1, 0/0 etc
        child (string): genotype of child
        idx (integer): 0 or 1. If 0 then the allele inherited from the parent
                       appears to the left of | symbol in the phased genotype,
                       and to the right of the | symbol if idx is 1
    Returns:
        The phased genotype of the child if it is unambiguously resolvable.
        Otherwise the genotype is returned as is, i.e. unphased

    >>> phase1('0/1', '0/1', 0)
    '0/1'
    >>> phase1('0/2', '1/2', 0)
    '2|1'
    >>> phase1('0/3', '1/3', 1)
    '1|3'
    >>> phase1('./.', '0/1', 1)
    '0/1'
    >>> phase1('./.', './.', 0)
    './.'
    """
    #ca, _, cb = child
    ca, cb = child.split('/')
    if ca == '.':
        sep = '/'
    elif ca == cb and ca in p1:
        sep = '|'
    elif ca in p1 and cb not in p1:
        sep = '|'
    elif cb in p1 and ca not in p1:
        cb, ca = ca, cb
        sep = '|'
    else:
        sep = '/'
    if sep == '|' and idx == 1:
        ca, cb = cb, ca
    return ca + sep + cb


def phaseX(mother, child):
    """Phase X chromosome genotype of a male child given the genotype of the
    mother.  This function only checks that the genotype of the child
    could have been inherited from the mother.

    Arguments:
        mother (string): genotype of mother expressed in the form 0/1, 0/0 etc
        child (string): genotype of child

    Returns:
        The phased genotype of the child if it is unambiguously resolvable.
        Otherwise the genotype is returned as is, i.e. unphased

    >>> phaseX('0/1', '0/1')
    '0/1'
    >>> phaseX('0/1', '1/1')
    '1|1'
    >>> phaseX('0/2', '1/1')
    '1/1'
    >>> phaseX('./.', '3/3')
    '3/3'
    """
    #ca, _, cb = child.split('/')
    ca, cb = child.split('/')
    if ca != cb:  # we have X chromosome heterozygosity
        return child
    if ca in mother:
        return '%s|%s' % (ca, ca)
    else:
        return child

# phayseY is entirely superfluous.  One could just call phaseX with the
# father's genome and get the right answer. Still ...
def phaseY(father, child):
    """Phase Y chromosome genotype of a male child given the genotype of the
    father. This function only checks that the genotype of the child
    could have been inherited from the father.

    Arguments:
        mother (string): genotype of father expressed in the form 0/1, 0/0 etc
        child (string): genotype of child

    Returns:
        The phased genotype of the child if it is unambiguously resolvable.
        Otherwise the genotype is returned as is, i.e. unphased

    >>> phaseY('0/1', '0/1')
    '0/1'
    >>> phaseY('0/1', '1/1')
    '1|1'
    >>> phaseY('0/2', '1/1')
    '1/1'
    >>> phaseY('./.', '3/3')
    '3/3'

    """
    #ca, _, cb = child
    ca, cb = child.split('/')
    if ca != cb:  # we have Y chromosome heterozygosity
        return child
    if ca in father:
        return '%s|%s' % (ca, ca)
    else:
        return child

set_gtypes()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
