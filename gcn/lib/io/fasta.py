"""
.. module:: bed
    :platform: Unix, Windows, MacOSX
    :synopsis: Read BED files

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Parser for BED (Browser Extensible Data)files. A detailed description of
the format can be found at the UCSC_ site.

.. _UCSC: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
"""
import os
import pysam
from gcn.data.complement import COMPLEMENT

def revcomplement(seq):
    """Reverse Complements a sequence"""
    revcomseq = ""
    revseq = seq[::-1]
    for neu in revseq:
        revcomseq += COMPLEMENT[neu]
    return revcomseq
      
class Fasta:
  def __init__(self,ref_fn,samtools_bin='samtools'):
    self.ref = ref_fn
    self.samtools_bin = samtools_bin
    self.refidx = self.create_index()
    self.fafp = None
    self.chromosomes = []
  
  def create_index(self):
     #print 'generating faidx...'
     cmd='%s faidx %s' % (self.samtools_bin,self.ref)
     refIdx=self.ref+'.fai'
     if not os.path.exists(refIdx):
       os.system(cmd)
     return refIdx
   
  def load(self):
    self.fafp = pysam.Fastafile(self.ref)
    self.chromosomes = self.fafp.references

  def get_bases(self,chrom,pos,bp_len=1):
    ''' note that pos should be 0-based'''
    
    if self.fafp:
      if chrom in self.chromosomes:
        return self.fafp.fetch(chrom,pos,pos+bp_len)
      else:
        return None
    else:
      self.load()
      if chrom in self.chromosomes:
        return fafp.fetch(chrom,pos,pos+bp_len)
      else:
        return None

# if __name__ == "__main__":
#     import sys
#     bedfile = sys.argv[1]
#     for line in anyopen.open(bedfile):
#         v = parseline(line)
# 
#     print(repr(v))
#     print(v)
