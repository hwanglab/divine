"""
.. module:: converter
    :platform: Unix, Windows, MacOSX
    :synopsis: hg coordinate version convertor

.. moduleauthor:: changjin.hong@gmail.com

This load ref genome and refseq tx for HGVS nomenclature.

import pyhgvs
import pyhgvs.utils
from pygr.seqdb import SequenceFileDB

refseq_acc='NM_000492'
chrom='chr7'
offset=117307123
ref='AGAG'
alt='A'

genome_fn = 'hg19.fa'
refseq_fn = 'refGene.txt'

genome = SequenceFileDB(genome_fn)
fp = open(refseq_fn,'r')
refseq = pyhgvs.utils.read_transcripts(fp)
fp.close()

transcript = refseq.get(refseq_acc)
hgvs_name = pyhgvs.format_hgvs_name(chrom, offset, ref, alt, genome, transcript)

"""

import os
import pyhgvs
import pyhgvs.utils
from gcn.etc import fileconfig
from pygr.seqdb import SequenceFileDB
from gcn.lib.databases.refgene import get_ucsc_chrom

CHROMOSOMES, ICHROMOSOMES = get_ucsc_chrom()

class Hgvs2:

    def __init__(self,refseq_fn=fileconfig.FILECONFIG['REFGENE']):
        """ register hg19 reference genome sequence and NCBI RefSeq transcript coordinates"""
        self.genome_fn = fileconfig.FILECONFIG['REFGENOME_UCSC']
        self.genome = None
        if not os.path.exists(self.genome_fn):
          raise IOError('Reference Genome Sequence (UCSC format) for %s is not found'%self.genome_fn)

        self.refseq_fn = refseq_fn
        self.refseq = None
        if not os.path.exists(self.refseq_fn):
          raise IOError('NCBI RefSeq transcript for %s is not found'%self.refseq_fn)
        
    def get_transcript(self, name):
        if self.refseq:
            return self.refseq.get(name)
        else:
            raise RuntimeError('first, load a transcript coordinate file![%s]'%self.refseq_fn)

    def load_resource(self):
        #load genome sequence
        print 'loading the genome sequence [%s] for HGVS...' % self.genome_fn
        self.genome = SequenceFileDB(self.genome_fn)
        print 'done.'
        
        #load refseq into dic
        print 'loading the refseq transcript [%s] for HGVS...' % self.refseq_fn
        fp = open(self.refseq_fn,'r')
        self.refseq = pyhgvs.utils.read_transcripts(fp)
        fp.close()
        print 'done.'
    
    def close_resource(self):
      self.genome.close()
      self.refseq = None
    
    def to_cDNA(self, chrom, offset, ref, alt, refseq_acc):
        """ convert to HGVS nomenclature """
        transcript = self.get_transcript(refseq_acc)
        
        if not chrom.startswith('chr'):
          chrom = 'chr%s'%chrom
          if chrom not in CHROMOSOMES:
            return None
          
        if not chrom in self.genome.keys():
          return None
          
        cdna = pyhgvs.format_hgvs_name(chrom, offset, ref, alt, self.genome, transcript)
        if cdna:
          itms = cdna.split(':')
          if len(itms)>1:
            cdna = itms[1]
          else:
            cdna = cdna
          return cdna
        else:
          return None
         
    def gdna_to_vcf(self,gdna):
        return pyhgvs.gdna_to_vcf(gdna,self.genome)
       
    def to_chrom_coordinate(self,cDNA):
        try:
            chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(cDNA, self.genome, \
                                       get_transcript = self.get_transcript)

            return chrom, offset, ref, alt
        except:
            print "[%s] cannot be coverted to chromosome coordinate"%cDNA
            return None, None, None, None


if __name__ in "__main__":
    hgvs2 = Hgvs2()
    
    chrom, offset, ref, alt, refseq_acc = ['chr1', 897325,'G','C','NM_198317']
    cdna = hgvs2.to_cDNA(chrom, offset, ref, alt, refseq_acc)
    
    if cdna:
      print cdna

    chrom, offset, ref, alt, refseq_acc = ['chr1', 236899899,'TC','T','NM_001278344']
    cdna = hgvs2.to_cDNA(chrom, offset, ref, alt, refseq_acc)
    
    if cdna:
      print cdna

   
