#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: anyopen
    :platform: Unix, Windows, MacOSX
    :synopsis: Transparent opening of compressed and uncompressed files

.. moduleauthor:: ; changjin.hong@gmail.com

Helper to open files of various types.  This module defines a function
*openfile* which is passed the name of the file and returns a file object
opened for reading.  The function can open files compressed using gzip,
bzip2, unix compress in addition to plain text uncompressed files.
"""

import os, sys, argparse
from gcn.lib.io import anyopen
from gcn.lib.io import fasta

class Hgmd_in_vcf:
  
  def __init__(self,hgmd_table,ref_file):
    self.tsv = hgmd_table
    self.ref = ref_file
    self.vcf = None
          
  def _write_head(self,fp2):
    fp2.write('##fileformat=VCFv4.2\n')
    fp2.write('##INFO=<ID=HGMD_ACC,Number=.,Type=String,Description="HGMD Accession Number",Source="HGMD",Version="03/01/2015">\n')
    fp2.write('##INFO=<ID=Gene,Number=.,Type=String,Description="HGNC Gene Symbol">\n')
    fp2.write('##INFO=<ID=VC,Number=7,Type=String,Description="HGMD Variant Class (DM=Disease causing (pathological) mutation, DM?=Likely disease causing (likely pathological) mutation, DFP=Disease associated polymorphism with additional supporting functional evidence, DP=Disease associated polymorphism, FP=Polymorphism affecting the structure, function or expression of a gene but with no disease association reported yet, FTV=Frameshift or truncating variant with no disease association reported yet, R=Retired entry)">\n')
    fp2.write('##INFO=<ID=PMID,Number=.,Type=String,Description="PubMed identification number">\n')
    fp2.write('##INFO=<ID=PHENO,Number=.,Type=String,Description="Associated Phenotype Description">\n')
    fp2.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    
  def to_vcf(self, out_vcf):
    self.vcf = out_vcf
    fp2= open(self.vcf,'w')
    
    #write vcf header
    self._write_head(fp2)
    
    fp = open(self.tsv,'r')
    cFa = fasta.Fasta(self.ref)
    cFa.load()
    
    qual = '100'
    filter = 'PASS'
    MISMATCH,INS,DEL = range(1,4)
    fp.next()
    for i in fp:
      
      items = i.rstrip().split('\t')
      #print items
      hgmd_acc,gene,strand,chrom,pos_start,pos_end,ref,alt,refseq_tx,hgvs,variant_class,rsid,pubmed,additional_pubmed = items
      
      if strand == '-':
        if ref != '.':
          ref = fasta.revcomplement(ref)
        if alt != '.':
          alt = fasta.revcomplement(alt)
      
      additional_pubmed = additional_pubmed.replace(" ", "")
      if pubmed != '.':
        if additional_pubmed != '.':
          pubmed+=',%s'%additional_pubmed
      elif additional_pubmed != '.':
        pubmed = additional_pubmed
      
      pos = int(pos_start)
      ref0 = cFa.get_bases(chrom,pos,1)
      if ref0 is None: continue
                        
      alts = alt.split(',')
      visited_ref_indels = False
      alts2 = []
      for alt in alts:
        if ref=='.':
          vtype = INS
        elif alt=='.':
          vtype = DEL
        elif len(ref) > len(alt):
          vtype = DEL
        elif len(ref) < len(alt):
          vtype = INS
        else:
          vtype = MISMATCH

        if vtype != MISMATCH:
          if not visited_ref_indels:
            visited_ref_indels = True
            if ref == '.':
              ref = ref0
            else:
              ref = ref0 + ref

          if alt == '.':
            alt = ref0
          else:
            alt = ref0 + alt
        alts2.append(alt)

      alt = ','.join(alts2)

      info = 'HGMD_ACC=%s'%hgmd_acc

      info += ';Gene=%s'%gene
      info += ';VC=%s'%variant_class
      if pubmed!='.':
        info += ';PMID=%s'%pubmed

      if vtype != MISMATCH:
        pos+=1
        
      fp2.write('%s\n' % '\t'.join([chrom,str(pos),rsid,ref,alt,qual,filter,info]))
    
    fp.close()
    fp2.close()
    cFa.fafp.close()

def main():
  
  desc = '-'
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-i', '--input', dest='hgmd_table', type=str, required=True, help='-')
  parser.add_argument('-r', '--ref', dest='ref_file', type=str, required=True, help='-')
  parser.add_argument('-o', '--output', dest='hgmd_vcf', type=str, default=True, help='-')
  args = parser.parse_args()
  
  cH = Hgmd_in_vcf(args.hgmd_table,args.ref_file)
  cH.to_vcf(args.hgmd_vcf)
  
if __name__ == "__main__":
  main()
