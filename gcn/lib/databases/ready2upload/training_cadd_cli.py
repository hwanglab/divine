#!/usr/bin/env python

import os,sys,argparse,re
import subprocess as sp
import gcn.lib.utils.lib_utils
import gcn.lib.io.vcf as vcf
from gcn.lib.databases.snpdb import ClinvarDB
from gcn.lib.databases.snpdb import HgmdDB
from gcn.lib.databases.snpdb import KGDB
from gcn.lib.io.vcfutils import normalize_variant
from scipy.stats import beta

import cPickle as pickle

def train_conservation_coeff(tr_varant_filt_vcf, include_hgmd):
  
  cadd_trset = {}
  mnp_cadd_trset = {}
  gerp_trset = {}

  v = vcf.VCFParser(tr_varant_filt_vcf)
  
  for rec in v:
    v.parseinfo(rec)
    
    found = False
    if include_hgmd:
      for cclass in rec.info.CLINSIG_CLASS:
        if not found:
          if 'HGMD' in cclass:
            found = True
            break
          
    if not found:
      for cclass in rec.info.CLINSIG_CLASS:
        if not found:
          if 'CLINVARDB' in cclass or '1kMAF' in cclass:
            found = True
            break

    if not found: continue
    
    varlist = normalize_variant(rec.chrom, rec.pos, rec.ref, rec.alt)
    mut_type = varlist[0][-1]
    
    if ':' in rec.id[0]:
      mut_type = 'mnp'
    
    #aaconv = './.'
    aaconv = '.'
    if rec.info.CADD_raw: 
      #to get CADD_raw (average)
      cadd_trset, mnp_cadd_trset, aaconv = vcf.get_CADD_scores_tr(mut_type, rec.info.CADD_aa, rec.info.CADD_raw, cadd_trset, mnp_cadd_trset)
    
    #to_get GERP score
    if rec.info.GerpConserve:
      gerp_trset = vcf.get_GERP_scores_tr(mut_type, aaconv, rec.info.GerpRSScore, gerp_trset)
      
  v.stream.close()
  return cadd_trset, mnp_cadd_trset, gerp_trset

def run_beta_fit(cadd_trset, mnp_cadd_trset, gerp_trset):
  '''
from scipy import stats  
import numpy as np  
import matplotlib.pylab as plt

# create some normal random noisy data
ser = 50*np.random.rand() * np.random.normal(10, 10, 100) + 20

# plot normed histogram
plt.hist(ser, normed=True)

# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]  
xmin, xmax = min(xt), max(xt)  
lnspc = np.linspace(xmin, xmax, len(ser))

ab,bb,cb,db = stats.beta.fit(ser)  
pdf_beta = stats.beta.pdf(lnspc, ab, bb,cb, db)  
plt.plot(lnspc, pdf_beta, label="Beta")

plt.show()
  '''
  
  cadd_trset_param = {}
  for aaconv in cadd_trset.keys():
     a,b,loc2,scale2 = beta.fit(cadd_trset[aaconv])
     mean2 = beta.mean(a,b,loc2,scale2)
     cadd_trset_param[aaconv] = [a,b,loc2,scale2,mean2]

  mnp_cadd_trset_param = {}
  for aaconv in mnp_cadd_trset.keys(): 
    a,b,loc2,scale2 = beta.fit(mnp_cadd_trset[aaconv])
    mean2 = beta.mean(a,b,loc2,scale2)
    mnp_cadd_trset_param[aaconv] = [a,b,loc2,scale2,mean2]

  gerp_trset_param = {}
  for aaconv in gerp_trset.keys():
    a,b,loc2,scale2 = beta.fit(gerp_trset[aaconv])
    mean2 = beta.mean(a,b,loc2,scale2)
    gerp_trset_param[aaconv] = [a,b,loc2,scale2,mean2]

  return cadd_trset_param, mnp_cadd_trset_param, gerp_trset_param

def main():
  parser = argparse.ArgumentParser(description="training cadd cli")
  parser.add_argument('--hgmd', action='store_const', dest = 'hgmd', required=False, default = False, const = True, help='want to include hgmd for training? It requires license. [False]')
  parser.add_argument('-r', action='store', dest='kg_sample_rate', required=False, default = .3, type=float, help='sampling rate for benign variants in 1kG')
  parser.add_argument('-o', action='store', dest='out_dir', required=True, help='output dir')
  parser.add_argument('--debug', action='store_const', dest='debug', required=False, default=False, const=True, help='debug?[False]')
  args=parser.parse_args()
  
  if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)

  #load HgmdDB (requires license!)
  if args.hgmd:
    try:
      hgmd = HgmdDB()
      training = hgmd.select_all_hgmd()
    except:
      print 'error: check if HGMD database is avail!'
      sys.exit(1)
	
  if False:
    #load ClinvarDB
    clnvar = ClinvarDB()
    training.extend(clnvar.select_all_clinvar())

    #load KGDB
    kgdb = KGDB()
    training.extend(kgdb.select_snps(0.1, 0.5, sample_rate=args.kg_sample_rate , snp_tag='benign_1kMAF'))
  
  #convert to vcf
  tr_vcfs = []

  clitags = ['benign','pathogenic']

  for c in range(2):
    tr_vcf = os.path.join(args.out_dir,'clin_%s_tr.vcf'%clitags[c])
    if not args.debug or not os.path.exists(tr_vcf):
      tr_vcf_body = tr_vcf + '.body' 
      fp2 = lib_utils.open2(tr_vcf_body,'w')
      printed = {}
      
      for chrom, pos, id, ref, alt, clisig in training:
        if clitags[c] in clisig:
          prim_key = '%s_%s_%s_%s' % (chrom,pos,ref,alt)
          if prim_key in printed: continue
          printed[prim_key]=clisig
          fp2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\tCLINSIG_CLASS=%s\n'%(chrom,pos,id,ref,alt,'100','PASS',clisig))
      fp2.close()
      
      #sort
      tr_vcf_body_so = tr_vcf_body + '.sorted'
      cmd = 'sort -k1,1 -k2,2n %s > %s' % (tr_vcf_body,tr_vcf_body_so)
      p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
      output, err = p.communicate()
      rc = p.returncode
          
      #append header
      tr_vcf_header = tr_vcf+'.head'
      fp2 = lib_utils.open2(tr_vcf_header,'w')
      fp2.write('##fileformat=VCFv4.2\n')
      fp2.write('##INFO=<ID=CLINSIG_CLASS,Number=7,Type=String,Description="benign_CLINVARDB,pathogenic_CLINVARDB,vus_CLINVARDB,benign_HGMDDB,pathogenic_HGMDDB,vus_HGMDDB,benign_1kMAF",Source="CLINVAR,HGMD",Version="03/01/2015">\n')
      fp2.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
      fp2.close()
      cmd = 'cat %s %s > %s' % (tr_vcf_header,tr_vcf_body_so,tr_vcf)
      p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
      output, err = p.communicate()
      rc = p.returncode
      os.unlink(tr_vcf_header)
      os.unlink(tr_vcf_body)
      os.unlink(tr_vcf_body_so)

    tr_vcfs.append(tr_vcf)
  
  training = None
  #run gcn
  gcn_dir = os.environ.get('GCN', None)
  if gcn_dir:
    annotpipe_bin = os.path.join(gcn_dir, 'gcn', 'bin', 'annotpipe.py')
  else:
    print 'error: cannot find annotpipe.py!'
    sys.exit(1)

  tr_varant_vcfs = []
  cadd_trset_params = []
  mnp_cadd_trset_params =[]
  gerp_trset_params =[]
  cadd_trsets = []
  mnp_cadd_trsets =[]
  gerp_trsets =[]
  
  filter_bin = os.path.join(gcn_dir,'gcn','lib','utils','filter_cj.py')   
  for c,tr_vcf in enumerate(tr_vcfs):
    tr_varant_vcf = os.path.join(args.out_dir,'clin_%s_tr_varant.vcf'%clitags[c])
    cmd = 'python %s -i %s -o %s'%(annotpipe_bin,tr_vcf,tr_varant_vcf)
    if not args.debug or not os.path.exists(tr_varant_vcf):
      print cmd
      p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
      output, err = p.communicate()
      rc = p.returncode
  
    #filter varant vcf
    tr_varant_filt_vcf = os.path.join(args.out_dir,'clin_%s_tr_varant_filt.vcf'%clitags[c])
    filterconf_tr = os.path.join(gcn_dir,'gcn','config','filter_tr_%s.conf'%clitags[c])
    
    cmd = 'python %s -i %s -o %s -f %s --no_genotype' %(filter_bin,tr_varant_vcf,tr_varant_filt_vcf,filterconf_tr)
    if not args.debug or not os.path.exists(tr_varant_filt_vcf):
      print cmd
      p = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
      output, err = p.communicate()
      rc = p.returncode
  
    cadd_trset, mnp_cadd_trset, gerp_trset = train_conservation_coeff(tr_varant_filt_vcf,args.hgmd)
    
    #to get parameters fitted in beta dist
    cadd_trset_param, mnp_cadd_trset_param, gerp_trset_param = run_beta_fit(cadd_trset, mnp_cadd_trset, gerp_trset)
    
    cadd_trset_params.append(cadd_trset_param)
    mnp_cadd_trset_params.append(mnp_cadd_trset_param)
    gerp_trset_params.append(gerp_trset_param)
    
    cadd_trsets.append(cadd_trset)
    mnp_cadd_trsets.append(mnp_cadd_trset)
    gerp_trsets.append(gerp_trset)
    
    print 'done.'

  pyv = os.path.join(args.out_dir,'clin_tr.pyv')
  fp2 = open(pyv,'wb')
  pickle.dump([cadd_trset_params, mnp_cadd_trset_params, gerp_trset_params], fp2)
  fp2.close()

if __name__ == "__main__":
  main()
