"""
.. module:: filter
    :platform: Unix, Windows, MacOSX
    :synopsis: Filters VCF record based on different fields

.. developed by Changjin Hong

Module to filter a vcf record based on different criteria:
 - Filter field (include or exclude)
 - Info field flag (include or exclude)
 - Info field value (check with operator)
 - Region based filter (include or exclude)

 
 e.g.:
[fltr]
incl=PASS
[infoflag]
excl=DB137
[infoval]
KGAF=lte:0.05 # TBD
[reg]
incl=exonic
[freq]
incl=lowgmaf

 """
import ConfigParser
# import gcn.lib.io.snpeff as snpeff
from gcn.lib.varann.vartype.varant import varant_parser as vp
from gcn.lib.databases.snpdb import classify_clnsig
from gcn.lib.databases.snpdb import ClinvarDB
import gcn.lib.io.vcf as vcf
import re
import os
from html5lib.treewalkers.base import UNKNOWN
from sre_constants import NEGATE

def str2bool(v):
  if v.lower() in ("yes", "true", "y"):
    return True
  elif v.lower() in ("no", "false", "n"):
    return False
  else:
    return v
  
def conv_str2format(val):
  try:
    val = str2bool(val)
  except ValueError:
    pass
    
  try:
    val = float(val)
  except ValueError:
    pass
  return val

def predict_gender_from_VCF(single_vcf,sample_id):

  from gcn.data import pseudoautosomal_genes
  UNKNOWN,MALE,FEMALE = range(3)

  f = Filter()
  f.geneincl = pseudoautosomal_genes.PSEUDO_AUTO_GENES

  v = vcf.VCFParser(single_vcf,sampleids=[sample_id])

  gender = UNKNOWN
  chrmXY = [0,0,0]

  msg = "predicting gender from the sample [%s]" %single_vcf
  print msg

  for rec in v:
    if rec['chrom'] == 'chrY' or rec['chrom'] == 'Y':
      v.parsegenotypes(rec)
      v.parseinfo(rec)

      if rec[v.samples[0]].GT != './.':
        chrmXY[1] += 1
      
      if not f.in_gene(rec,f.geneincl):
        chrmXY[2] += 1
    elif rec['chrom'] == 'chrX' or rec['chrom'] == 'X':
      chrmXY[0] += 1
  
  if chrmXY[0] > 0:
    chrY2X_rate = 1.*(chrmXY[1]+chrmXY[2])/chrmXY[0]
    if chrY2X_rate > 0.01:
      gender = MALE
    else:
      gender = FEMALE
  elif chrmXY[2] > 0:
    gender = MALE
  else:
    gender = UNKNOWN

  v.stream.close()
  msg = "gender identified [%d], Done."%gender
  print msg
  
  return gender
   
class Filter:
    """Class to filter a VCF record"""

    def __init__(self, filterconf=None):
        self.iflgincl = []
        self.iflgexcl = []
        self.invalincl = []
        self.invalexcl = []
        self.fltrincl = []
        self.fltrexcl = []
        self.regincl = []
        self.regexcl = []
        self.freqincl = []
        self.cli_freqincl = []
        self.geneincl = []
        self.geneexcl = []
        self.vtincl = []
        self.vtexcl = []
        self.dconf = {}
        self.clinvar_genes = None
        self.NO_HGMD,self.HGMD_GENE,self.HGMD_POS,self.HGMD_ALT = range(4)
        
        #print "filterconf",  filterconf
        if filterconf and os.path.exists(filterconf):
          self.parseconf(filterconf)

    def parseconf(self, filterconf):
        """Parse filter configuration file to filter conditions"""
        config = ConfigParser.RawConfigParser()
        config.read(filterconf)
        MIN_MAFCLI = 0.05
        MIN_MAF = 0.01
        if config.has_section('fltr'):
            self.fltrincl, self.fltrexcl = self.conf_flags(config.items('fltr'))
        if config.has_section('reg'):
            self.regincl, self.regexcl = self.conf_flags(config.items('reg'))
        if config.has_section('freq'):
            self.freqincl, _ = self.conf_flags(config.items('freq'))
            if not self.freqincl:
              self.freqincl.append(MIN_MAF)
        if config.has_section('freq_cli'):
            self.cli_freqincl, _ = self.conf_flags(config.items('freq_cli'))
            if not self.cli_freqincl:
              self.cli_freqincl.append(MIN_MAFCLI)            
        if config.has_section('vartype'):
            self.vtincl, self.vtexcl = self.conf_flags(config.items('vartype'))
        if config.has_section('infoflag'):
            self.iflgincl, self.iflgexcl = self.conf_flags(config.items('infoflag'))

        if config.has_section('infoval'):
            self.dconf = self.conf_items(config.items('infoval'))
          

    def store_genelist(self, select_mode, genelist=None):
        """Parse filter configuration file to filter conditions"""

        if genelist is None: return []
        
        fp = open(genelist,'r')
        in_block = False
        for i in fp:
          if i.startswith('#'):
            if in_block: break
            if i.startswith('#%s'%select_mode):
              in_block = True
              gene_sel=[]
          else:
            gene_sel.append(i.rstrip())
        fp.close()
        
        if not in_block:
          raise ValueError("the gene list file [%s] does not contain a section, '%s'"%(genelist,select_mode))
        
        return gene_sel
    
    def conf_items(self, flagitems):
      dConfOpt = {}
      for opt, val in flagitems:
        if not opt in dConfOpt:
          dConfOpt[opt] = conv_str2format(val) 
      return dConfOpt
      
    def conf_flags(self, flagitems):
        incl = []
        excl = []
        for opt, val in flagitems:
            if opt == 'incl':
                vals = val.split(':')
                for v in vals:
                  incl.append(conv_str2format(v))
            else:
                vals = val.split(':')
                for v in vals:
                  excl.append(conv_str2format(v))
                    
        return incl, excl
        
    def infoval(self, sample):
        """Get genotype alleles"""
        pass
      
    def infoflag(self, sample):
        pass
      
    def issnp(self, rec):
        if len(rec.ref) == 1:
            if len(rec.alt) == 1 and len(rec.alt[0]) == 1:
                return True
            else:
                if len(rec.alt[0]) == 1:
                    return True
        return False
    
    def get_clinvar_pathogenes(self):
      clnvar = ClinvarDB()
      self.clinvar_genes = clnvar.get_clinvar_genes()
    
    def get_clinvar_class(self, rec):
      
      UNKNOWN = 0
      clinvar_sigs = []
      if rec.info.CLNSIG:
        clnsigs = rec.info.CLNSIG
        #check if it was normalized
        if '__' in clnsigs[0]:
          clnsigs = clnsigs[0].split('__')

        for sigs in clnsigs:
          sigs = sigs.split('|')
          for sig in sigs:
            if sig=='.': continue
            clinvar_sigs.append(sig)
        return classify_clnsig(clinvar_sigs)
      else:
        return UNKNOWN

    def get_hgmd_class(self, rec):
      FP,UNKNOWN,GENE,PREDICTED,LIKELY_DM,DM = range(-1,5)
      hgmd_pos = ['DM']
      hgmd_likely_pos = ['DM?']
      hgmd_predicted_pos = ['DFP','DP']
      hgmd_neg = ['FP','R']
      vc_found = UNKNOWN
      hgmd_match_types = [self.HGMD_POS,self.HGMD_ALT]
      
      if rec.info.HGMDDB:
        '''
        there are 6 different hgmd class
        DFP (disease-associated polymorphisms with additional supporting functional evidence. These are reported to be in significant association with disease (p<0.05) and to have evidence for being of direct functional importance (e.g. as a consequence of altered expression, mRNA studies etc))
        DM (disease-causing mutations, pathological mutations reported to be disease causing in the original literature report. The other three tags are used for polymorphisms. )
        DM? (likely DM)
        DP (disease-associated polymorphisms. These are reported to be in significant association with disease (p<0.05) and are assumed to be functional (e.g. as a consequence of location, evolutionary conservation, replication studies etc), although there may be as yet no direct evidence (e.g. from an expression study) of function)
        FP (in vitro/laboratory or in vivo functional polymorphisms. These are reported to affect the structure, function or expression of the gene (or gene product), but with no disease association reported as yet)
        R (Retired entry)
        '''
        hgmd_vars = rec.info.HGMD_DESC
        for var in hgmd_vars:
          mObj = None
          mObj = re.search(r'(\d+)\((.+)\)',var)
          if mObj:
            hgmd_match = int(mObj.group(1))
            if hgmd_match in hgmd_match_types:
              vc = mObj.group(2).split('|')[1]
              if vc in hgmd_pos:
                vc_found = vc
                break
              elif vc in hgmd_likely_pos:
                if vc_found<LIKELY_DM:
                  vc_found = LIKELY_DM
              elif vc in hgmd_predicted_pos:
                if vc_found<PREDICTED:
                  vc_found = PREDICTED
              elif vc in hgmd_neg:
                if vc_found==UNKNOWN:
                  vc_found = FP
          else:
            print 'check HGMD_DESC[%s]'%var
      
      return vc_found

    def in_pathogene(self, rec):
      
      if self.in_gene(rec, self.clinvar_genes):#check first clinvar
        return True
      elif 'HGMDDB' in rec.info:
        '''
        there are 6 different hgmd class
        DFP (disease-associated polymorphisms with additional supporting functional evidence. These are reported to be in significant association with disease (p<0.05) and to have evidence for being of direct functional importance (e.g. as a consequence of altered expression, mRNA studies etc))
        DM (disease-causing mutations, pathological mutations reported to be disease causing in the original literature report. The other three tags are used for polymorphisms. )
        DM? (likely DM)
        DP (disease-associated polymorphisms. These are reported to be in significant association with disease (p<0.05) and are assumed to be functional (e.g. as a consequence of location, evolutionary conservation, replication studies etc), although there may be as yet no direct evidence (e.g. from an expression study) of function)
        FP (in vitro/laboratory or in vivo functional polymorphisms. These are reported to affect the structure, function or expression of the gene (or gene product), but with no disease association reported as yet)
        R (Retired entry)
        '''
        for var in rec.info.HGMD_DESC:
          #print var #debug
          mObj = re.search(r'(\d+)\((.+)\)',var)
          if mObj:
            hgmd_match = int(mObj.group(1))
            if hgmd_match in [self.HGMD_GENE]:
              return True
          else:
            print 'check HGMD_DESC[%s]'%var
            
      return False
      
    def in_gene(self, rec, genes):
        vpop = vp.parse(rec.info)
        for altnum, val in vpop.items():
            for gene, gd in val.items():
                if gene in genes:
                    return True
        
    def vtexonic(self, rec):
        """
        Check if variant is in the exonic region (varant annotated)
        TODO: apply a max distance of donor/acceptor distance; a cascade filter  
        """
        hpm = 'StopGain StopLoss StartLoss NonSyn FrameShiftInsert FrameShiftDelete NonFrameShiftInsert NonFrameShiftDelete'.split()
        splc = ['SpliceDonor', 'SpliceAcceptor']
        splc_coding = ['ESS', 'ESE']
        warn = ['CDS_NOT_MULTIPLE_OF_3'] 
        intronic = ['CodingIntronic','NonCodingIntronic']
        exonic = ['CodingExonic','NonCodingExonic']
        utr = ['UTR5','UTR3']

        vpop = vp.parse(rec.info)
        regions = []

        for altnum, val in vpop.items():
            for gene, gd in val.items():
                if gd:
                    #to extract recessive/dominant
                    inherited = 'as_recessive'
                    for mim in gd['MIM_PHENS']:
                        if 'AUTOSOMAL_DOMINANT' in mim: #TODO:https://ghr.nlm.nih.gov/handbook/inheritance/inheritancepatterns
                          inherited = 'as_dominant'
                          
                    for t in gd['TRANSCRIPTS']:
                        if self.regincl=='all' or t.region in self.regincl:
                          region = [t.trans_id, t.region, t.protein_len, inherited]
                          if t.region in exonic: #CodingExonic,NonCodingExonic
                            if t.splice in splc_coding:
                              region_tag = 'splc_coding'
                              if t.mutation in hpm:
                                region_tag += ';hpm'
                              region.append(region_tag)
                              regions.append(region)
                            elif t.mutation in hpm:
                              region.append('hpm')
                              regions.append(region)
                            elif t.warning in warn:
                              region.append('warn')
                              regions.append(region)
                          elif self.dconf['splice_dist']>0:
                            dists = []
                            tcdna = t.cdna.split('_')
                            T = len(tcdna)
                            if t.region in intronic:
                              mObj = re.search(r'c\.(.+)[+-](\d+)',tcdna[0])
                              if mObj:
                                dists.append(int(mObj.group(2)))
                              if T>1:
                                mObj = re.search(r'(.+)[+-](\d+)',tcdna[1])
                                if mObj:
                                  dists.append(int(mObj.group(2)))
                            elif t.region in utr:
                              mObj = re.search(r'c\.[-\*](\d+)',tcdna[0])
                              if mObj:
                                dists.append(int(mObj.group(1)))
                              if T>1:
                                mObj = re.search(r'[-\*](\d+)',tcdna[1])
                                if mObj:
                                  dists.append(int(mObj.group(1)))
                            
                            if dists:
                              if min(dists) <= self.dconf['splice_dist']:
                                region.append('splc_ext_intron')
                                regions.append(region)
        
        if regions:
          return regions
        elif self.dconf['regulome']:
          for altnum, val in vpop.items():
            for gene, gd in val.items():
              if gd:
                inherited = 'as_recessive'
                for mim in gd['MIM_PHENS']:
                    if 'AUTOSOMAL_DOMINANT' in mim: #TODO:https://ghr.nlm.nih.gov/handbook/inheritance/inheritancepatterns
                      inherited = 'as_dominant'
                      break
                for t in gd['TRANSCRIPTS']:
                  if rec.info.RegulomeScore:
                    return [t.trans_id, t.region, t.protein_len, inherited, 'regulome']
                  else:
                    return [t.trans_id, None, None, inherited, None]

    def novel(self, rec):
        """check if variant is novel"""
        if not rec.info.KGDB:
            return True
        else:
            return False
    
    def lowfreq(self, vargmaf, freq=0.05):
        for el in vargmaf:
            args = el.split(':')
            for varfreq in args:
                if float(varfreq) <= freq:
                    return True
        return False

    def hifreq(self, vargmaf, freq=0.05):
        for el in vargmaf:
            args = el.split(':')
            for varfreq in args:
                if float(varfreq) > freq:
                    return True
        return False
          
    def gmaf(self, rec, freq=0.05):
        """Check if variant occurs at a lower frequency than specified
        this def is obsolete and it is replaced by gmaf_stringent()
        """
        
        low_afkg = False
        low_afesp = False
        low_afexac = False
        # This check is for FB annotations
        if rec.info.KGAF:
            for el in rec.info.KGAF:
                if float(el) <= freq:
                    low_afkg = True
                    break
            if low_afkg:
                return True
        if rec.info.ESPAF:
            for el in rec.info.ESPAF:
                if float(el) <= freq:
                    low_afesp = True
                    break
            if low_afesp:
                return True
        if rec.info.EXACAF:
            for el in rec.info.EXACAF:
                if float(el) <= freq:
                    low_afexac = True
                    break
            if low_afexac:
                return True

        # For varant annotations
        if ((rec.info.KGDB and self.lowfreq(rec.info.KGAF, freq)) or (not rec.info.KGDB)) \
            or ((rec.info.ESPDB and self.lowfreq(rec.info.ESPAF, freq)) or (not rec.info.ESPDB)) \
            or ((rec.info.EXACDB and self.lowfreq(rec.info.EXACAF, freq)) or (not rec.info.EXACDB)):
            return True
        else:
            return False
        return False
 
    def _get_maf_varant(self, info_maf):
        maf = []
        for el in info_maf:
          if ':' in el:
            min_af = 1.
            for af in el.split(':'):
              af = float(af)
              if af <= min_af:
                min_af = af
            maf.append(min_af)
          else:
            maf.append(float(el))
        return maf
        
    def gmaf_stringent(self, rec, freq=0.05):
        """this is more strict version of gmaf().
        if at least one maf is higher than freq, then it returns fail
        if all maf (as long as it exists) is lower than freq, then it return true
        if any maf is not available, then it return true.  
        """
        
        mafs = []
        if 'kgaf' in self.dconf and rec.info.KGDB:
          mafs.append(self._get_maf_varant(rec.info.KGAF))
        if 'espaf' in self.dconf and rec.info.ESPDB:
          mafs.append(self._get_maf_varant(rec.info.ESPAF))
        if 'exacaf' in self.dconf and rec.info.EXACDB:
          mafs.append(self._get_maf_varant(rec.info.EXACAF))
        
        if not mafs: return 'novel'
        
        M = len(mafs)
        N = len(mafs[0])
        conditions = ['rare']*N
        for n in range(N):
          for m in range(M):
            if mafs[m][n] > freq:
              conditions = 'frequent'
              break
        
        if 'rare' in conditions: return 'rare'
        else: return 'frequent'
                      
    def prociflag(self, rec):
        incl = False
        excl = False
        if self.iflgincl:
            for el in self.iflgincl:
                if el in rec.info:
                    incl = True
        else:
            incl = True
        
        if self.iflgexcl:
            for el in self.iflgexcl:
                if el in rec.info:
                    excl = True
        if incl and not excl:
            return True
        else:
            return False

    def retain(self, rec):
        """Function to check if variant satisfies the filter conditions"""
        
        retain = fltrretain = regretain = freqretain = generetain = vtretain = iflgretain = False
        class_tag = ''

        if self.fltrincl:
            for f in rec.filter:
                if f in self.fltrincl:
                    fltrretain = True
        else:
            fltrretain = True

        if self.fltrexcl:
            for f in rec.filter:
                if f in self.fltrexcl:
                    return False, None, '1'

        if self.geneincl:
           if self.in_gene(rec,self.geneincl):
              generetain = True
        else:
            generetain = True

        if self.vtincl:
            if 'snp' in self.vtincl and self.issnp(rec):
                vtretain = True
            elif 'indel' in self.vtincl and not self.issnp(rec):
                vtretain = True
        else:
            vtretain = True
            
        if self.prociflag(rec):
            iflgretain = True

        clinvar_dm = self.get_clinvar_class(rec)
        if clinvar_dm == -1:
          return retain, None, '1c'

        hgmd_vc = self.get_hgmd_class(rec)
        if hgmd_vc == -1:
          return retain, None, '1c'

        if clinvar_dm == 1 or hgmd_vc > 2:
          maf_tag = self.gmaf_stringent(rec, self.cli_freqincl[0])
          if maf_tag == 'frequent':
            return retain, None, '2c'
          else:
            class_tag += 'c'
        else:
          maf_tag = self.gmaf_stringent(rec, self.freqincl[0])
          if maf_tag == 'frequent':
            return retain, None, '2'
        
        freqretain = True
        if maf_tag == 'rare':
          class_tag += '3'
        elif maf_tag == 'novel':
          class_tag += '4'
                  
        region_score = {'CodingExonic':5,'CodingIntronic':4,'NonCodingExonic':3,'UTR5':2,'UTR3':2,'NonCodingIntronic':1}
        mut_type = {'splc_coding;hpm':5,'hpm':4,'warn':3,'splc_coding':2,'splc_ext_intron':1}
        non_codings = ['NonCodingExonic','NonCodingIntronic']
        
        #print rec #debug
        region_annots = self.vtexonic(rec) #[gene,region,protein_len,inherited,(splc_coding,hpm,warn,splc_ext_intron,regulome)]

        reg_score = []
        type_score = []
        txes = {}
        
        if not region_annots:
          return retain, txes, '2'
        
        regretain = True
        for trans_id,region,protein_len,inherited,mut_annot in region_annots:
          if region is None:
            if 'c' in class_tag:
              class_tag += 'I'
              txes[trans_id] = inherited
              return True, txes, class_tag
            else:
              return retain, txes, '2'
          else:
            reg_score.append(region_score[region])
            type_score.append(mut_type[mut_annot])
        
        #choose the most damaged region ---------------------->
        max_regsc = max(reg_score)
        maxi = [i for i,sc in enumerate(reg_score) if sc == max_regsc]

        max_typesc = 0
        for j in maxi:
          if type_score[j]>=max_typesc:
            max_typesc = type_score[j]

        maxj = [j for j in maxi if type_score[j] == max_typesc]
        #<-----------finally, we choose which variant is the most significant!
        
        for j in maxj:
          trans_id,region,protein_len,inherited,mut_annot = region_annots[j]
          if trans_id not in txes:
            txes[trans_id]=inherited

        if mut_annot == 'splc_coding;hpm':
          class_tag += 'S'
        elif mut_annot == 'hpm':
          class_tag += 'e'
        elif mut_annot == 'warn':
          class_tag += 'w'
        elif mut_annot == 'splc_coding':
          class_tag += 's'
        elif mut_annot == 'splc_ext_intron':
          class_tag += 'i'

        if region in non_codings:
          class_tag += 'n'
        
        if 'c' not in class_tag:
          if self.in_pathogene(rec):
            class_tag += 'g'

        if rec.info.LCR: #low complexity
          if 'c' not in class_tag: 
            return retain, txes, class_tag
          
        if fltrretain and regretain and freqretain and generetain and vtretain and iflgretain:
          retain = True

        return retain, txes, class_tag 

if __name__ == "__main__":
    import argparse
    from gcn.lib.utils.filter_cj import *
    import gcn.lib.utils.fileutils as fileutils
    import os
    import gcn.lib.io.vcf as vcf
    
    parser = argparse.ArgumentParser(description='Filter variants in the vcf')
    parser.add_argument('-i', dest='infile', help='input vcf file')
    parser.add_argument('-s', dest='sample_id', default=None, help='sample_id to work on the input vcf [None]')
    parser.add_argument('-f', dest='filterconf', default=None, help='filterfile')
    parser.add_argument('-o', dest='outfile', default=None, help='output file')
    parser.add_argument('-l', dest="genelist",default=None, help='Comma separated input list to include')
    parser.add_argument('--no_genotype', action='store_const', dest='skip_parse_genotype', required=False, default=False, const=True, help='want to skip parsing genotypes [False]')
    
    options = parser.parse_args()
    
    #create filter class
    f = Filter(options.filterconf)
    
    #read gene list file to tell which genes are to be included/excluded
    f.geneincl = f.store_genelist('incl',options.genelist)
    f.geneexcl = f.store_genelist('excl',options.genelist)
    
    min_depth = 0
    if 'min_depth' in f.dconf:
      min_depth = f.dconf['min_depth']
    
    f.get_clinvar_pathogenes()

    if options.sample_id:
      v = vcf.VCFParser(options.infile,sampleids = [options.sample_id])
    else:
      v = vcf.VCFParser(options.infile)

    if options.skip_parse_genotype:
      gt_sample = False
    else:
      if v._sampleids:
        sidx = v._sampleids[0]
        gt_sample = True
      else:
        gt_sample = False

    ostream = open(options.outfile, 'w')
    v.add_meta_info("CLASS_TAG", "1", "String","1:confirmed_benign, b:inherited from healthy parents, 2:frequent, 3:rare, 4:novel, s:splice, i:intronic_coding, e:exonic, n:noncoding_exonic, c:previously_known_pathogenic_vairant, g:previously_known_pathogenic_gene, I:previously_known_pathogenic_variant_in_non_exonic, h:compound_het, H:homozygous, D:autosomal_dominant, d:de-novo mutation from parents, w:warning(e.g.,CDS_not_multiple_3)")
    v.writeheader(ostream)
    
    for rec in v:
      v.parseinfo(rec)

      if gt_sample:
        v.parsegenotypes(rec)
        
        if min_depth>0 and int(rec[v.samples[sidx]].DP) < min_depth:
          continue
      
        if rec[v.samples[sidx]].GT == './.':
          continue

      retain, dTxes, rec.info.CLASS_TAG = f.retain(rec)
      if retain:
        v.write(ostream, rec)

    ostream.close()
    v.stream.close()
