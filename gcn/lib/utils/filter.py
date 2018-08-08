"""
.. module:: filter
    :platform: Unix, Windows, MacOSX
    :synopsis: Filters VCF record based on different fields

.. moduleauthor:: TCS

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
import os

class Filter:
    """Class to filter a VCF record"""

    def __init__(self, filterconf):
        self.iflgincl = []
        self.iflgexcl = []
        self.invalincl = []
        self.invalexcl = []
        self.fltrincl = []
        self.fltrexcl = []
        self.regincl = []
        self.regexcl = []
        self.freqincl = []
        self.freqexcl = []
        self.geneincl = []
        self.geneexcl = []
        self.vtincl = []
        self.vtexcl = []
        
        #print "filterconf",  filterconf
        if isinstance(filterconf, unicode) and os.path.exists(filterconf):
            self.parseconf(filterconf)
        else:
            self.parsefl(filterconf)
    def parseconf(self, filterconf):
        """Parse filter configuration file to filter conditions"""
        config = ConfigParser.RawConfigParser()
        config.read(filterconf)
        if config.has_section('fltr'):
            self.fltrincl, self.fltrexcl = self.conf_flags(config.items('fltr'))
        if config.has_section('reg'):
            self.regincl, self.regexcl = self.conf_flags(config.items('reg'))
        if config.has_section('freq'):
            self.freqincl, self.freqexcl = self.conf_flags(config.items('freq'))
        if config.has_section('vartype'):
            self.vtincl, self.vtexcl = self.conf_flags(config.items('vartype'))
        if config.has_section('infoflag'):
            self.iflgincl, self.iflgexcl = self.conf_flags(config.items('infoflag'))
    
    def parsefl(self, filterlist):
        """Parse filter configuration file to filter conditions"""
        self.fltrincl, self.fltrexcl = filterlist[0],  []
        self.regincl, self.regexcl = filterlist[1],  []
        self.freqincl, self.freqexcl = filterlist[2],  []
        self.geneincl, self.geneexcl = filterlist[3],  []
        
    def conf_flags(self, flagitems):
        incl = []
        excl = []
        for opt, val in flagitems:
            if opt == 'incl':
                incl.extend(val.split(':'))
            else:
                excl.extend(val.split(':'))
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
        
    def exonic(self, rec):
        """Check if variant is in the exonic region"""
        FLTR_EFF = ["SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", "START_LOST"\
        , "EXON_DELETED", "FRAME_SHIFT", "STOP_GAINED", "STOP_LOST", "RARE_AMINO_ACID", "NON_SYNONYMOUS_CODING"\
        , "CODON_CHANGE", "CODON_INSERTION", "CODON_CHANGE_PLUS_CODON_INSERTION"\
        , "CODON_DELETION", "CODON_CHANGE_PLUS_CODON_DELETION" \
        , "NON_SYNONYMOUS_START", "START_GAINED"]
        
        for eff in rec.info.EFF:
          s = 1
            #if snpeff.parse_snpeff(eff)['Effect'] in FLTR_EFF:
            #    return True
        return False

    def vtexonic(self, rec):
        """Check if variant is in the exonic region (varant annotated)"""
        hpm = 'StopGain StopLoss StartLoss NonSyn FrameShiftInsert FrameShiftDelete NonFrameShiftInsert NonFrameShiftDelete'.split()
        splc = ['SpliceDonor', 'SpliceAcceptor']
        warn = ['CDS_NOT_MULTIPLE_OF_3'] 
        
        vpop = vp.parse(rec.info)
        for altnum, val in vpop.items():
            for gene, gd in val.items():
                if gd:
                    for t in gd['TRANSCRIPTS']:
                        if t.mutation in hpm:
                            return True
                        if t.splice in splc:
                            return True
                        if t.warning in warn:
                            return True
        return False

    def gene(self, rec,  gl):
        """Check if variant is in the gene list"""
        for eff in rec.info.EFF:
          s = 1
            #if snpeff.parse_snpeff(eff)['Gene_Name'] in gl:
            #    return True
        return False

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
    
    def gmaf(self, rec, freq=0.05):
        """Check if variant occurs at a lower frequency than specified"""
        low_afkg = False
        low_afesp = False
        low_afexac = False
        # This check is for FB annotations
        if rec.info.AFKG:
            for el in rec.info.AFKG:
                if el <= freq * 100:
                    low_afkg = True
                    break
            if low_afkg:
                return True
        if rec.info.AFESP:
            for el in rec.info.AFESP:
                if el <= freq * 100:
                    low_afesp = True
                    break
            if low_afesp:
                return True
        if rec.info.AFEXAC:
            for el in rec.info.AFEXAC:
                if el <= freq * 100:
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
	#print self.__dict__
	#print rec
        retain = fltrretain = regretain = freqretain = generetain = vtretain = iflgretain = False
        if self.fltrincl:
            for f in rec.filter:
                if f in self.fltrincl:
                    fltrretain = True
        else:
            fltrretain = True
        if self.regincl:
           if 'exonic' in self.regincl and self.exonic(rec):
                regretain = True
           elif 'vthpv' in self.regincl and self.vtexonic(rec):
                regretain = True
        else:
            regretain = True
        
        if self.freqincl:
           if self.gmaf(rec,  float(self.freqincl[0])):
               freqretain = True
        else:
            freqretain = True
        if self.geneincl:
           if self.gene(rec, self.geneincl):
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
            
        if fltrretain and regretain and freqretain and generetain and vtretain and iflgretain:
            retain = True
        return retain

if __name__ == "__main__":
    import argparse
    from gcn.lib.utils.filter import *
    import gcn.lib.utils.fileutils as fileutils
    #import gcn.lib.io.snpeff as snpeff
    import os
    import gcn.lib.io.vcf as vcf
    
    parser = argparse.ArgumentParser(description='Filter variants in the vcf')
    parser.add_argument('-i', dest='infile', help='input vcf file')
    parser.add_argument('-f', dest='filterconf', default="", help='filterfile')
    parser.add_argument('-o', dest='outfile', default=None, help='output file')
    parser.add_argument('-l', dest="filterlist",default="", help='Comma separated input list to include')
    
    options = parser.parse_args()
    fl = []
    if options.filterlist:
        for el in open(options.filterlist,  'rU'):
            fl.append(el.strip().split(','))
        f = Filter(fl)
    else:
        f = Filter(options.filterconf)
    v = vcf.VCFParser(options.infile)
    ostream = open(options.outfile,  'w')
    v.writeheader(ostream)
    for rec in v:
        v.parseinfo(rec)
        v.parsegenotypes(rec)
        if f.retain(rec):
            v.write(ostream, rec)
