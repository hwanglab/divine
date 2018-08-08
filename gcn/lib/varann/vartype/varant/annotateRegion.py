#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan (rajgopal.srinivasan@tcs.com); modified by changjin.hong@gmail.com
#
"""
.. module:: annotateRegion
    :platform: Unix, Windows, MacOSX
    :synopsis: Region Annotation

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

For a given genome coordinate this module identifies its region i.e.
Intergenic, Intronic, Exonic, 5'UTR or 3'UTR.
"""
import os
import numpy as np
from gcn.etc import dbconfig, fileconfig
from gcn.lib.utils import lib_utils

class RefGeneUcscTB():
    
    def __init__(self,work_dir='/tmp',refgene_fn=fileconfig.FILECONFIG['REFGENE'],logger=None):

        self.refGene_fn = refgene_fn
        if not os.path.exists(self.refGene_fn):
            raise IOError('refseq file [%s] not exist!'%self.refGene_fn)

        if not os.path.exists(work_dir):
            raise IOError('work directory you provided [%s] not exist!'%work_dir)

        self.work_dir = work_dir
        self.cds_stats = []
        self.boundary = {}
        self.logger = logger
        self.bed_fn = None

    def collapse_bed(self,tmp_bed,job_name,ext_bp):
        msg = 'sorting bed file ... @ %s' % job_name
        lib_utils.msgout('notice', msg)
        if self.logger: self.logger.info(msg)

        tmp_so_bed = os.path.join(self.work_dir, 'refGene_e%d_so.bed' % ext_bp)
        # sort
        lib_utils.sort_tsv_by_col2(tmp_bed, [1, 2, 3], ['V', 'n', 'n'], True, tmp_so_bed)

        msg = 'merging exon coordinates overlapped each other... @ %s' % job_name
        lib_utils.msgout('notice', msg)
        if self.logger: self.logger.info(msg)

        # merge boundaries if any overlapped
        fp = open(tmp_so_bed, 'r')
        fp2 = open(self.bed_fn, 'w')

        chromp, e1p, e2p, annotp = fp.next().rstrip().split('\t')
        e1p = int(e1p)
        e2p = int(e2p)

        wrapup = 1;
        merge = 2
        fp.seek(0)
        for i in fp:
            chrom, e1, e2, annot = i.rstrip().split('\t')
            e1 = int(e1)
            e2 = int(e2)
            if chrom == chromp:
                if e2p < e1:
                    action = wrapup
                else:
                    action = merge
            else:
                action = wrapup

            if action == wrapup:
                fp2.write('%s\t%d\t%d\t%s\n' % (chromp, e1p, e2p, annotp))
                chromp, e1p, e2p, annotp = chrom, e1, e2, annot
            elif action == merge:
                if e2p < e2:
                    e2p = e2
                    annotp += '|%s' % annot
        fp2.write('%s\t%d\t%d\t%s\n' % (chromp, e1p, e2p, annotp))
        fp.close()
        fp2.close()

        os.unlink(tmp_so_bed)

        msg = 'done. @ %s' % job_name
        lib_utils.msgout('notice', msg)
        if self.logger: self.logger.info(msg)


    def create_bed(self, ext_bp=0, reuse=False):

        job_name = 'RefGeneUcscTB.create_bed'
        
        self.bed_fn = os.path.join(self.work_dir,'refGene_e%d_so_merged.bed'%ext_bp)
        
        msg = 'creating a bed file[%s] containing RefGene coding region (cmpl/incmpl/unk) @ %s'%(self.bed_fn,job_name)
        
        lib_utils.msgout('notice',msg)
        if self.logger: self.logger.info(msg)
        
        if reuse and lib_utils.check_if_file_valid(self.bed_fn):
            msg = 'reuse bed file [%s] generated previously @ %s'%(self.bed_fn,job_name)
            lib_utils.msgout('notice',msg)
            if self.logger: self.logger.info(msg)
            return self.bed_fn

        #to get a working directory
        tmp_bed = os.path.join(self.work_dir,'refGene_e%d.bed'%ext_bp)
        
        fp = open(self.refGene_fn,'r')
        fp2= open(tmp_bed,'w')
        for i in fp:
            j=i.rstrip().split('\t')
            chrom = j[2]
            
            for e1,e2 in zip(j[9].split(',')[:-1],j[10].split(',')[:-1]):
                e1_ext=int(e1)-ext_bp
                e2_ext=int(e2)+ext_bp
                fp2.write('%s\t%d\t%d\t%s;%s\n'%(chrom,e1_ext,e2_ext,j[12],j[1]))
        fp2.close()
        fp.close()
        
        self.collapse_bed(tmp_bed,job_name,ext_bp)
        os.unlink(tmp_bed)

        return self.bed_fn
    
    def get_boundary(self,cds_stats=['cmpl','incmpl','unk','none'],ext_bp=0):
        
        job_name = 'RefGeneUcscTB.get_boundary'
        if self.bed_fn is None:
            raise RuntimeError('Bed file should be set first!')

        msg = 'storing coding region boundaries from [%s] @ %s'%(self.bed_fn,job_name)
        lib_utils.msgout('notice',msg)
        if self.logger: self.logger.info(msg)

        maxNumExon = int(1e6)
        fp = open(self.bed_fn,'r')
        chromp, e1, e2, _ = fp.next().rstrip().split('\t')
        j = 0
        fp.seek(0)
        for i in fp:
            chrom, e1, e2, _ = i.rstrip().split('\t')
            if chrom not in self.boundary:
                if j>0:
                    self.boundary[chromp] = np.delete(self.boundary[chromp],range(j,maxNumExon),0)
                self.boundary[chrom] = np.zeros((maxNumExon,2),dtype=int)
                chromp = chrom
                j = 0
                
            self.boundary[chrom][j,0]=int(e1)
            self.boundary[chrom][j,1]=int(e2)
            j += 1
        if j>0:
            self.boundary[chromp] = np.delete(self.boundary[chromp],range(j,maxNumExon),0)
        fp.close()

        msg = 'done. @ %s'%job_name
        lib_utils.msgout('notice',msg)
        if self.logger: self.logger.info(msg)
        
        return self.boundary
    
class RegionAnnotation():

    def __init__(self, spos, epos, ref, alt, refgene_entry):
        self.ref = ref
        self.alt = alt
        self.frames = refgene_entry.frames
        self.strand = refgene_entry.strand
        self.refseq_acc = refgene_entry.refseq_acc
        self.start_trns_pos = int(refgene_entry.str_transcript)
        self.end_trns_pos = int(refgene_entry.end_transcript)
        self.cds_start = int(refgene_entry.str_cds)
        self.cds_end = int(refgene_entry.end_cds)
        self.exon_str = [int(e) for e in refgene_entry.str_exon_pos.split(',')]
        self.exon_end = [int(e) for e in refgene_entry.end_exon_pos.split(',')]

        if self.strand == '+':
            self.spos = spos - 1
            self.epos = epos - 1
        else:
            self.spos = spos
            self.epos = epos

        # For SNP, INDELs
        self.poslist = range(self.spos, self.epos + 1)

        # For MNPs
        if len(self.ref) > 1 and len(self.alt) > 1:
            self.poslist = [ele[0] for ele in zip(self.poslist,
                                                  self.ref, self.alt)
                            if ele[1] != ele[2]]

        self.splice_site_range = 2
        self.region = None
        self.splice_site = None
        self.exons = None
        self.cdna_poslist = []

    def annotate_positive_strand_cds(self, position):
        region = None
        exon_num = 0
        cdna_pos = 0
        cdna_flag = False
        for x, y in zip(self.exon_str, self.exon_end):
            exon_num += 1
            if y < self.cds_start:
                exon_flag = False
                if x <= position < y:
                    region = 'UTR5'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
            elif x <= self.cds_start < y:
                exon_flag = True
                if x <= position < self.cds_start:
                    region = 'UTR5'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
                elif self.cds_start <= position < self.cds_end \
                                                    and position < y:
                    region = 'CodingExonic'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
                elif self.cds_end <= position < y:
                    region = 'UTR3'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
            elif x <= self.cds_end < y:
                exon_flag = False
                if self.cds_end <= position < y:
                    region = 'UTR3'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
                elif x <= position < self.cds_end:
                    region = 'CodingExonic'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
            elif x > self.cds_end:
                if x <= position < y:
                    region = 'UTR3'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
            elif exon_flag is True:
                if x <= position < y:
                    region = 'CodingExonic'
                    cdna_pos += position - x + 1
                    cdna_flag = True
                    break
            cdna_pos += y - x
        if cdna_flag is False:
            cdna_pos = None
        return region, exon_num, cdna_pos

    def annotate_negative_strand_cds(self, position):
        region = None
        self.exon_str.reverse()
        self.exon_end.reverse()
        exon_num = 0
        cdna_pos = 0
        cdna_flag = False
        for x, y in zip(self.exon_end, self.exon_str):
            exon_num += 1
            if y > self.cds_end:
                exon_flag = False
                if x >= position > y:
                    region = 'UTR5'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
            elif x >= self.cds_end > y:
                exon_flag = True
                if x >= position > self.cds_end:
                    region = 'UTR5'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
                elif self.cds_end >= position > self.cds_start \
                                                    and position > y:
                    region = 'CodingExonic'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
                elif self.cds_start >= position > y:
                    region = 'UTR3'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
            elif x >= self.cds_start > y:
                exon_flag = False
                if self.cds_start >= position > y:
                    region = 'UTR3'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
                elif x >= position > self.cds_start:
                    region = 'CodingExonic'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
            elif x < self.cds_start:
                if x >= position > y:
                    region = 'UTR3'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
            elif exon_flag is True:
                if x >= position > y:
                    region = 'CodingExonic'
                    cdna_pos += x - position + 1
                    cdna_flag = True
                    break
            cdna_pos += x - y
        if cdna_flag is False:
            cdna_pos = None
        self.exon_str.reverse()
        self.exon_end.reverse()
        #exon_num = len(self.exon_str) - (exon_num - 1)
        return region, exon_num, cdna_pos

    def annotate_splice_site(self):
        if self.region[-8:] == 'boundary':
            reg1, reg2, _ = self.region.split('_', 2)
            if reg1 in ['UTR5', 'CodingExonic', 'UTR3', 'NonCodingExonic'] \
                     and reg2 in ['CodingIntronic', 'NonCodingIntronic']:
                if self.strand == '+':
                    self.splice_site = 'SpliceDonor'
                elif self.strand == '-':
                    self.splice_site = 'SpliceAcceptor'
            elif reg1 in ['CodingIntronic', 'NonCodingIntronic'] and reg2 in \
                        ['UTR5', 'CodingExonic', 'UTR3', 'NonCodingExonic']:
                if self.strand == '+':
                    self.splice_site = 'SpliceAcceptor'
                elif self.strand == '-':
                    self.splice_site = 'SpliceDonor'

        elif self.region[-8:] == 'Intronic':
            exon_end = self.exon_end
            exon_str = self.exon_str
            intron_cnt = 0
            del exon_str[0]
            del exon_end[-1]
            for intr_start, intr_end in zip(exon_end, exon_str):
                intron_cnt += 1
                if self.strand == '+':
                    if intr_start <= self.poslist[0] < intr_end and \
                                    intr_start <= self.poslist[-1] < intr_end:
                        dnr_dist = self.poslist[0] - intr_start + 1
                        acpt_dist = intr_end - self.poslist[-1]
                        if dnr_dist <= self.splice_site_range:
                            dnr_dist = 0
#                         else:
#                             dnr_dist -= 2
                        if acpt_dist <= self.splice_site_range:
                            acpt_dist = 0
#                         else:
#                             acpt_dist -= 2
                        break
                elif self.strand == '-':
                    if intr_start < self.poslist[0] <= intr_end and \
                                    intr_start < self.poslist[-1] <= intr_end:
                        acpt_dist = self.poslist[0] - intr_start
                        dnr_dist = intr_end - self.poslist[-1] + 1
                        if dnr_dist <= self.splice_site_range:
                            dnr_dist = 0
#                         else:
#                             dnr_dist -= 2
                        if acpt_dist <= self.splice_site_range:
                            acpt_dist = 0
#                         else:
#                             acpt_dist -= 2
                            break
            if dnr_dist == 0:
                self.splice_site = 'SpliceDonor'
            elif acpt_dist == 0:
                self.splice_site = 'SpliceAcceptor'
            else:
                self.splice_site = \
                            'Exon%s_SpliceDonor_%s__Exon%s_SpliceAcceptor_%s' \
                            % (intron_cnt, dnr_dist, intron_cnt + 1, acpt_dist)

    def annotate_nc(self, position, strand):
        region = None
        exon_num = 0
        cdna_pos = 0
        if strand == '+':
            for x, y in zip(self.exon_end, self.exon_str):
                exon_num += 1
                if y <= position < x:
                    region = 'Exonic'
                    cdna_pos += position - y + 1
                    break
                cdna_pos += x - y
        elif strand == '-':
            self.exon_str.reverse()
            self.exon_end.reverse()
            for x, y in zip(self.exon_end, self.exon_str):
                exon_num += 1
                if y < position <= x:
                    region = 'Exonic'
                    cdna_pos += x - position + 1
                    break
                cdna_pos += x - y
            self.exon_str.reverse()
            self.exon_end.reverse()
            exon_num = len(self.exon_end) - (exon_num - 1)
        if region is None:
            region = 'Intronic'
            exon_num = None
            cdna_pos = 0
        return region, exon_num, cdna_pos

    def annotate(self):
        reglist = []
        exons = []
        for position in self.poslist:
            reg = None
            if self.strand == '+':
                if position < self.start_trns_pos or \
                                position >= self.end_trns_pos:
                    reg = 'Intergenic'
                    exon_num = None
                    cdna_pos = 0
                else:
                    if self.cds_start != self.cds_end: #check if start/stop codon is known/compl?
                        reg, exon_num, cdna_pos = \
                                    self.annotate_positive_strand_cds(position)
            elif self.strand == '-':
                if position - 1 < self.start_trns_pos or \
                            position - 1 >= self.end_trns_pos:
                        reg = 'Intergenic'
                        exon_num = None
                        cdna_pos = 0
                else:
                    if self.cds_start != self.cds_end: #check if start/stop codon is known/compl?
                        reg, exon_num, cdna_pos = \
                                    self.annotate_negative_strand_cds(position)

            if reg is None and (self.refseq_acc.startswith('NR_') or
                            (len(list(set(self.frames.split(',')))) == 1 and
                                list(set(self.frames.split(',')))[0] == '-1')): #check if this tx is uncharacterized (start/stop coding is unknown or non-coding RNA)
                t, exon_num, cdna_pos = self.annotate_nc(position, self.strand)
                reg = 'NonCoding' + t

            if reg is None:
                reg = 'CodingIntronic'
                exon_num = None
                cdna_pos = 0

            if reg not in reglist:
                reglist.append(reg)
            if exon_num and exon_num not in exons:
                exons.append(exon_num)

            if cdna_pos:
                self.cdna_poslist.append(cdna_pos)

        self.cdna_poslist.sort()
        if len(reglist) > 1:
            self.region = "_".join(reglist) + '_boundary'
        else:
            self.region = reglist[0]
            
        if exons:
            exons = [str(e) for e in exons]
            self.exons = '__'.join(exons)
            
        self.annotate_splice_site()
