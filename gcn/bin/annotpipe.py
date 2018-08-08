"""
.. module:: annotpipe
	:platform: Unix, Windows, MacOSX
	:synopsis: A wraper to call VARANT

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is a wrapper to call VARANT. The inputs are -
1. Unannotated VCF file path
2. Path to create annotated vcf file (Option)
"""
import os
import argparse
import time
import datetime
from gcn.config import lib_config
from gcn.etc.fileconfig import getlogger, FILECONFIG
from gcn.lib.varann.vartype.varant import annotator

class AnnotPline:

	def __init__(self, invcf, max_exac_maf, capture_kit_name, probe_ext_bp, \
					outvcf=None, cosmic_on=False, hgmd_on=False, \
					dblink=False, log_dir=None):
		
		self.invcf = invcf
		ts = datetime.datetime.fromtimestamp(
							time.time()).strftime('%Y%m%d_%H%M')
		
		self.capture_kit_name =  capture_kit_name
		self.probe_ext_bp = probe_ext_bp
		self.cosmic_on = cosmic_on
		self.hgmd_on = hgmd_on
		self.dblink = dblink
		self.max_exac_af = max_exac_maf

		if outvcf:
			self.outvcf = outvcf
		else:
			self.outvcf = os.path.splitext(invcf)[0] + '.varant_%s.vcf' % ts

		if log_dir:
			FILECONFIG['LOGFILE'] = os.path.join(log_dir,'varant_%s.log' % ts)
		elif outvcf:
			FILECONFIG['LOGFILE'] = os.path.splitext(self.outvcf)[0] + \
									'.varant_%s.log' % ts
		else:
			FILECONFIG['LOGFILE'] = os.path.splitext(invcf)[0] +\
									'.varant.%s.log' % ts

		self.logger = getlogger()
		status = self._check_env()  # Check Environment variable setting
		if status == 1:
			d = 'There seems to be problem with setting the '\
				'environment variables..'
			self.logger.error(d)
			print d
		else:
			self.logger.info('Environment variable check successful..')

	def _check_env(self):
		flag = 0
		envs = ['GCN', 'GCN_DATA_DIR', 'GCN_DB_DIR']
		for env in envs:
			loc = lib_config.gcn_path(env)
			if loc:
				if not os.path.exists(loc):
					self.logger.error('Path set to %s environment does not '
										'exist..' % env)
					flag = 1
			else:
				self.logger.error('%s environment variable is not set..' % env)
				flag = 1
		return flag

	def annotate_varant(self):
		'Annotate the VCF using VARANT'
		self.logger.info('Input file = %s, Output file = %s' %
						 (self.invcf, self.outvcf))
		
		options = [self.capture_kit_name, self.probe_ext_bp, self.cosmic_on, self.hgmd_on, self.max_exac_af, self.dblink]
		
		annotator.main(self.invcf, self.outvcf, self.logger, options)
		
		d = 'Annotation complete [%s;%s].'%(self.invcf,self.outvcf)
		
		print d
		self.logger.info(d)

def main():
	"""Main script to annotate VCF file"""

	desc = 'VCF Annotator'
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('-i', '--inputVCF', dest='invcf', type=str,
						help='Input VCF file path')
	parser.add_argument('-c', '--captureKitName', dest='capture_kit_name', type=str, required = False, default = None,
						help='Input Capture kit name [None]')
	parser.add_argument('-e', '--probeFlankingB', dest='probe_flanking_bp', type=int, required = False, default = 50,
						help='Input Capture kit probe extension bp [50]')
	parser.add_argument('-f', '--maxExacMAF', dest='max_exac_maf', type=float, required = False, default = 0.05,
						help='max ExAC allele frequency [0.05]')
	parser.add_argument('--cosmic', action='store_const', dest='cosmic', required=False, default=False, const=True, help='enable COSMIC database [False]')
	parser.add_argument('--hgmd', action='store_const', dest='hgmd', required=False, default=False, const=True, help='enable HGMD (requires a license)[False]')
	parser.add_argument('--dblink', action='store_const', dest='dblink', required=False, default=False, const=True, help='upload variants to database[False]')
	parser.add_argument('-o', '--outVCF', dest='outvcf', type=str,
						default=None, help='Output VCF file path')
	parser.add_argument('-l', '--logDir', dest='logdir', type=str,
						default=None, help='Log directory file path')
	
	args = parser.parse_args()
	
	ap = AnnotPline(args.invcf, args.max_exac_maf, args.capture_kit_name, args.probe_flanking_bp, args.outvcf, args.cosmic, args.hgmd, args.dblink, args.logdir)
	
	ap.annotate_varant()

if __name__ == "__main__":
	main()
