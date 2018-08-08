#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan and modified by changjin.hong@gmail.com
#
"""
.. module:: preprocess_clinvar
		:platform: Unix, Windows, MacOSX
		:synopsis: Transparent opening of compressed and uncompressed files

.. moduleauthor:: ; changjin.hong@gmail.com
"""
import os, re

NA, GENE, POS, ALT, OVR = range(5)
BENIGN, VUS, PATHOGENIC = range(-1, 2)
clitag = {BENIGN: 'benign', VUS: 'vus', PATHOGENIC: 'pathogenic'}

def classify_clnsig(clinvar_sigs):
	if '2' in clinvar_sigs:
		if '5' in clinvar_sigs:
			return VUS
		else:
			return BENIGN
	elif '5' in clinvar_sigs:
		return PATHOGENIC
	elif '3' in clinvar_sigs:
		if '4' in clinvar_sigs:
			return VUS
		else:
			return BENIGN
	elif '4' in clinvar_sigs:
		return PATHOGENIC
	else:
		return VUS

def get_clinvar_class(clnsigs):
	"""
	:param clnsigs: list of clinvar significance (e.g, annotator.ClinvarReport) 
	:return: voting for its clinical significance
	"""
	clinvar_sigs = []
	if clnsigs:
		for sigs in clnsigs:
			sigs = sigs.split('|')
			for sig in sigs:
				if sig == '.': continue
				clinvar_sigs.append(sig)
		return classify_clnsig(set(clinvar_sigs))
	else:
		return VUS
