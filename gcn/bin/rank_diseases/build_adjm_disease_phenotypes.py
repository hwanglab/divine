#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''
import os, sys, argparse, re, math, datetime, time, logging
from gcn.lib.utils import lib_hposim, lib_utils

def main():
	parser = argparse.ArgumentParser(description="build_adjm_disease_phenotypes")
	parser.add_argument('-o','--work_dir', dest='work_dir', required=True, help='output directory')
	
	args = parser.parse_args()
	
	cHpo = lib_hposim.HpoSim()
	
	query_dir = '%s/query'%args.work_dir
	lib_utils.ensure_dir(query_dir)
	
	query_hpos = cHpo.disease_hpo_to_file(query_dir)
	
	result_dir = '%s/results'%args.work_dir
	lib_utils.ensure_dir(result_dir)
	
	for disId, query_hpo_fn in query_hpos:
		out_fn = '%s/%s.tsv'%(result_dir,disId)
		if not os.path.exists(out_fn):
			cHpo.run(query_hpo_fn,out_fn)

if __name__ == '__main__':
	main()