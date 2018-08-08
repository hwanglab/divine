#!/usr/bin/env python
'''
author: changjin.hong@nationwidechildren
date: 10/20/2014
objective: to evaluate mapping/variant call quality based on GIAB data
TODO: 1. to implement // comp
'''

import os, sys, argparse, shutil, pickle, time
from multiprocessing import Queue, Process
from Queue import Empty
import gcn.lib.utils.lib_utils as lib_utils
import lib_geneontology

def split_file(file2,ncpu,min_num_lines=100):
	print 'splitting [%s] into %d pieces'%(file2,ncpu)
	D,fname,fbase,fext=lib_utils.separateDirFn2(file2)
	out_prefix='%s/%s_parts'%(D,fbase)
	
	#compute the number of lines per a splitted file
	N = lib_utils.count_num_lines(file2)
	L = round(1.*N/ncpu)+1
	if L<=min_num_lines:
		file2b = '%s_00'%out_prefix
		print 'the input file size is too small to split.'
		shutil.copyfile(file2,file2b)
		return [file2b]
	else:
		cmd = 'split -d -l %d %s %s'%(L,file2,out_prefix)
		lib_utils.gen_msg_time('cmd',cmd,'split_file') #debug
		upair_fns = []
		for j in range(ncpu):
			part_fn='%s%02d'%(out_prefix,j)
			if lib_utils.check_if_file_valid(part_fn):
				upair_fns.append(part_fn)
			else:
				lib_utils.gen_msg_time('error','file[%s] does not exist'%part_fn,'split_file')
		print 'done.'
		return upair_fns

#=======================
def run_par_defs(function,q):
	while True:
		try:
			x=q.get(block=False)
			function(x)
		except Empty:
			break
		
# ================
def run_par_defs_wait(function,work_queue,ncpu):
	processes=[Process(target=run_par_defs, args=(function,work_queue,)) for j in range(ncpu)]
	for t in processes:
		t.start()
	for t in processes:
		t.join()
	time.sleep(0.5)
	
def par_go_sim(input_vec):
	[binD,go_dir,workD,query_fn,sim_method,out_pyv,debug] = input_vec
	#setup an class to run gene-ontology enrichment
	print 'computing GO sim on [%s] ...'%query_fn
	cGO = lib_geneontology.Fastsemsim(go_dir,work_dir=workD,bin_dir=binD,debug=debug)
	
	cGO.query_fn = query_fn
	cGO.map_uniprot_to_refgene()
	cGO.run_sim(workD,sim_method)

	funSim_fn = "%s/%s_funSim3.tsv" % (workD,cGO.query_fn.split('/')[-1])
	_ = cGO.get_funsim_mat(funSim_fn)
	
	fp2 = open(out_pyv,'wb')
	print 'registering output file into pyv[%s]'%out_pyv
	pickle.dump([funSim_fn,cGO.gosim_fns],fp2)
	fp2.close()
	print 'done.'

def main():
	parser = argparse.ArgumentParser(description="build GO using fastsemsim [changjin.hong@nationwidechildrens.org]")
	parser.add_argument('--version', action='version', version='%(prog)s 0.2')
	parser.add_argument('-b', action='store', dest='fastsemsim_bin', required=True, help='')
	parser.add_argument('-G', dest='resource_dir', required=True, help='gene ontology resource dir')
	parser.add_argument('-s', action='store', dest='sim_method', required=False, default='SimRel', help='')
	parser.add_argument('-S', action='store', dest='method_id', required=False, default=1, type = int, help='assign an integer number for sim_method[1]')
	
	parser.add_argument('--prep', action='store_const', dest='prepare_input', required=False, default=False, const=True, help='prepare pairs of uniprot query files[False]')
	
	parser.add_argument('-q', action='store', dest='upair_fn', required=True, help='')
	parser.add_argument('-o', action='store', dest='out_prefix', required=True, help='')
	parser.add_argument('-n', action='store', dest='ncpu', required=False, type = int, default = 1, help='specify the number of cpus to utilize')
	
	args=parser.parse_args()
	
	workD = '%s_work'%args.out_prefix
	if not os.path.exists(workD):
		os.makedirs(workD)
		
	if args.prepare_input:
		cGO = lib_geneontology.Fastsemsim(args.resource_dir,workD,args.fastsemsim_bin,False)
		if lib_utils.check_if_file_valid(args.upair_fn):
			sys.exit(1)
		cGO.gen_all_pairs_file(args.upair_fn)
		upair_fns = split_file(args.upair_fn,args.ncpu)
		print 'check %s'%upair_fns
		print 'then, run this program again to cal go sim on the splitted input files ...' 
		sys.exit(0)
	else:
		#split
		upair_fns = split_file(args.upair_fn,args.ncpu)
		ncpu = len(upair_fns)

	if ncpu >1 :
		work_queue = Queue()

	#start for loop to dispatch jobs //
	out_pyvs = []
	for upair_fn in upair_fns:
		print 'dispatching fastsemsim on [%s] ...'%upair_fn
		out_pyv = '%s_out.pyv'%upair_fn
		input_vec = [args.gosim_dir,args.resource_dir,workD,upair_fn,args.sim_method,out_pyv,True]
		if ncpu == 1:
			par_go_sim(input_vec)
		else:
			work_queue.put(input_vec)
		out_pyvs.append(out_pyv)
		print 'done'
	
	print 'job dispatch done. now wait...'
	
	if ncpu > 1:
		run_par_defs_wait(par_go_sim,work_queue,ncpu)
	
	#collect outputs via pyv
	merged_funSim_fn,merged_goSim_fn = collect_outputs(args.out_prefix,args.sim_method,args.method_id,out_pyvs)
	
	if True:
		print 'cleanup work_dir files...'
		lib_utils.unlink_fns(upair_fns)
		lib_utils.unlink_fns(out_pyvs)
		#shutil.rmtree(workD)
	print 'done.'

def reformat_go_sim_fns(go_sim_fns,out_fn,method_id=1): #method_id (1) means SimRel
	
	suflabs=['BP','MF','CC']
	fp2=lib_utils.open2(out_fn,'w')
	v = 0
	for key,go_sim_fn in go_sim_fns.iteritems(): #BP,MF,CC
		print 'appending root node at the end of [%s]'%go_sim_fn
		suflab = suflabs[v]
		fp = lib_utils.open2(go_sim_fn,'r')
		go_sim_fn2 = lib_utils.file_tag2(go_sim_fn,'category',None)
		fp.next() #strip off head
		for i in fp:
			j=i.rstrip().split('\t')
			if len(j)==2:
				j.append('-1.')
			fp2.write('%s\t%s\n'%('\t'.join(j), suflab)) #uniprot1,uniprot2,score,BP
		fp.close()
		print 'done.'
		v += 1
	fp2.close()
	
	print 'sorting...'
	#to get temporary file to sort
	out_fn2 = lib_utils.file_tag2(out_fn,'sort',None)
	temp_sort_dir,_,_,_ = lib_utils.separateDirFn2(out_fn)
	lib_utils.sort_tsv_by_col2(out_fn,[1,2,4],['V','V','V'],False,out_fn2,temp_dir=temp_sort_dir)
	os.rename(out_fn2,out_fn)
	print 'done.'
	
	#groupping
	
	print 'collapsing GO sim scores to make the format easier to import SQL [%s] ...'%out_fn
	out_fn2 = lib_utils.file_tag2(out_fn,'dense',None)
	fp2 = lib_utils.open2(out_fn2,'w')
	#heads = '#uniprot1\tuniprot2\tscore_mode\tBP\tMF\tCC\tmethod_id'
	#fp2.write('%s\n'%heads)

	fp=lib_utils.open2(out_fn,'r')

	visit1 = True
	idx = {'BP':0,'MF':1,'CC':2}
	prev_key = None
	gosim_holder = ['-1','-1','-1'] #-1 means N/A
	
	for i in fp:
		prot1,prot2,score,go_class = i.rstrip().split('\t')
		key = '%s\t%s'%(prot1,prot2)
		if key != prev_key:
			if visit1:
				visit1 = False
			else: #wrap up
				fp2.write('%s\t%s\t%d\n'%(prev_key,'\t'.join(gosim_holder),method_id))
				gosim_holder = ['-1','-1','-1']
			gosim_holder[idx[go_class]] = score
			prev_key = key
		else: #keep storing values
			gosim_holder[idx[go_class]] = score
	fp.close()
	
	#don't forget the last entry
	fp2.write('%s\t%s\t%d\n'%(prev_key,'\t'.join(gosim_holder),method_id))
	fp2.close()
	os.rename(out_fn2,out_fn)
	
	print 'done.'

def collect_outputs(out_prefix,sim_method,method_id,out_pyvs):
	funSim_fns = []
	gosim_reformat_fns = []
	
	print 'collecting job outputs...'
	for pyv in out_pyvs:
		fp = open(pyv,"rb")
		funSim_fn,gosim_fns = pickle.load(fp)
		funSim_fns.append(funSim_fn)
		fp.close()

		gosim_reformat_fn = lib_utils.file_tag2(gosim_fns['biological_process'],'merged',None)
		reformat_go_sim_fns(gosim_fns,gosim_reformat_fn,method_id)
		gosim_reformat_fns.append(gosim_reformat_fn)
	
	print 'merging splitted outputs %s,%s...'%(funSim_fns,gosim_reformat_fns)
	merged_funSim_fn = '%s_funSim.tsv.gz'%out_prefix
	lib_utils.merge_multiple_gz_in_order(funSim_fns,False,merged_funSim_fn)
	
	merged_goSim_fn = '%s_%s.tsv.gz'%(out_prefix,sim_method)
	lib_utils.merge_multiple_gz_in_order(gosim_reformat_fns,False,merged_goSim_fn)
	print 'check [%s,%s]'%(merged_funSim_fn,merged_goSim_fn)
	print 'done.'
	
	return merged_funSim_fn,merged_goSim_fn

if __name__ == '__main__':
	main()
