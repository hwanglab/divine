#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''
import os
from gcn.lib.utils import lib_utils
from ConfigParser import SafeConfigParser

def gcn_path(entry,section='config'):
	
	sparser = SafeConfigParser()
	
	divine_root_dir = os.environ.get('DIVINE')
	if not divine_root_dir:
		raise EnvironmentError("set DIVINE variable properly!")

	config_fn = os.path.join(divine_root_dir,'gcn','config','divine.conf')
	if not lib_utils.check_if_file_valid(config_fn):
		raise IOError("check if the configuration file[%s] is valid!" % config_fn)
	
	sparser.read(config_fn)
	
	try:
		path2 = sparser.get(section, entry)
		if not path2.startswith('/'):
			path2 = os.path.join(divine_root_dir,path2)
		return path2
	except:
		print 'WARNING:The config file [%s] does not contain an entry [%s] in the section [%s]'%(config_fn,entry,section)
		return None
