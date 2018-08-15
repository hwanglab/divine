#!/usr/bin/env python
'''
author: changjin.hong@gmail.com
'''

import os, sys, urllib, argparse, shutil, imp, re
import subprocess as sp

VERSION = '0.1.2'
author_email = 'hongc2@ccf.org'

import hashlib

def md5Checksum(filePath):
	with open(filePath, 'rb') as fh:
		m = hashlib.md5()
		while True:
			data = fh.read(8192)
			if not data:
				break
			m.update(data)
		return m.hexdigest()

def cmd_exists(cmd):
	return sp.call("type " + cmd, shell=True,\
		stdout=sp.PIPE, stderr=sp.PIPE) == 0

def has_module(module):
	try:
		imp.find_module(module)
		print '[%s] is in path or installed'%module
		found = True
	except ImportError:
		print '[%s] is not in path or installed'%module
		found = False
	return found
	
def syscmd(cmd):
	subproc = sp.Popen(cmd, stdout=sp.PIPE, shell=True)
	output,error = subproc.communicate()
	if output:
		output=output.strip()
	if error:
		error=error.strip()
	exit_code = subproc.wait()
	if exit_code>0:
		raise RuntimeError('command[%s] failed;%s'%(cmd,error))
	return output, error
	
def check_py_ver():
	verStr = '%d.%d'%(sys.version_info[0],sys.version_info[1])
	return verStr

def check_network_on():
	try :
		stri = "https://www.google.com"
		data = urllib.urlopen(stri)
		print "Connected"
		return True
	except:
		print "not connected"
		return False

def get_pylib_prefix():
	return 'python%s'%check_py_ver()

class Setup():
	def __init__(self,url_fn=None):
		#to get the current directory
		self.install_dir = os.path.dirname(os.path.abspath(__file__))
		self.py_libs_dir = os.path.join(self.install_dir,'python_libs','pkg')
		self.py_libs_prefix = ['fastsemsim-0.9.4','hgvs','hpo_similarity','genmod']
		self.py_libs = ['fastsemsim','pyhgvs','hpo_similarity','genmod']
		
		self.py_modules = ['ConfigParser','backports-abc','html5lib==0.999999999',\
			'backports.ssl-match-hostname','certifi','decorator',\
				'matplotlib','networkx','nose','numpy','pandas','pygr',\
					'pyparsing','pysam','python-dateutil','pytz','scipy',\
						'scikit-learn','singledispatch','six','tornado','xlwt','dill']

		self.pylib_prefix = 'python_libs'
		self.pylib = os.path.join(self.install_dir,self.pylib_prefix)
		self.pylib_build_prefix = 'site-packages'
		self.pylib_build = os.path.join(self.pylib,'lib',get_pylib_prefix(),self.pylib_build_prefix)
		
		if not os.path.exists(self.pylib_build):
			os.makedirs(self.pylib_build)

		self.gcn = 'gcn'
		self.gcndata = 'gcndata'
		self.gcndb = 'gcndb'
		self.gcnlog = 'gcn/logs/gcn.log'
		log_dir = os.path.join(self.install_dir,'gcn','logs')
		if os.path.exists(log_dir):
			os.makedirs(log_dir)

		self.user_has_pythonpath=os.environ.get('PYTHONPATH',None)
		
		if url_fn is None:
			self.url_fn = os.path.join(self.install_dir,'resource','support_pkg_urls.txt')
		else:
			self.url_fn = url_fn

		if not os.path.exists(self.url_fn):
			raise RuntimeError('a config file[%s] for resource packages does not exist!'%self.url_fn)

	def install_module_by_pip(self,pkg_name):
		print 'installing %s ...'%pkg_name
		cmd = "pip install --user %s"%pkg_name
		syscmd(cmd)
	
	def download_data(self,ucategory=None,keep_download_files=True):
		
		print "downloading Divine resource files [%s] ..."%ucategory
		
		if not cmd_exists('wget'):
			raise RuntimeError('wget should be in path!')
		
		#to access a file containing url file
		dn_dir = os.path.join(self.install_dir,'resource','v%s'%VERSION)
		dn_fns = []
		
		if not os.path.exists(dn_dir):
			os.makedirs(dn_dir)
			
		fp = open(self.url_fn)
		for i in fp:
			if i.startswith('#'):continue
			category,res_name,res_url,org_md5 = i.strip().split('\t')
			
			if ucategory and (ucategory != category):
				continue
			dn_fn = os.path.join(dn_dir,res_name)

			cmd = 'wget -O %s -c %s'%(dn_fn,res_url)
			dn_fns.append(dn_fn)
			syscmd(cmd)

		fp.close()
		print "Done."

		if (False):
			print "Check if md5 is matched ..."
			your_md5 = md5Checksum(dn_fn)
			if your_md5 != org_md5:
				raise RuntimeError('%s: Expected MD5 [%s] is not same as the one you generated[%s]'%(dn_fn,org_md5,your_md5))

		print "extracting resource files (it may take a couple of hours to complete the installation and be patient) ...\nAfter the installation, add some environment variables to your bash shell configuration file (~/.bashrc or ~/bash_profile).\nFollow the instruction at the end of the installation"

		extracted = []
		#to backup data
		#tar cvJf - {gcndb,gcn/bin/prioritize/examples} | split -b 4GB - divine-0.1.2_data.tar.xz.
		#to extract multi-volume tar.bz2
		for fn in dn_fns:
			if 'tar.xz' in fn:
				mObj = re.search(r'(\S+)\.tar\.xz\.\S+',fn)
				if mObj:
					fbase='%s.tar.xz'%mObj.group(1)
					arch_path = os.path.join(dn_dir,fbase)
					cmd = "cat %s.* | tar -xvJf - -C %s"%(arch_path,self.install_dir)
				else:
					arch_path = os.path.join(dn_dir,fn)
					cmd = "tar -xvJf %s -C %s"%(arch_path,self.install_dir)

				if arch_path not in extracted:
					extracted.append(arch_path)
					print "extracting [%s] ..."%arch_path
					syscmd(cmd)
		print "Done."
		
		if not keep_download_files:
			for fn in extracted:
				syscmd('rm -rf %s*'%fn)
	
	def uninstall_resource(self):
		shutil.rmtree(os.path.join(self.install_dir,'gcndb'))
		shutil.rmtree(os.path.join(self.install_dir,'gcndata'))
		
	def install_python_libs(self):
		#to make sure that PYTHONPATH includes the directory to install python modules
		
		for module in self.py_modules:
			if not has_module(module):
				print 'installing [%s] ...'%module
				self.install_module_by_pip(module)
				print 'done.'
		
		#download python_lib zip file
		self.download_data(ucategory='MODULES')

		for module in self.py_libs_prefix:
			module_path = os.path.join(self.py_libs_dir,module)
			setup_py = os.path.join(module_path,'setup.py')
			if os.path.exists(module_path) and os.path.exists(setup_py):
				cmd = "export PYTHONPATH=%s:%s"%(self.pylib,self.pylib_build)
				if self.user_has_pythonpath:
					cmd += ":$PYTHONPATH"
					
				cmd += ";cd %s;python ./setup.py install --prefix=%s"%(module_path,self.pylib)
				print 'installing [%s] ...'%module
				print cmd
				syscmd(cmd)
				print 'done.'
			else:
				raise IOError('check if [%s] contains all necessary files;enable the option, --update_db'%module_path)

		print 'python modules are installed successfully in [%s]'%self.pylib
	
	def uninstall_python_libs(self):

		for module in self.py_modules:
			if has_module(module):
				cmd = "pip uninstall -y %s"%module
				try:
					print 'uninstalling [%s]'%module
					syscmd(cmd)
					print 'done.'
				except:
					print "%s may not be installed previously."%module
			else:
				print "%s may not be installed previously."%module

		path2del = os.path.join(self.pylib,'lib')
		if os.path.exists(path2del):
			print 'cleaning up python module lib [%s] ...'%path2del
			shutil.rmtree(path2del)
			print 'done.'
			
		path2del = os.path.join(self.pylib,'bin')
		if os.path.exists(path2del):
			print 'cleaning up python module lib [%s] ...'%path2del
			shutil.rmtree(path2del)
			print 'done.'

		print 'python modules are uninstalled successfully in [%s]'%self.pylib
	
	def msg_config(self,setup_mode='install'):
		
		#detect the user default shell type
		shell_bin,err_msg = syscmd("echo $0")
		if shell_bin and 'csh' in shell_bin:
			set_keyword = 'setenv'
			set_link = ' '
			shell_cnf_fn = '$HOME/.cshrc'
		else:
			set_keyword = 'export'
			set_link = '='
			shell_cnf_fn= '$HOME/.profile, $HOME/.bash_profile, or $HOME/.bashrc'

		hline = "-------------------------------"

		if setup_mode == 'install':
			desc1='add'
			desc2='to'
		
		elif setup_mode == 'remove':
			desc1='remove'
			desc2='from'
		
		print "\n"
		print '-Open %s your shell script (e.g., %s).'%(desc1,desc2,shell_cnf_fn)
		print "-Then, %s the following variables %s PATH and PYTHONPATH in your shell script if not set(e.g., %s)."%(desc1,desc2,shell_cnf_fn)
		print hline

		msg = "%s DIVINE%s%s\n" % (set_keyword, set_link, self.install_dir)
		msg += "%s PATH=$HOME/.local/bin:$PATH\n"%(set_keyword)
		msg += "%s PATH%s$DIVINE/gcn/bin/prioritize:$PATH\n"%(set_keyword,set_link)
		msg += "%s PYTHONPATH%s$DIVINE"%(set_keyword,set_link)
		
		if self.user_has_pythonpath: msg += ":$PYTHONPATH\n"
		else: msg += "\n"
		
		# msg += "%s PYTHONPATH%s$DIVINE/%s:$PYTHONPATH\n"%(set_keyword,set_link,self.pylib_prefix)
		msg += "%s PYTHONPATH%s$DIVINE/%s/%s/%s/%s:$PYTHONPATH\n"%\
			(set_keyword,set_link,self.pylib_prefix,'lib',get_pylib_prefix(),self.pylib_build_prefix)
		print msg

def main():
	parser = argparse.ArgumentParser(description="Divine setup (v%s) [author:%s]"%(VERSION,author_email))
	parser.add_argument('--install', action='store_const', dest='install', required=False, default=False, const=True, help='install Divine')
	parser.add_argument('--update_db', action='store_const', dest='update_db', required=False, default=False, const=True, help='update_db [False]')
	parser.add_argument('-c', action='store', dest='url_fn', required=False, default=None, help='a file containing URL to download divine resource files')
	parser.add_argument('--uninstall', action='store_const', dest='uninstall', required=False, default=False, const=True, help='uninstall Divine')
	parser.add_argument('--remove_db', action='store_const', dest='remove_db', required=False, default=False, const=True, help='remove_db [False]')

	args = parser.parse_args()

	if (args.install and args.uninstall) or (not args.install and not args.uninstall):
		raise RuntimeError('choose --install or --uninstall')
	
	#to install python libraries
	cs = Setup(url_fn=args.url_fn)
	if args.install:
		if not check_network_on():
			raise RuntimeError('it requires network connection to setup Divine!')
		cs.install_python_libs()
		cs.download_data(ucategory='EXAMPLES')

		if args.update_db:
			cs.download_data(ucategory='DATA')
			cs.download_data(ucategory='DATA_DELTA')
		cs.msg_config('install')
	elif args.uninstall:
		cs.uninstall_python_libs()
		if args.remove_db:
			cs.delete_resource()
		cs.msg_config('remove')
	
if __name__ == '__main__':
	main()
