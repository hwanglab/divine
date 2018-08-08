#!/usr/bin/env python
import os, sys, shutil
import argparse
import codecs, json
import urllib, socket

class Hsa:
	
	def __init__(self,hsa,desc):
		self.hsa = hsa
		self.desc = desc
		self.genes = None
		self.disIds = {} #a value indicates coverage ratio of disease genes

class Kegg:
	
	def __init__(self,temp_dir=None,hsa_fn=None):
		
		self.cHsa = {}
		self.hsa_fn = hsa_fn
		
		self.temp_dir = temp_dir
		if temp_dir:
			if not os.path.exists(temp_dir):
				os.makedirs(temp_dir)
		
		self.kegg_url = "http://rest.kegg.jp/list/pathway/hsa"
		self.hsa_url = "http://togows.dbcls.jp/entry/pathway/%s/genes.json"
		'''
		to get disease list from pathway
				http://togows.org/entry/kegg-pathway/hsa04610/diseases.json
		to get omim from each disease
				http://www.genome.jp/dbget-bin/www_bget?ds:H01434
		to get an edge direction in pathway
				http://rest.kegg.jp/get/hsa04610/kgml
		'''
		
	def get_hsa(self):
		
		self.cHsa = {}
		fp = open(self.hsa_fn,'r')
		fp.next()
		for i in fp:
			hsa,desc,geneStr = i[:-1].split('\t')
			kG = Hsa(hsa,desc)
			kG.genes = geneStr.split(',')
			self.cHsa[hsa] = kG
		fp.close()
	
	def get_hsa_genes(self,hsa,kG):
		try:
			gene_json_fn = '%s/%s.json'%(self.temp_dir,hsa)
			url2 = self.hsa_url%hsa
			print url2 #debug
			urllib.urlretrieve(url2,gene_json_fn)
		except:
			raise RuntimeError('cannot access to url [%s]'%url2)
		
		fp2=codecs.open(gene_json_fn,'r','utf-8')
		hsa_genes = json.load(fp2)
		fp2.close()

		genes = []
		if not hsa_genes:
			print 'kegg pathway [%s] is empty, check url [%s]'%(hsa,url2)
			return genes

		for kgid,desc in hsa_genes[0].iteritems():
			gene = desc.split(';')[0].strip()
			genes.append(gene)
		return genes
	
	def print_pathways(self,out_fn):
		
		timeout = 10
		socket.setdefaulttimeout(timeout)
		
		#1. to get hsa list
		kG2 = []

		hsa_list_fn = '%s/hsa_list.tsv'%self.temp_dir
		print hsa_list_fn #debug

		try:
			urllib.urlretrieve(self.kegg_url,hsa_list_fn)
		except:
			raise RuntimeError('cannot access to url [%s]'%self.kegg_url)
		
		fp = open(hsa_list_fn,'r')
		fpw= open(out_fn,'w')
		fpw.write('#hsa:id\tpathway_desc\tgenes\n')
		for i in fp:
			i = i.strip()
			if i:
				hsa,desc = i.split('\t')
				hsa = hsa.split(':')[1].strip()
				desc = desc.split('Homo sapiens')[0][:-2].strip()

				kG = Hsa(hsa,desc)
				
				#2. to access 
				kG.genes = self.get_hsa_genes(hsa,kG)
				kG2.append(kG)
				
				#3. to print
				fpw.write('%s\t%s\t%s\n'%(kG.hsa,kG.desc,','.join(kG.genes)))
		fp.close()
		fpw.close()
		shutil.rmtree(self.temp_dir) #debug
		return kG2

	def load_pathways(self, kegg_genes_fn):
		fp = open(kegg_genes_fn, 'r')
		heads = fp.next()[1:-1].split('\t')
		hsa_to_genes = []
		for i in fp:
			hsa, desc, genes = i[:-1].split('\t')
			for gene in genes.split(','):
				hsa_to_genes.append([hsa, gene])
		fp.close()
		return heads, hsa_to_genes

def list2dataFrame(myList):
	import pandas as pd
	heads_lite = ['hsa_id', 'gene']
	df = pd.DataFrame.from_records(myList, columns=heads_lite)
	df = df.drop_duplicates()
	df['flag'] = 1
	df = df.pivot(index='gene', columns='hsa_id', values='flag')
	return df

def main():
	parser = argparse.ArgumentParser(description="collect genes per each KEGG pathway [changjin.hong@gmail.com]")
	parser.add_argument('-o', action='store', dest='out_fn', required=True, help='output file')
	
	args=parser.parse_args()
	
	cK = Kegg(temp_dir='/tmp/kegg_pathway')
	cK.print_pathways(args.out_fn)

if __name__ == '__main__':
	main()
