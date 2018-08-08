#!/usr/bin/env python

from gcn.lib.io import anyopen, db
from gcn.lib.utils import lib_utils
import decimal
from collections import namedtuple

def close(cnx):
	if cnx and cnx.open:
		cursor=cnx.cursor()
		cursor.close()
		cnx.close()
		print 'close db'
		
def browse(con,sqlCmd):
	if not con:
		raise RuntimeError('database connection is not established!')
	else:
		#print 'running sql query [%s]'%sqlCmd
		cursor = con.cursor()
		cursor.execute(sqlCmd)
		rows = cursor.fetchall()
		#print 'Done.'
	return rows

def to_file(rows,Header,out_fn,fmode='wb'):
	#check if out_fn can be writable
	fp2 = anyopen.openfile(out_fn,fmode)
	if isinstance(Header,basestring):
		headStr = Header
	else:
		headStr = lib_utils.joined(Header,'\t')
	fp2.write('#%s\n'%headStr)
	
	if len(rows)>0:
		decimal_idx = get_decimal_idx(rows[0])
		fix_record = False
		if len(decimal_idx)>0:
			fix_record = True
			
		for i, r in enumerate(rows):
			r = list(r)
			if fix_record:
				r = reformat_fields(r, decimal_idx)
			fp2.write('%s\n'%lib_utils.joined(r,'\t'))
	fp2.close()
	
def get_decimal_idx(Sqlrow):
	decimal_idx = []
	for i,r in enumerate(Sqlrow):
		if type(r) is decimal.Decimal:
			decimal_idx.append(i)
	return decimal_idx

def reformat_fields(r, decimal_idx):
	for i in decimal_idx:
		r[i] = float(r[i])
	return r

def get_col_names(curex):
	if curex:
		col_names = [tuple[0] for tuple in curex.description]
	else:
		raise RuntimeError('the cursor user provides does not exist!')
	return col_names

def iter_tupfetchall(curex):
	"Return all rows from a cursor as a tuple"

	fields = get_col_names(curex)
	seltup = namedtuple('seltup', fields)
	for row in curex.fetchall():
		yield seltup._make(row)

def dictfetchall(curex):
	"Return all rows from a cursor as a dict"

	fields = get_col_names(curex)
	try:
		rows = curex.fetchall()
		if rows: #to check if it has at least one record
			return [
				dict(zip(curex, row))
				for row in rows
				]
		else:#otherwise, return empty list
			return []
	except:#otherwise, return empty list
		return []
