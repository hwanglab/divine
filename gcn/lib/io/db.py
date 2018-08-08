#
# COPYRIGHT (C) 2002-2011 Rajgopal Srinivasan
#
"""
.. module:: db
	:platform: Unix, Windows, MacOSX
	:synopsis: Class for accessing sqlite databases

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com); modified by changjin.hong@gmail.com

Class for creating and accessing sqlite databases
"""

from __future__ import absolute_import
import sqlite3
from collections import namedtuple
from gcn.etc.dbconfig import DBCONFIG
from gcn.etc import fileconfig
from gcn.lib.io import fileutils
import os

Column = namedtuple('Column', 'name type notnull default pk')

class DB(object):

	__slots__ = ['db', 'conn']

	"""Class to work with sqlite3 databases

	An instance of this class has 2 attributes
		- db   (string): Name of the database
		- conn (object): A connection to the database
	"""

	def __init__(self):
		"""Class initialization"""
		self.db = self.conn = None

	def load(self, db=None, name=None):
		"""Load a database by providing either the sqlite database filename or
		the name of the resource. When a resource name is provided the
		corresponding file name is looked up from the `DBCONFIG` dictionary
		in the `gcn.etc.dbconfig` module.  It is an error to specify both the
		filename and the resource name.

		Args:
			db (str): File name of sqlite database
			name (str): Name of resource for which database is to be loaded
		"""
		if self.db:
			return

		if db and name:
			raise ValueError('Only one of *db* or *name* should be specified')

		if name:
			db = DBCONFIG[name.upper()]['name']
		self.db = db
		self.conn = sqlite3.connect(self.db, check_same_thread=False)
		self.conn.text_factory = str

	def set_version(self, dbname, file_date, version, entry_cnt):
		"""Set the version information in release_info table of DB.

		Args:
			dbname (str):    Database name.
							 Eg. dbsnp, kgdb, espdb, exacdb, refgene etc
			file_date (str):    Last modified date of the raw file.
								Format-YYYY-MM-DD
			version (str):    Version of Database.
								Note - For some databases there may not be
								exact version number so for those version
								is represented by 'vYYYY_MM' i.e. year and month
								when the raw file was last modified.
								Eg. for refgene version = v2013_11
									for omim version = v2013_11
									for dbSNP version = v137
									etc.
			entry_cnt (int):    Number of entries in DB
		"""
		if not self.db:
			raise AttributeError('Database had not yet been loaded')
		cur = self.conn.cursor()
		if 'release_info' not in self.tablenames():
			stmt1 = """CREATE TABLE release_info (dbname TEXT,
												file_date TEXT,
												version TEXT,
												entry_count INTEGER)"""
			cur.execute(stmt1)
		else:
			stmt2 = "DELETE FROM release_info where dbname='%s'" % dbname
			cur.execute(stmt2)
		stmt3 = "INSERT INTO release_info VALUES ('%s', '%s', '%s', %d)" % \
						(dbname, file_date, version, entry_cnt)
		cur.execute(stmt3)
		cur.close()

	def get_version(self):
		"""Returns the version information from release_info table of DB

		Returns:
			version_info (list):    List of named tuples.
									This is list because some databases are
									compiled from more than one database.
									For Eg. clnsnp database is compiled from
									NCBI-GAD and NHGRI-GWAS. So this database
									will hold the version details of these
									two databases as well.
		"""

		Version = namedtuple('Version', ['dbname', 'file_date',
										 'version', 'entry_count'])
		if not self.db:
			raise AttributeError('Database had not yet been loaded')
		stmt = 'SELECT * FROM release_info'
		version_info = []
		for rec in self.execute(stmt):
			version_info.append(Version._make(rec))
		return version_info

	def tablenames(self):
		"""List of all tables in the database"""

		if not self.db:
			raise AttributeError('Database had not yet been loaded')

		cur = self.conn.cursor()
		cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
		result = cur.fetchall()
		cur.close()
		return [el[0] for el in result]

	def hastable(self, name):
		"""Check if database has a table with the given name

		Args:
			name (str): Table name
		Returns:
			True if table is present and False otherwise
		"""
		return name in self.tablenames()

	def schema(self, tablename):
		"""Retrieve schema for a table

		Args:
			tablename (str): Name of the table

		Returns:
			List of columns in the table. Each column is a namedtuple with the
			following fields
				- name (str): Column name
				- type (str): Value type
				- notnull (bool): Whether the column can have Null values
				- default : Default value for the column
				- pk : Whether column is or part-of primary key
		"""
		cur = self.conn.cursor()
		cur.execute('pragma table_info(%s)' % (tablename,))
		result = cur.fetchall()
		cur.close()
		return [Column._make(r[1:]) for r in result]

	def createtable(self, schema, drop=False):
		"""Create a table.

		Args:
			schema (str): Schema for the table
			drop (bool):  Drop the table if it already exists
		"""
		args = ' '.join(schema.strip().splitlines(0)).split()
		name = args[2].strip()
		if self.hastable(name):
			if not drop:
				return
			cur = self.conn.cursor()
			cur.execute('drop table %s' % (name,))
			self.conn.commit()
			cur.close()
		cur = self.conn.cursor()
		cur.execute(schema)
		cur.close()

	def createdb(self, inname, outname, force):
		"""Create the database

		Args:
			inname (str): Name of the input file from which the database
							will be created
			outname (str): Name of sqlite3 database
			force (bool): If True overwrite existing database even if it
							is newer than `inname`
		"""
		if os.path.exists(outname):
			newer = fileutils.file_newer(inname, outname)
		else:
			newer = True

		self.logger = fileconfig.getlogger()

		if not newer and not force:
			self.logger.info('Not Updating. %s database already uptodate' \
								 % inname)
		else:
			self.inname = inname
			self.outname = outname
			self._makedb()
		return

	def make_namedtuple(self, tablename, name):
		"""Create a nameduple to represent a row in a table.

		Args:
			tablename (str): Name of table for which namedtuple is to be
							 created
			name (str): Name of the namedtuple class

		Returns:
			A namedtuple class to represent a row in the table `tablename`
		"""
		cur = self.conn.cursor()
		cur.execute('pragma table_info(%s)' % (tablename,))
		result = cur.fetchall()
		cur.close()
		columns = [r[1].replace('info_', '') for r in result]
		return namedtuple(name, columns)

	def execute(self, stmt, args=None):
		"""Execute an SQL query against the database

		Arguments:
			stmt (string): The sql statement
			args (list/tuple): Arguments for the SQL statement
		Returns:
			An iterator on the results
		"""
		try:
			cursor = self.conn.cursor()
		except AttributeError:
			raise AttributeError('DB has not been loaded')
		else:
			if args is None:
				return cursor.execute(stmt)
			else:
				return cursor.execute(stmt, args)

	def executemany(self, stmt, args):
		"""Execute many SQL queries against the database

		Arguments:
			stmt (string): The sql statement
			args (list/tuple): Arguments for the SQL statement
		Returns:
			An iterator on the results
		"""
		try:
			cursor = self.conn.cursor()
		except AttributeError:
			raise AttributeError('DB has not been loaded')
		else:
			return cursor.executemany(stmt, args)
