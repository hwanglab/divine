#
# COPYRIGHT (C) 2016 Changjin Hong 
#
"""
.. module:: refgene
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing refgene information from an sqlite database

.. moduleauthor:: Changjin Hong (changjin.hong@gmail.com)

Class for accessing gene ontology information from an sqlite database. A geneontology
is an instance of the `Geneontology` namedtuple
"""

from gcn.lib.io import db
from gcn.etc.dbconfig import DBCONFIG
import os

class Geneontology(db.DB):
  """Class to retrieve gene ontology information from Geneontology database"""
  
  def __init__(self):
    """Class Initialization
    Argument
      name (string): Name of the database or filename for the database
    """
    name = "GENEONTOLOGY"
    super(Geneontology,self).__init__()
    if name in DBCONFIG:
      self.load(name=name)
    elif os.path.exists(name):#accepting a database file name directly
      self.load(db=name)
    else:
      raise ValueError('No such database %s' % name)

  def iterate_funsim(self, genes1, genes2, min_score=0.90,min_denominator=2):
    
    geneStr1 = "'"+"','".join(genes1)+"'"
    geneStr2 = "'"+"','".join(genes2)+"'"
    max_limit_entries = 1e6
    
    stmt = "select gene1, gene2, score from funsim \
      where ((gene1 in (%s) and gene2 in (%s)) \
      or (gene1 in (%s) and gene2 in (%s))) \
      and score >= %g \
      and denominator >= %d \
      order by score desc limit %d" % \
      (geneStr1,geneStr2,\
       geneStr2,geneStr1,\
       min_score,\
       min_denominator,\
       max_limit_entries)

    results = self.execute(stmt)
    for row in results:
      yield (row[0], row[1], float(row[2]))
  
  def get_funsim(self,queries,targets,min_score=0.90,min_denominator=2):
    target_enriched = {}
    for g1,g2,sc in self.iterate_funsim(\
                                   sorted(queries),\
                                   sorted(targets),\
                                   min_score=min_score,\
                                   min_denominator=min_denominator):
      if g1 in queries:
        if g2 not in target_enriched:
          target_enriched[(g1,g2)] = sc
      elif g1 not in target_enriched:
        target_enriched[(g2,g1)] = sc
    return target_enriched