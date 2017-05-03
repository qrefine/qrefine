from __future__ import division

import os,sys,time
import libtbx.load_env
from cluster_wrapper import Clusterer,Result
from test_reg_00_base import test_base
qrefine_path = libtbx.env.find_in_repositories("qrefine")

qr_reg_data = os.path.join(qrefine_path, "regression/datasets/cluster")

class test_chunk(test_base):

  def __init__(self):
    self.pdbs = self.pdbs(qr_reg_data)
    self.func = Clusterer()

  def check_assertions(self,result):
    """
    Exercise to test clustering indices, chunk indices, and the number of atoms in each chunk. 
    """
    print "checking result",result
    if self.db.old.find_one({"pdb_code": result.pdb_code}) is None:
      self.insert(result)
    else:
      result_db = self.db.old.find_one({"pdb_code": result.pdb_code})
      assert result_db['chunk']       == result.clusters
      assert result_db['chunks']      == result.chunks
      assert result_db['chunk_sizes'] == result.chunk_sizes

  def insert(self,result):
    self.db.old.insert_one({"pdb_code"      : result.pdb_code,
                     "chunk"           : result.clusters,
                     "chunks"          : result.chunks,
                     "chunk_sizes"     : result.chunk_sizes,
                     "num_of_chunks"   : result.num_of_chunks,
                     "num_of_clusters" : result.num_of_clusters})