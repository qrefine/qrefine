from __future__ import division
import os
import sys
import time
import os.path
import iotbx.pdb
import libtbx.load_env
from Cluster_mp import Clusterer
from libtbx.easy_mp import parallel_map
from qrefine.core.restraints import from_qm
from pymongo import MongoClient


db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
utils_path= os.path.join(qr_path,"utils") 
qr_reg_data = os.path.join(qrefine_path, "regression/datasets/cluster")


def check_assertions(result):
  """
  Exercise to test clustering indices, chunk indices, and the number of atoms in each chunk. 
  """
  if db.old.find_one({"pdb_code": result.pdb_code}) is None:
    insert(result)
  else:
    result_db = db.old.find_one({"pdb_code": result.pdb_code})
    assert result_db['chunk']    == result.clusters
    assert result_db['chunks']      == result.chunks
    assert result_db['chunk_sizes'] == result.chunk_sizes

def insert(result):
  db.old.insert_one({"pdb_code"        : result.pdb_code,
                     "chunk"        : result.clusters,
                     "chunks"          : result.chunks,
                     "chunk_sizes"     : result.chunk_sizes,
                     "num_of_chunks"   : result.num_of_chunks,
                     "num_of_clusters" : result.num_of_clusters})

class Result(object):
  def __init__(self, pdb_code, clusters, chunks, chunk_sizes):
    self.pdb_code        = pdb_code
    self.clusters        = clusters
    self.chunks          = chunks
    self.chunk_sizes     = chunk_sizes
    self.num_of_clusters = len(clusters)
    self.num_of_chunks   = len(chunks)

def run(parallel=True):
  pdbs=[]
  for pdb_file in os.listdir(qr_reg_data):
    pdbs.append(os.path.join(qr_reg_data,pdb_file))
  if(parallel):
    qsub_command  = 'qsub  -N reg_test_cluster -m ae -q qr  -l nodes=1:ppn=4'
    test_results =parallel_map(
      func=Clusterer(),
      iterable=pdbs,
      method='pbs',
      preserve_exception_message=True,
      processes=len(pdbs),
      qsub_command=qsub_command,
      use_manager=True)

    for test_result in test_results:
      check_assertions(test_results)
  else:
    for pdb_file in os.listdir(pdb_path):
      clusterer = Clusterer(os.path.join(pdb_path, pdb_file),12)
      check_assertions(clusterer.process())

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run()
  print "Total time (all regression): %6.2f"%(time.time()-t0)
