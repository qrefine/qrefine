from __future__ import division
import os
import sys
import time
import os.path
import pymongo
import iotbx.pdb
from qrefine.core.plugin import pyoink
from qrefine.core.restraints import from_qm
from pymongo import MongoClient
from qrefine.core import completion
try:
  from jpype import startJVM
except ImportError, e:
  raise Sorry(str(e))

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

class Result(object):
    def __init__(self,pdb_code,clusters,chunks,chunk_sizes):
      self.pdb_code        = pdb_code
      self.clusters        = clusters
      self.chunks          = chunks
      self.chunk_sizes     = chunk_sizes
      self.num_of_clusters = len(clusters)
      self.num_of_chunks   = len(chunks)

def check_assertions(cr):
    result_db = db.old.find_one({"pdb_code": cr.pdb_code})
    assert result_db['clusters']    == cr.clusters
    assert result_db['chunks']      == cr.chunks
    assert result_db['chunk_sizes'] == cr.chunk_sizes

def insert(result):
    db.old.insert_one({"pdb_code": result.pdb_code,
                       "clusters": result.clusters,
                       "chunks": result.chunks,
                       "num_of_clusters":result.num_of_clusters,
                       "num_of_chunks":result.num_of_chunks,
                       "chunk_sizes": result.chunk_sizes})

def process(fq,pdb_file):
  chunks=[]
  chunk_sizes=[]
  for chunk in fq.fragments.qm_pdb_hierarchies:
    res_in_chunk=[]
    atom_tot_per_residue = 0
    for chain in chunk.only_model().chains():
      for residue_group in chain.residue_groups():
        res_in_chunk.append(residue_group.resid())
        for atom_group in residue_group.atom_groups():
          atom_tot_per_residue += atom_group.atoms_size()
          chunk_sizes.append(atom_tot_per_residue)
    chunks.append(res_in_chunk)
  return Result(pdb_file,
                fq.fragments.clusters,
                chunks,
                chunk_sizes)

def run(pdb_file, maxnum_residues_in_cluster=15):
  """
  Exercise clustering indices, and chunk indices, number of atoms in each chunk. 
  """
#  completion.run(pdb_file)
#  pdb_file = os.path.basename(pdb_file)[:-4] + "_complete.pdb"
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  sites_cart = ph.atoms().extract_xyz()
  fq = from_qm(
    pdb_hierarchy       = ph,
    crystal_symmetry    = cs,
    use_cluster_qm    = True,
    maxnum_residues_in_cluster = int(maxnum_residues_in_cluster),
    yoink_jar_path    = qr_path + "/./core/plugin/yoink/Yoink-0.0.1.jar",
    yoink_dat_path    = qr_path + "/./core/plugin/yoink/dat")

  result = process(fq,pdb_file)
  if db.old.find_one({"pdb_code": result.pdb_code}) is None:
     insert(result)
  else:
     check_assertions(result)

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(file)
  print "Total time (all tests): %6.2f"%(time.time()-t0)

