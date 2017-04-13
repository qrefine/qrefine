from __future__ import division
import sys
import os
import time
import os.path
import iotbx.pdb
#from iotbx.pdb import hierarchy
from qrefine.core.plugin import pyoink
from qrefine.core.restraints import from_qm
import pymongo
from pymongo import MongoClient
from qrefine.core import completion
try:
  from jpype import startJVM
except ImportError, e:
  raise Sorry(str(e))


db = MongoClient('localhost', 27017).pyoink
qr_path = "/Applications/phenix-dev-2733/modules/qrefine"


pdbs = [
"1il5.pdb",
"1va7.pdb",
"1w3b.pdb",
"2bx9.pdb",
"2dwq.pdb",
"2ec2.pdb",
"2gd5.pdb",
"2ghj.pdb",
"2oeq.pdb",
"2oy0.pdb",
"2pro.pdb",
"2x10.pdb",
"2yht.pdb",
"3en7.pdb",
"3i05.pdb",
"3i3u.pdb",
"3kyi.pdb",
"3lrj.pdb",
"3tz7.pdb",
"3tz9.pdb",
"4c0m.pdb",
"4cc1.pdb",
"4ctd.pdb",
"4drw.pdb",
"4i21.pdb",
"4k2r.pdb",
"4lgh.pdb",
"4rnf.pdb",
"4xa1.pdb",
"4xg9.pdb",
"4ynv.pdb",
"4ynw.pdb",
"5arj.pdb",
"5d12.pdb",
"5eta.pdb",
"5fv7.pdb",
"5ghv.pdb",
"5ldj.pdb",
"5lwe.pdb",
"5t68.pdb",
]


class cluster_result(object):
    def __init__(self,pdb_code,clusters,chunks,chunk_sizes):
      self.pdb_code = pdb_code
      self.clusters = clusters
      self.num_of_clusters = len(clusters)
      self.chunks = chunks
      self.num_of_chunks   = len(chunks)
      self.chunk_sizes=chunk_sizes


def check_assertions(cr):
    result_db = db.old.find_one({"pdb_code": cr.pdb_code})
    assert result_db['clusters'] == cr.clusters
    assert result_db['chunks'] == cr.chunks
    assert result_db['chunk_sizes'] == cr.chunk_sizes

def insert(cr):
    db.old.insert_one({"pdb_code": cr.pdb_code,
                       "clusters": cr.clusters,
                       "chunks": cr.chunks,
                       "num_of_clusters":cr.num_of_clusters,
                       "num_of_chunks":cr.num_of_chunks,
                       "chunk_sizes": cr.chunk_sizes})


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
  return cluster_result(pdb_file,
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

  cr = process(fq,pdb_file)
  if db.old.find_one({"pdb_code": cr.pdb_code}) is None:
     insert(cr)
  else:
     check_assertions(cr)

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  for file in pdbs:
     run(file)
  print "Total time (all tests): %6.2f"%(time.time()-t0)

