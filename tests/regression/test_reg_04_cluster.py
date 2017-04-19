from __future__ import division
import sys
import time
import os.path
import iotbx.pdb
import libtbx
from libtbx.easy_mp import parallel_map
from pymongo import MongoClient
from qrefine.core import completion
from qrefine.core.plugin.yoink import pyoink
from qrefine.core.restraints import from_qm
try:
  from jpype import startJVM
except ImportError, e:
  raise Sorry(str(e))

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
utils_path= os.path.join(qr_path,"utils") 
pdb_path= os.path.join(qrefine_path,"tests/regression/data/p1")

class Result(object):
    def __init__(self, pdb_code, clusters, chunks, chunk_sizes):
        self.pdb_code = pdb_code
        self.clusters = clusters
        self.chunks = chunks
        self.chunk_sizes = chunk_sizes
        self.num_of_clusters = len(clusters)
        self.num_of_chunks = len(chunks)

class Clusterer(object):
   def __init__(self,pdb,maxnum_residues_in_cluster):
      self.pdb = pdb
      self.pdb_inp = iotbx.pdb.input(self.pdb)
      self.ph = self.pdb_inp.construct_hierarchy()
      self.cs = self.pdb_inp.crystal_symmetry()
      self.sites_cart = self.ph.atoms().extract_xyz()
      self.fq = from_qm(
             use_cluster_qm=True,
             pdb_hierarchy       = self.ph,
             crystal_symmetry    = self.cs,
             maxnum_residues_in_cluster = int(maxnum_residues_in_cluster),
             yoink_jar_path    = qr_path + "/plugin/yoink/Yoink-0.0.1.jar",
             yoink_dat_path    = qr_path + "/plugin/yoink/dat")

   def process(self):
       chunks = []
       chunk_sizes = []
       for chunk in self.fq.fragments.qm_pdb_hierarchies:
           res_in_chunk = []
           atom_tot_per_residue = 0
           for chain in chunk.only_model().chains():
               for residue_group in chain.residue_groups():
                   res_in_chunk.append(residue_group.resid())
                   for atom_group in residue_group.atom_groups():
                       atom_tot_per_residue += atom_group.atoms_size()
                       chunk_sizes.append(atom_tot_per_residue)
           chunks.append(res_in_chunk)
       return Result(self.pdb_file,
                     self.fq.fragments.clusters,
                     chunks,
                     chunk_sizes)

   def __call__(self):
       return self.process()

class Pbs_Cluster():
    def __init__(self):
      self.qsub = None

    def qsub(self):
        qsub_file = os.join(utils_path, "/qsub.pbs")
        qsub_command = """ -N %s -v arguments=" %s  "  """
        qsub_command = qsub_command + "  > /dev/null"
        os.system(qsub_command)


def check_all(results):
  for result in results():
    if db.old.find_one({"pdb_code": result.pdb_code}) is None:
        insert(result)
    else:
        check_assertions(result)
    check_assertions(result)

def check_assertions(result):
    """
      Exercise to test clustering indices, chunk indices, and the number of atoms in each chunk. 
    """
    result_db = db.old.find_one({"pdb_code": result.pdb_code})
    assert result_db['clusters']    == result.clusters
    assert result_db['chunks']      == result.chunks
    assert result_db['chunk_sizes'] == result.chunk_sizes

def insert(result):
    db.old.insert_one({"pdb_code"        : result.pdb_code,
                       "clusters"        : result.clusters,
                       "chunks"          : result.chunks,
                       "chunk_sizes"     : result.chunk_sizes,
                       "num_of_chunks"   : result.num_of_chunks,
                       "num_of_clusters" : result.num_of_clusters,
                       })

def serial_run(pdbs):
  maxnum_residues_in_cluster = 15
  test_results = []
  for pdb_file in pdbs:
    clusterer = Clusterer(pdb_file, maxnum_residues_in_cluster)
    test_results.append(clusterer.process())
  check_all(test_results)

def parallel_run(pdbs):      
  """ first attempt to run in parallel on cluster"""
  q = Pbs_Cluster()
  test_results = parallel_map(
        func=Clusterer(),
        iterable=pdbs,
        qsub_command= q.qsub_command,
        processes=len(pdbs),
        method='pbs',
        preserve_exception_message=True,
        use_manager=True)
  check_all(test_results)
  
def pdbs():
  pdbs=[]
  for file in os.listdir(pdb_path): 
     print "process:", (os.path.join(pdb_path,file)) 
     pdbs.append((os.path.join(pdb_path,file)))
  return pdbs 
  
if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  if(1):
    parallel_run(pdbs())
  if(1):
    serial_run(pdbs())
  print "Total time (all tests): %6.2f"%(time.time()-t0)
