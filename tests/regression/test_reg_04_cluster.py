from __future__ import division
import sys
import time
import os.path
import iotbx.pdb
from libtbx.easy_mp import parallel_map
from pymongo import MongoClient
from qrefine.core import completion
from qrefine.core.plugin import pyoink
from qrefine.core.restraints import from_qm
try:
  from jpype import startJVM
except ImportError, e:
  raise Sorry(str(e))

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
pdb_path= os.path.join(qrefine_path,"tests/regression/data/p1") 
utils_path= os.path.join(qr_path,"utils") 

class Cluster_Pdb(object):
   def __init__(self,pdb,fq):
      self.pdb = pdb
      self.fq = fq
      
   def __call__(self,pdb):
    return process(fq,self.pdb):

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
  Exercise to test clustering indices, chunk indices, and the number of atoms in each chunk. 
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
    yoink_jar_path    = qr_path + "/plugin/yoink/Yoink-0.0.1.jar",
    yoink_dat_path    = qr_path + "/plugin/yoink/dat")

  result = process(fq,pdb_file)
  if db.old.find_one({"pdb_code": result.pdb_code}) is None:
     insert(result)
  else:
     check_assertions(result)

def qsub():
  qsub_file = os.join(utils_path,"/qsub.pbs")
  qsub_command = """ -N %s -v arguments=" %s  "  """
  qsub_command = qsub_command + "  > /dev/null"
  os.system(qsub_command)      
      
def parallel_run(pdbs):      
  """ first attempt to run in parallel on cluster"""
  test_results = parallel_map(
        func=Cluster_Pdb(),
        iterable=pdbs,
        qsub_command= self.qsub_command,
        processes=len(pdbs),
        method='pbs',
        preserve_exception_message=True,
        use_manager=True)        

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  pdbs=[]
  for file in os.listdir(pdb_path): 
     print "process:", (os.path.join(pdb_path,file)) 
     pdbs.append((os.path.join(pdb_path,file))) 
  parallel_run(pdbs)  
  print "Total time (all tests): %6.2f"%(time.time()-t0)

