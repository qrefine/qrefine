from __future__ import division
import sys
import time
import os.path
import iotbx.pdb
import libtbx.load_env
from libtbx.easy_mp import parallel_map
from qrefine.core.restraints import from_qm,from_cctbx
from pymongo import MongoClient

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

def check_assertions(result):
  """
  Exercise to test restraints are being calculated consistently. 
  """
  if db.old.find_one({"pdb_code": result.pdb_code}) is None:
      insert(result)
  else:
    result_db = db.old.find_one({"pdb_code": result.pdb_code})
    assert result_db['energy']    == result.chunks
    assert result_db['gradients'] == result.chunk_sizes

def insert(result):
  db.old.insert_one({"pdb_code"  : result.pdb_code,
                     "energy"    : result.energy,
                     "gradients" : result.gradients})

class Result(object):
  def __init__(self, pdb_code, energy, gradients):
    self.pdb_code = pdb_code
    self.energy = energy
    self.gradients = gradients

class Restraints(object):
  def __init__(self,pdb,manager):
    self.pdb = pdb
    self.pdb_inp = iotbx.pdb.input(self.pdb)
    self.ph = self.pdb_inp.construct_hierarchy()
    self.cs = self.pdb_inp.crystal_symmetry()
    self.sites_cart = self.ph.atoms().extract_xyz()
    self.manager = manager

  def process(self):
    energy, gradients = self.manager.target_and_gradients(self.sites_cart)
    return Result(self.pdb_code,energy, gradients)

  def __call__(self):
    return self.process()

def calculators(pdb_file):
   return[
       "moapc",
       "pyscf",
       "terachem",
       "turbomole"
   ]

def run(pdb_path, manager, parallel):
  pdbs=[]
  for pdb_file in os.listdir(pdb_path):
    pdbs.append(pdb_file)
  if(parallel):
    test_results = parallel_map(
          func=Restraints(),
          iterable=pdbs,
          qsub_command=qsub_command,
          processes=len(pdbs),
          method='pbs',
          preserve_exception_message=True,
          use_manager=True)
    for test_result in test_results:
      check_assertions(test_result)
  else:
   for calculator in calculators(pdb_file):
     for pdb_file in pdbs:
      restraint = Restraints(os.path.join(pdb_path, pdb_file), calculator)
      check_assertions(restraint.process())

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  if(1):
    run(args[0],parallel=False)
  if (0):
    run(args[0],parallel=True)
  print "Total time (all regression): %6.2f"%(time.time()-t0)
