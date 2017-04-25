from __future__ import division
import sys
import time
import os.path
import iotbx.pdb
import libtbx.load_env
from libtbx.easy_mp import parallel_map
from qrefine.core.restraints import from_qm,from_cctbx
import restraint_wrapper
from pymongo import MongoClient

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
utils_path= os.path.join(qr_path,"utils")
qr_reg_data = os.path.join(qrefine_path, "regression/datasets/cluster")

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

def calculators(pdb_file):
   return[
       "moapc",
       "pyscf",
       "terachem",
       "turbomole"
   ]

def run():
  pdbs=[]
  for pdb_file in os.listdir(qr_reg_data):
    pdbs.append(os.path.join(qr_reg_data,pdb_file))
  if(1):
    print "running parallel test"
    test_results = parallel_map(
          func=restraint_wrapper(manager),
          iterable=pdbs,
          qsub_command=qsub_command,
          processes=len(pdbs),
          method='pbs',
          preserve_exception_message=True,
          use_manager=True)
    for test_result in test_results:
      check_assertions(test_result)
  else:
   print "running serial"
   for calculator in calculators(pdb_file):
     for pdb_file in pdbs:
      restraint = Restraints(os.path.join(qr_reg_data, pdb_file), calculator)
      check_assertions(restraint.process())

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run()
  print "Total time (all regression): %6.2f"%(time.time()-t0)
