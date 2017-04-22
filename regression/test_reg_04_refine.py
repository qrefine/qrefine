from __future__ import division
import sys
import time
import os.path
import iotbx.pdb
import libtbx.load_env
from libtbx.easy_mp import parallel_map
from pymongo import MongoClient

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
utils_path = os.path.join(qr_path,"utils")

def check_assertions(result):
  """
  Exercise to refinements are carried out successfully. 
  """
  if db.old.find_one({"pdb_code": result.pdb_code}) is None:
    insert(result)
  else:
    past_result = db.old.find_one({"pdb_code": result.pdb_code})
    assert past_result['pdb_code']    == result.pdb_code
    assert past_result['rsmd_inital'] == result.rsmd_inital
    assert past_result['rsmd_final']  == result.rsmd_final
    assert past_result['r_start']     == result.r_start
    assert past_result['r_work']      == result.r_work

def insert(result):
    db.old.insert_one({"pdb_code"    : result.pdb_code,
                       "rsmd_inital" : result.rsmd_inital,
                       "rsmd_final"  : result.rsmd_final,
                       "r_start"     : result.r_start,
                       "r_work"      : result.r_work})

class Result(object):
    def __init__(self,
                 pdb_code,
                 rsmd_inital,
                 rsmd_final,
                 r_start,
                 r_work):
        self.pdb_code    = pdb_code
        self.rsmd_inital = rsmd_inital
        self.rsmd_final  = rsmd_final
        self.r_start     = r_start
        self.r_work      = r_work

class Driver(object):
    def __init__(self, pdb, driver):
        self.pdb = pdb
        self.pdb_inp = iotbx.pdb.input(self.pdb)
        self.ph = self.pdb_inp.construct_hierarchy()
        self.cs = self.pdb_inp.crystal_symmetry()
        self.sites_cart = self.ph.atoms().extract_xyz()
        self.driver = driver

    def refine(self):
        self.driver(self.sites_cart)
        return Result(self.pdb_code,
                      self.rsmd_inital,
                      self.rsmd_final,
                      self.r_start,
                      self.r_work)

    def __call__(self):
        return self.refine()

def run(pdb_path,parallel=False):
  pdbs=[]
  for pdb_file in os.listdir(pdb_path):
      pdbs.append(os.path.join(pdb_path))
  if(parallel):
    parallel_map(
      func=Driver(),
      iterable=pdbs,
      qsub_command= qsub_command,
      processes=len(pdbs),
      method='pbs',
      preserve_exception_message=True,
      use_manager=True)
  else:
    for pdb_file in pdbs:
      driver = Driver(os.path.join(pdb_path, pdb_file))
      check_assertions(driver.refine())

if(__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  if(1):
      run(args[0])
  if(0):
      run(args[0],parallel=True)
  print "Total time (all regression): %6.2f"%(time.time()-t0)
