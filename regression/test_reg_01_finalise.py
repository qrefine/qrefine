from __future__ import division
import os
import shutil
import time
import libtbx.load_env
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from command_line import finalise
from regression.results import expected

from pymongo import MongoClient

db = MongoClient('localhost', 27017).pyoink

qrefine_path = libtbx.env.find_in_repositories("qrefine")
pdb_dir_p1 = os.path.join(qrefine_path, "regression/datasets/p1/")
pdb_dir_cluster = os.path.join(qrefine_path, "regression/datasets/cluster/")

def check_assertions(result):
  """
  Exercise structure completion by finalise.py
  """
  result_db = db.old.find_one({"pdb_code": result.pdb_code})
  assert result_db['charge']     == result.charge
  assert result_db['super_cell'] == result.super_cell
  assert result_db['completion'] == result.completion
  assert result_db['finalise']   == result.capping
  assert result_db['failing']    == result.failing

def run(prefix = "tst_reg_01_finalise"):
  complete_pdbs(expected.cluster, pdb_dir_cluster)
  complete_pdbs(expected.p1, pdb_dir_p1)

def complete_pdbs(expected_list, pdb_dir):
  shutil.rmtree("./tmp",ignore_errors=True)
  os.makedirs("./tmp")
  no_error_list = finalise.run(pdb_dir, nproc=4, only_code='1il5', )
  assert 0
  shutil.rmtree("./tmp")
  assert approx_equal(expected_list.sort(), no_error_list.sort()), \
         '%s has different pdbs passing finalise.py'%pdb_dir

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
