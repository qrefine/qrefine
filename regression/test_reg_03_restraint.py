from __future__ import division
import sys
import time
import os.path
import iotbx.pdb
import libtbx.load_env
from qrefine.core.restraints import from_qm,from_cctbx
from restraint_wrapper import Restraints, Result
from test_reg_00_base import test_base


qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
utils_path= os.path.join(qr_path,"utils")
qr_reg_data = os.path.join(qrefine_path, "regression/datasets/cluster")

class test_restraint(test_base):
  def __init__(self):
    self.pdbs = self.pdbs(qr_reg_data)
    self.maxnum_residues_in_cluster = 12
    self.func = Restraints()

  def check_assertions(result):
    """
    Exercise to test restraints are being calculated consistently. 
    """
    if self.db.old.find_one({"pdb_code": result.pdb_code}) is None:
        insert(result)
    else:
      result_db = self.db.old.find_one({"pdb_code": result.pdb_code})
      assert result_db['energy']    == result.chunks
      assert result_db['gradients'] == result.chunk_sizes

  def insert(result):
    self.db.old.insert_one({"pdb_code"  : result.pdb_code,
                     "energy"    : result.energy,
                     "gradients" : result.gradients})


  def calculators(pdb_file):
     return[
         "moapc",
         "pyscf",
         "terachem",
         "turbomole"
     ]

