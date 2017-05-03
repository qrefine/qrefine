from __future__ import division
import os.path
import libtbx.load_env
from test_reg_00_base import test_base

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")
qr_reg_data = os.path.join(qrefine_path, "regression/datasets/cluster")

class test_refine(test_base):

  def check_assertions(self,result):
    """
    Exercise to refinements are carried out successfully. 
    """
    if self.db.old.find_one({"pdb_code": result.pdb_code}) is None:
      self.insert(result)
    else:
      past_result = self.db.old.find_one({"pdb_code": result.pdb_code})
      assert past_result['pdb_code']    == result.pdb_code
      assert past_result['rsmd_inital'] == result.rsmd_inital
      assert past_result['rsmd_final']  == result.rsmd_final
      assert past_result['r_start']     == result.r_start
      assert past_result['r_work']      == result.r_work

  def insert(self,result):
    self.db.old.insert_one({"pdb_code"  : result.pdb_code,
                     "rsmd_inital" : result.rsmd_inital,
                     "rsmd_final"  : result.rsmd_final,
                     "r_start"     : result.r_start,
                     "r_work"      : result.r_work})


