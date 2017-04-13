from __future__ import division

import os
import shutil
import time

import libtbx.env
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry

import expected
from command_line import finalise

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

def run(prefix = "tst_reg_00"):
  """
  Exercise structure completion by finalise.py
  """
  pdb_dir_p1 =os.path.join(qrefine_path,"tests/regression/data/p1/")
  pdb_dir_cluster=os.path.join(qrefine_path,"tests/regression/data/cluster/")
  complete_pdbs(expected.cluster, pdb_dir_cluster)
  complete_pdbs(expected.p1, pdb_dir_p1)

def complete_pdbs(expected_list, pdb_dir):
  if(not (os.path.isdir(pdb_dir))):
    raise Sorry(pdb_dir + " not exist.  Please get its repository on GitHub")
  input_var = str(
    raw_input("run finalise.py for all pdbs in %s will take quite a while (30~60 minutes), continue Y/N : " % pdb_dir))
  if(input_var == "Y"):
    shutil.rmtree("./tmp",ignore_errors=True)
    os.makedirs("./tmp")
    no_error_list = finalise.run(pdb_dir, nproc=4, only_code='1il5', )
    assert 0
    shutil.rmtree("./tmp")
    expected_list.sort()
    no_error_list.sort()
    assert approx_equal(expected_list, no_error_list),'%s has different pdbs passing finalise.py'%pdb_dir
  else:
    print "skip"

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
