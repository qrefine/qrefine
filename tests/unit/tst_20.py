from __future__ import division
import iotbx.pdb
import os
import time
import libtbx.load_env
from qrefine.super_cell import expand
import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix = "tst_20"):
  """
  Make sure expand works for this particular file.
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  pdb_name = os.path.join(qr_unit_tests_data,
    "1bdw_ala_refine_001_complete_minimized.pdb_modified.pdb")
  pdb_inp = iotbx.pdb.input(pdb_name)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  super_cell = expand(
    pdb_hierarchy=ph,
    crystal_symmetry=cs)

if __name__ == '__main__':
  t0 = time.time()
  prefix = "tst_20"
  if(1):
    run(prefix)
    print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
  else:
    print prefix + ":  Skipped    "
  run_tests.clean_up(prefix)
