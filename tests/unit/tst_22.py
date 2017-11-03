from __future__ import division
import os
import time
import libtbx.load_env
import run_tests
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Make sure 'qr.charges tst_22.pdb' runs without errors (finishes successfully).
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  pdb_name = os.path.join(qr_unit_tests_data, "tst_22.pdb")
  cmd = "qr.charges %s"%pdb_name
  if(0): print cmd
  easy_run.call(cmd)

if __name__ == '__main__':
  t0 = time.time()
  prefix = "tst_22"
  if(1):
    run(prefix)
    print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
  else:
    print prefix + ":  Skipped    "
  run_tests.clean_up(prefix)
