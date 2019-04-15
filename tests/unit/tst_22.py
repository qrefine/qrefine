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
  pdb_name = os.path.join(qr_unit_tests_data, "tst_22.pdb")
  cmd = "qr.charges %s verbose=False"%pdb_name
  if(0): print cmd
  r = easy_run.go(cmd)
  # Make sure no
  assert len(r.stderr_lines)==0, r.stderr_lines
  assert len(r.stdout_lines)==0, r.stdout_lines

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
