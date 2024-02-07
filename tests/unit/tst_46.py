from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import time
import libtbx.load_env
from qrefine.tests.unit import run_tests
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Make sure real-space-refinement works (runs).
  """
  pdb_name = os.path.join(qr_unit_tests_data, "tst_46.pdb")
  mtz_name = os.path.join(qr_unit_tests_data, "tst_46.map")
  cmd = "qr.refine %s %s restraints=cctbx clustering=false"%(
    pdb_name, mtz_name)
  if(0): print(cmd)
  r = easy_run.go(cmd)
  # Make sure no
  assert len(r.stderr_lines)==0, r.stderr_lines

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
