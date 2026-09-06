from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import libtbx.load_env
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Make sure 'qr.charges tst_22.pdb' runs without errors (finishes successfully).
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  pdb_name = os.path.join(qr_unit_tests_data, "tst_22.pdb")
  cmd = "qr.charges %s verbose=False"%pdb_name
  r = easy_run.go(cmd)
  # Make sure no
  assert len(r.stderr_lines)==0, r.stderr_lines
  assert len(r.stdout_lines)==0, r.stdout_lines

if(__name__ == '__main__'):
  run()
