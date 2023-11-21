from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import time
import os
from qrefine.tests.unit import run_tests
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")
qr_unit_tests_data = os.path.join(qr_unit_tests,"data_files")


f_info = """~  # clusters  :  4
~  list of atoms per cluster:
~    [29, 23, 21, 21]
~  list of atoms per fragment:
~    [29, 23, 21, 21]"""

def run(prefix):
  """
  test execution of mode=gtest

  XXX TEST FAILS: some args are not recognized (g_scan, g_mode). No result checks.
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  pdb = os.path.join(qr_unit_tests_data,'helix.pdb')
  cmd = f"qr.gtest {pdb} restraints=cctbx g_scan=3 g_mode=1"
  assert easy_run.call(cmd)==0
  files = ["1-3.npy",'1-3/0_cluster.pdb','1-3/0_frag.pdb','1-3/fragment_info.txt']
  for f in files:
    assert os.path.isfile(f), f
  with open('1-3/fragment_info.txt','r') as f:
    lines = f.read().splitlines() 
    ref = f_info.splitlines()
    for i,l in enumerate(lines):
      print(l,ref[i])
      assert l == ref[i]


if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
