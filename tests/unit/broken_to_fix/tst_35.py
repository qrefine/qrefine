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


def run(prefix):
  """
  test principle execution of mode=gtest. Check if files are present.
  ToDo: check other modes.

  """
  pdb = os.path.join(qr_unit_tests_data,'helix.pdb')
  cmd = f"qr.gtest {pdb} restraints=cctbx g_scan=4 g_mode=1"
  assert easy_run.call(cmd)==0
  files = ["1-4.npy",'1-4/0_cluster.pdb','1-4/3_cluster.pdb','1-4/0_frag.pdb','1-4/3_frag.pdb','1-4/fragment_info.txt']
  for f in files:
    assert os.path.isfile(f), f


if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
