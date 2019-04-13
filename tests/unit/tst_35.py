from __future__ import division

import time
import os
import run_tests
import iotbx.pdb
import libtbx.load_env


qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")


def run(prefix = "tst_35"):
  """
  test execution of mode=gtest
  """

  # run_tests.assert_folder_is_empty(prefix=prefix)
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  # pdb_inp = os.path.join(qr_unit_tests,"data_files","2lvr.pdb")
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx","mode=gtest","g_scan=20","g_mode=1"],
                    pdb_name = 'm00_poor.pdb', mtz_name='')
  assert os.path.isfile('1-20.npy')

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_35"
  run(prefix)
  print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
  run_tests.clean_up(prefix)