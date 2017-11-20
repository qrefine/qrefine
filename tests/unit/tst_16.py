from __future__ import division
import iotbx.pdb
import os
from scitbx.array_family import flex
from libtbx import easy_pickle
import time
import run_tests
from libtbx.test_utils import approx_equal
import libtbx.load_env

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix = "tst_16"):
  """
  Exercise gradients match:
    - small vs large box:
      -- using clustering vs not using clustering.
  Non-P1 case (P212121)
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  data_file_prefix = "2olx"
  common_args = ["restraints=cctbx", "mode=opt", "nproc=1"]
  r = run_tests.run_cmd(prefix,
    args     = common_args+["clustering=true",
                            "dump_gradients=cluster_true.pkl"],
    pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix),
    mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix))
  r = run_tests.run_cmd(prefix,
    args     = common_args+["clustering=false",
                           "dump_gradients=cluster_false.pkl"],
    pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix),
    mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix))
  #
  g1 = flex.vec3_double(easy_pickle.load("cluster_false.pkl"))
  g2 = flex.vec3_double(easy_pickle.load("cluster_true.pkl"))
  assert g1.size() == g2.size()
  diff = g1-g2
  if(0):
    for i, diff_i in enumerate(diff):
      print i, diff_i#, g1[i], g2[i]
    print
  assert approx_equal(diff.max(), [0,0,0])

if __name__ == '__main__':
  rc = run_tests.runner(function=run, prefix="tst_16", disable=False)
  assert not rc, 'tst_00 rc: %s' % rc
