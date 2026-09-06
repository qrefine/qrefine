from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import easy_pickle
import time
from qrefine.tests.unit import run_tests
import libtbx.load_env
from libtbx.test_utils import approx_equal

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Exercise gradients match:
    - small vs large box:
      -- using clustering vs not using clustering.
  CCTBX only
  """
  for data_file_prefix in ["2ona_box_L", "2ona_box_S"]:
    common_args = ["restraints=cctbx",
                   "mode=opt"]
    args = common_args+["clustering=false","dump_gradients=cluster_false.pkl"]
    r = run_tests.run_cmd(prefix,
      args     = args,
      pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix),
      mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix))

    args = common_args+["clustering=true", "dump_gradients=cluster_true.pkl"]
    r = run_tests.run_cmd(prefix,
      args     = args,
      pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix),
      mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix))

    #
    g1 = easy_pickle.load("cluster_false.pkl")
    g2 = easy_pickle.load("cluster_true.pkl")
    g1 = g1.as_double()
    g2 = g2.as_double()
    assert g1.size() == g2.size()
    diff = g1-g2
    if(0):
      for i, diff_i in enumerate(diff):
        print(i+1, diff_i, g1[i], g2[i])
      print()
    assert approx_equal(flex.max(diff), 0, 1.0E-4)

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
