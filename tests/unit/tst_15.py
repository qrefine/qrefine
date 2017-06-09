from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import easy_pickle
import time
import run_tests
from libtbx.test_utils import approx_equal

def run(prefix = "tst_15"):
  """
  Exercise gradients match:
    - small vs large box:
      -- using clustering vs not using clustering.
  """
  for data_file_prefix in ["2ona_box_L", "2ona_box_S"]:
    common_args = ["restraints=cctbx", "mode=opt", "nproc=1"]
    r = run_tests.run_cmd(prefix,
      args     = common_args+["clustering=true",
                              "dump_gradients=cluster_true.pkl"],
      pdb_name = "data_files/%s.pdb"%data_file_prefix,
      mtz_name = "data_files/%s.mtz"%data_file_prefix)
    r = run_tests.run_cmd(prefix,
      args     = common_args+["clustering=false", 
                             "dump_gradients=cluster_false.pkl"],
      pdb_name = "data_files/%s.pdb"%data_file_prefix,
      mtz_name = "data_files/%s.mtz"%data_file_prefix)
    #
    g1 = flex.vec3_double(easy_pickle.load("cluster_false.pkl"))
    g2 = flex.vec3_double(easy_pickle.load("cluster_true.pkl"))
    assert g1.size() == g2.size()
    diff = g1-g2
    #for i, diff_i in enumerate(diff):
    #  print i, diff_i#, g1[i], g2[i]
    #print 
    assert approx_equal(diff.max(), [0,0,0])

if __name__ == '__main__':
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"


