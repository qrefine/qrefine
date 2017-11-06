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

def run(prefix = "tst_18"):
  """
  Exercise gradients match:
  -- pdbs with altlocs
      -- using clustering with less clusters vs not using clustering.
      -- using clustering with more clusters vs not using clustering.
  """
  import multiprocessing
  nproc = str(multiprocessing.cpu_count())
  for data_file_prefix in [ "h_altconf_complete", "h_altconf_2_complete"]:
    for maxnum in ["15", "2"]:
      common_args = ["restraints=cctbx", "mode=opt", "nproc="+nproc] +\
                ["altloc_method=subtract","maxnum_residues_in_cluster="+maxnum]
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
      assert approx_equal(diff.max(), [0,0,0], [1.0E-3,1.0E-3,1.0E-3])

if __name__ == '__main__':
  run_tests.runner(function=run, prefix="tst_18", disable=False)
