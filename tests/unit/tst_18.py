from __future__ import division
from __future__ import absolute_import
import iotbx.pdb
import os
from scitbx.array_family import flex
from libtbx import easy_pickle
import time
from qrefine.tests.unit import run_tests
from libtbx.test_utils import approx_equal
import libtbx.load_env

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Exercise gradients match:
  -- pdbs with altlocs
      -- using clustering with less clusters vs not using clustering.
      -- using clustering with more clusters vs not using clustering.
  """
  import multiprocessing
  nproc = str(multiprocessing.cpu_count())
  for data_file_prefix in [ "h_altconf_complete", "h_altconf_2_complete"]:
    pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix)
    mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix)
    for maxnum in ["15", "2"]:
      common_args = [
        "restraints=cctbx",
        "mode=opt",
        "altloc_method=subtract",
        "maxnum_residues_in_cluster="+maxnum]
      r = run_tests.run_cmd(prefix,
        args = common_args+["clustering=true", "dump_gradients=cluster_true.pkl"],
        pdb_name = pdb_name,
        mtz_name = mtz_name)
      r = run_tests.run_cmd(prefix,
        args = common_args+["clustering=false", "dump_gradients=cluster_false.pkl"],
        pdb_name = pdb_name,
        mtz_name = mtz_name)
      #
      g1 = easy_pickle.load("cluster_false.pkl").as_double()
      g2 = easy_pickle.load("cluster_true.pkl").as_double()
      assert g1.size() == g2.size()
      diff = flex.abs(g1-g2)
      assert approx_equal(flex.max(diff), 0, 1.e-6)

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=True)
