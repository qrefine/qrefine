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

def run(prefix):
  """
  Exercise refinement  match:
  -- pdbs with altlocs
      -- using subtract vs using average
  """
  import multiprocessing
  nproc = str(multiprocessing.cpu_count())
  for data_file_prefix in ["h_altconf_complete", "h_altconf_2_complete"]:
    for maxnum in ["15"]:
      common_args = ["stpmax=0.5", "restraints=cctbx", "mode=refine", "parallel.nproc="+nproc, "clustering=true"] +\
        ["gradient_only=true", "maxnum_residues_in_cluster="+maxnum]
      r = run_tests.run_cmd(prefix,
        args     = common_args+["altloc_method=subtract", "output_file_name_prefix=subtract-"+data_file_prefix],
        pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix),
        mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix))
      r = run_tests.run_cmd(prefix,
        args     = common_args+["altloc_method=average", "output_file_name_prefix=average-"+data_file_prefix],
        pdb_name = os.path.join(qr_unit_tests_data,"%s.pdb"%data_file_prefix),
        mtz_name = os.path.join(qr_unit_tests_data,"%s.mtz"%data_file_prefix))
      sites_cart_average = iotbx.pdb.input("%s/average-%s_refined.pdb"%(prefix,data_file_prefix)).\
                       construct_hierarchy().extract_xray_structure().sites_cart()
      sites_cart_subtract = iotbx.pdb.input("%s/subtract-%s_refined.pdb"%(prefix,data_file_prefix)).\
                       construct_hierarchy().extract_xray_structure().sites_cart()
      rmsd_diff = sites_cart_average.rms_difference(sites_cart_subtract)
      assert approx_equal(rmsd_diff, 0, 0.2), 'rmsd diff between subtract and average is too large %0.3f' % rmsd_diff

if __name__ == '__main__':
  rc = run_tests.runner(function=run, prefix="tst_19", disable=False)
  assert not rc, 'tst_19 rc: %s' % rc
