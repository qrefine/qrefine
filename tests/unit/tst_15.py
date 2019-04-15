from __future__ import division
import os
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import easy_pickle
import time
import run_tests
import libtbx.load_env
from libtbx.test_utils import approx_equal

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Exercise gradients match:
    - small vs large box:
      -- using clustering vs not using clustering.
  """
  for restraints in ["cctbx","qm"]:
  #for restraints in ["cctbx",]:
  # XXX qm option is not supposed to work fulfull the test with 2ona_box_S
  # XXX qm option is currently suspected to be broken for 2ona_box_L
    for data_file_prefix in ["2ona_box_L", "2ona_box_S"]:
      if(restraints is "qm" and data_file_prefix is "2ona_box_S"): continue
      common_args = ["restraints=%s"%restraints,
                     "mode=opt",
                     "parallel.nproc=1",
                     "quantum.engine_name=mopac",
                     "two_buffers=true"]
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
      if(restraints is "cctbx"):
        ##tight comparison
        ## clustering cctbx should match the atomic gradients
        ## at x,y,z directions
        g1 = g1.as_double()
        g2 = g2.as_double()
        assert g1.size() == g2.size()
        diff = g1-g2
        if(0):
          for i, diff_i in enumerate(diff):
            print i+1, diff_i, g1[i], g2[i]
          print
        assert approx_equal(g1, g2, 1.0E-4)
      else:
        ## loose comparison
        ## clustering qm just checks the norm of gradients from
        ## x, y, z directions
        assert approx_equal(g1.norm(), g2.norm(), g1.norm()*0.05)
        g1_norms = flex.sqrt(g1.dot())
        g2_norms = flex.sqrt(g2.dot())
        for i in range(g1_norms.size()):
          #print i+1, g1_norms[i], g2_norms[i], g1_norms[i]*0.2
          assert approx_equal(g1_norms[i], g2_norms[i], g1_norms[i]*0.2)

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
