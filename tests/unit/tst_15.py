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
from qrefine.command_line import granalyse

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Exercise gradients match:
    - small vs large box:
      -- using clustering vs not using clustering.
  QM only.
  """
  for restraints in ["qm"]:
  # XXX qm option is not supposed to work fulfull the test with 2ona_box_S
  # XXX qm option is currently suspected to be broken for 2ona_box_L
    for data_file_prefix in ["2ona_box_L", "2ona_box_S"]:
      if(restraints == "qm" and data_file_prefix == "2ona_box_S"): continue
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
      if 1:
        g1 = g1.as_double()
        g2 = g2.as_double()
        assert g1.size() == g2.size()
        diff = g1-g2
        if(1):
          rs = flex.double()
          for a, b in zip(g1.as_double(), g2.as_double()):
            r = abs(abs(a)-abs(b))/(abs(a)+abs(b))*2.*100
            rs.append(r)
        print("min/max/mean:", rs.min_max_mean().as_tuple())

        d = max(granalyse.get_grad_wdelta(ref=g1, g=g2))
        print("get_grad_wdelta:", d)
        assert d < 6.

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
