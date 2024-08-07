from __future__ import division
from __future__ import absolute_import

import os
import time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from qrefine.tests.unit import run_tests
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

def run(prefix):
  """
  Exercise standard (cctbx-based restraints) optimization (no data required).
  """
  args = ["restraints=cctbx mode=opt use_convergence_test=False",
          "number_of_micro_cycles=3",
          "max_iterations_refine=100"]
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  run_tests.run_cmd(
    prefix   = prefix,
    args     = args,
    mtz_name = "")
  # Check results
  pdb_inp = iotbx.pdb.input(
    file_name = os.path.join(prefix,"m00_poor_refined.pdb"))
  model_1 = mmtbx.model.manager(model_input = pdb_inp, log=null_out())
  model_1.process(make_restraints=True)
  s1 = model_1.geometry_statistics().result()
  assert s1.bond.mean < 0.005

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
