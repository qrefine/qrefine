from __future__ import division

import os
import time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
import run_tests
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

def run(prefix):
  """
  Exercise standard (cctbx-based restraints) optimization (no data required).
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  run_tests.run_cmd(
    prefix   = prefix, 
    args     = ["restraints=cctbx mode=opt use_convergence_test=False"], 
    mtz_name = "")
  # Check results
  pdb_inp = iotbx.pdb.input(
    file_name = os.path.join(prefix,"m00_poor_refined.pdb"))
  model_1 = mmtbx.model.manager(model_input = pdb_inp, build_grm=True, 
    log=null_out())
  s1 = model_1.geometry_statistics().result()
  assert s1.bond.mean < 0.005
  ##########################
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  run_tests.run_cmd(
    prefix   = prefix, 
    args     = ["restraints=cctbx mode=opt use_convergence_test=False"])
  # Check results
  pdb_inp = iotbx.pdb.input(
    file_name = os.path.join(prefix,"m00_poor_refined.pdb"))
  model_2 = mmtbx.model.manager(model_input = pdb_inp, build_grm=True, 
    log=null_out())
  s2 = model_2.geometry_statistics().result()
  assert s2.bond.mean < 0.005
  ###########################
  assert approx_equal(s1.bond.mean, s2.bond.mean, 0.0001)

if(__name__ == "__main__"):
  prefix="tst_30"
  rc = run_tests.runner(function=run, prefix=prefix, disable=False)
  assert not rc, '%s rc: %s' % (prefix, rc)
