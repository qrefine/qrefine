from __future__ import division
from __future__ import absolute_import

import os
import time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from qrefine.tests.unit import run_tests

def run(prefix):
  """
  Exercise standard (cctbx-based restraints) refinement with all defaults.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  # run_tests.run_cmd(prefix, args = ["restraints=cctbx stpmax=0.000000001"]) # force to fail for testing
  run_tests.run_cmd(prefix, args = ["restraints=cctbx"])
  # Check results
  xrs_refined = iotbx.pdb.input(
    file_name = os.path.join(prefix,"m00_poor_refined.pdb")).xray_structure_simple()
  d = flex.sqrt((xrs_good.sites_cart() - xrs_poor.sites_cart()).dot())
  assert flex.mean(d) > 0.10
  d = flex.sqrt((xrs_good.sites_cart() - xrs_refined.sites_cart()).dot())
  assert flex.mean(d) < 0.05
  # Check R-factors
  r_start, r_final = None,None
  ofo = open("%s.log"%prefix,"r")
  for l in ofo.readlines():
    if(l.strip().endswith("n_fev: 0")):
      r_start = float(l.split()[2])
    if(l.strip().startswith("Best r_work:")):
      r_final = float(l.split()[2])
  assert r_start > 0.1
  assert r_final < 0.04
  # Make sure output model actually corresponds to reported R-factor
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs_refined)
  fmodel.update_all_scales()
  assert fmodel.r_work() < 0.04
  assert abs(r_final-fmodel.r_work())<0.0005, abs(r_final-fmodel.r_work())

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
