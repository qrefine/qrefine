from __future__ import division

import os
import time
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
import run_tests

def run(prefix):
  """
  Exercise standard (cctbx-based restraints) refinement with all defaults.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  run_tests.run_cmd(prefix, args = ["restraints=cctbx"])
  # Check results
  xrs_refined = iotbx.pdb.input(
    file_name = os.path.join(prefix,"m00_poor_refined.pdb")).xray_structure_simple()
  d = flex.sqrt((xrs_good.sites_cart() - xrs_poor.sites_cart()).dot())
  assert flex.mean(d) > 0.10
  d = flex.sqrt((xrs_good.sites_cart() - xrs_refined.sites_cart()).dot())
  assert flex.mean(d) < 0.02, 'Mean difference %0.3f is greater than 0.02' % (
    flex.mean(d))
  # Check R-factors
  r_start, r_final = None,None
  ofo = open("%s.log"%prefix,"r")
  for l in ofo.readlines():
    if(l.strip().startswith("0 Rw:")):
      r_start = float(l.split()[2])
    if(l.strip().startswith("Best r_work:")):
      r_final = float(l.split()[2])
  assert r_start > 0.1, r_start
  assert r_final < 0.006
  # Make sure output model actually corresponds to reported R-factor
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs_refined)
  fmodel.update_all_scales()
  assert fmodel.r_work() < 0.006
  run_tests.clean_up(prefix)

if(__name__ == "__main__"):
  prefix = "tst_00"
  t0 = time.time()
  run(prefix)
  print prefix +":  OK  " + "Time: %6.2f (s)"%(time.time()-t0)
