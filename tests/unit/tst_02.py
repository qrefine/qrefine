from __future__ import division

import time
import run_tests

def run(prefix = "tst_02"):
  """
  Assert if use_convergence_test=False refinement runs number_of_macro_cycles.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  run_tests.run_cmd(prefix,args = ["restraints=cctbx",
                                   "number_of_macro_cycles=5",
                                   "number_of_micro_cycles=10",
                                   "use_convergence_test=False"])
  ofo = open("%s.log"%prefix,"r")
  found = False
  not_found = False
  for l in ofo.readlines():
    if(l.strip().startswith("8 Rw:")):
      found = True
    if(l.strip().startswith("12 Rw:")):
      not_found = True
  assert found
  assert not not_found
  run_tests.clean_up(prefix)

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_02"
  run(prefix)
  print prefix +":  OK  " + "Time: %6.2f (s)"%(time.time()-t0)

