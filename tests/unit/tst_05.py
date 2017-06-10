from __future__ import division

import time
import run_tests

def run(prefix = "tst_05"):
  """
  Exercise max_atoms.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx","max_atoms=3"])
  assert r.stdout_lines == ['Sorry: Too many atoms.']
  run_tests.clean_up(prefix)

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_05"
  run(prefix)
  print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
