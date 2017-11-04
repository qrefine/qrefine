from __future__ import division

import time
import run_tests

def run(prefix):
  """
  Exercise max_atoms.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx","max_atoms=3"])
  assert r.stdout_lines == ['Sorry: Too many atoms.']

if(__name__ == "__main__"):
  run_tests.runner(function=run, prefix="tst_05", disable=False)
