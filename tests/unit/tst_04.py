from __future__ import division

import time
import run_tests

def run(prefix = "tst_04"):
  """
  Exercise stop if not P1.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example(
                                          pdb_name = "p212121.pdb",
                                          mtz_name =  "p212121.mtz")
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx"],
                    pdb_name = "p212121.pdb",
                    mtz_name = "p212121.mtz")
  assert r.stdout_lines == ['Sorry: Only P1 is supported.']

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
