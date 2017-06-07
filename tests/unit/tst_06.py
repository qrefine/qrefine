from __future__ import division

import time
import run_tests
import iotbx.pdb

def run(prefix = "tst_06"):
  """
  Exercise stop if altlocs.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example(
                                          pdb_name = "data_files/altlocs.pdb",
                                          mtz_name = "data_files/altlocs.mtz")
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx"],
                    pdb_name = "data_files/altlocs.pdb",
                    mtz_name = "data_files/altlocs.mtz")
  assert r.stdout_lines == \
    ['Sorry: Alternative conformations are not supported.']
  run_tests.clean_up(prefix,mtz_name = "data_files/altlocs.mtz")

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
