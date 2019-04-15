from __future__ import division

import time, os
import run_tests
import iotbx.pdb

def run(prefix):
  """
  Exercise altlocs.
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example(
    pdb_name = "altlocs.pdb",
    mtz_name = "altlocs.mtz")
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx"],
                    pdb_name = "altlocs.pdb",
                    mtz_name = "altlocs.mtz")
  run_tests.clean_up(prefix,mtz_name = "altlocs.mtz")

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = os.path.basename(__file__).replace(".py","")
  run(prefix)
  run_tests.clean_up(prefix)
