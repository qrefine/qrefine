from __future__ import division
from __future__ import absolute_import

import time, os
from qrefine.tests.unit import run_tests
import iotbx.pdb

def run(prefix):
  """
  Exercise altlocs: just makes sure all runs.

  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example(
    pdb_name = "altlocs.pdb",
    mtz_name = "altlocs.mtz")
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx"],
                    pdb_name = "altlocs.pdb",
                    mtz_name = "altlocs.mtz")

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
