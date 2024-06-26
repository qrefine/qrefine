from __future__ import division
from __future__ import absolute_import

import time, os
from qrefine.tests.unit import run_tests

def run(prefix):
  """
  Exercise max_atoms.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  r = run_tests.run_cmd(prefix,
                    args = ["restraints=cctbx","max_atoms=3"])
  found = False
  with open("%s.log"%prefix,"r") as fo:
    for l in fo.readlines():
      if "Too many atoms. Can take forever or crash." in l:
        found = True
        break
  assert found

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
