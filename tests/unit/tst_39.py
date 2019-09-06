from __future__ import division

import os
import run_tests
from mmtbx.pair_interaction import tst_01

def run(prefix):
  """
  Exercise standard test from pair_interactions.
  """
  tst_01.run()

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
