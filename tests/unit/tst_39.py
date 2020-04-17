from __future__ import division
from __future__ import absolute_import

import os
from qrefine.tests.unit import run_tests
from mmtbx.pair_interaction import tst_01

def run(prefix):
  """
  Exercise standard test from pair_interactions.
  """
  tst_02.run()

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
