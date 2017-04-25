from __future__ import division

import os
import time
import libtbx.load_env
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_reg_tests = os.path.join(qrefine, "regression")
qr_reg_data = os.path.join(qrefine, "regression/datasets")

def run():
  regression_tests = [
    #"test_reg_01_finalise.py",
    "test_reg_02_chunk.py",
    #"test_reg_03_restraint.py",
    #"test_reg_04_refine.py",
  ]

  for regression_test in regression_tests:
        print "Running: {} ".format(regression_test)
        easy_run.call("cctbx.python {} ".format(  
           os.path.join(qr_reg_tests,regression_test)))

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Total time (all regression): %6.2f"%(time.time()-t0)
  print "OK"
