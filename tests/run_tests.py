from __future__ import division

import time

from libtbx import easy_run

def run_regression_tests():
  regression_tests = [
    "test_reg_00.py", #nigel
    "test_reg_01.py", #nigel
    "test_reg_02.py"
  ]
  for file_name in regression_tests:
    print "Running regression test:", file_name
    easy_run.call("cctbx.python %s"%file_name)


if(__name__ == "__main__"):
  t0 = time.time()
  run_regression_tests()
  print "Total time (all tests): %6.2f"%(time.time()-t0)
  print "OK"
