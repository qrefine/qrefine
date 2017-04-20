from __future__ import division

import os
import time
import libtbx
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_reg_tests = os.path.join(qrefine, "tests/regression")
qr_reg_data = os.path.join(qrefine, "tests/regression/data")

pdb_path= os.path.join(qr_reg_data,"p1")
cluster_path = os.path.join(qr_reg_data,"cluster")
babel_path = os.path.join(qr_reg_data,"babel")
charmm_path = os.path.join(qr_reg_data,"charmm")

def datasets():
  return [pdb_path,
          cluster_path,
          babel_path,
          charmm_path]

def run():
  regression_tests = [
    #"test_reg_00_charge.py",
    #"test_reg_01_finalise.py",
    #"test_reg_02_finalise.py",
    #"test_reg_03_finalise.py",
    "test_reg_04_cluster.py",
  ]

  for regression_test in regression_tests:
      for dataset in datasets():
        regression_test = os.path.join(qr_reg_tests,regression_test)
        print "Running regression test: {}  on dataset: {} ".format(regression_test,dataset)
        easy_run.call("cctbx.python %s %s"%(regression_test,dataset))

if(__name__ == "__main__"):
  t0 = time.time()
  # TODO: tests can be run on an external dataset, pass in arg for filepath.
  run()
  print "Total time (all tests): %6.2f"%(time.time()-t0)
  print "OK"
