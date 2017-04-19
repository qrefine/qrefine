from __future__ import division

import os
import time
import requests
#from libtbx import easy_run

def build_qrefine():
  print "Building qrefine"
  import bootstrap
  bootstrap.run()
  os.remove('bootstrap.py')

def get_latest_qrefine():
  print "Downloading bootstrap.py"
  f = open('bootstrap.py', 'w')
  response = requests.get('https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py')
  f.write(response.content)
  f.close()
  build_qrefine()

def run():
  get_latest_qrefine()
  tests = [
    "qr.test --unit"
    "qr.test --regression",]

  for test in tests:
    print "Running test:", test
    easy_run.call(test)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Total time (all tests): %6.2f"%(time.time()-t0)
  print "OK"


