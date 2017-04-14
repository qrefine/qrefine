from __future__ import division

import os
import time
import iotbx.pdb
from libtbx import easy_run,easy_qsub

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_reg_tests = os.path.join(qrefine, "tests/regression")
phenix_source = os.path.join(qrefine,"../../build/setpaths_all.csh")

PBS_LIMIT=69000
#TODO fix paths
pdb_mirror = "/home/pdb/mirror/pub/pdb/data/structures/divided/pdb"
#work directory pbs
where = "/home/mark/working/save"

def submit(commands):
  easy_qsub.run(
    phenix_source=phenix_source,
    where=where,
    size_of_chunks=50,
    qsub_cmd="qsub -q qr",
    commands=commands)

def run_batch():
  commands = []
  for dirname, dirs, files in os.walk(pdb_mirror):
    for file in files:
      pdb_file_name =  os.path.join(dirname, file)
      assert os.path.isfile(pdb_file_name)
      cmd = "phenix.pdb.hierarchy %s >& %s.log" % (pdb_file_name,file)
      #cmd = "qr.run_tests --test_reg_04_cluster.py  %s >& %s.log " % (pdb_file_name, pdb_code)
      commands.append(cmd)
  print "Total:", len(commands)
  submit(commands)

def run_regression_tests():
  regression_tests = [
      "test_reg_00_charge.py",
      "test_reg_01_capping.py",
      "test_reg_02_completion.py",
      "test_reg_03_finalise.py",
      "test_reg_04_cluster.py",
      ]

  for file_name in regression_tests:
    file_name = os.path.join(qr_reg_tests,file_name)
    print "Running regression test:", file_name
    easy_run.call("cctbx.python %s"%file_name)

if(__name__ == "__main__"):
  t0 = time.time()
  run_regression_tests()
  print "Total time (all tests): %6.2f"%(time.time()-t0)
  print "OK"
