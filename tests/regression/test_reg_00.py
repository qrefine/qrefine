from __future__ import division

import os
import shutil
import time

from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry

import batch_run_finalise


def run(prefix = "tst_reg_00"):
  """
  Exercise qr-test-cluster and qr-test-p1 structure completion by finalise.py
  """
  expected_list_cluster = ["5diz",
                           "5d12",
                           "4xa1",
                           "4rnf",
                           "4k2r",
                           "4fsx",
                           "4drw",
                           "4ctd",
                           "3uj4",
                           "3uds",
                           "3tz9",
                           "3dtj",
                           "2x10",
                           "2pro",
                           "2oy0",
                           "2oeq",
                           "2ghj",
                           "1y1l",
                           "1va7",
                           "1ok9",
                           "1il5"]
  pdb_dir_cluster="../../qr-tests-cluster/01"
  expected_list_p1 =  ["5e61",
                       "5e5v",
                       "5cgo",
                       "4xfo",
                       "4wxt",
                       "4w71",
                       "4w67",
                       "4w5y",
                       "4uiw",
                       "4uiv",
                       "4uiu",
                       "4uit",
                       "4uby",
                       "4rp6",
                       "4onk",
                       "4lzt",
                       "4lzl",
                       "4itk",
                       "3t4f",
                       "3pzz",
                       "3ovj",
                       "3osm",
                       "3o2h",
                       "2y3k",
                       "2y3j",
                       "2w9r",
                       "2ona",
                       "2omq",
                       "2ol9",
                       "2i1u",
                       "2f30",
                       "2f2n",
                       "2akf",
                       "1vfy",
                       "1v7s",
                       "1rfs",
                       "1opd",
                       "1lzt",
                       "1ly2",
                       "1i07",
                       "1a7y"]
  pdb_dir_p1 = "../../qr-tests-p1/02_pdb_selected/" 
  complete_pdbs(expected_list_cluster, pdb_dir_cluster)
  complete_pdbs(expected_list_p1, pdb_dir_p1)

def complete_pdbs(expected_list, pdb_dir):
  if(not (os.path.isdir(pdb_dir))):
    raise Sorry(pdb_dir + " not exist.  Please get its repository on GitHub")
  input_var = str(
    raw_input("run finalise.py for all pdbs in %s will take quite a while (30~60 minutes), continue Y/N : " % pdb_dir))
  if(input_var == "Y"):
    shutil.rmtree("./tmp",ignore_errors=True)
    os.makedirs("./tmp")
    no_error_list = batch_run_finalise.run(pdb_dir,nproc=4,only_code='1il5',)
    assert 0
    shutil.rmtree("./tmp")
    expected_list.sort()
    no_error_list.sort()
    assert approx_equal(expected_list, no_error_list),'%s has different pdbs passing finalise.py'%pdb_dir
  else:
    print "skip"

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
