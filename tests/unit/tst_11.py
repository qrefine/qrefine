from __future__ import division

import time
import shutil
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
import os
import libtbx.load_env
import batch_run_finalise
from qrefine.skip import skip

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix = "tst_11"):
  """
  Exercise structure completion by finalise.py
  """
  pdb_dir_cluster=os.path.join(qr_unit_tests,"babel_pdbs","clusters")
  expected_list_cluster = []
  for filename in os.listdir(pdb_dir_cluster):
    if filename.endswith('.pdb'):
      filename = filename.replace('_refine_001','')
      if filename not in skip:
        expected_list_cluster.append(os.path.splitext(filename)[0])
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
  pdb_dir_p1 = os.path.join(qrefine,"regression","datasets","p1")
  complete_pdbs(expected_list_cluster, pdb_dir_cluster)
  complete_pdbs(expected_list_p1, pdb_dir_p1)

def complete_pdbs(expected_list, pdb_dir):
  if(not (os.path.isdir(pdb_dir))):
    raise Sorry(pdb_dir + " not exist.  Please get its repository on GitHub")
  #input_var = str(
  #  raw_input("run finalise.py for all pdbs in %s will take quite a while (30~60 minutes), continue Y/N : " % pdb_dir))
  if 1: #(input_var == "Y") or 1:
    print '"%s"' % pdb_dir
    no_error_list = batch_run_finalise.run(pdb_dir,
                                           nproc=6,
                                           #only_code='1il5',
    )
    print no_error_list
    #shutil.rmtree("./tmp")
    expected_list.sort()
    no_error_list.sort()
    print expected_list
    print no_error_list
    assert approx_equal(expected_list, no_error_list),'%s has different pdbs passing finalise.py'%pdb_dir
  else:
    print "skip"

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_11"
  if(0):
    run(prefix)
    print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
  else:
    print prefix + ":  Skipped    "
