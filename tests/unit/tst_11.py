from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import time, os
import shutil
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
import os
import libtbx.load_env
from qrefine.tests.unit import batch_run_finalise
from qrefine.tests.unit import skip
from qrefine.tests.unit import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix, verbose=False, only_code='3'):
  """
  Exercise structure completion by finalise.py
  Runs only one example out of long list.
  """
  pdb_dir_cluster=os.path.join(qr_unit_tests,"babel_pdbs","clusters")
  expected_list_cluster = []
  for filename in os.listdir(pdb_dir_cluster):
    if filename.endswith('.pdb'):
      filename = filename.replace('_refine_001','')
      if filename not in skip.skip:
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
                       "1a7y"
                       ]
  pdb_dir_p1 = os.path.join(qrefine,"tests", "regression","datasets","p1")
  complete_pdbs(expected_list_cluster, pdb_dir_cluster, only_code, verbose)
  complete_pdbs(expected_list_p1, pdb_dir_p1, only_code, verbose)

def complete_pdbs(expected_list, pdb_dir, only_code, verbose):
  if(not (os.path.isdir(pdb_dir))):
    raise Sorry(pdb_dir + " not exist.  Please get its repository on GitHub")
  if 1:
    if verbose: print('"%s"' % pdb_dir)
    no_error_list = batch_run_finalise.run(
      pdb_dir, nproc=1, only_code=only_code, verbose=verbose)
    if verbose: print(no_error_list)
    if only_code is not None:
      assert only_code in no_error_list, only_code
    else:
      expected_list.sort()
      no_error_list.sort()
      assert approx_equal(expected_list, no_error_list),'%s has different pdbs passing finalise.py'%pdb_dir
    if verbose: print("expected",expected_list)
    if verbose: print("no error list",no_error_list)
  else:
    if verbose: print("skip")

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
