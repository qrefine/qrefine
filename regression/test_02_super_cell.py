import os, sys
import time
from StringIO import StringIO
import iotbx
from iotbx import pdb
import libtbx.load_env
from libtbx import easy_run
from libtbx import easy_mp
import test_data
import qrefine.core.finalise
import qrefine.core.charges as charges
import qrefine.core.completion
from qrefine.core.utils import hierarchy_utils

qr_repo_parent = libtbx.env.find_in_repositories("qrefine")

#unit test

def run(prefix = "tst_reg_02_super_cell"):
  """
  Exercise structure preparation including charge, capping, completion
  """
  tests = [

   # [test_helix, 80],
    ]
  def get_test(i, j):
    func = tests[i][0]
    if j is not None: return func(j)
    else: return func()
  ####
  argss = []
  for i, (func, j) in enumerate(tests):
    if j==1: argss.append((i,None))
    else:
      for p in range(j):
        argss.append((i,p))
  # for testing the regression
  if 0:
    argss = [
      [9, 60],
      ]
  #
  passed=0
  failed=0
  print dir(easy_mp)
  for args, res, errstr in easy_mp.multi_core_run(get_test, argss, 6):
    if errstr:
      print '-'*80
      print args
      print 'RESULT - ERROR   : %s %s' % (tests[args[0]][0].func_name, args)
      print errstr
      print '-'*80
      failed+=1
    else:
      print 'RESULT - SUCCESS : %s %s' % (tests[args[0]][0].func_name, args)
      passed+=1
  print '\n\tpassed : %d\n\tfailed : %s' % (passed, failed)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
