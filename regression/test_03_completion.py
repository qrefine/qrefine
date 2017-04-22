import os, sys
import time
import iotbx
import test_data
from iotbx import pdb
import libtbx.load_env
from libtbx import easy_run
from libtbx import easy_mp
from StringIO import StringIO
import qrefine.core.completion
import qrefine.core.finalise
import qrefine.core.charges
from qrefine.core.utils import hierarchy_utils

qr_repo_parent = libtbx.env.find_in_repositories("qrefine")
pdb_dir = '%s/qr-core/regression/babel' % qr_repo_parent
babel_dir = os.path.join(pdb_dir, 'finalise')
cluster_dir = os.path.join(pdb_dir, 'chunk')
cluster_files = os.listdir(cluster_dir)

def test_capping_of_cluster_complete(only_i=None):
  for i, cluster_file in enumerate(cluster_files):
    if only_i is not None and i!=only_i: continue
    if cluster_file.find('finalise')>-1: continue
    if cluster_file.endswith(".pdb") and ("temp" not in cluster_file):
      cluster_file_path = os.path.join(cluster_dir, cluster_file)
      if not os.path.exists(cluster_file):
        os.symlink(cluster_file_path, cluster_file)
      cluster_file_path = cluster_file
      cmd = "phenix.python %s/qr-core/completion.py %s model_completion=False" % (
        qr_repo_parent,
        cluster_file_path,
        )
      print cmd
      easy_run.call(cmd)    
      result_file = cluster_file_path[:-4] + "_capping.pdb" 
      babel_file = os.path.join(babel_dir, cluster_file[:-4] + "_babel.pdb")
      result_size = len(pdb.input(result_file).atoms())    
      babel_size =  len(pdb.input(babel_file).atoms())
      assert result_size ==  babel_size,\
        '%s atom size after babel finalise: %d, after run_cluster_complete: %d' %(cluster_file, babel_size, result_size)

def run(prefix = "tst_reg_03_completion"):
  """
  Exercise structure preparation including charge, finalise, completion
  """
  tests = [

    [test_capping_of_cluster_complete, 65],
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
