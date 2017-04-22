import os, sys
import time
import iotbx
import test_data
from iotbx import pdb
import libtbx.load_env
from libtbx import easy_mp,env
from libtbx import easy_run
from StringIO import StringIO
from qrefine.core import charges
from qrefine.core.utils import hierarchy_utils

qr_repo_parent = libtbx.env.find_in_repositories("qrefine")
pdb_dir  = os.path.join(qr_repo_parent,"regression/datasets")

#unit test

c_terminal_capping = '''
CRYST1   16.291   18.744   30.715  90.00  90.00  90.00 P 1
ATOM     65  N   GLY A  96       8.848   9.743  18.892  1.00 80.00           N
ATOM     66  CA  GLY A  96       8.100   9.974  20.115  1.00 80.00           C
ATOM     67  C   GLY A  96       7.693   8.670  20.785  1.00 80.00           C
ATOM     68  O   GLY A  96       7.804   8.505  21.994  1.00 80.00           O
ATOM     69  H   GLY A  96       8.430  10.144  18.053  1.00 80.00           H
ATOM     70  HA2 GLY A  96       8.710  10.548  20.813  1.00 80.00           H
ATOM     71  HA3 GLY A  96       7.199  10.545  19.890  1.00 80.00           H
ATOM     72  N   GLY A  97       7.185   7.734  19.986  1.00 80.00           N
ATOM     73  CA  GLY A  97       6.758   6.466  20.521  1.00 80.00           C
ATOM     74  C   GLY A  97       7.928   5.704  21.152  1.00 80.00           C
ATOM     75  O   GLY A  97       7.822   5.143  22.238  1.00 80.00           O
ATOM     76  H   GLY A  97       7.062   7.832  18.978  1.00 80.00           H
ATOM     77  HA2 GLY A  97       5.995   6.627  21.282  1.00 80.00           H
ATOM     78  HA3 GLY A  97       6.333   5.856  19.724  1.00 80.00           H
ATOM     79  N   GLY A  98       9.077   5.703  20.473  1.00 80.00           N
ATOM     80  CA  GLY A  98      10.259   5.038  20.999  1.00 80.00           C
ATOM     81  C   GLY A  98      10.776   5.677  22.259  1.00 80.00           C
ATOM     82  O   GLY A  98      11.190   5.000  23.207  1.00 80.00           O
ATOM     83  H   GLY A  98       9.215   6.150  19.567  1.00 80.00           H
ATOM     84  HA2 GLY A  98      10.022   3.996  21.214  1.00 80.00           H
ATOM     85  HA3 GLY A  98      11.051   5.063  20.251  1.00 80.00           H
ATOM     86  N   GLY A  99      10.812   6.998  22.323  1.00 80.00           N
ATOM     87  CA  GLY A  99      11.291   7.723  23.477  1.00 80.00           C
ATOM     88  C   GLY A  99      10.302   7.784  24.656  1.00 80.00           C
ATOM     89  O   GLY A  99      10.656   8.287  25.715  1.00 80.00           O
ATOM     90  OXT GLY A  99       9.077   7.282  24.507  1.00 80.00           O
ATOM     91  H   GLY A  99      10.505   7.607  21.565  1.00 80.00           H
ATOM     92  HA2 GLY A  99      12.210   7.256  23.832  1.00 80.00           H
ATOM     93  HA3 GLY A  99      11.525   8.746  23.181  1.00 80.00           H
ATOM     94  HXT GLY A  99       8.983   6.901  23.568  1.00 80.00           H
''',

def test_capping_of_C_terminal():
  tf = 'c_terminal_capping.pdb'
  f=file(tf,'wb')
  f.write(c_terminal_capping)
  f.close()
  #we can switch to qr.finalise
  cmd = 'qr.finalise model_completion=False %s' % (tf)
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_capping.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  must_find = ['OXT']
  for atom in hierarchy.atoms():
    if atom.name.strip() in must_find:
      must_find.remove(atom.name.strip())
  assert not must_find

def test_capping_of_cluster_complete(only_i=None):
  babel_dir = os.path.join(pdb_dir, 'finalise')
  cluster_dir = os.path.join(pdb_dir, 'chunk')

  cluster_files = os.listdir(cluster_dir)
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

def run(prefix = "tst_reg_01"):
  """
  Exercise structure preparation including charge, finalise, completion
  """
  tests = [
    [test_capping_of_C_terminal, 1],
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
