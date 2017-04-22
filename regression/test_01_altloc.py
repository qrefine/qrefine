import os, sys
import time
from iotbx import pdb
import libtbx.load_env
from libtbx import easy_run
from libtbx import easy_mp
import qrefine.core.charges as charges
import test_data

qr_repo_parent = libtbx.env.find_in_repositories("qrefine")

#unit test?

def test_terminal_and_alt_loc(residue):
  tf = '%s_terminal.pdb' % residue
  f=file(tf, "wb")
  f.write(test_data.pdbs["%s_terminal" % residue])
  f.close()
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  cmd = 'iotbx.python %s/qr-core/finalise.py %s' % (qr_repo_parent, tf)
  print cmd
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  must_find = ['H2', 'H3', 'OXT']
  if residue!='PRO': must_find.append('H1')
  for atom in hierarchy.atoms():
    if residue=='PRO': assert atom.name.strip()!='H1'
    if atom.name.strip() in must_find:
      must_find.remove(atom.name.strip())
  assert not must_find, 'must find %s is not empty' % must_find

def test_PRO_terminal_and_alt_loc():
  test_terminal_and_alt_loc('PRO')

def test_GLY_terminal_and_alt_loc():
  test_terminal_and_alt_loc('GLY')


def run(prefix = "tst_reg_01_altloc"):
  """
  Exercise structure preparation including charge, finalise, completion
  """
  tests = [
    [test_GLY_terminal_and_alt_loc, 1],
    [test_PRO_terminal_and_alt_loc, 1],
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
