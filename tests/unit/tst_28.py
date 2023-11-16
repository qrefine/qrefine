from __future__ import print_function
from __future__ import absolute_import
import os, sys
from qrefine.tests.unit import run_tests
from libtbx import easy_run
import libtbx.load_env

qrefine_path = libtbx.env.find_in_repositories("qrefine")

pdb_lines = '''
CRYST1   72.470   66.336   68.552  90.00  90.00  90.00 P 1
ATOM    156  N   CYS A  12      58.225  95.502  96.618  1.00100.41           N
ATOM    157  CA  CYS A  12      57.161  94.968  95.784  1.00104.22           C
ATOM    158  C   CYS A  12      57.817  94.476  94.511  1.00 99.49           C
ATOM    159  O   CYS A  12      57.818  95.168  93.493  1.00102.24           O
ATOM    160  CB  CYS A  12      56.104  96.001  95.493  1.00105.99           C
ATOM    161  SG  CYS A  12      55.020  95.413  94.223  1.00102.06           S
ATOM    162  H   CYS A  12      58.768  96.021  96.199  0.00100.41           H
ATOM    163  HA  CYS A  12      56.701  94.246  96.240  0.00104.22           H
ATOM    164  HB2 CYS A  12      55.598  96.193  96.298  0.00105.99           H
ATOM    165  HB3 CYS A  12      56.520  96.832  95.215  0.00105.99           H
ATOM    166  HG  CYS A  12      55.661  95.132  93.248  0.00102.06           H
ATOM   1338  N   CYSSS  33      50.685  96.618  95.502  1.00100.41           N
ATOM   1339  CA  CYSSS  33      51.749  95.784  94.968  1.00104.22           C
ATOM   1340  C   CYSSS  33      51.093  94.511  94.476  1.00 99.49           C
ATOM   1341  O   CYSSS  33      51.092  93.493  95.168  1.00102.24           O
ATOM   1342  CB  CYSSS  33      52.806  95.493  96.001  1.00105.99           C
ATOM   1343  SG  CYSSS  33      53.890  94.223  95.413  1.00102.06           S
ATOM   1344  H   CYSSS  33      50.142  96.199  96.021  0.00100.41           H
ATOM   1345  HA  CYSSS  33      52.209  96.240  94.246  0.00104.22           H
ATOM   1346  HB2 CYSSS  33      53.312  96.298  96.193  0.00105.99           H
ATOM   1347  HB3 CYSSS  33      52.390  95.215  96.832  0.00105.99           H
ATOM   1348  HG  CYSSS  33      53.249  93.248  95.132  0.00102.06           H
'''

def run(prefix):
  """
  Exercise qr.finalise + capping

  XXX TEST FAILS (qr.finalise + capping). Nigel?

  """
  fn='test_cys_cys_sym.pdb'
  f=open(fn, 'w')
  f.write(pdb_lines)
  f.close()
  cmd = 'qr.finalise %s action="capping"' % (fn)
  if 0: print(cmd)
  rc = easy_run.go(cmd)
  os.remove(fn)
  fnc = '%s_capping.pdb' % fn.replace('.pdb','')
  f=open(fnc, 'r')
  lines=f.read()
  f.close()
  assert ' HG  CYS A  12' not in lines
  assert ' HG  CYSSS  33' not in lines
  cmd = 'qr.charges %s verbose=1' % (fnc)
  if 0: print(cmd)
  rc = easy_run.go(cmd)
  assert 'Charge: 0' in rc.stdout_lines
  os.remove(fnc)
  return rc

if(__name__=='__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=True)
