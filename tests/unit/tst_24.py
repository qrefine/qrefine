from __future__ import print_function
from __future__ import absolute_import
import os, sys
from qrefine.tests.unit import run_tests
from libtbx import easy_run
import libtbx.load_env

qrefine_path = libtbx.env.find_in_repositories("qrefine")

pdb_lines = '''
CRYST1   81.798   80.030   79.804  90.00  90.00  90.00 P 1
ATOM   1163  N   CYS A  78      11.141  13.078  17.954  1.00  8.15           N
ANISOU 1163  N   CYS A  78      368   2009    718    -21     61    -57       N
ATOM   1164  CA  CYS A  78      10.932  12.280  19.167  1.00  8.14           C
ANISOU 1164  CA  CYS A  78      335   2060    698     37     66   -116       C
ATOM   1165  C   CYS A  78      12.067  11.265  19.234  1.00  7.84           C
ANISOU 1165  C   CYS A  78      482   1942    556    -72    125   -132       C
ATOM   1166  O   CYS A  78      12.293  10.604  18.236  1.00 10.09           O
ANISOU 1166  O   CYS A  78      760   2315    759     85    -15   -185       O
ATOM   1167  CB  CYS A  78       9.596  11.643  19.135  1.00  9.23           C
ANISOU 1167  CB  CYS A  78      350   2159    997     -5     45    -59       C
ATOM   1168  SG  CYS A  78       9.430  10.394  20.434  1.00  9.81           S
ANISOU 1168  SG  CYS A  78      707   2210    809   -188     50    -25       S
ATOM   1169  H   CYS A  78      10.663  12.707  17.133  1.00  8.15           H
ATOM   1170  HA  CYS A  78      10.991  12.924  20.044  1.00  8.14           H
ATOM   1171  HB2 CYS A  78       8.830  12.403  19.291  1.00  9.23           H
ATOM   1172  HB3 CYS A  78       9.453  11.155  18.172  1.00  9.23           H
HETATM 1832 CU    CU A 201       7.545   9.241  20.314  0.40  8.15          Cu
ANISOU 1832 CU    CU A 201      503   1903    690    -55    119      0      Cu
'''

def run(prefix):
  fn='test_cu_cys.pdb'
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
  assert ' HG  CYS A  78' not in lines
  cmd = 'qr.charges %s verbose=1' % (fnc)
  if 0: print(cmd)
  rc = easy_run.go(cmd)
  assert 'Charge: 0' in rc.stdout_lines
  os.remove(fnc)
  return rc

if(__name__=='__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
