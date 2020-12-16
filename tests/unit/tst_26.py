from __future__ import print_function
from __future__ import absolute_import
import os, sys
from qrefine.tests.unit import run_tests
from libtbx import easy_run
import libtbx.load_env

qrefine_path = libtbx.env.find_in_repositories("qrefine")

pdb_lines = '''
CRYST1   72.470   66.336   68.552  90.00  90.00  90.00 P 1
ATOM    387  N   HIS A  30      62.619  25.986  37.359  1.00 66.84           N
ATOM    388  CA  HIS A  30      63.258  26.030  36.050  1.00 70.57           C
ATOM    389  C   HIS A  30      64.699  26.498  36.196  1.00 70.51           C
ATOM    390  O   HIS A  30      64.980  27.444  36.921  1.00 73.92           O
ATOM    391  CB  HIS A  30      62.568  26.958  35.058  1.00 70.79           C
ATOM    392  CG  HIS A  30      61.106  26.715  34.861  1.00 68.99           C
ATOM    393  ND1 HIS A  30      60.132  27.545  35.365  1.00 77.35           N
ATOM    394  CD2 HIS A  30      60.459  25.708  34.234  1.00 70.51           C
ATOM    395  CE1 HIS A  30      58.941  27.084  35.013  1.00 79.15           C
ATOM    396  NE2 HIS A  30      59.114  25.973  34.318  1.00 70.69           N
ATOM    397  H   HIS A  30      61.945  26.509  37.464  0.00 66.84           H
ATOM    398  HA  HIS A  30      63.202  25.127  35.700  0.00 70.57           H
ATOM    399  HB2 HIS A  30      62.691  27.873  35.355  0.00 70.79           H
ATOM    400  HB3 HIS A  30      63.012  26.877  34.200  0.00 70.79           H
ATOM    401  HD2 HIS A  30      60.851  24.972  33.822  0.00 70.51           H
ATOM    402  HE1 HIS A  30      58.123  27.475  35.220  0.00 79.15           H
ATOM    403  HE2 HIS A  30      58.487  25.495  33.975  0.00 70.69           H
HETATM  541 ZN    ZN A 101      60.278  29.235  36.302  1.00 76.89          ZN
TER
'''

def run(prefix):
  fn='test_zn_his_charge.pdb'
  f=open(fn, 'wb')
  f.write(bytes(pdb_lines,encoding='utf8'))
  f.close()
  cmd = 'qr.charges %s verbose=1' % (fn)
  if 0: print(cmd)
  rc = easy_run.go(cmd)
  assert 'Charge: 0' in rc.stdout_lines
  os.remove(fn)
  return rc

if(__name__=='__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
