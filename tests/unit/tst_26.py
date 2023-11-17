from __future__ import print_function
from __future__ import absolute_import
import os, sys
from qrefine.tests.unit import run_tests
from libtbx import easy_run
import libtbx.load_env

qrefine_path = libtbx.env.find_in_repositories("qrefine")

pdb_lines = '''
CRYST1   72.470   66.336   68.552  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A  30      62.619  25.986  37.359  1.00 66.84           N
ATOM      2  CA  HIS A  30      63.258  26.030  36.050  1.00 70.57           C
ATOM      3  C   HIS A  30      64.699  26.498  36.196  1.00 70.51           C
ATOM      4  O   HIS A  30      64.980  27.444  36.921  1.00 73.92           O
ATOM      5  CB  HIS A  30      62.568  26.958  35.058  1.00 70.79           C
ATOM      6  CG  HIS A  30      61.106  26.715  34.861  1.00 68.99           C
ATOM      7  ND1 HIS A  30      60.132  27.545  35.365  1.00 77.35           N
ATOM      8  CD2 HIS A  30      60.459  25.708  34.234  1.00 70.51           C
ATOM      9  CE1 HIS A  30      58.941  27.084  35.013  1.00 79.15           C
ATOM     10  NE2 HIS A  30      59.114  25.973  34.318  1.00 70.69           N
ATOM     11  OXT HIS A  30      65.598  25.924  35.581  1.00 70.51           O
ATOM     12  H   HIS A  30      61.823  26.617  37.454  0.00 66.84           H
ATOM     13  H2  HIS A  30      62.284  25.060  37.536  0.00 66.84           H
ATOM     14  H3  HIS A  30      63.282  26.241  38.063  0.00 66.84           H
ATOM     15  HA  HIS A  30      63.199  25.023  35.638  0.00 70.57           H
ATOM     16  HB2 HIS A  30      62.682  27.983  35.411  0.00 70.79           H
ATOM     17  HB3 HIS A  30      63.050  26.842  34.087  0.00 70.79           H
ATOM     18  HD2 HIS A  30      60.915  24.854  33.756  0.00 70.51           H
ATOM     19  HE1 HIS A  30      57.991  27.538  35.253  0.00 79.15           H
ATOM     20  HE2 HIS A  30      58.371  25.406  33.911  0.00 70.69           H
HETATM   21 ZN    ZN A 101      60.278  29.235  36.302  1.00 76.89          Zn
'''

def run(prefix):
  """
  Exercise qr.charges

  XXX TEST FAILS (qr.charges). Nigel?

  """
  fn='test_zn_his_charge.pdb'
  f=open(fn, 'w')
  f.write(pdb_lines)
  f.close()
  cmd = 'qr.charges %s verbose=1' % (fn)
  if 1: print(cmd)
  rc = easy_run.go(cmd)
  assert 'Charge: 2' in rc.stdout_lines
  os.remove(fn)
  return rc

if(__name__=='__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
