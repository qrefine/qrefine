from __future__ import print_function
from __future__ import absolute_import
import os, sys
from qrefine.tests.unit import run_tests
from libtbx import easy_run
import libtbx.load_env

qrefine_path = libtbx.env.find_in_repositories("qrefine")

pdb_lines = '''
CRYST1   70.834   80.043   78.266  90.00  90.00  90.00 P 1
ATOM     44  N   GLY A 334      10.917  28.755   9.158  1.00 23.45           N
ATOM     45  CA  GLY A 334      10.037  28.371  10.245  1.00 27.12           C
ATOM     46  C   GLY A 334       9.876  26.874  10.373  1.00 28.22           C
ATOM     47  O   GLY A 334       9.588  26.191   9.391  1.00 23.33           O
ATOM     48  H   GLY A 334      11.754  28.624   9.309  0.00 23.45           H
ATOM     49  HA2 GLY A 334      10.391  28.713  11.081  0.00 27.12           H
ATOM     50  HA3 GLY A 334       9.161  28.761  10.102  0.00 27.12           H
ATOM     51  N   ARG A 335      10.062  26.354  11.587  1.00 27.32           N
ATOM     52  CA  ARG A 335       9.963  24.915  11.809  1.00 35.88           C
ATOM     53  C   ARG A 335       8.575  24.393  11.467  1.00 28.28           C
ATOM     54  O   ARG A 335       8.441  23.379  10.774  1.00 26.15           O
ATOM     55  H   ARG A 335      10.245  26.811  12.292  0.00 27.32           H
ATOM     56  HA  ARG A 335      10.599  24.466  11.230  0.00 35.88           H
ATOM     57  CB AARG A 335      10.318  24.590  13.262  0.38 45.61           C
ATOM     58  CG AARG A 335      10.137  23.131  13.642  0.38 51.31           C
ATOM     59  CD AARG A 335      11.476  22.461  13.896  0.38 62.60           C
ATOM     60  NE AARG A 335      11.949  21.723  12.728  0.38 83.94           N
ATOM     61  CZ AARG A 335      12.241  20.427  12.732  0.38 75.12           C
ATOM     62  NH1AARG A 335      12.107  19.719  13.845  0.38 71.36           N
ATOM     63  NH2AARG A 335      12.665  19.836  11.624  0.38 70.05           N
ATOM     64  HB2AARG A 335      11.248  24.820  13.413  0.00 45.60           H
ATOM     65  HB3AARG A 335       9.751  25.119  13.846  0.00 45.60           H
ATOM     66  HG2AARG A 335       9.608  23.073  14.453  0.00 51.51           H
ATOM     67  HG3AARG A 335       9.694  22.663  12.917  0.00 51.51           H
ATOM     68  HD2AARG A 335      12.135  23.139  14.113  0.00 57.35           H
ATOM     69  HD3AARG A 335      11.385  21.836  14.632  0.00 57.35           H
ATOM     70  HE AARG A 335      12.066  22.160  11.996  0.00 77.91           H
ATOM     71 HH11AARG A 335      11.832  20.099  14.566  0.00 64.92           H
ATOM     72 HH12AARG A 335      12.296  18.880  13.846  0.00 64.92           H
ATOM     73 HH21AARG A 335      12.753  20.292  10.900  0.00 71.58           H
ATOM     74 HH22AARG A 335      12.853  18.997  11.629  0.00 71.58           H
ATOM     75  N   GLU A 336       7.526  25.071  11.943  1.00 27.79           N
ATOM     76  CA  GLU A 336       6.174  24.576  11.703  1.00 28.50           C
ATOM     77  C   GLU A 336       5.869  24.541  10.213  1.00 24.18           C
ATOM     78  O   GLU A 336       5.267  23.583   9.714  1.00 24.02           O
ATOM     79  CB  GLU A 336       5.148  25.428  12.454  1.00 32.66           C
ATOM     80  CG  GLU A 336       5.358  25.460  13.962  1.00 56.63           C
ATOM     81  CD  GLU A 336       4.578  24.382  14.691  1.00 98.70           C
ATOM     82  OE1 GLU A 336       3.732  24.730  15.542  1.00110.33           O
ATOM     83  OE2 GLU A 336       4.815  23.186  14.419  1.00 95.45           O
ATOM     84  H   GLU A 336       7.570  25.801  12.396  0.00 27.79           H
ATOM     85  HA  GLU A 336       6.109  23.669  12.040  0.00 28.50           H
ATOM     86  HB2 GLU A 336       5.201  26.340  12.128  0.00 32.66           H
ATOM     87  HB3 GLU A 336       4.262  25.070  12.286  0.00 32.66           H
ATOM     88  HG2 GLU A 336       6.300  25.328  14.152  0.00 56.63           H
ATOM     89  HG3 GLU A 336       5.069  26.321  14.303  0.00 56.63           H
'''

def run(prefix):
  fn='test_one_alt_loc.pdb'
  f=open(fn, 'w')
  f.write(pdb_lines)
  f.close()
  cmd = 'qr.charges %s verbose=1' % (fn)
  if 0: print(cmd)
  rc = easy_run.go(cmd)
  assert 'Charge: -2' in rc.stdout_lines
  os.remove(fn)
  return rc

if(__name__=='__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
