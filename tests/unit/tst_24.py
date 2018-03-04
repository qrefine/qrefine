import os, sys
import run_tests
from libtbx import easy_run

pdb_lines = '''
CRYST1  133.085  137.481   99.277  90.00  90.00  90.00 P 1
SCALE1      0.007514  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007274  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010073        0.00000
HETATM 2128  C1  MTNZO 101       7.968  -1.961   4.199  1.00 65.67           C
HETATM 2129  C2  MTNZO 101       9.061  -2.465   3.517  1.00 58.19           C
HETATM 2130  C3  MTNZO 101       8.677  -3.663   2.717  1.00 55.96           C
HETATM 2131  C4  MTNZO 101       9.757  -4.506   2.044  1.00 59.97           C
HETATM 2132  C5  MTNZO 101       7.568  -4.266   3.282  1.00 63.60           C
HETATM 2133  C6  MTNZO 101       8.010  -5.326   4.292  1.00 60.07           C
HETATM 2134  C7  MTNZO 101       6.664  -4.873   2.205  1.00 53.07           C
HETATM 2135  C8  MTNZO 101       7.527  -0.643   3.566  1.00 57.50           C
HETATM 2136  C9  MTNZO 101       8.310  -1.762   5.679  1.00 64.29           C
HETATM 2137  N1  MTNZO 101       6.826  -3.062   4.018  1.00 49.95           N
HETATM 2138  O1  MTNZO 101       6.340  -3.478   5.267  1.00 61.58           O
HETATM 2139  S1  MTNZO 101      11.367  -3.684   2.227  1.00 49.68           S
HETATM 2140  H2  MTNZO 101       9.720  -1.888   3.203  0.00 58.19           H
HETATM 2141  H11 MTNZO 101       7.060  -3.692   5.630  0.00 30.00           H
HETATM 2142  H41 MTNZO 101       9.556  -4.602   1.132  0.00 59.97           H
HETATM 2143  H42 MTNZO 101       9.790  -5.353   2.448  0.00 59.97           H
HETATM 2144  H61 MTNZO 101       8.594  -5.932   3.875  0.00 60.07           H
HETATM 2145  H62 MTNZO 101       7.257  -5.787   4.612  0.00 60.07           H
HETATM 2146  H63 MTNZO 101       8.450  -4.910   5.010  0.00 60.07           H
HETATM 2147  H71 MTNZO 101       7.174  -5.429   1.645  0.00 53.07           H
HETATM 2148  H72 MTNZO 101       5.986  -5.377   2.614  0.00 53.07           H
HETATM 2149  H73 MTNZO 101       6.279  -4.185   1.696  0.00 53.07           H
HETATM 2150  H81 MTNZO 101       8.242  -0.035   3.578  0.00 57.50           H
HETATM 2151  H82 MTNZO 101       7.261  -0.796   2.678  0.00 57.50           H
HETATM 2152  H83 MTNZO 101       6.806  -0.288   4.052  0.00 57.50           H
HETATM 2153  H91 MTNZO 101       8.687  -2.550   6.022  0.00 64.29           H
HETATM 2154  H92 MTNZO 101       8.924  -1.057   5.766  0.00 64.29           H
HETATM 2155  H93 MTNZO 101       7.528  -1.556   6.156  0.00 64.29           H
'''
def run(prefix):
  fn='test_mtn_1.pdb'
  f=file(fn, 'wb')
  f.write(pdb_lines)
  f.close()
  cmd = 'qr.finalise action=capping %s' % fn
  if 0: print cmd
  easy_run.go(cmd)
  fnc = '%s_capping.pdb' % fn.replace('.pdb','')
  f=file(fnc, 'rb')
  lines = f.read()
  f.close()
  assert lines.find(' H1 ')>-1
  os.remove(fn)
  os.remove(fnc)
  return 0

if __name__=='__main__':
  rc = run_tests.runner(function=run, prefix="tst_24", disable=False)
  assert not rc, 'tst_24 rc: %s' % rc
