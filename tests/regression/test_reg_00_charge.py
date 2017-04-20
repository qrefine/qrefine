import os, sys
import time
from StringIO import StringIO
import iotbx
from iotbx import pdb
import libtbx.load_env
from libtbx import easy_run
from libtbx import easy_mp,env

import finalise
import charges
import completion
from utils import hierarchy_utils

qr_repo_parent = libtbx.env.find_in_repositories("qrefine")

### TODO: after this test works, split it into charge|completion|finalise


pdbs = {"PRO_terminal" : """
CRYST1   42.664   49.718   66.065 109.25  94.96  99.24 P 1           4
ATOM    874  N   PRO A 115      -9.167  -7.159   4.783  1.00 30.39           N
ATOM    875  CA  PRO A 115      -9.350  -8.630   5.081  1.00 30.17           C
ATOM    876  C   PRO A 115      -9.382  -9.587   3.862  1.00 32.77           C
ATOM    877  O   PRO A 115     -10.069  -9.304   2.880  1.00 33.78           O
ATOM    878  CB  PRO A 115     -10.669  -8.662   5.838  1.00 28.12           C
ATOM    879  CG  PRO A 115     -10.665  -7.351   6.592  1.00 29.06           C
ATOM    880  CD  PRO A 115     -10.013  -6.329   5.660  1.00 30.28           C
ATOM    881  N  AGLY A 116      -8.745 -10.760   4.020  0.50 31.15           N
ATOM    882  CA AGLY A 116      -8.770 -11.833   3.012  0.50 27.99           C
ATOM    883  C  AGLY A 116      -7.959 -13.081   3.370  0.50 26.73           C
ATOM    884  O  AGLY A 116      -7.545 -13.840   2.487  0.50 26.23           O
ATOM    885  N  BGLY A 116      -8.621 -10.686   3.937  0.50 35.25           N
ATOM    886  CA BGLY A 116      -8.079 -11.393   2.744  0.50 36.56           C
ATOM    887  C  BGLY A 116      -8.989 -12.091   1.734  0.50 36.99           C
ATOM    888  O  BGLY A 116      -9.994 -11.547   1.293  0.50 37.57           O
ATOM    889  N  AGLY A 117      -7.748 -13.313   4.661  0.50 25.75           N
ATOM    890  CA AGLY A 117      -6.947 -14.457   5.099  0.50 24.75           C
ATOM    891  C  AGLY A 117      -5.660 -14.520   4.307  0.50 24.36           C
ATOM    892  O  AGLY A 117      -5.284 -13.541   3.660  0.50 23.17           O
ATOM    893  N  BGLY A 117      -8.582 -13.285   1.315  0.50 36.89           N
ATOM    894  CA BGLY A 117      -9.277 -14.011   0.264  0.50 35.40           C
ATOM    895  C  BGLY A 117      -8.318 -14.993  -0.358  0.50 37.53           C
ATOM    896  O  BGLY A 117      -8.109 -14.988  -1.568  0.50 39.46           O
ATOM    897  N  ASER A 118      -4.970 -15.660   4.366  0.50 25.06           N
ATOM    898  CA ASER A 118      -3.753 -15.855   3.578  0.50 25.99           C
ATOM    899  C  ASER A 118      -3.903 -15.046   2.311  0.50 28.77           C
ATOM    900  O  ASER A 118      -4.893 -15.203   1.584  0.50 30.82           O
ATOM    901  CB ASER A 118      -3.569 -17.326   3.226  0.50 24.41           C
ATOM    902  OG ASER A 118      -4.540 -17.742   2.278  0.50 25.05           O
ATOM    903  N  BSER A 118      -7.724 -15.839   0.481  0.50 38.64           N
ATOM    904  CA BSER A 118      -6.635 -16.713   0.055  0.50 38.40           C
ATOM    905  C  BSER A 118      -5.353 -15.909  -0.093  0.50 39.49           C
ATOM    906  O  BSER A 118      -4.263 -16.472  -0.176  0.50 39.26           O
ATOM    907  CB BSER A 118      -6.980 -17.441  -1.256  0.50 37.71           C
ATOM    908  OG BSER A 118      -8.020 -16.788  -1.972  0.50 34.61           O
""",
  'GLY_terminal' : '''
CRYST1   16.614   15.265   16.902  90.00  90.00  90.00 P 1
HETATM    1  N   GLY A   1      -2.476  -2.684  -3.356  1.00 20.00      A    N+1
HETATM    5  CA  GLY A   1      -2.620  -1.390  -2.713  1.00 20.00      A    C  
HETATM    8  C   GLY A   1      -1.777  -1.333  -1.440  1.00 20.00      A    C  
HETATM    9  O   GLY A   1      -1.349  -2.331  -0.968  1.00 20.00      A    O  
HETATM   10  N   ALA A   2      -1.518  -0.057  -0.799  1.00 20.00      A    N  
HETATM   12  CA  ALA A   2      -0.808  -0.011   0.467  1.00 20.00      A    C  
HETATM   14  CB  ALA A   2      -1.723   0.085   1.686  1.00 20.00      A    C  
HETATM   18  C   ALA A   2       0.356   0.977   0.488  1.00 20.00      A    C  
HETATM   19  O   ALA A   2       0.271   2.005  -0.094  1.00 20.00      A    O  
HETATM   20  N   ALA A   3       1.566   0.657   1.223  1.00 20.00      A    N  
HETATM   22  CA  ALA A   3       2.684   1.584   1.220  1.00 20.00      A    C  
HETATM   24  CB  ALA A   3       3.994   0.806   1.318  1.00 20.00      A    C  
HETATM   28  C   ALA A   3       2.525   2.581   2.367  1.00 20.00      A    C  
HETATM   29  O   ALA A   3       2.383   2.166   3.546  1.00 20.00      A    O  
''',
  'helix' : """
CRYST1   16.291   18.744   30.715  90.00  90.00  90.00 P 1
SCALE1      0.061384  0.000000  0.000000        0.00000
SCALE2      0.000000  0.053350  0.000000        0.00000
SCALE3      0.000000  0.000000  0.032557        0.00000
ATOM      1  N   GLY A  87       6.046  11.922   5.000  1.00 80.00           N
ATOM      2  CA  GLY A  87       6.777  11.037   5.889  1.00 80.00           C
ATOM      3  C   GLY A  87       7.208  11.727   7.197  1.00 80.00           C
ATOM      4  O   GLY A  87       7.552  11.016   8.143  1.00 80.00           O
ATOM      5  HA2 GLY A  87       6.151  10.181   6.142  1.00 80.00           H
ATOM      6  HA3 GLY A  87       7.670  10.673   5.382  1.00 80.00           H
ATOM      7  HT1 GLY A  87       5.481  11.402   4.424  1.00 80.00           H
ATOM      8  HT2 GLY A  87       6.362  12.695   4.526  1.00 80.00           H
ATOM      9  N   GLY A  88       7.248  13.042   7.248  1.00 80.00           N
ATOM     10  CA  GLY A  88       7.861  13.744   8.387  1.00 80.00           C
ATOM     11  C   GLY A  88       7.082  13.437   9.673  1.00 80.00           C
ATOM     12  O   GLY A  88       7.692  13.174  10.712  1.00 80.00           O
ATOM     13  H   GLY A  88       6.872  13.660   6.529  1.00 80.00           H
ATOM     14  HA2 GLY A  88       8.894  13.419   8.512  1.00 80.00           H
ATOM     15  HA3 GLY A  88       7.847  14.820   8.212  1.00 80.00           H
ATOM     16  N   GLY A  89       5.753  13.529   9.640  1.00 80.00           N
ATOM     17  CA  GLY A  89       5.000  13.335  10.878  1.00 80.00           C
ATOM     18  C   GLY A  89       5.204  11.952  11.421  1.00 80.00           C
ATOM     19  O   GLY A  89       5.421  11.767  12.629  1.00 80.00           O
ATOM     20  H   GLY A  89       5.191  13.727   8.812  1.00 80.00           H
ATOM     21  HA2 GLY A  89       5.328  14.057  11.626  1.00 80.00           H
ATOM     22  HA3 GLY A  89       3.937  13.484  10.690  1.00 80.00           H
ATOM     23  N   GLY A  90       5.138  10.940  10.551  1.00 80.00           N
ATOM     24  CA  GLY A  90       5.317   9.584  11.016  1.00 80.00           C
ATOM     25  C   GLY A  90       6.702   9.325  11.564  1.00 80.00           C
ATOM     26  O   GLY A  90       6.869   8.577  12.510  1.00 80.00           O
ATOM     27  H   GLY A  90       4.967  11.033   9.550  1.00 80.00           H
ATOM     28  HA2 GLY A  90       4.592   9.374  11.802  1.00 80.00           H
ATOM     29  HA3 GLY A  90       5.140   8.893  10.192  1.00 80.00           H
ATOM     30  N   GLY A  91       7.717   9.953  10.967  1.00 80.00           N
ATOM     31  CA  GLY A  91       9.078   9.825  11.482  1.00 80.00           C
ATOM     32  C   GLY A  91       9.167  10.416  12.883  1.00 80.00           C
ATOM     33  O   GLY A  91       9.768   9.836  13.795  1.00 80.00           O
ATOM     34  H   GLY A  91       7.631  10.546  10.141  1.00 80.00           H
ATOM     35  HA2 GLY A  91       9.362   8.774  11.523  1.00 80.00           H
ATOM     36  HA3 GLY A  91       9.773  10.354  10.830  1.00 80.00           H
ATOM     37  N   GLY A  92       8.569  11.589  13.066  1.00 80.00           N
ATOM     38  CA  GLY A  92       8.571  12.212  14.379  1.00 80.00           C
ATOM     39  C   GLY A  92       7.807  11.339  15.374  1.00 80.00           C
ATOM     40  O   GLY A  92       8.274  11.127  16.495  1.00 80.00           O
ATOM     41  H   GLY A  92       8.086  12.120  12.341  1.00 80.00           H
ATOM     42  HA2 GLY A  92       9.595  12.336  14.730  1.00 80.00           H
ATOM     43  HA3 GLY A  92       8.094  13.191  14.327  1.00 80.00           H
ATOM     44  N   GLY A  93       6.648  10.831  14.954  1.00 80.00           N
ATOM     45  CA  GLY A  93       5.852  10.007  15.846  1.00 80.00           C
ATOM     46  C   GLY A  93       6.543   8.707  16.220  1.00 80.00           C
ATOM     47  O   GLY A  93       6.368   8.241  17.352  1.00 80.00           O
ATOM     48  H   GLY A  93       6.247  10.970  14.026  1.00 80.00           H
ATOM     49  HA2 GLY A  93       5.644  10.561  16.761  1.00 80.00           H
ATOM     50  HA3 GLY A  93       4.903   9.766  15.366  1.00 80.00           H
ATOM     51  N   GLY A  94       7.316   8.108  15.306  1.00 80.00           N
ATOM     52  CA  GLY A  94       8.069   6.911  15.672  1.00 80.00           C
ATOM     53  C   GLY A  94       9.095   7.229  16.749  1.00 80.00           C
ATOM     54  O   GLY A  94       9.242   6.502  17.717  1.00 80.00           O
ATOM     55  H   GLY A  94       7.436   8.416  14.341  1.00 80.00           H
ATOM     56  HA2 GLY A  94       7.389   6.147  16.049  1.00 80.00           H
ATOM     57  HA3 GLY A  94       8.588   6.519  14.797  1.00 80.00           H
ATOM     58  N   GLY A  95       9.799   8.362  16.572  1.00 80.00           N
ATOM     59  CA  GLY A  95      10.718   8.799  17.602  1.00 80.00           C
ATOM     60  C   GLY A  95       9.986   9.053  18.925  1.00 80.00           C
ATOM     61  O   GLY A  95      10.478   8.633  19.982  1.00 80.00           O
ATOM     62  H   GLY A  95       9.748   8.967  15.753  1.00 80.00           H
ATOM     63  HA2 GLY A  95      11.479   8.036  17.764  1.00 80.00           H
ATOM     64  HA3 GLY A  95      11.208   9.721  17.290  1.00 80.00           H
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
TER
END
""",
  'water' : '''
HETATM    1  O   HOH A   1      -0.210   0.000  -0.296  1.00 20.00      A    O  
HETATM    2  H1  HOH A   1       0.733   0.000  -0.296  1.00 20.00      A    H  
HETATM    3  H2  HOH A   1      -0.524   0.000   0.593  1.00 20.00      A    H  
''',
  'pro' : '''
ATOM      1  N   PRO A   2     -25.598  25.222  75.213  1.00101.82           N
ATOM      2  CA  PRO A   2     -25.442  24.229  76.278  1.00101.04           C
ATOM      3  C   PRO A   2     -24.538  24.735  77.395  1.00102.01           C
ATOM      4  O   PRO A   2     -24.789  25.799  77.962  1.00106.79           O
ATOM      5  CB  PRO A   2     -24.800  23.031  75.558  1.00103.62           C
ATOM      6  CG  PRO A   2     -24.532  23.488  74.125  1.00 98.65           C
ATOM      7  CD  PRO A   2     -24.661  24.975  74.108  1.00 98.90           C
ATOM      8  H2  PRO A   2     -26.468  25.439  75.129  1.00101.82           H
ATOM      9  H3  PRO A   2     -25.124  25.959  75.418  1.00101.82           H
ATOM     10  HA  PRO A   2     -26.304  23.972  76.641  1.00101.04           H
ATOM     11  HB2 PRO A   2     -23.977  22.703  75.922  1.00103.62           H
ATOM     12  HB3 PRO A   2     -25.448  22.309  75.544  1.00103.62           H
ATOM     13  HG2 PRO A   2     -23.635  23.224  73.866  1.00 98.65           H
ATOM     14  HG3 PRO A   2     -25.186  23.086  73.532  1.00 98.65           H
ATOM     15  HD2 PRO A   2     -23.805  25.394  74.287  1.00 98.90           H
ATOM     16  HD3 PRO A   2     -25.039  25.275  73.267  1.00 98.90           H
''',
  'gly' : '''
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  H1  GLY A   1      -9.778   5.000   6.325  1.00 16.77           H
ATOM      6  H2  GLY A   1      -8.884   3.890   6.608  1.00 16.77           H
ATOM      7  H3  GLY A   1      -8.340   5.184   6.231  1.00 16.77           H
ATOM      8  HA2 GLY A   1      -9.928   3.856   4.426  1.00 16.57           H
ATOM      9  HA3 GLY A   1      -8.858   4.970   4.084  1.00 16.57           H
''',
  'c_terminal_capping' : '''
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
        }

def test_qxyz_non_zero():
  def _check_non_zero_charge(filename):
    for line in open(filename, 'rb').readlines():
      tmp = line.split()
      if len(tmp)<2: continue
      assert float(tmp[0])!=0, 'no partial charge %s' % line
  for residue in ['pro',
                  'gly',
                  ]:
    tf='%s.pdb' % residue
    f=file(tf, "wb")
    f.write(pdbs[residue])
    f.close()
    pdb_inp = pdb.input(tf)
    hierarchy = pdb_inp.construct_hierarchy()
    charges.write_pdb_hierarchy_qxyz_file(hierarchy,
                                               'test_%s.dat' % residue,
                                             )
    _check_non_zero_charge('test_%s.dat' % residue)

def test_qxyz_xyzq():
  tf='water.pdb'
  f=file(tf, "wb")
  f.write(pdbs["water"])
  f.close()
  pdb_inp = pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  if  os.path.exists('test_water.dat'): os.remove('test_water.dat')
  charges.write_pdb_hierarchy_qxyz_file(hierarchy,
                                             'test_water.dat',
                                             )
  assert not os.path.exists('test_water.dat')
  charges.write_pdb_hierarchy_qxyz_file(hierarchy,
                                             'test_water.dat',
                                             exclude_water=False,
                                             )
  tst_str = '''\
3  
  
-0.408  -0.21  0.0  -0.296    
0.204  0.733  0.0  -0.296    
0.204  -0.524  0.0  0.593'''
  lines = open('test_water.dat', 'rb').read()
  assert lines.strip()==tst_str, '%s %s' % (tst_str, lines)
  os.remove('test_water.dat')
  charges.write_pdb_hierarchy_xyzq_file(hierarchy,
                                             'test_water.dat',
                                             exclude_water=False,
                                             )
  tst_str = '''\
-0.21  0.0  -0.296  -0.408    
0.733  0.0  -0.296  0.204    
-0.524  0.0  0.593  0.204'''
  lines = open('test_water.dat', 'rb').read()
  assert lines.strip()==tst_str, '%s'% (lines)
  os.remove('test_water.dat')

def test_terminal_and_alt_loc(residue):
  tf = '%s_terminal.pdb' % residue
  f=file(tf, "wb")
  f.write(pdbs["%s_terminal" % residue])
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

def test_1yjp_charge():
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  tf='%s/qr-core/tests/datasets/1yjp.pdb' % qr_repo_parent
  try: os.symlink(tf, os.path.basename(tf))
  except: pass
  tf = os.path.basename(tf)
  pdb_inp = pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  try:
    charge = charges.calculate_pdb_hierarchy_charge(hierarchy)
    assert 0
  except Exception, e:
    assert e.message.find('no hydrogens')>-1
  cmd = 'iotbx.python %s/qr-core/finalise.py %s' % (qr_repo_parent, tf)
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  charge = charges.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  assert charge==0, 'charge of 1yjp should be zero not %s' % charge

def test_terminal_charge(residue, charge=0):
  # must run after other PRO
  tf = '%s_terminal_complete.pdb' % residue
  ppf = hierarchy_utils .get_processed_pdb(pdb_filename=tf)
  inter_residue_bonds = charges.get_inter_residue_bonds(ppf)
  # should the hierarchy come from ppf???
  pdb_inp = iotbx.pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  hetero_charges = charges.get_hetero_charges(pdb_inp)
  if not hetero_charges:
    # some defaults
    hetero_charges = charges.default_ion_charges
  total_charge = charges.calculate_pdb_hierarchy_charge(
    hierarchy,
    hetero_charges=hetero_charges,
    inter_residue_bonds=inter_residue_bonds,
    verbose=True,
  )
  assert total_charge==charge, "total_charge: %d, charge:%d"%(total_charge,charge)

def test_PRO_terminal_charge():
  test_terminal_charge('PRO')

def test_GLY_terminal_charge():
  test_terminal_charge('GLY')

def test_capping_of_C_terminal():
  tf = 'c_terminal_capping.pdb'
  f=file(tf,'wb')
  f.write(pdbs['c_terminal_capping'])
  f.close()
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  cmd = 'iotbx.python %s/qr-core/finalise.py model_completion=False %s' % (
    qr_repo_parent,
    tf,
    )
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_capping.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  must_find = ['OXT']
  for atom in hierarchy.atoms():
    if atom.name.strip() in must_find:
      must_find.remove(atom.name.strip())
  assert not must_find

def test_helix():
  tf = 'helix.pdb'
  f=file(tf, "wb")
  f.write(pdbs["helix"])
  f.close()
  pdb_inp=pdb.input(tf)
  hierarchy = pdb_inp.construct_hierarchy()
  charge = charges.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  assert charge==0, 'charge of helix should be zero not %s' % charge
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  cmd = 'iotbx.python %s/qr-core/finalise.py %s' % (qr_repo_parent, tf)
  easy_run.call(cmd)
  pdb_inp = pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  must_find = ['H1', 'H2', 'H3', 'OXT']
  for atom in hierarchy.atoms():
    if atom.name.strip() in must_find:
      must_find.remove(atom.name.strip())
  assert not must_find
  pdb_inp=pdb.input(tf.replace('.pdb', '_complete.pdb'))
  hierarchy = pdb_inp.construct_hierarchy()
  charge = charges.calculate_pdb_hierarchy_charge(hierarchy, verbose=1)
  assert charge==1, 'charge of helix should be one not %s' % charge

def test_charge_for_charmm_pdbs(only_i=None):
  from charges import calculate_pdb_hierarchy_charge
  charge_dict = {'3kyi': -12,
                 '2oy0': 16,
                 '1y1l': -8,
                 '3dtj': 0,
                 '3tz9': -10,
                 '4rnf': -13,
                 '2jee': -32,
                 '4k2r': -1,
                 '5d12': -11,
                 '1il5': 0,
                 '1va7': -8,
                 '2oeq': -5,
                 '4drw': -16,
                 '4xa1': -45,
                 '2ghj': 46,
                 '1ok9': -14, 
		 '3nak': 0, 
		 '2x10': -22, 
		 '3oe9': 27, 
		 '4fsx': -42, 
		 '4p7h': -19,
		 '5diz': -25, 
		 '3uds': -8, 
		 '3uj4': 0, 
		 '4ctd': -44,
                 '1byz': -4, 
                 '1lzt': 8, 
		 '1vfy': -1, 
		 '1m24': -2, 
		 '1i07': 4, 
		 '1opd': -6, 
  		 '1vbw': 8, 
		 '1rfs': -4, 
		 '1ly2': 7, 
		 '1a7y': 0, 
		 '1v7s': 8,
	         '2wpz': 0, 
		 '2akf': -3,
		 '4lzt': 8,
		 '3ovj': 4,
		 '4lzl': -5,
		 '4w71': 0,
		 '4uiv': 4,
		 '4rp6': -1,
		 '3u29': 0,
		 '2omq': -4,
		 '2ol9': 0,
		 '4xfo': 0,
		 '2xmu': -10,
		 '2xmt': -10,
		 '4uit': 4,
		 '4uiu': 4,
		 '4onk': 0,
		 '2f2n': 8,
		 '4uiw': 4,
		 '5cgo': 6,
		 '2y3j': 0,
		 '2i1u': -1,
		 '4w67': 0, 
		 '3osm': 9,
		 '4wxt': -4,
		 '3o2h': -1,
		 '2f30': 8,
		 '4w5y': 0,
		 '5e5v': 0,
		 '2ona': 0, 
		 '4itk': -9,
		 '5e61': 0,
		 '4kdw': -17,
		 '4z0w': -2,
		 '4uby': 0,
		 '2w9r': -1,
		 '3t4f': 0,
		 '3pzz': 0,
		 '2y3k': 0}
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  pdb_dir = '%s/qr-core/tests/charmm' % qr_repo_parent
  pdb_files = os.listdir(pdb_dir)
  for i, pdb_file in enumerate(pdb_files):
    if only_i is not None and i!=only_i: continue
    if pdb_file[:-4] not in charge_dict: continue
    if pdb_file.endswith(".pdb"):
      if pdb_file in ['2i1u.pdb',
                      '4wxt.pdb',
      ]:
        # disuphide bridge = in 2i1u 2.94 Phenix says yes, Charmm says no
        continue
      pdb_file_path = os.path.join(pdb_dir, pdb_file)
      charge = charges.get_total_charge_from_pdb(pdb_file_path)
      assert charge==charge_dict[pdb_file[:-4]], \
        '%s charge is %d, charmm charge is %d,  no matchy matchy' % (
          pdb_file[:-4],
          charge,
          charge_dict[pdb_file[:-4]],
        )

def test_charge_of_neutral_terminal():
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  charge = 0
  tf = "%s/qr-core/tests/babel/clusters/neutral_nterminal.pdb" % qr_repo_parent
  charge_neutral_nterminal = charges.get_total_charge_from_pdb(tf)
  assert charge == charge_neutral_nterminal, 'no match %s %s' % (
    charge,
    charge_neutral_nterminal,
    )
  tf = "%s/qr-core/tests/babel/clusters/neutral_cterminal.pdb" % qr_repo_parent
  charge_neutral_cterminal = charges.get_total_charge_from_pdb(tf)
  assert charge == charge_neutral_cterminal, 'no match %s %s' % (
    charge,
    charge_neutral_cterminal,
    )

def test_capping_of_cluster_complete(only_i=None):
  assert  qr_repo_parent, 'Set environmental variable %s' % qr_repo_parent_env
  pdb_dir = '%s/qr-core/tests/babel' % qr_repo_parent
  babel_dir = os.path.join(pdb_dir, 'capping')
  cluster_dir = os.path.join(pdb_dir, 'clusters')
  cluster_files = os.listdir(cluster_dir)
  for i, cluster_file in enumerate(cluster_files):
    if only_i is not None and i!=only_i: continue
    if cluster_file.find('capping')>-1: continue
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
        '%s atom size after babel capping: %d, after run_cluster_complete: %d' %(cluster_file, babel_size, result_size)

def run(prefix = "tst_reg_01"):
  """
  Exercise structure preparation including charge, capping, completion
  """
  tests = [
    [test_GLY_terminal_and_alt_loc, 1],
    [test_PRO_terminal_and_alt_loc, 1],
    [test_capping_of_C_terminal, 1],
    [test_charge_of_neutral_terminal, 1],
    [test_qxyz_non_zero, 1],
    [test_helix, 1],
    [test_qxyz_xyzq, 1],
    [test_1yjp_charge, 1],
    [test_charge_for_charmm_pdbs, 80],
    [test_capping_of_cluster_complete, 65],
    # need files from alt loc test
    [test_GLY_terminal_charge, 1],
    [test_PRO_terminal_charge, 1],
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
  # for testing the tests
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
