from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import easy_pickle
import time, sys
from qrefine.tests.unit import run_tests
import libtbx.load_env
from libtbx.test_utils import approx_equal
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

pdb_str_in1 = """
REMARK Residue 94 has altloc
CRYST1   18.273   20.824   31.291  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  87       7.109  12.926   5.576  1.00 80.00           N
ATOM      2  CA  GLY A  87       7.840  12.041   6.465  1.00 80.00           C
ATOM      3  C   GLY A  87       8.271  12.731   7.773  1.00 80.00           C
ATOM      4  O   GLY A  87       8.615  12.020   8.719  1.00 80.00           O
ATOM         H1  GLY A  87       7.550  12.998   4.795  1.00 80.00           H
ATOM         H2  GLY A  87       6.288  12.590   5.426  1.00 80.00           H
ATOM         H3  GLY A  87       7.036  13.739   5.955  1.00 80.00           H
ATOM      5  HA2 GLY A  87       7.287  11.273   6.677  1.00 80.00           H
ATOM      6  HA3 GLY A  87       8.627  11.706   6.006  1.00 80.00           H
ATOM      7  N   GLY A  88       8.311  14.046   7.824  1.00 80.00           N
ATOM      8  CA  GLY A  88       8.924  14.748   8.963  1.00 80.00           C
ATOM      9  C   GLY A  88       8.145  14.441  10.249  1.00 80.00           C
ATOM     10  O   GLY A  88       8.755  14.178  11.288  1.00 80.00           O
ATOM     11  H   GLY A  88       7.992  14.563   7.215  1.00 80.00           H
ATOM     12  HA2 GLY A  88       9.849  14.473   9.063  1.00 80.00           H
ATOM     13  HA3 GLY A  88       8.929  15.704   8.799  1.00 80.00           H
ATOM     14  N   GLY A  89       6.816  14.533  10.216  1.00 80.00           N
ATOM     15  CA  GLY A  89       6.063  14.339  11.454  1.00 80.00           C
ATOM     16  C   GLY A  89       6.267  12.956  11.997  1.00 80.00           C
ATOM     17  O   GLY A  89       6.484  12.771  13.205  1.00 80.00           O
ATOM     18  H   GLY A  89       6.346  14.700   9.516  1.00 80.00           H
ATOM     19  HA2 GLY A  89       6.344  14.993  12.113  1.00 80.00           H
ATOM     20  HA3 GLY A  89       5.119  14.489  11.289  1.00 80.00           H
ATOM     21  N   GLY A  90       6.201  11.944  11.127  1.00 80.00           N
ATOM     22  CA  GLY A  90       6.380  10.588  11.592  1.00 80.00           C
ATOM     23  C   GLY A  90       7.765  10.329  12.140  1.00 80.00           C
ATOM     24  O   GLY A  90       7.932   9.581  13.086  1.00 80.00           O
ATOM     25  H   GLY A  90       6.056  12.028  10.283  1.00 80.00           H
ATOM     26  HA2 GLY A  90       5.724  10.399  12.281  1.00 80.00           H
ATOM     27  HA3 GLY A  90       6.207   9.975  10.860  1.00 80.00           H
ATOM     28  N   GLY A  91       8.780  10.957  11.543  1.00 80.00           N
ATOM     29  CA  GLY A  91      10.141  10.829  12.058  1.00 80.00           C
ATOM     30  C   GLY A  91      10.230  11.420  13.459  1.00 80.00           C
ATOM     31  O   GLY A  91      10.831  10.840  14.371  1.00 80.00           O
ATOM     32  H   GLY A  91       8.702  11.456  10.847  1.00 80.00           H
ATOM     33  HA2 GLY A  91      10.400   9.895  12.077  1.00 80.00           H
ATOM     34  HA3 GLY A  91      10.762  11.284  11.468  1.00 80.00           H
ATOM     35  N   GLY A  92       9.632  12.593  13.642  1.00 80.00           N
ATOM     36  CA  GLY A  92       9.634  13.216  14.955  1.00 80.00           C
ATOM     37  C   GLY A  92       8.870  12.343  15.950  1.00 80.00           C
ATOM     38  O   GLY A  92       9.337  12.131  17.071  1.00 80.00           O
ATOM     39  H   GLY A  92       9.226  13.037  13.028  1.00 80.00           H
ATOM     40  HA2 GLY A  92      10.546  13.343  15.259  1.00 80.00           H
ATOM     41  HA3 GLY A  92       9.226  14.095  14.905  1.00 80.00           H
ATOM     42  N   GLY A  93       7.711  11.835  15.530  1.00 80.00           N
ATOM     43  CA  GLY A  93       6.915  11.011  16.422  1.00 80.00           C
ATOM     44  C   GLY A  93       7.606   9.711  16.796  1.00 80.00           C
ATOM     45  O   GLY A  93       7.431   9.245  17.928  1.00 80.00           O
ATOM     46  H   GLY A  93       7.377  11.955  14.747  1.00 80.00           H
ATOM     47  HA2 GLY A  93       6.719  11.512  17.229  1.00 80.00           H
ATOM     48  HA3 GLY A  93       6.066  10.811  15.999  1.00 80.00           H
ATOM     49  N  ATYR A  94       8.379   9.112  15.882  0.50 80.00           N
ATOM     50  CA ATYR A  94       9.132   7.915  16.248  0.50 80.00           C
ATOM     51  C  ATYR A  94      10.158   8.233  17.325  0.50 80.00           C
ATOM     52  O  ATYR A  94      10.305   7.506  18.293  0.50 80.00           O
ATOM     53  CB ATYR A  94       9.821   7.315  15.022  0.50 30.00           C
ATOM     54  CG ATYR A  94       8.864   6.790  13.975  0.50 30.00           C
ATOM     55  CD1ATYR A  94       8.372   5.492  14.045  0.50 30.00           C
ATOM     56  CD2ATYR A  94       8.453   7.592  12.920  0.50 30.00           C
ATOM     57  CE1ATYR A  94       7.498   5.009  13.091  0.50 30.00           C
ATOM     58  CE2ATYR A  94       7.577   7.117  11.960  0.50 30.00           C
ATOM     59  CZ ATYR A  94       7.103   5.826  12.051  0.50 30.00           C
ATOM     60  OH ATYR A  94       6.234   5.346  11.098  0.50 30.00           O
ATOM     61  H  ATYR A  94       8.477   9.376  15.069  0.50 80.00           H
ATOM     62  HA ATYR A  94       8.507   7.263  16.601  0.50 80.00           H
ATOM     63  HB2ATYR A  94      10.389   7.990  14.619  0.50 30.00           H
ATOM     64  HB3ATYR A  94      10.400   6.592  15.309  0.50 30.00           H
ATOM     65  HD1ATYR A  94       8.636   4.940  14.745  0.50 30.00           H
ATOM     66  HD2ATYR A  94       8.772   8.463  12.857  0.50 30.00           H
ATOM     67  HE1ATYR A  94       7.177   4.138  13.149  0.50 30.00           H
ATOM     68  HE2ATYR A  94       7.310   7.665  11.258  0.50 30.00           H
ATOM     69  HH ATYR A  94       6.081   5.946  10.530  0.50 30.00           H
ATOM     70  N  BTYR A  94       8.479   9.112  15.882  0.50 80.00           N
ATOM     71  CA BTYR A  94       9.232   7.915  16.248  0.50 80.00           C
ATOM     72  C  BTYR A  94      10.258   8.233  17.325  0.50 80.00           C
ATOM     73  O  BTYR A  94      10.405   7.506  18.293  0.50 80.00           O
ATOM     74  CB BTYR A  94       9.921   7.315  15.022  0.50 30.00           C
ATOM     75  CG BTYR A  94      10.801   6.124  15.334  0.50 30.00           C
ATOM     76  CD1BTYR A  94      10.247   4.873  15.580  0.50 30.00           C
ATOM     77  CD2BTYR A  94      12.182   6.252  15.381  0.50 30.00           C
ATOM     78  CE1BTYR A  94      11.047   3.784  15.865  0.50 30.00           C
ATOM     79  CE2BTYR A  94      12.990   5.166  15.666  0.50 30.00           C
ATOM     80  CZ BTYR A  94      12.418   3.936  15.907  0.50 30.00           C
ATOM     81  OH BTYR A  94      13.216   2.852  16.193  0.50 30.00           O
ATOM     82  H  BTYR A  94       8.620   9.399  15.084  0.50 80.00           H
ATOM     83  HA BTYR A  94       8.607   7.263  16.601  0.50 80.00           H
ATOM     84  HB2BTYR A  94       9.244   7.047  14.381  0.50 30.00           H
ATOM     85  HB3BTYR A  94      10.459   8.001  14.596  0.50 30.00           H
ATOM     86  HD1BTYR A  94       9.323   4.768  15.552  0.50 30.00           H
ATOM     87  HD2BTYR A  94      12.570   7.081  15.218  0.50 30.00           H
ATOM     88  HE1BTYR A  94      10.664   2.952  16.028  0.50 30.00           H
ATOM     89  HE2BTYR A  94      13.914   5.266  15.695  0.50 30.00           H
ATOM     90  HH BTYR A  94      14.022   3.087  16.184  0.50 30.00           H
ATOM     91  N   GLY A  95      10.862   9.366  17.148  1.00 80.00           N
ATOM     92  CA  GLY A  95      11.781   9.803  18.178  1.00 80.00           C
ATOM     93  C   GLY A  95      11.049  10.057  19.501  1.00 80.00           C
ATOM     94  O   GLY A  95      11.541   9.637  20.558  1.00 80.00           O
ATOM     95  HA2 GLY A  95      12.468   9.131  18.309  1.00 80.00           H
ATOM     96  HA3 GLY A  95      12.229  10.614  17.892  1.00 80.00           H
ATOM     97  H  AGLY A  95      10.815   9.873  16.455  0.50 80.00           H
ATOM     98  H  BGLY A  95      10.769   9.889  16.472  0.50 80.00           H
ATOM     99  N   GLY A  96       9.911  10.747  19.468  1.00 80.00           N
ATOM    100  CA  GLY A  96       9.163  10.978  20.691  1.00 80.00           C
ATOM    101  C   GLY A  96       8.756   9.674  21.361  1.00 80.00           C
ATOM    102  O   GLY A  96       8.867   9.509  22.570  1.00 80.00           O
ATOM    103  H   GLY A  96       9.562  11.084  18.757  1.00 80.00           H
ATOM    104  HA2 GLY A  96       9.701  11.503  21.305  1.00 80.00           H
ATOM    105  HA3 GLY A  96       8.370  11.500  20.491  1.00 80.00           H
ATOM    106  N   GLY A  97       8.248   8.738  20.562  1.00 80.00           N
ATOM    107  CA  GLY A  97       7.821   7.470  21.097  1.00 80.00           C
ATOM    108  C   GLY A  97       8.991   6.708  21.728  1.00 80.00           C
ATOM    109  O   GLY A  97       8.885   6.147  22.814  1.00 80.00           O
ATOM    110  H   GLY A  97       8.146   8.826  19.713  1.00 80.00           H
ATOM    111  HA2 GLY A  97       7.130   7.614  21.762  1.00 80.00           H
ATOM    112  HA3 GLY A  97       7.427   6.935  20.390  1.00 80.00           H
ATOM    113  N   GLY A  98      10.140   6.707  21.049  1.00 80.00           N
ATOM    114  CA  GLY A  98      11.322   6.042  21.575  1.00 80.00           C
ATOM    115  C   GLY A  98      11.839   6.681  22.835  1.00 80.00           C
ATOM    116  O   GLY A  98      12.253   6.004  23.783  1.00 80.00           O
ATOM    117  H   GLY A  98      10.251   7.086  20.285  1.00 80.00           H
ATOM    118  HA2 GLY A  98      11.113   5.112  21.752  1.00 80.00           H
ATOM    119  HA3 GLY A  98      12.020   6.053  20.902  1.00 80.00           H
ATOM    120  N   GLY A  99      11.875   8.002  22.899  1.00 80.00           N
ATOM    121  CA  GLY A  99      12.354   8.727  24.053  1.00 80.00           C
ATOM    122  C   GLY A  99      11.365   8.788  25.232  1.00 80.00           C
ATOM    123  O   GLY A  99      11.719   9.291  26.291  1.00 80.00           O
ATOM    124  OXT GLY A  99      10.140   8.286  25.083  1.00 80.00           O
ATOM    125  H   GLY A  99      11.615   8.511  22.256  1.00 80.00           H
ATOM    126  HA2 GLY A  99      13.178   8.316  24.359  1.00 80.00           H
ATOM    127  HA3 GLY A  99      12.572   9.633  23.783  1.00 80.00           H
TER
"""

pdb_str_in2 = """
REMARK as as pdb_str_in1 but no altlocs
CRYST1   18.273   20.824   31.291  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  87       7.109  12.926   5.576  1.00 80.00           N
ATOM      2  CA  GLY A  87       7.840  12.041   6.465  1.00 80.00           C
ATOM      3  C   GLY A  87       8.271  12.731   7.773  1.00 80.00           C
ATOM      4  O   GLY A  87       8.615  12.020   8.719  1.00 80.00           O
ATOM         H1  GLY A  87       7.550  12.998   4.795  1.00 80.00           H
ATOM         H2  GLY A  87       6.288  12.590   5.426  1.00 80.00           H
ATOM         H3  GLY A  87       7.036  13.739   5.955  1.00 80.00           H
ATOM      5  HA2 GLY A  87       7.287  11.273   6.677  1.00 80.00           H
ATOM      6  HA3 GLY A  87       8.627  11.706   6.006  1.00 80.00           H
ATOM      7  N   GLY A  88       8.311  14.046   7.824  1.00 80.00           N
ATOM      8  CA  GLY A  88       8.924  14.748   8.963  1.00 80.00           C
ATOM      9  C   GLY A  88       8.145  14.441  10.249  1.00 80.00           C
ATOM     10  O   GLY A  88       8.755  14.178  11.288  1.00 80.00           O
ATOM     11  H   GLY A  88       7.992  14.563   7.215  1.00 80.00           H
ATOM     12  HA2 GLY A  88       9.849  14.473   9.063  1.00 80.00           H
ATOM     13  HA3 GLY A  88       8.929  15.704   8.799  1.00 80.00           H
ATOM     14  N   GLY A  89       6.816  14.533  10.216  1.00 80.00           N
ATOM     15  CA  GLY A  89       6.063  14.339  11.454  1.00 80.00           C
ATOM     16  C   GLY A  89       6.267  12.956  11.997  1.00 80.00           C
ATOM     17  O   GLY A  89       6.484  12.771  13.205  1.00 80.00           O
ATOM     18  H   GLY A  89       6.346  14.700   9.516  1.00 80.00           H
ATOM     19  HA2 GLY A  89       6.344  14.993  12.113  1.00 80.00           H
ATOM     20  HA3 GLY A  89       5.119  14.489  11.289  1.00 80.00           H
ATOM     21  N   GLY A  90       6.201  11.944  11.127  1.00 80.00           N
ATOM     22  CA  GLY A  90       6.380  10.588  11.592  1.00 80.00           C
ATOM     23  C   GLY A  90       7.765  10.329  12.140  1.00 80.00           C
ATOM     24  O   GLY A  90       7.932   9.581  13.086  1.00 80.00           O
ATOM     25  H   GLY A  90       6.056  12.028  10.283  1.00 80.00           H
ATOM     26  HA2 GLY A  90       5.724  10.399  12.281  1.00 80.00           H
ATOM     27  HA3 GLY A  90       6.207   9.975  10.860  1.00 80.00           H
ATOM     28  N   GLY A  91       8.780  10.957  11.543  1.00 80.00           N
ATOM     29  CA  GLY A  91      10.141  10.829  12.058  1.00 80.00           C
ATOM     30  C   GLY A  91      10.230  11.420  13.459  1.00 80.00           C
ATOM     31  O   GLY A  91      10.831  10.840  14.371  1.00 80.00           O
ATOM     32  H   GLY A  91       8.702  11.456  10.847  1.00 80.00           H
ATOM     33  HA2 GLY A  91      10.400   9.895  12.077  1.00 80.00           H
ATOM     34  HA3 GLY A  91      10.762  11.284  11.468  1.00 80.00           H
ATOM     35  N   GLY A  92       9.632  12.593  13.642  1.00 80.00           N
ATOM     36  CA  GLY A  92       9.634  13.216  14.955  1.00 80.00           C
ATOM     37  C   GLY A  92       8.870  12.343  15.950  1.00 80.00           C
ATOM     38  O   GLY A  92       9.337  12.131  17.071  1.00 80.00           O
ATOM     39  H   GLY A  92       9.226  13.037  13.028  1.00 80.00           H
ATOM     40  HA2 GLY A  92      10.546  13.343  15.259  1.00 80.00           H
ATOM     41  HA3 GLY A  92       9.226  14.095  14.905  1.00 80.00           H
ATOM     42  N   GLY A  93       7.711  11.835  15.530  1.00 80.00           N
ATOM     43  CA  GLY A  93       6.915  11.011  16.422  1.00 80.00           C
ATOM     44  C   GLY A  93       7.606   9.711  16.796  1.00 80.00           C
ATOM     45  O   GLY A  93       7.431   9.245  17.928  1.00 80.00           O
ATOM     46  H   GLY A  93       7.377  11.955  14.747  1.00 80.00           H
ATOM     47  HA2 GLY A  93       6.719  11.512  17.229  1.00 80.00           H
ATOM     48  HA3 GLY A  93       6.066  10.811  15.999  1.00 80.00           H
ATOM     70  N   TYR A  94       8.479   9.112  15.882  1.00 80.00           N
ATOM     71  CA  TYR A  94       9.232   7.915  16.248  1.00 80.00           C
ATOM     72  C   TYR A  94      10.258   8.233  17.325  1.00 80.00           C
ATOM     73  O   TYR A  94      10.405   7.506  18.293  1.00 80.00           O
ATOM     74  CB  TYR A  94       9.921   7.315  15.022  1.00 30.00           C
ATOM     75  CG  TYR A  94      10.801   6.124  15.334  1.00 30.00           C
ATOM     76  CD1 TYR A  94      10.247   4.873  15.580  1.00 30.00           C
ATOM     77  CD2 TYR A  94      12.182   6.252  15.381  1.00 30.00           C
ATOM     78  CE1 TYR A  94      11.047   3.784  15.865  1.00 30.00           C
ATOM     79  CE2 TYR A  94      12.990   5.166  15.666  1.00 30.00           C
ATOM     80  CZ  TYR A  94      12.418   3.936  15.907  1.00 30.00           C
ATOM     81  OH  TYR A  94      13.216   2.852  16.193  1.00 30.00           O
ATOM     82  H   TYR A  94       8.620   9.399  15.084  1.00 80.00           H
ATOM     83  HA  TYR A  94       8.607   7.263  16.601  1.00 80.00           H
ATOM     84  HB2 TYR A  94       9.244   7.047  14.381  1.00 30.00           H
ATOM     85  HB3 TYR A  94      10.459   8.001  14.596  1.00 30.00           H
ATOM     86  HD1 TYR A  94       9.323   4.768  15.552  1.00 30.00           H
ATOM     87  HD2 TYR A  94      12.570   7.081  15.218  1.00 30.00           H
ATOM     88  HE1 TYR A  94      10.664   2.952  16.028  1.00 30.00           H
ATOM     89  HE2 TYR A  94      13.914   5.266  15.695  1.00 30.00           H
ATOM     90  HH  TYR A  94      14.022   3.087  16.184  1.00 30.00           H
ATOM     91  N   GLY A  95      10.862   9.366  17.148  1.00 80.00           N
ATOM     92  CA  GLY A  95      11.781   9.803  18.178  1.00 80.00           C
ATOM     93  C   GLY A  95      11.049  10.057  19.501  1.00 80.00           C
ATOM     94  O   GLY A  95      11.541   9.637  20.558  1.00 80.00           O
ATOM     95  HA2 GLY A  95      12.468   9.131  18.309  1.00 80.00           H
ATOM     96  HA3 GLY A  95      12.229  10.614  17.892  1.00 80.00           H
ATOM     98  H   GLY A  95      10.769   9.889  16.472  1.00 80.00           H
ATOM     99  N   GLY A  96       9.911  10.747  19.468  1.00 80.00           N
ATOM    100  CA  GLY A  96       9.163  10.978  20.691  1.00 80.00           C
ATOM    101  C   GLY A  96       8.756   9.674  21.361  1.00 80.00           C
ATOM    102  O   GLY A  96       8.867   9.509  22.570  1.00 80.00           O
ATOM    103  H   GLY A  96       9.562  11.084  18.757  1.00 80.00           H
ATOM    104  HA2 GLY A  96       9.701  11.503  21.305  1.00 80.00           H
ATOM    105  HA3 GLY A  96       8.370  11.500  20.491  1.00 80.00           H
ATOM    106  N   GLY A  97       8.248   8.738  20.562  1.00 80.00           N
ATOM    107  CA  GLY A  97       7.821   7.470  21.097  1.00 80.00           C
ATOM    108  C   GLY A  97       8.991   6.708  21.728  1.00 80.00           C
ATOM    109  O   GLY A  97       8.885   6.147  22.814  1.00 80.00           O
ATOM    110  H   GLY A  97       8.146   8.826  19.713  1.00 80.00           H
ATOM    111  HA2 GLY A  97       7.130   7.614  21.762  1.00 80.00           H
ATOM    112  HA3 GLY A  97       7.427   6.935  20.390  1.00 80.00           H
ATOM    113  N   GLY A  98      10.140   6.707  21.049  1.00 80.00           N
ATOM    114  CA  GLY A  98      11.322   6.042  21.575  1.00 80.00           C
ATOM    115  C   GLY A  98      11.839   6.681  22.835  1.00 80.00           C
ATOM    116  O   GLY A  98      12.253   6.004  23.783  1.00 80.00           O
ATOM    117  H   GLY A  98      10.251   7.086  20.285  1.00 80.00           H
ATOM    118  HA2 GLY A  98      11.113   5.112  21.752  1.00 80.00           H
ATOM    119  HA3 GLY A  98      12.020   6.053  20.902  1.00 80.00           H
ATOM    120  N   GLY A  99      11.875   8.002  22.899  1.00 80.00           N
ATOM    121  CA  GLY A  99      12.354   8.727  24.053  1.00 80.00           C
ATOM    122  C   GLY A  99      11.365   8.788  25.232  1.00 80.00           C
ATOM    123  O   GLY A  99      11.719   9.291  26.291  1.00 80.00           O
ATOM    124  OXT GLY A  99      10.140   8.286  25.083  1.00 80.00           O
ATOM    125  H   GLY A  99      11.615   8.511  22.256  1.00 80.00           H
ATOM    126  HA2 GLY A  99      13.178   8.316  24.359  1.00 80.00           H
ATOM    127  HA3 GLY A  99      12.572   9.633  23.783  1.00 80.00           H
TER
"""

def run(prefix):
  """
  Exercise gradients match: using clustering vs not using clustering.
  Altlocs.
  
  XXX TEST FAILS DUE TO MOPAC ERROR (no results check)
  XXX Same test with xtb (properly set) is OK
  
  """
  for i, pdb_str_in in enumerate([pdb_str_in1, pdb_str_in2]):
    if(i==0): print("Altlocs present", "-"*30)
    else:     print("No altlocs", "-"*30)
    pdb_in = "%s.pdb"%prefix
    open(pdb_in, "w").write(pdb_str_in)
    #
    # for fast_interaction in [True, False]:
    for fast_interaction in [True]:
      print("fast_interaction:", fast_interaction)
      for restraints in ["cctbx", "qm"]:
        print("  restraints:", restraints)
        for two_buffers in [False, True]:
          print("    two_buffers=", two_buffers)
          for clustering in ["true", "false"]:
            print("      clustering=", clustering)

            if(not clustering): expansion=True
            else:               expansion=False

            cmd = " ".join([
              "qr.refine",
              pdb_in,

              "expansion=%s"%str(expansion),

              "mode=opt",
              "altloc_method=subtract",
              "fast_interaction=%s"%fast_interaction,
              "stpmax=0.2",
              "gradient_only=true",
              "clustering=%s"%clustering,
              "dump_gradients=cluster_%s.pkl"%clustering,
              "restraints=%s"%restraints,
              "quantum.engine_name=mopac",
              "number_of_micro_cycles=1",
              "max_iterations_refine=5",
              "two_buffers=%s"%str(two_buffers),
              "> %s.log"%prefix])
            #print cmd
            assert easy_run.call(cmd)==0
          g1 = easy_pickle.load("cluster_false.pkl")
          g2 = easy_pickle.load("cluster_true.pkl")
          g1 = g1.as_double()
          g2 = g2.as_double()
          diff = flex.abs(g1-g2)
          print("        min/max/mean of (gradient1 - gradient2):", \
              diff.min_max_mean().as_tuple())
          os.remove("cluster_false.pkl")
          os.remove("cluster_true.pkl")

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=True)
