from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from qrefine import qr
from qrefine.tests.unit import run_qrefine_in, run_fmodel

pdb_str_good = """
CRYST1   18.273   20.824   31.291  90.00  90.00  90.00 P 1
SCALE1      0.054726  0.000000  0.000000        0.00000
SCALE2      0.000000  0.048022  0.000000        0.00000
SCALE3      0.000000  0.000000  0.031958        0.00000
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
ATOM     75  CG BTYR A  94       9.047   6.368  14.228  0.50 30.00           C
ATOM     76  CD1BTYR A  94       8.057   6.853  13.380  0.50 30.00           C
ATOM     77  CD2BTYR A  94       9.212   4.994  14.328  0.50 30.00           C
ATOM     78  CE1BTYR A  94       7.258   5.992  12.655  0.50 30.00           C
ATOM     79  CE2BTYR A  94       8.415   4.125  13.604  0.50 30.00           C
ATOM     80  CZ BTYR A  94       7.441   4.629  12.771  0.50 30.00           C
ATOM     81  OH BTYR A  94       6.643   3.770  12.050  0.50 30.00           O
ATOM     82  H  BTYR A  94       8.620   9.399  15.084  0.50 80.00           H
ATOM     83  HA BTYR A  94       8.607   7.263  16.601  0.50 80.00           H
ATOM     84  HB2BTYR A  94      10.213   8.035  14.441  0.50 30.00           H
ATOM     85  HB3BTYR A  94      10.718   6.842  15.309  0.50 30.00           H
ATOM     86  HD1BTYR A  94       7.931   7.771  13.300  0.50 30.00           H
ATOM     87  HD2BTYR A  94       9.869   4.652  14.891  0.50 30.00           H
ATOM     88  HE1BTYR A  94       6.600   6.329  12.091  0.50 30.00           H
ATOM     89  HE2BTYR A  94       8.536   3.206  13.679  0.50 30.00           H
ATOM     90  HH BTYR A  94       6.864   2.977  12.216  0.50 30.00           H
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

pdb_str_poor = """
CRYST1   18.273   20.824   31.291  90.00  90.00  90.00 P 1
SCALE1      0.054726  0.000000  0.000000        0.00000
SCALE2      0.000000  0.048022  0.000000        0.00000
SCALE3      0.000000  0.000000  0.031958        0.00000
ATOM      1  N   GLY A  87       7.070  12.940   5.537  1.00 80.00           N
ATOM      2  CA  GLY A  87       7.820  12.016   6.447  1.00 80.00           C
ATOM      3  C   GLY A  87       8.313  12.732   7.747  1.00 80.00           C
ATOM      4  O   GLY A  87       8.517  11.999   8.604  1.00 80.00           O
ATOM      5  H1  GLY A  87       7.681  13.079   4.693  1.00 80.00           H
ATOM      6  H2  GLY A  87       6.128  12.494   5.303  1.00 80.00           H
ATOM      7  H3  GLY A  87       7.036  13.862   6.147  1.00 80.00           H
ATOM      8  HA2 GLY A  87       7.450  11.079   6.816  1.00 80.00           H
ATOM      9  HA3 GLY A  87       8.613  11.855   6.070  1.00 80.00           H
ATOM     10  N   GLY A  88       8.146  13.910   7.994  1.00 80.00           N
ATOM     11  CA  GLY A  88       9.099  14.930   9.023  1.00 80.00           C
ATOM     12  C   GLY A  88       8.265  14.454  10.301  1.00 80.00           C
ATOM     13  O   GLY A  88       8.603  14.264  11.332  1.00 80.00           O
ATOM     14  H   GLY A  88       8.032  14.520   7.358  1.00 80.00           H
ATOM     15  HA2 GLY A  88       9.971  14.369   9.005  1.00 80.00           H
ATOM     16  HA3 GLY A  88       8.899  15.829   8.758  1.00 80.00           H
ATOM     17  N   GLY A  89       6.810  14.537  10.189  1.00 80.00           N
ATOM     18  CA  GLY A  89       6.177  14.367  11.483  1.00 80.00           C
ATOM     19  C   GLY A  89       6.150  13.113  11.953  1.00 80.00           C
ATOM     20  O   GLY A  89       6.325  12.670  13.259  1.00 80.00           O
ATOM     21  H   GLY A  89       6.189  14.623   9.568  1.00 80.00           H
ATOM     22  HA2 GLY A  89       6.211  14.972  12.177  1.00 80.00           H
ATOM     23  HA3 GLY A  89       5.032  14.502  11.130  1.00 80.00           H
ATOM     24  N   GLY A  90       6.020  12.063  11.289  1.00 80.00           N
ATOM     25  CA  GLY A  90       6.361  10.754  11.471  1.00 80.00           C
ATOM     26  C   GLY A  90       7.873  10.384  12.173  1.00 80.00           C
ATOM     27  O   GLY A  90       8.078   9.702  13.263  1.00 80.00           O
ATOM     28  H   GLY A  90       6.149  11.951  10.278  1.00 80.00           H
ATOM     29  HA2 GLY A  90       5.704  10.304  12.403  1.00 80.00           H
ATOM     30  HA3 GLY A  90       6.101   9.973  10.748  1.00 80.00           H
ATOM     31  N   GLY A  91       8.901  10.912  11.462  1.00 80.00           N
ATOM     32  CA  GLY A  91      10.296  10.883  11.915  1.00 80.00           C
ATOM     33  C   GLY A  91      10.188  11.529  13.468  1.00 80.00           C
ATOM     34  O   GLY A  91      10.885  10.694  14.466  1.00 80.00           O
ATOM     35  H   GLY A  91       8.512  11.267  10.823  1.00 80.00           H
ATOM     36  HA2 GLY A  91      10.489  10.050  12.226  1.00 80.00           H
ATOM     37  HA3 GLY A  91      10.768  11.136  11.386  1.00 80.00           H
ATOM     38  N   GLY A  92       9.544  12.486  13.759  1.00 80.00           N
ATOM     39  CA  GLY A  92       9.585  13.049  15.078  1.00 80.00           C
ATOM     40  C   GLY A  92       8.732  12.440  15.963  1.00 80.00           C
ATOM     41  O   GLY A  92       9.276  12.281  16.919  1.00 80.00           O
ATOM     42  H   GLY A  92       9.287  13.030  12.832  1.00 80.00           H
ATOM     43  HA2 GLY A  92      10.586  13.291  15.293  1.00 80.00           H
ATOM     44  HA3 GLY A  92       9.047  14.212  14.955  1.00 80.00           H
ATOM     45  N   GLY A  93       7.752  11.794  15.512  1.00 80.00           N
ATOM     46  CA  GLY A  93       7.008  11.103  16.375  1.00 80.00           C
ATOM     47  C   GLY A  93       7.699   9.588  16.888  1.00 80.00           C
ATOM     48  O   GLY A  93       7.498   9.106  17.797  1.00 80.00           O
ATOM     49  H   GLY A  93       7.389  11.814  14.738  1.00 80.00           H
ATOM     50  HA2 GLY A  93       6.854  11.313  17.106  1.00 80.00           H
ATOM     51  HA3 GLY A  93       5.981  10.755  16.093  1.00 80.00           H
ATOM     52  N  ATYR A  94       8.512   9.270  15.798  0.50 80.00           N
ATOM     53  CA ATYR A  94       9.112   7.782  16.119  0.50 80.00           C
ATOM     54  C  ATYR A  94       9.959   8.302  17.167  0.50 80.00           C
ATOM     55  O  ATYR A  94      10.456   7.375  18.149  0.50 80.00           O
ATOM     56  CB ATYR A  94       9.821   7.232  14.848  0.50 30.00           C
ATOM     57  CG ATYR A  94       9.040   6.599  14.173  0.50 30.00           C
ATOM     58  CD1ATYR A  94       8.186   5.352  13.892  0.50 30.00           C
ATOM     59  CD2ATYR A  94       8.446   7.625  12.930  0.50 30.00           C
ATOM     60  CE1ATYR A  94       7.653   4.927  12.909  0.50 30.00           C
ATOM     61  CE2ATYR A  94       7.534   7.009  12.127  0.50 30.00           C
ATOM     62  CZ ATYR A  94       7.134   6.005  12.236  0.50 30.00           C
ATOM     63  OH ATYR A  94       6.420   5.413  10.932  0.50 30.00           O
ATOM     64  H  ATYR A  94       8.674   9.551  15.139  0.50 80.00           H
ATOM     65  HA ATYR A  94       8.626   7.114  16.670  0.50 80.00           H
ATOM     66  HB2ATYR A  94      10.383   8.116  14.516  0.50 30.00           H
ATOM     67  HB3ATYR A  94      10.292   6.573  15.174  0.50 30.00           H
ATOM     68  HD1ATYR A  94       8.596   4.753  14.680  0.50 30.00           H
ATOM     69  HD2ATYR A  94       8.908   8.313  12.754  0.50 30.00           H
ATOM     70  HE1ATYR A  94       7.048   4.163  12.995  0.50 30.00           H
ATOM     71  HE2ATYR A  94       7.279   7.747  11.344  0.50 30.00           H
ATOM     72  HH ATYR A  94       5.940   6.072  10.633  0.50 30.00           H
ATOM     73  N  BTYR A  94       8.361   8.983  15.778  0.50 80.00           N
ATOM     74  CA BTYR A  94       9.138   8.033  16.094  0.50 80.00           C
ATOM     75  C  BTYR A  94      10.234   8.034  17.132  0.50 80.00           C
ATOM     76  O  BTYR A  94      10.399   7.457  18.179  0.50 80.00           O
ATOM     77  CB BTYR A  94      10.088   7.277  14.837  0.50 30.00           C
ATOM     78  CG BTYR A  94       8.983   6.508  14.074  0.50 30.00           C
ATOM     79  CD1BTYR A  94       7.987   6.721  13.321  0.50 30.00           C
ATOM     80  CD2BTYR A  94       9.378   4.804  14.229  0.50 30.00           C
ATOM     81  CE1BTYR A  94       7.082   6.179  12.756  0.50 30.00           C
ATOM     82  CE2BTYR A  94       8.383   4.090  13.415  0.50 30.00           C
ATOM     83  CZ BTYR A  94       7.507   4.484  12.736  0.50 30.00           C
ATOM     84  OH BTYR A  94       6.725   3.902  11.886  0.50 30.00           O
ATOM     85  H  BTYR A  94       8.687   9.407  15.229  0.50 80.00           H
ATOM     86  HA BTYR A  94       8.790   7.298  16.404  0.50 80.00           H
ATOM     87  HB2BTYR A  94      10.018   8.130  14.535  0.50 30.00           H
ATOM     88  HB3BTYR A  94      10.792   6.947  15.390  0.50 30.00           H
ATOM     89  HD1BTYR A  94       7.913   7.673  13.348  0.50 30.00           H
ATOM     90  HD2BTYR A  94       9.698   4.755  15.027  0.50 30.00           H
ATOM     91  HE1BTYR A  94       6.557   6.500  12.023  0.50 30.00           H
ATOM     92  HE2BTYR A  94       8.565   3.196  13.574  0.50 30.00           H
ATOM     93  HH BTYR A  94       6.971   3.003  12.199  0.50 30.00           H
ATOM     94  N   GLY A  95      10.804   9.377  17.240  1.00 80.00           N
ATOM     95  CA  GLY A  95      11.726   9.813  18.235  1.00 80.00           C
ATOM     96  C   GLY A  95      11.082  10.096  19.404  1.00 80.00           C
ATOM     97  O   GLY A  95      11.694   9.459  20.382  1.00 80.00           O
ATOM     98  HA2 GLY A  95      12.277   9.205  18.203  1.00 80.00           H
ATOM     99  HA3 GLY A  95      12.286  10.505  17.915  1.00 80.00           H
ATOM    100  H  AGLY A  95      10.662   9.690  16.471  0.50 80.00           H
ATOM    101  H  BGLY A  95      10.717   9.816  16.549  0.50 80.00           H
ATOM    102  N   GLY A  96      10.071  10.597  19.566  1.00 80.00           N
ATOM    103  CA  GLY A  96       9.314  10.797  20.815  1.00 80.00           C
ATOM    104  C   GLY A  96       8.588   9.849  21.555  1.00 80.00           C
ATOM    105  O   GLY A  96       8.823   9.454  22.401  1.00 80.00           O
ATOM    106  H   GLY A  96       9.671  11.146  18.712  1.00 80.00           H
ATOM    107  HA2 GLY A  96       9.731  11.502  21.155  1.00 80.00           H
ATOM    108  HA3 GLY A  96       8.189  11.500  20.688  1.00 80.00           H
ATOM    109  N   GLY A  97       8.096   8.821  20.365  1.00 80.00           N
ATOM    110  CA  GLY A  97       7.952   7.496  21.140  1.00 80.00           C
ATOM    111  C   GLY A  97       9.103   6.540  21.629  1.00 80.00           C
ATOM    112  O   GLY A  97       8.956   6.240  22.790  1.00 80.00           O
ATOM    113  H   GLY A  97       8.010   8.865  19.609  1.00 80.00           H
ATOM    114  HA2 GLY A  97       7.027   7.484  21.673  1.00 80.00           H
ATOM    115  HA3 GLY A  97       7.581   6.740  20.200  1.00 80.00           H
ATOM    116  N   GLY A  98       9.969   6.582  21.108  1.00 80.00           N
ATOM    117  CA  GLY A  98      11.428   5.922  21.457  1.00 80.00           C
ATOM    118  C   GLY A  98      11.910   6.578  22.661  1.00 80.00           C
ATOM    119  O   GLY A  98      12.449   5.990  23.911  1.00 80.00           O
ATOM    120  H   GLY A  98      10.390   7.019  20.436  1.00 80.00           H
ATOM    121  HA2 GLY A  98      11.147   5.108  21.786  1.00 80.00           H
ATOM    122  HA3 GLY A  98      12.070   6.155  20.810  1.00 80.00           H
ATOM    123  N   GLY A  99      11.999   7.920  22.991  1.00 80.00           N
ATOM    124  CA  GLY A  99      12.526   8.869  23.893  1.00 80.00           C
ATOM    125  C   GLY A  99      11.308   8.718  25.274  1.00 80.00           C
ATOM    126  O   GLY A  99      11.594   9.320  26.400  1.00 80.00           O
ATOM    127  OXT GLY A  99      10.098   8.159  25.260  1.00 80.00           O
ATOM    128  H   GLY A  99      11.784   8.666  22.189  1.00 80.00           H
ATOM    129  HA2 GLY A  99      13.014   8.377  24.215  1.00 80.00           H
ATOM    130  HA3 GLY A  99      12.731   9.761  23.809  1.00 80.00           H
TER
END
"""

# Keeps each qr.refine run short: the cost is dominated by the per-step model
# statistics (rama_z, clashscore), not by the number of function evaluations.
refine_fast = ["number_of_weight_search_cycles=1", "number_of_refine_cycles=1"]

def get_model(pdb_file):
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(file_name = pdb_file), log = null_out())
  model.setup_scattering_dictionaries(
    scattering_table = qr.get_default_params().scattering_table)
  return model

def r_work_of(pdb_file, r_fm):
  fmodel = mmtbx.f_model.manager(
    f_obs          = r_fm.f_obs,
    r_free_flags   = r_fm.flags,
    xray_structure = get_model(pdb_file).get_xray_structure())
  fmodel.update_all_scales()
  return fmodel.r_work()

def write_inputs(prefix):
  """
  Write good and poor models into the current folder; return absolute paths.
  """
  pdb_good = os.path.abspath("%s_good.pdb"%prefix)
  with open(pdb_good, "w") as fo:
    fo.write(pdb_str_good)
  pdb_poor = os.path.abspath("%s_poor.pdb"%prefix)
  with open(pdb_poor, "w") as fo:
    fo.write(pdb_str_poor)
  return pdb_good, pdb_poor

def check_refined(r_qr, pdb_poor, r_fm, r_start):
  """
  Altlocs and occupancies preserved, R factor dropped.
  """
  poor_occupancies = get_model(pdb_poor).get_hierarchy().atoms().extract_occ()
  h = r_qr.model.get_hierarchy()
  assert sorted(h.altloc_indices().keys()) == ["", "A", "B"]
  assert approx_equal(h.atoms().extract_occ(), poor_occupancies)
  r_end = r_work_of(r_qr.pdb, r_fm)
  print("r_start: %6.4f r_end: %6.4f" % (r_start, r_end))
  assert r_end < 0.5 * r_start, (r_start, r_end)
  return r_end

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine with alternative conformations, cctbx restraints,
  altlocs.handling=conformers; native with QM is refused. The other routes
  are in tst_39 (expansion, exclude) and tst_43 (clustering, native).
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb_good, pdb_poor = write_inputs(prefix)
  r_fm = run_fmodel(prefix=prefix, args=[pdb_good, "high_res=2"])
  mtz = os.path.abspath(r_fm.mtz)
  r_start = r_work_of(pdb_poor, r_fm)
  r_qr = run_qrefine_in(folder="conformers", prefix=prefix+"_conformers",
    args=[pdb_poor, mtz, "restraints=cctbx", "clustering=False",
      "altlocs.handling=conformers"] + refine_fast)
  check_refined(r_qr=r_qr, pdb_poor=pdb_poor, r_fm=r_fm, r_start=r_start)
  assert r_qr.std_out_text.find(
    "alternative conformations require minimizer=lbfgs") > -1
  r_qr = run_qrefine_in(folder="sorry", prefix=prefix+"_sorry",
    sorry_expected=True, args=[pdb_poor, mtz, "restraints=qm",
      "quantum.engine_name=mopac", "altlocs.handling=native", "dry_run=True"])
  assert r_qr.std_out_text.find(
    "altlocs.handling=native needs restraints=cctbx") > -1, r_qr.std_out_text

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
