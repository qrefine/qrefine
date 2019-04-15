from __future__ import division
import os
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import easy_pickle
import time, sys
import run_tests
import libtbx.load_env
from libtbx.test_utils import approx_equal
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

pdb_str_in = """
CRYST1   70.364   66.358   73.268  90.00  90.00  90.00 P 1
SCALE1      0.014212  0.000000  0.000000        0.00000
SCALE2      0.000000  0.015070  0.000000        0.00000
SCALE3      0.000000  0.000000  0.013649        0.00000
ATOM      1  N   SER A   1      24.469  33.098  30.908  1.00 88.66           N
ATOM      2  CA  SER A   1      25.006  32.552  32.148  1.00 91.92           C
ATOM      3  C   SER A   1      25.343  33.653  33.133  1.00 91.25           C
ATOM      4  O   SER A   1      25.876  34.689  32.740  1.00 95.01           O
ATOM      5  CB  SER A   1      26.261  31.704  31.877  1.00 89.62           C
ATOM      6  OG  SER A   1      25.914  30.347  31.651  1.00 94.99           O
ATOM      7  H1  SER A   1      24.691  32.567  30.229  1.00 88.66           H
ATOM      8  H2  SER A   1      23.582  33.151  30.966  1.00 88.66           H
ATOM      9  H3  SER A   1      24.806  33.910  30.771  1.00 88.66           H
ATOM     10  HA  SER A   1      24.319  31.987  32.535  1.00 91.92           H
ATOM     11  HB2 SER A   1      26.733  32.054  31.105  1.00 89.62           H
ATOM     12  HB3 SER A   1      26.867  31.766  32.632  1.00 89.62           H
ATOM     13  HG  SER A   1      26.610  29.900  31.504  1.00 94.99           H
ATOM     14  N   VAL A   2      25.042  33.416  34.412  1.00 84.96           N
ATOM     15  CA  VAL A   2      25.447  34.295  35.505  1.00 87.97           C
ATOM     16  C   VAL A   2      26.185  33.479  36.560  1.00 84.17           C
ATOM     17  O   VAL A   2      25.784  32.350  36.874  1.00 81.72           O
ATOM     18  CB  VAL A   2      24.239  35.047  36.116  1.00 89.60           C
ATOM     19  CG1 VAL A   2      23.109  34.090  36.377  1.00 91.48           C
ATOM     20  CG2 VAL A   2      24.635  35.758  37.418  1.00 85.72           C
ATOM     21  H   VAL A   2      24.591  32.731  34.669  0.00 84.96           H
ATOM     22  HA  VAL A   2      26.044  34.974  35.153  0.00 87.97           H
ATOM     23  HB  VAL A   2      23.948  35.717  35.478  0.00 89.60           H
ATOM     24 HG11 VAL A   2      22.359  34.571  36.759  0.00 91.48           H
ATOM     25 HG12 VAL A   2      22.836  33.676  35.543  0.00 91.48           H
ATOM     26 HG13 VAL A   2      23.402  33.404  36.997  0.00 91.48           H
ATOM     27 HG21 VAL A   2      23.864  36.221  37.781  0.00 85.72           H
ATOM     28 HG22 VAL A   2      24.951  35.104  38.061  0.00 85.72           H
ATOM     29 HG23 VAL A   2      25.340  36.399  37.236  0.00 85.72           H
ATOM     30  N   CYS A   3      27.274  34.055  37.102  1.00 76.61           N
ATOM     31  CA  CYS A   3      28.072  33.388  38.109  1.00 78.13           C
ATOM     32  C   CYS A   3      27.919  34.113  39.444  1.00 83.42           C
ATOM     33  O   CYS A   3      26.986  34.896  39.643  1.00 87.64           O
ATOM     34  CB  CYS A   3      29.522  33.281  37.635  1.00 78.38           C
ATOM     35  SG  CYS A   3      30.585  34.726  37.772  1.00 63.83           S
ATOM     36  H   CYS A   3      27.557  34.838  36.888  0.00 76.61           H
ATOM     37  HA  CYS A   3      27.761  32.480  38.247  0.00 78.13           H
ATOM     38  HB2 CYS A   3      29.941  32.560  38.130  0.00 78.38           H
ATOM     39  HB3 CYS A   3      29.508  33.014  36.703  0.00 78.38           H
ATOM     40  N   ALA A   4      28.823  33.829  40.379  1.00 76.53           N
ATOM     41  CA  ALA A   4      28.694  34.351  41.728  1.00 78.90           C
ATOM     42  C   ALA A   4      29.464  35.642  41.972  1.00 92.51           C
ATOM     43  O   ALA A   4      29.208  36.308  42.982  1.00102.42           O
ATOM     44  CB  ALA A   4      29.135  33.286  42.736  1.00 86.89           C
ATOM     45  H   ALA A   4      29.516  33.336  40.249  0.00 76.53           H
ATOM     46  HA  ALA A   4      27.757  34.571  41.845  0.00 78.90           H
ATOM     47  HB1 ALA A   4      29.048  33.637  43.636  0.00 86.89           H
ATOM     48  HB2 ALA A   4      28.576  32.499  42.641  0.00 86.89           H
ATOM     49  HB3 ALA A   4      30.060  33.047  42.570  0.00 86.89           H
ATOM     50  N   ALA A   5      30.418  35.988  41.105  1.00 97.63           N
ATOM     51  CA  ALA A   5      31.097  37.276  41.179  1.00 98.34           C
ATOM     52  C   ALA A   5      30.084  38.420  41.164  1.00100.45           C
ATOM     53  O   ALA A   5      29.041  38.334  40.505  1.00 89.07           O
ATOM     54  CB  ALA A   5      32.069  37.411  40.008  1.00 85.11           C
ATOM     55  H   ALA A   5      30.686  35.483  40.462  0.00 97.63           H
ATOM     56  HA  ALA A   5      31.593  37.323  42.011  0.00 98.34           H
ATOM     57  HB1 ALA A   5      32.520  38.268  40.059  0.00 85.11           H
ATOM     58  HB2 ALA A   5      32.725  36.698  40.048  0.00 85.11           H
ATOM     59  HB3 ALA A   5      31.580  37.352  39.173  0.00 85.11           H
ATOM     60  N   ALA A   6      30.410  39.501  41.893  1.00 99.80           N
ATOM     61  CA  ALA A   6      29.436  40.572  42.125  1.00103.38           C
ATOM     62  C   ALA A   6      29.089  41.290  40.832  1.00105.41           C
ATOM     63  O   ALA A   6      27.920  41.618  40.585  1.00 98.78           O
ATOM     64  CB  ALA A   6      29.964  41.569  43.165  1.00114.34           C
ATOM     65  H   ALA A   6      31.180  39.628  42.255  0.00 99.80           H
ATOM     66  HA  ALA A   6      28.626  40.165  42.469  0.00103.38           H
ATOM     67  HB1 ALA A   6      29.307  42.269  43.305  1.00114.34           H
ATOM     68  HB2 ALA A   6      30.905  41.731  42.994  0.00114.34           H
ATOM     69  HB3 ALA A   6      29.504  42.415  43.050  0.00114.34           H
ATOM     70  N   ASN A   7      30.084  41.535  39.989  1.00104.87           N
ATOM     71  CA  ASN A   7      29.841  42.020  38.640  1.00106.00           C
ATOM     72  C   ASN A   7      30.121  40.895  37.646  1.00102.72           C
ATOM     73  O   ASN A   7      31.156  40.852  36.959  1.00 99.10           O
ATOM     74  CB  ASN A   7      30.668  43.274  38.359  1.00112.46           C
ATOM     75  CG  ASN A   7      29.845  44.553  38.481  1.00116.30           C
ATOM     76  OD1 ASN A   7      28.816  44.578  39.158  1.00116.27           O
ATOM     77  ND2 ASN A   7      30.291  45.617  37.813  1.00122.68           N
ATOM     78  H   ASN A   7      30.915  41.425  40.182  0.00104.87           H
ATOM     79  HA  ASN A   7      28.912  42.280  38.543  0.00106.00           H
ATOM     80  HB2 ASN A   7      31.414  43.312  38.978  0.00112.46           H
ATOM     81  HB3 ASN A   7      31.043  43.218  37.466  0.00112.46           H
ATOM     82 HD21 ASN A   7      29.856  46.358  37.845  0.00122.68           H
ATOM     83 HD22 ASN A   7      31.014  45.562  37.350  0.00122.68           H
ATOM     84  N   CYS A   8      29.186  39.944  37.606  1.00 96.46           N
ATOM     85  CA  CYS A   8      29.223  38.939  36.561  1.00 84.79           C
ATOM     86  C   CYS A   8      28.879  39.632  35.275  1.00 88.83           C
ATOM     87  O   CYS A   8      27.815  40.244  35.160  1.00 93.29           O
ATOM     88  CB  CYS A   8      28.254  37.786  36.819  1.00 85.00           C
ATOM     89  SG  CYS A   8      28.015  36.669  35.373  1.00 78.50           S
ATOM     90  H   CYS A   8      28.537  39.869  38.165  0.00 96.46           H
ATOM     91  HA  CYS A   8      30.107  38.541  36.529  0.00 84.79           H
ATOM     92  HB2 CYS A   8      28.580  37.264  37.569  0.00 85.00           H
ATOM     93  HB3 CYS A   8      27.394  38.151  37.080  0.00 85.00           H
ATOM     94  N   GLN A   9      29.798  39.547  34.323  1.00 87.34           N
ATOM     95  CA  GLN A   9      29.725  40.249  33.057  1.00 89.41           C
ATOM     96  C   GLN A   9      28.900  39.482  32.039  1.00 90.15           C
ATOM     97  O   GLN A   9      29.121  39.591  30.833  1.00 97.42           O
ATOM     98  CB  GLN A   9      31.148  40.497  32.602  1.00 92.86           C
ATOM     99  CG  GLN A   9      31.986  40.818  33.822  1.00 99.03           C
ATOM    100  CD  GLN A   9      33.456  40.784  33.534  1.00101.99           C
ATOM    101  OE1 GLN A   9      33.870  40.259  32.504  1.00 91.19           O
ATOM    102  NE2 GLN A   9      34.268  41.312  34.459  1.00102.88           N
ATOM    103  H   GLN A   9      30.504  39.062  34.401  0.00 87.34           H
ATOM    104  HA  GLN A   9      29.266  41.098  33.157  0.00 89.41           H
ATOM    105  HB2 GLN A   9      31.499  39.715  32.147  0.00 92.86           H
ATOM    106  HB3 GLN A   9      31.178  41.231  31.969  0.00 92.86           H
ATOM    107  HG2 GLN A   9      31.745  41.697  34.154  0.00 99.03           H
ATOM    108  HG3 GLN A   9      31.783  40.183  34.526  0.00 99.03           H
ATOM    109 HE21 GLN A   9      33.937  41.671  35.167  0.00102.88           H
ATOM    110 HE22 GLN A   9      35.120  41.292  34.345  0.00102.88           H
ATOM    111  N   ARG A  10      27.944  38.730  32.541  1.00 85.67           N
ATOM    112  CA  ARG A  10      27.035  37.835  31.849  1.00 92.75           C
ATOM    113  C   ARG A  10      27.612  37.294  30.548  1.00 94.63           C
ATOM    114  O   ARG A  10      27.067  37.591  29.482  1.00101.80           O
ATOM    115  CB  ARG A  10      25.714  38.562  31.592  1.00105.14           C
ATOM    116  CG  ARG A  10      24.735  38.497  32.783  1.00106.96           C
ATOM    117  CD  ARG A  10      23.401  39.137  32.456  1.00105.52           C
ATOM    118  NE  ARG A  10      22.308  38.506  33.192  1.00118.28           N
ATOM    119  CZ  ARG A  10      21.496  37.576  32.680  1.00123.87           C
ATOM    120  NH1 ARG A  10      21.651  37.155  31.417  1.00111.56           N
ATOM    121  NH2 ARG A  10      20.525  37.059  33.436  1.00118.06           N
ATOM    122  H   ARG A  10      27.793  38.729  33.388  0.00 85.67           H
ATOM    123  HA  ARG A  10      26.887  37.065  32.420  0.00 92.75           H
ATOM    124  HB2 ARG A  10      25.899  39.491  31.385  0.00105.14           H
ATOM    125  HB3 ARG A  10      25.288  38.177  30.810  0.00105.14           H
ATOM    126  HG2 ARG A  10      24.595  37.571  33.036  0.00106.96           H
ATOM    127  HG3 ARG A  10      25.128  38.944  33.549  0.00106.96           H
ATOM    128  HD2 ARG A  10      23.433  40.082  32.671  0.00105.52           H
ATOM    129  HD3 ARG A  10      23.233  39.068  31.503  0.00105.52           H
ATOM    130  HE  ARG A  10      22.179  38.749  34.007  0.00118.28           H
ATOM    131 HH11 ARG A  10      22.277  37.484  30.928  0.00111.56           H
ATOM    132 HH12 ARG A  10      21.124  36.556  31.096  0.00111.56           H
ATOM    133 HH21 ARG A  10      20.425  37.325  34.248  0.00118.06           H
ATOM    134 HH22 ARG A  10      20.000  36.460  33.111  0.00118.06           H
ATOM    135  N   PRO A  11      28.700  36.510  30.587  1.00 97.76           N
ATOM    136  CA  PRO A  11      29.236  35.939  29.345  1.00 94.92           C
ATOM    137  C   PRO A  11      28.159  35.169  28.610  1.00100.25           C
ATOM    138  O   PRO A  11      27.426  34.364  29.196  1.00100.17           O
ATOM    139  CB  PRO A  11      30.356  35.012  29.826  1.00 87.60           C
ATOM    140  CG  PRO A  11      30.720  35.506  31.124  1.00 88.53           C
ATOM    141  CD  PRO A  11      29.430  35.986  31.749  1.00 95.29           C
ATOM         OXT PRO A  11      27.998  35.338  27.401  1.00100.25           O
ATOM    142  HA  PRO A  11      29.555  36.610  28.721  0.00 94.92           H
ATOM    143  HB2 PRO A  11      30.053  34.092  29.873  0.00 87.60           H
ATOM    144  HB3 PRO A  11      31.112  35.031  29.219  0.00 87.60           H
ATOM    145  HG2 PRO A  11      31.130  34.810  31.660  0.00 88.53           H
ATOM    146  HG3 PRO A  11      31.365  36.228  31.055  0.00 88.53           H
ATOM    147  HD2 PRO A  11      28.948  35.266  32.185  0.00 95.29           H
ATOM    148  HD3 PRO A  11      29.584  36.670  32.419  0.00 95.29           H
TER
ATOM    149  N   PHE A  29      32.397  30.508  39.724  1.00 77.91           N
ATOM    150  CA  PHE A  29      32.556  31.745  38.984  1.00 76.50           C
ATOM    151  C   PHE A  29      32.926  31.449  37.545  1.00 69.81           C
ATOM    152  O   PHE A  29      33.617  30.458  37.279  1.00 73.98           O
ATOM    153  CB  PHE A  29      33.612  32.568  39.637  1.00 81.67           C
ATOM    154  CG  PHE A  29      33.423  32.684  41.095  1.00 85.76           C
ATOM    155  CD1 PHE A  29      34.099  31.852  41.957  1.00 91.51           C
ATOM    156  CD2 PHE A  29      32.529  33.616  41.613  1.00 87.86           C
ATOM    157  CE1 PHE A  29      33.895  31.964  43.331  1.00 97.51           C
ATOM    158  CE2 PHE A  29      32.332  33.734  42.967  1.00 83.24           C
ATOM    159  CZ  PHE A  29      33.003  32.909  43.829  1.00 87.57           C
ATOM    160  H   PHE A  29      33.125  30.063  39.832  0.00 77.91           H
ATOM         H2  PHE A  29      32.075  30.688  40.545  1.00 77.91           H
ATOM         H3  PHE A  29      31.822  29.969  39.289  1.00 77.91           H
ATOM    161  HA  PHE A  29      31.720  32.237  38.985  0.00 76.50           H
ATOM    162  HB2 PHE A  29      34.480  32.175  39.458  0.00 81.67           H
ATOM    163  HB3 PHE A  29      33.615  33.455  39.244  0.00 81.67           H
ATOM    164  HD1 PHE A  29      34.691  31.216  41.624  0.00 91.51           H
ATOM    165  HD2 PHE A  29      32.057  34.168  41.033  0.00 87.86           H
ATOM    166  HE1 PHE A  29      34.356  31.407  43.916  0.00 97.51           H
ATOM    167  HE2 PHE A  29      31.743  34.374  43.298  0.00 83.24           H
ATOM    168  HZ  PHE A  29      32.863  32.980  44.746  0.00 87.57           H
ATOM    169  N   HIS A  30      32.435  32.280  36.614  1.00 66.84           N
ATOM    170  CA  HIS A  30      33.074  32.324  35.305  1.00 70.57           C
ATOM    171  C   HIS A  30      34.515  32.792  35.451  1.00 70.51           C
ATOM    172  O   HIS A  30      34.796  33.738  36.176  1.00 73.92           O
ATOM    173  CB  HIS A  30      32.384  33.252  34.313  1.00 70.79           C
ATOM    174  CG  HIS A  30      30.922  33.009  34.116  1.00 68.99           C
ATOM    175  ND1 HIS A  30      29.948  33.839  34.620  1.00 77.35           N
ATOM    176  CD2 HIS A  30      30.275  32.002  33.489  1.00 70.51           C
ATOM    177  CE1 HIS A  30      28.757  33.378  34.268  1.00 79.15           C
ATOM    178  NE2 HIS A  30      28.930  32.267  33.573  1.00 70.69           N
ATOM    179  H   HIS A  30      31.761  32.803  36.719  0.00 66.84           H
ATOM    180  HA  HIS A  30      33.018  31.421  34.955  0.00 70.57           H
ATOM    181  HB2 HIS A  30      32.507  34.167  34.610  0.00 70.79           H
ATOM    182  HB3 HIS A  30      32.828  33.171  33.455  0.00 70.79           H
ATOM    183  HD2 HIS A  30      30.667  31.266  33.077  0.00 70.51           H
ATOM    184  HE1 HIS A  30      27.939  33.769  34.475  0.00 79.15           H
ATOM    185  HE2 HIS A  30      28.303  31.789  33.230  0.00 70.69           H
ATOM    186  N   GLN A  31      35.432  32.140  34.737  1.00 76.35           N
ATOM    187  CA  GLN A  31      36.835  32.546  34.836  1.00 78.58           C
ATOM    188  C   GLN A  31      37.014  34.005  34.435  1.00 83.34           C
ATOM    189  O   GLN A  31      37.648  34.778  35.164  1.00 81.98           O
ATOM    190  CB  GLN A  31      37.723  31.641  33.986  1.00 74.26           C
ATOM    191  CG  GLN A  31      37.968  30.296  34.632  1.00 73.57           C
ATOM    192  CD  GLN A  31      38.731  29.321  33.741  1.00 79.39           C
ATOM    193  OE1 GLN A  31      39.337  29.704  32.736  1.00 86.62           O
ATOM    194  NE2 GLN A  31      38.689  28.046  34.105  1.00 77.37           N
ATOM    195  H   GLN A  31      35.272  31.482  34.207  0.00 76.35           H
ATOM    196  HA  GLN A  31      37.106  32.455  35.763  0.00 78.58           H
ATOM    197  HB2 GLN A  31      37.309  31.509  33.119  0.00 74.26           H
ATOM    198  HB3 GLN A  31      38.573  32.081  33.830  0.00 74.26           H
ATOM    199  HG2 GLN A  31      38.464  30.427  35.455  0.00 73.57           H
ATOM    200  HG3 GLN A  31      37.115  29.902  34.874  0.00 73.57           H
ATOM    201 HE21 GLN A  31      38.257  27.815  34.812  0.00 77.37           H
ATOM    202 HE22 GLN A  31      39.094  27.451  33.634  0.00 77.37           H
ATOM    203  N   VAL A  32      36.415  34.396  33.302  1.00 78.93           N
ATOM    204  CA  VAL A  32      36.457  35.775  32.817  1.00 74.01           C
ATOM    205  C   VAL A  32      36.097  36.760  33.927  1.00 81.89           C
ATOM    206  O   VAL A  32      36.871  37.674  34.238  1.00 84.99           O
ATOM    207  CB  VAL A  32      35.521  35.923  31.602  1.00 86.69           C
ATOM    208  CG1 VAL A  32      35.127  37.354  31.404  1.00 90.42           C
ATOM    209  CG2 VAL A  32      36.183  35.383  30.331  1.00 87.13           C
ATOM    210  H   VAL A  32      35.972  33.862  32.794  0.00 78.93           H
ATOM    211  HA  VAL A  32      37.362  35.984  32.537  0.00 74.01           H
ATOM    212  HB  VAL A  32      34.723  35.402  31.781  0.00 86.69           H
ATOM    213 HG11 VAL A  32      34.539  37.425  30.636  0.00 90.42           H
ATOM    214 HG12 VAL A  32      34.665  37.675  32.194  0.00 90.42           H
ATOM    215 HG13 VAL A  32      35.921  37.891  31.253  0.00 90.42           H
ATOM    216 HG21 VAL A  32      35.577  35.486  29.581  0.00 87.13           H
ATOM    217 HG22 VAL A  32      36.999  35.877  30.156  0.00 87.13           H
ATOM    218 HG23 VAL A  32      36.393  34.443  30.449  0.00 87.13           H
ATOM    219  N   CYS A  33      34.922  36.576  34.548  1.00 77.84           N
ATOM    220  CA  CYS A  33      34.433  37.442  35.624  1.00 77.49           C
ATOM    221  C   CYS A  33      35.387  37.543  36.807  1.00 82.56           C
ATOM    222  O   CYS A  33      35.242  38.441  37.643  1.00 86.85           O
ATOM    223  CB  CYS A  33      33.110  36.917  36.157  1.00 74.32           C
ATOM    224  SG  CYS A  33      31.860  36.747  34.937  1.00 77.49           S
ATOM    225  H   CYS A  33      34.383  35.936  34.351  0.00 77.84           H
ATOM    226  HA  CYS A  33      34.342  38.322  35.226  0.00 77.49           H
ATOM    227  HB2 CYS A  33      33.260  36.054  36.574  0.00 74.32           H
ATOM    228  HB3 CYS A  33      32.792  37.516  36.851  0.00 74.32           H
ATOM    229  N   VAL A  34      36.334  36.638  36.924  1.00 85.00           N
ATOM    230  CA  VAL A  34      37.235  36.606  38.061  1.00 82.82           C
ATOM    231  C   VAL A  34      38.677  36.825  37.617  1.00 85.41           C
ATOM    232  O   VAL A  34      39.576  37.027  38.446  1.00 75.97           O
ATOM    233  CB  VAL A  34      37.026  35.291  38.831  1.00 87.08           C
ATOM    234  CG1 VAL A  34      38.087  35.044  39.776  1.00103.47           C
ATOM    235  CG2 VAL A  34      35.730  35.343  39.595  1.00 93.75           C
ATOM    236  H   VAL A  34      36.477  36.019  36.344  0.00 85.00           H
ATOM    237  HA  VAL A  34      37.035  37.335  38.669  0.00 82.82           H
ATOM    238  HB  VAL A  34      37.015  34.578  38.174  0.00 87.08           H
ATOM    239 HG11 VAL A  34      37.920  34.209  40.240  0.00103.47           H
ATOM    240 HG12 VAL A  34      38.933  34.988  39.305  0.00103.47           H
ATOM    241 HG13 VAL A  34      38.123  35.769  40.419  0.00103.47           H
ATOM    242 HG21 VAL A  34      35.605  34.511  40.078  0.00 93.75           H
ATOM    243 HG22 VAL A  34      35.754  36.081  40.224  0.00 93.75           H
ATOM    244 HG23 VAL A  34      34.994  35.471  38.976  0.00 93.75           H
ATOM    245  N   GLY A  35      38.905  36.864  36.311  1.00 89.72           N
ATOM    246  CA  GLY A  35      40.172  37.323  35.790  1.00 89.74           C
ATOM    247  C   GLY A  35      41.276  36.361  36.111  1.00 87.06           C
ATOM    248  O   GLY A  35      42.243  36.702  36.787  1.00 94.95           O
ATOM         OXT GLY A  35      41.224  35.204  35.695  1.00 87.06           O
ATOM    249  H   GLY A  35      38.335  36.627  35.712  0.00 89.72           H
ATOM    250  HA2 GLY A  35      40.107  37.434  34.829  0.00 89.74           H
ATOM    251  HA3 GLY A  35      40.381  38.194  36.163  0.00 89.74           H
TER
HETATM  252 ZN    ZN A 101      30.094  35.529  35.557  1.00 76.89          Zn
"""

def run(prefix):
  """
  Exercise gradients match: using clustering vs not using clustering.
  This is the case of several CYS coordinating ZN. Gradients match if two
  buffers are used, and they don't if one buffer is used.
  """
  pdb_in = "%s.pdb"%prefix
  open(pdb_in, "w").write(pdb_str_in)
  #
  for two_buffers in [False, True]:
    if(0): print "two_buffers=", two_buffers
    for clustering in ["true", "false"]:
      cmd = " ".join([
        "qr.refine",
        pdb_in,
        "mode=opt",
        "stpmax=0.2",
        "gradient_only=true",
        "clustering=%s"%clustering,
        "dump_gradients=cluster_%s.pkl"%clustering,
        "quantum.engine_name=mopac",
        "number_of_micro_cycles=1",
        "max_iterations_refine=5",
        "two_buffers=%s"%str(two_buffers),
        "> %s.log"%prefix])
      assert easy_run.call(cmd)==0
    g1 = easy_pickle.load("cluster_false.pkl")
    g2 = easy_pickle.load("cluster_true.pkl")
    g1 = g1.as_double()
    g2 = g2.as_double()
    diff = flex.abs(g1-g2)
    if(0):
      print "  min/max/mean of (gradient1 - gradient2):", \
        diff.min_max_mean().as_tuple()

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
