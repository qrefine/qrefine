from __future__ import division
from __future__ import absolute_import

import os
import iotbx.pdb
import mmtbx.f_model
from qrefine.tests.unit import run_qrefine, run_fmodel

pdb_str_good = """
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -8.595   3.702   6.627  1.00 16.77           N
ATOM      2  CA  GLY A   1      -8.864   3.218   5.243  1.00 16.57           C
ATOM      3  C   GLY A   1      -7.821   2.232   4.754  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.626   1.174   5.351  1.00 16.78           O
ATOM      5  H1  GLY A   1      -9.293   4.429   6.882  1.00 16.77           H
ATOM      6  H2  GLY A   1      -8.666   2.905   7.292  1.00 16.77           H
ATOM      7  H3  GLY A   1      -7.639   4.107   6.674  1.00 16.77           H
ATOM      8  HA2 GLY A   1      -8.883   4.067   4.559  1.00 16.57           H
ATOM      9  HA3 GLY A   1      -9.838   2.731   5.212  1.00 16.57           H
ATOM     10  N   ASN A   2      -7.148   2.583   3.660  1.00 15.02           N
ATOM     11  CA  ASN A   2      -6.114   1.739   3.061  1.00 14.10           C
ATOM     12  C   ASN A   2      -4.760   2.257   3.534  1.00 13.13           C
ATOM     13  O   ASN A   2      -4.191   3.185   2.958  1.00 11.91           O
ATOM     14  CB  ASN A   2      -6.224   1.743   1.541  1.00 15.38           C
ATOM     15  CG  ASN A   2      -7.495   1.081   1.046  1.00 14.08           C
ATOM     16  OD1 ASN A   2      -7.736  -0.097   1.307  1.00 17.46           O
ATOM     17  ND2 ASN A   2      -8.315   1.839   0.328  1.00 11.72           N
ATOM     18  H   ASN A   2      -7.298   3.458   3.158  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.233   0.714   3.412  1.00 14.10           H
ATOM     20  HB2 ASN A   2      -6.220   2.774   1.187  1.00 15.38           H
ATOM     21  HB3 ASN A   2      -5.376   1.202   1.122  1.00 15.38           H
ATOM     22 HD21 ASN A   2      -9.186   1.448  -0.032  1.00 11.72           H
ATOM     23 HD22 ASN A   2      -8.075   2.811   0.135  1.00 11.72           H
ATOM     24  N   ASN A   3      -4.242   1.647   4.597  1.00 12.26           N
ATOM     25  CA  ASN A   3      -2.956   2.044   5.148  1.00 11.74           C
ATOM     26  C   ASN A   3      -1.819   1.412   4.355  1.00 11.10           C
ATOM     27  O   ASN A   3      -1.905   0.257   3.929  1.00 10.42           O
ATOM     28  CB  ASN A   3      -2.857   1.636   6.618  1.00 12.15           C
ATOM     29  CG  ASN A   3      -1.574   2.113   7.270  1.00 12.82           C
ATOM     30  OD1 ASN A   3      -1.500   3.233   7.775  1.00 15.05           O
ATOM     31  ND2 ASN A   3      -0.554   1.262   7.263  1.00 13.48           N
ATOM     32  H   ASN A   3      -4.690   0.877   5.094  1.00 12.26           H
ATOM     33  HA  ASN A   3      -2.856   3.128   5.086  1.00 11.74           H
ATOM     34  HB2 ASN A   3      -2.887   0.548   6.688  1.00 12.15           H
ATOM     35  HB3 ASN A   3      -3.695   2.068   7.165  1.00 12.15           H
ATOM     36 HD21 ASN A   3       0.335   1.527   7.687  1.00 13.48           H
ATOM     37 HD22 ASN A   3      -0.660   0.343   6.833  1.00 13.48           H
ATOM     38  N   GLN A   4      -0.746   2.183   4.157  1.00 10.29           N
ATOM     39  CA  GLN A   4       0.428   1.733   3.419  1.00 10.53           C
ATOM     40  C   GLN A   4       1.631   2.542   3.914  1.00 10.24           C
ATOM     41  O   GLN A   4       2.130   3.457   3.260  1.00  8.86           O
ATOM     42  CB  GLN A   4       0.236   1.881   1.908  1.00  9.80           C
ATOM     43  CG  GLN A   4       1.318   1.206   1.080  1.00 10.25           C
ATOM     44  CD  GLN A   4       2.019   2.168   0.141  1.00 12.43           C
ATOM     45  OE1 GLN A   4       2.371   3.283   0.525  1.00 14.62           O
ATOM     46  NE2 GLN A   4       2.228   1.739  -1.098  1.00  9.05           N
ATOM     47  H   GLN A   4      -0.665   3.139   4.504  1.00 10.29           H
ATOM     48  HA  GLN A   4       0.610   0.681   3.639  1.00 10.53           H
ATOM     49  HB2 GLN A   4       0.237   2.941   1.657  1.00  9.80           H
ATOM     50  HB3 GLN A   4      -0.720   1.437   1.632  1.00  9.80           H
ATOM     51  HG2 GLN A   4       2.066   0.781   1.749  1.00 10.25           H
ATOM     52  HG3 GLN A   4       0.865   0.418   0.478  1.00 10.25           H
ATOM     53 HE21 GLN A   4       2.695   2.342  -1.775  1.00  9.05           H
ATOM     54 HE22 GLN A   4       1.920   0.806  -1.374  1.00  9.05           H
ATOM     55  N   GLN A   5       2.111   2.194   5.106  1.00 10.38           N
ATOM     56  CA  GLN A   5       3.251   2.860   5.726  1.00 11.39           C
ATOM     57  C   GLN A   5       4.502   2.033   5.451  1.00 11.52           C
ATOM     58  O   GLN A   5       4.672   0.946   6.014  1.00 12.05           O
ATOM     59  CB  GLN A   5       3.027   3.037   7.226  1.00 11.96           C
ATOM     60  CG  GLN A   5       1.810   3.879   7.574  1.00 10.81           C
ATOM     61  CD  GLN A   5       1.466   3.824   9.050  1.00 13.10           C
ATOM     62  OE1 GLN A   5       0.954   2.819   9.543  1.00 10.65           O
ATOM     63  NE2 GLN A   5       1.746   4.908   9.763  1.00 12.30           N
ATOM     64  H   GLN A   5       1.723   1.441   5.675  1.00 10.38           H
ATOM     65  HA  GLN A   5       3.385   3.846   5.279  1.00 11.39           H
ATOM     66  HB2 GLN A   5       3.902   3.523   7.656  1.00 11.96           H
ATOM     67  HB3 GLN A   5       2.891   2.055   7.678  1.00 11.96           H
ATOM     68  HG2 GLN A   5       2.009   4.918   7.312  1.00 10.81           H
ATOM     69  HG3 GLN A   5       0.950   3.511   7.014  1.00 10.81           H
ATOM     70 HE21 GLN A   5       1.537   4.930  10.762  1.00 12.30           H
ATOM     71 HE22 GLN A   5       2.171   5.719   9.314  1.00 12.30           H
ATOM     72  N   ASN A   6       5.373   2.547   4.586  1.00 11.99           N
ATOM     73  CA  ASN A   6       6.611   1.860   4.219  1.00 12.30           C
ATOM     74  C   ASN A   6       7.684   2.235   5.234  1.00 13.40           C
ATOM     75  O   ASN A   6       8.241   3.335   5.191  1.00 13.92           O
ATOM     76  CB  ASN A   6       7.031   2.229   2.801  1.00 12.13           C
ATOM     77  CG  ASN A   6       6.137   1.604   1.748  1.00 12.77           C
ATOM     78  OD1 ASN A   6       6.297   0.436   1.395  1.00 14.27           O
ATOM     79  ND2 ASN A   6       5.189   2.382   1.238  1.00 10.07           N
ATOM     80  H   ASN A   6       5.250   3.445   4.119  1.00 11.99           H
ATOM     81  HA  ASN A   6       6.455   0.782   4.263  1.00 12.30           H
ATOM     82  HB2 ASN A   6       6.982   3.312   2.686  1.00 12.13           H
ATOM     83  HB3 ASN A   6       8.049   1.881   2.630  1.00 12.13           H
ATOM     84 HD21 ASN A   6       4.558   2.016   0.526  1.00 10.07           H
ATOM     85 HD22 ASN A   6       5.092   3.345   1.560  1.00 10.07           H
ATOM     86  N   TYR A   7       7.977   1.317   6.149  1.00 14.70           N
ATOM     87  CA  TYR A   7       8.992   1.549   7.170  1.00 15.18           C
ATOM     88  C   TYR A   7      10.366   1.113   6.674  1.00 15.91           C
ATOM     89  O   TYR A   7      10.482   0.359   5.708  1.00 15.76           O
ATOM     90  CB  TYR A   7       8.637   0.803   8.458  1.00 15.35           C
ATOM     91  CG  TYR A   7       7.294   1.186   9.037  1.00 14.45           C
ATOM     92  CD1 TYR A   7       7.182   2.217   9.959  1.00 14.80           C
ATOM     93  CD2 TYR A   7       6.138   0.515   8.661  1.00 15.68           C
ATOM     94  CE1 TYR A   7       5.956   2.570  10.491  1.00 14.33           C
ATOM     95  CE2 TYR A   7       4.908   0.860   9.187  1.00 13.46           C
ATOM     96  CZ  TYR A   7       4.823   1.889  10.101  1.00 15.09           C
ATOM     97  OH  TYR A   7       3.600   2.237  10.628  1.00 14.39           O
ATOM     98  OXT TYR A   7      11.393   1.505   7.228  1.00 17.49           O
ATOM     99  H   TYR A   7       7.528   0.403   6.211  1.00 14.70           H
ATOM    100  HA  TYR A   7       9.035   2.614   7.396  1.00 15.18           H
ATOM    101  HB2 TYR A   7       9.398   1.020   9.208  1.00 15.35           H
ATOM    102  HB3 TYR A   7       8.618  -0.266   8.250  1.00 15.35           H
ATOM    103  HD1 TYR A   7       8.068   2.753  10.266  1.00 14.80           H
ATOM    104  HD2 TYR A   7       6.201  -0.291   7.945  1.00 15.68           H
ATOM    105  HE1 TYR A   7       5.886   3.375  11.207  1.00 14.33           H
ATOM    106  HE2 TYR A   7       4.019   0.328   8.884  1.00 13.46           H
ATOM    107  HH  TYR A   7       3.706   2.980  11.258  1.00 14.39           H
TER
HETATM  108  O   HOH A   8      -5.873   4.909   8.094  1.00 22.62           O
HETATM  109  O   HOH A   9      10.611   1.751   2.782  1.00 19.71           O
HETATM  110  O   HOH A  10     -11.354   1.724  -1.488  1.00 17.08           O
HETATM  111  O   HOH A  11      12.626   4.850   9.243  1.00 23.99           O
HETATM  112  O   HOH A  12      14.481   2.076   9.189  1.00 26.17           O
HETATM  113  O   HOH A  13      -2.729   3.257  10.671  1.00 39.15           O
HETATM  114  O   HOH A  14      -1.435   0.824  11.162  1.00 43.49           O
END
"""

pdb_str_poor = """
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -8.432   3.651   6.614  1.00 16.77           N
ATOM      2  CA  GLY A   1      -8.771   3.427   5.460  1.00 16.57           C
ATOM      3  C   GLY A   1      -7.899   2.524   4.543  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.817   1.432   5.238  1.00 16.78           O
ATOM      5  H1  GLY A   1      -9.043   4.291   7.085  1.00 16.77           H
ATOM      6  H2  GLY A   1      -8.524   2.687   7.413  1.00 16.77           H
ATOM      7  H3  GLY A   1      -7.812   3.811   6.485  1.00 16.77           H
ATOM      8  HA2 GLY A   1      -9.109   3.851   4.708  1.00 16.57           H
ATOM      9  HA3 GLY A   1      -9.881   2.952   4.970  1.00 16.57           H
ATOM     10  N   ASN A   2      -7.380   2.819   3.800  1.00 15.02           N
ATOM     11  CA  ASN A   2      -6.405   1.880   3.160  1.00 14.10           C
ATOM     12  C   ASN A   2      -4.915   2.484   3.406  1.00 13.13           C
ATOM     13  O   ASN A   2      -4.451   3.363   2.741  1.00 11.91           O
ATOM     14  CB  ASN A   2      -6.134   1.905   1.345  1.00 15.38           C
ATOM     15  CG  ASN A   2      -7.714   1.227   1.338  1.00 14.08           C
ATOM     16  OD1 ASN A   2      -7.553  -0.076   1.137  1.00 17.46           O
ATOM     17  ND2 ASN A   2      -8.544   1.560   0.075  1.00 11.72           N
ATOM     18  H   ASN A   2      -7.354   3.430   3.150  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.402   0.811   3.222  1.00 14.10           H
ATOM     20  HB2 ASN A   2      -6.029   3.042   1.347  1.00 15.38           H
ATOM     21  HB3 ASN A   2      -5.334   1.326   1.404  1.00 15.38           H
ATOM     22 HD21 ASN A   2      -9.324   1.384   0.029  1.00 11.72           H
ATOM     23 HD22 ASN A   2      -8.201   2.904   0.371  1.00 11.72           H
ATOM     24  N   ASN A   3      -3.949   1.424   4.504  1.00 12.26           N
ATOM     25  CA  ASN A   3      -2.781   1.767   5.152  1.00 11.74           C
ATOM     26  C   ASN A   3      -1.984   1.262   4.332  1.00 11.10           C
ATOM     27  O   ASN A   3      -1.950   0.169   3.961  1.00 10.42           O
ATOM     28  CB  ASN A   3      -2.696   1.480   6.752  1.00 12.15           C
ATOM     29  CG  ASN A   3      -1.470   2.168   7.174  1.00 12.82           C
ATOM     30  OD1 ASN A   3      -1.797   3.508   7.598  1.00 15.05           O
ATOM     31  ND2 ASN A   3      -0.776   1.330   7.292  1.00 13.48           N
ATOM     32  H   ASN A   3      -4.727   0.842   5.246  1.00 12.26           H
ATOM     33  HA  ASN A   3      -2.813   3.338   5.166  1.00 11.74           H
ATOM     34  HB2 ASN A   3      -2.797   0.765   6.992  1.00 12.15           H
ATOM     35  HB3 ASN A   3      -3.935   2.167   7.002  1.00 12.15           H
ATOM     36 HD21 ASN A   3       0.170   1.702   7.707  1.00 13.48           H
ATOM     37 HD22 ASN A   3      -0.515   0.496   6.846  1.00 13.48           H
ATOM     38  N   GLN A   4      -0.840   2.423   3.917  1.00 10.29           N
ATOM     39  CA  GLN A   4       0.128   1.485   3.593  1.00 10.53           C
ATOM     40  C   GLN A   4       1.724   2.610   3.942  1.00 10.24           C
ATOM     41  O   GLN A   4       2.416   3.455   3.506  1.00  8.86           O
ATOM     42  CB  GLN A   4       0.542   1.704   1.695  1.00  9.80           C
ATOM     43  CG  GLN A   4       1.540   1.055   1.102  1.00 10.25           C
ATOM     44  CD  GLN A   4       2.322   2.342   0.156  1.00 12.43           C
ATOM     45  OE1 GLN A   4       2.189   3.185   0.825  1.00 14.62           O
ATOM     46  NE2 GLN A   4       2.233   1.774  -0.804  1.00  9.05           N
ATOM     47  H   GLN A   4      -0.682   3.287   4.776  1.00 10.29           H
ATOM     48  HA  GLN A   4       0.319   0.613   3.840  1.00 10.53           H
ATOM     49  HB2 GLN A   4       0.364   3.032   1.879  1.00  9.80           H
ATOM     50  HB3 GLN A   4      -0.948   1.454   1.753  1.00  9.80           H
ATOM     51  HG2 GLN A   4       1.863   0.902   1.823  1.00 10.25           H
ATOM     52  HG3 GLN A   4       0.917   0.721   0.546  1.00 10.25           H
ATOM     53 HE21 GLN A   4       2.559   2.335  -1.710  1.00  9.05           H
ATOM     54 HE22 GLN A   4       2.070   0.947  -1.559  1.00  9.05           H
ATOM     55  N   GLN A   5       1.909   2.123   5.364  1.00 10.38           N
ATOM     56  CA  GLN A   5       3.000   2.979   5.796  1.00 11.39           C
ATOM     57  C   GLN A   5       4.426   1.882   5.232  1.00 11.52           C
ATOM     58  O   GLN A   5       4.644   1.004   6.243  1.00 12.05           O
ATOM     59  CB  GLN A   5       3.252   2.749   7.175  1.00 11.96           C
ATOM     60  CG  GLN A   5       1.655   3.933   7.447  1.00 10.81           C
ATOM     61  CD  GLN A   5       1.757   3.729   8.983  1.00 13.10           C
ATOM     62  OE1 GLN A   5       1.110   2.616   9.343  1.00 10.65           O
ATOM     63  NE2 GLN A   5       1.678   4.777   9.976  1.00 12.30           N
ATOM     64  H   GLN A   5       1.723   1.669   5.977  1.00 10.38           H
ATOM     65  HA  GLN A   5       3.514   3.610   5.238  1.00 11.39           H
ATOM     66  HB2 GLN A   5       4.167   3.469   7.534  1.00 11.96           H
ATOM     67  HB3 GLN A   5       3.014   2.128   7.619  1.00 11.96           H
ATOM     68  HG2 GLN A   5       1.734   4.900   7.402  1.00 10.81           H
ATOM     69  HG3 GLN A   5       1.059   3.393   6.924  1.00 10.81           H
ATOM     70 HE21 GLN A   5       1.484   4.914  10.570  1.00 12.30           H
ATOM     71 HE22 GLN A   5       2.385   6.020   9.089  1.00 12.30           H
ATOM     72  N   ASN A   6       5.263   2.342   4.520  1.00 11.99           N
ATOM     73  CA  ASN A   6       6.817   2.058   3.959  1.00 12.30           C
ATOM     74  C   ASN A   6       7.389   2.175   5.515  1.00 13.40           C
ATOM     75  O   ASN A   6       8.101   3.610   5.272  1.00 13.92           O
ATOM     76  CB  ASN A   6       7.321   2.105   3.101  1.00 12.13           C
ATOM     77  CG  ASN A   6       6.021   1.447   1.472  1.00 12.77           C
ATOM     78  OD1 ASN A   6       6.472   0.206   1.286  1.00 14.27           O
ATOM     79  ND2 ASN A   6       5.050   2.473   1.152  1.00 10.07           N
ATOM     80  H   ASN A   6       5.070   3.481   4.224  1.00 11.99           H
ATOM     81  HA  ASN A   6       6.500   0.696   4.268  1.00 12.30           H
ATOM     82  HB2 ASN A   6       7.241   3.175   2.829  1.00 12.13           H
ATOM     83  HB3 ASN A   6       8.057   2.132   2.479  1.00 12.13           H
ATOM     84 HD21 ASN A   6       4.352   1.910   0.574  1.00 10.07           H
ATOM     85 HD22 ASN A   6       4.963   3.097   1.793  1.00 10.07           H
ATOM     86  N   TYR A   7       7.701   1.194   5.946  1.00 14.70           N
ATOM     87  CA  TYR A   7       9.073   1.599   7.175  1.00 15.18           C
ATOM     88  C   TYR A   7      10.531   1.060   6.740  1.00 15.91           C
ATOM     89  O   TYR A   7      10.733   0.328   5.701  1.00 15.76           O
ATOM     90  CB  TYR A   7       8.474   0.641   8.233  1.00 15.35           C
ATOM     91  CG  TYR A   7       7.329   1.353   9.195  1.00 14.45           C
ATOM     92  CD1 TYR A   7       7.350   2.420  10.059  1.00 14.80           C
ATOM     93  CD2 TYR A   7       6.230   0.305   8.569  1.00 15.68           C
ATOM     94  CE1 TYR A   7       5.940   2.625  10.498  1.00 14.33           C
ATOM     95  CE2 TYR A   7       4.670   0.962   8.935  1.00 13.46           C
ATOM     96  CZ  TYR A   7       5.046   1.861   9.911  1.00 15.09           C
ATOM     97  OH  TYR A   7       3.442   2.474  10.691  1.00 14.39           O
ATOM     98  OXT TYR A   7      11.303   1.236   7.020  1.00 17.49           O
ATOM     99  H   TYR A   7       7.368   0.673   6.116  1.00 14.70           H
ATOM    100  HA  TYR A   7       8.767   2.543   7.488  1.00 15.18           H
ATOM    101  HB2 TYR A   7       9.332   0.902   8.958  1.00 15.35           H
ATOM    102  HB3 TYR A   7       8.726   0.019   8.399  1.00 15.35           H
ATOM    103  HD1 TYR A   7       8.126   2.754  10.062  1.00 14.80           H
ATOM    104  HD2 TYR A   7       6.335   0.012   7.894  1.00 15.68           H
ATOM    105  HE1 TYR A   7       5.746   3.127  11.297  1.00 14.33           H
ATOM    106  HE2 TYR A   7       3.916   0.063   8.971  1.00 13.46           H
ATOM    107  HH  TYR A   7       3.972   3.067  11.542  1.00 14.39           H
TER
HETATM  108  O   HOH A   8      -5.844   4.636   8.176  1.00 22.62           O
HETATM  109  O   HOH A   9      10.818   1.772   2.609  1.00 19.71           O
HETATM  110  O   HOH A  10     -11.403   2.028  -1.447  1.00 17.08           O
HETATM  111  O   HOH A  11      12.482   4.990   9.478  1.00 23.99           O
HETATM  112  O   HOH A  12      14.531   1.828   9.286  1.00 26.17           O
HETATM  113  O   HOH A  13      -2.757   3.192  10.832  1.00 39.15           O
HETATM  114  O   HOH A  14      -1.566   0.866  10.971  1.00 43.49           O
END
"""

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Refinement with CCTBX restraints, clustering is True
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  pdb_good = "%s_good.pdb"%prefix
  with open(pdb_good, "w") as fo:
    fo.write(pdb_str_good)
  #
  pdb_poor = "%s_poor.pdb"%prefix
  with open(pdb_poor, "w") as fo:
    fo.write(pdb_str_poor)
  #
  r_fm = run_fmodel(prefix = prefix, args=[pdb_good, "high_res=2"])
  r_qr = run_qrefine(prefix = prefix, args=[pdb_poor, r_fm.mtz, 
    "restraints=cctbx", "stpmax=1", "clustering=True"])
  fmodel = mmtbx.f_model.manager(
    f_obs          = r_fm.f_obs,
    r_free_flags   = r_fm.flags,
    xray_structure = r_qr.model.get_xray_structure())
  fmodel.update_all_scales()
  assert fmodel.r_work() < 0.005, fmodel.r_work()

if(__name__ == "__main__"):
  run()