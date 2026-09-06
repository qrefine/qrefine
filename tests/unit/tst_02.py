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
ATOM      1  N   GLY A   1      -9.342   3.680   5.865  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.670   3.065   4.593  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.491   2.366   3.943  1.00 16.16           C
ATOM      4  O   GLY A   1      -8.164   1.231   4.288  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.851   3.051   2.996  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.696   2.499   2.286  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.409   2.955   2.974  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.645   3.784   2.475  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.728   2.910   0.819  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.923   2.340   0.080  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -7.873   1.224  -0.437  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -9.007   3.106   0.028  1.00 11.72           N
ATOM     13  N   ASN A   3      -5.180   2.385   4.153  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.998   2.711   4.935  1.00 11.74           C
ATOM     15  C   ASN A   3      -2.776   1.979   4.394  1.00 11.10           C
ATOM     16  O   ASN A   3      -2.864   0.848   3.909  1.00 10.42           O
ATOM     17  CB  ASN A   3      -4.211   2.352   6.406  1.00 12.15           C
ATOM     18  CG  ASN A   3      -3.248   3.075   7.328  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -3.519   4.187   7.779  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -2.115   2.443   7.612  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.625   2.643   4.481  1.00 10.29           N
ATOM     22  CA  GLN A   4      -0.368   2.069   4.005  1.00 10.53           C
ATOM     23  C   GLN A   4       0.768   2.698   4.797  1.00 10.24           C
ATOM     24  O   GLN A   4       1.014   3.902   4.680  1.00  8.86           O
ATOM     25  CB  GLN A   4      -0.189   2.301   2.507  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.062   1.662   1.926  1.00 10.25           C
ATOM     27  CD  GLN A   4       1.519   2.334   0.646  1.00 12.43           C
ATOM     28  OE1 GLN A   4       1.290   3.526   0.442  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.171   1.571  -0.223  1.00  9.05           N
ATOM     30  N   GLN A   5       1.454   1.888   5.598  1.00 10.38           N
ATOM     31  CA  GLN A   5       2.574   2.339   6.420  1.00 11.39           C
ATOM     32  C   GLN A   5       3.842   1.672   5.896  1.00 11.52           C
ATOM     33  O   GLN A   5       4.119   0.510   6.209  1.00 12.05           O
ATOM     34  CB  GLN A   5       2.337   2.013   7.892  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.231   2.831   8.539  1.00 10.81           C
ATOM     36  CD  GLN A   5       0.811   2.281   9.887  1.00 13.10           C
ATOM     37  OE1 GLN A   5       0.121   1.265   9.968  1.00 10.65           O
ATOM     38  NE2 GLN A   5       1.227   2.952  10.955  1.00 12.30           N
ATOM     39  N   ASN A   6       4.610   2.412   5.098  1.00 11.99           N
ATOM     40  CA  ASN A   6       5.854   1.904   4.519  1.00 12.30           C
ATOM     41  C   ASN A   6       6.985   2.186   5.499  1.00 13.40           C
ATOM     42  O   ASN A   6       7.609   3.248   5.476  1.00 13.92           O
ATOM     43  CB  ASN A   6       6.114   2.538   3.158  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.165   2.032   2.089  1.00 12.77           C
ATOM     45  OD1 ASN A   6       4.673   0.907   2.162  1.00 14.27           O
ATOM     46  ND2 ASN A   6       4.905   2.864   1.087  1.00 10.07           N
ATOM     47  N   TYR A   7       7.252   1.220   6.372  1.00 14.70           N
ATOM     48  CA  TYR A   7       8.309   1.356   7.366  1.00 15.18           C
ATOM     49  C   TYR A   7       9.687   1.246   6.721  1.00 15.91           C
ATOM     50  O   TYR A   7       9.831   0.712   5.621  1.00 15.76           O
ATOM     51  CB  TYR A   7       8.155   0.297   8.461  1.00 15.35           C
ATOM     52  CG  TYR A   7       6.835   0.362   9.195  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.646   1.248  10.247  1.00 14.80           C
ATOM     54  CD2 TYR A   7       5.778  -0.464   8.836  1.00 15.68           C
ATOM     55  CE1 TYR A   7       5.441   1.310  10.922  1.00 14.33           C
ATOM     56  CE2 TYR A   7       4.569  -0.409   9.504  1.00 13.46           C
ATOM     57  CZ  TYR A   7       4.406   0.480  10.546  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.205   0.538  11.214  1.00 14.39           O
ATOM     59  OXT TYR A   7      10.688   1.689   7.284  1.00 17.49           O
TER
HETATM   60  O   HOH A   8      -6.584   4.782   6.726  1.00 22.62           O
HETATM   61  O   HOH A   9       9.506   2.325   2.651  1.00 19.71           O
HETATM   62  O   HOH A  10     -11.771   2.883  -1.893  1.00 17.08           O
HETATM   63  O   HOH A  11      12.458   4.657   8.804  1.00 23.99           O
HETATM   64  O   HOH A  12      14.153   2.233   8.233  1.00 26.17           O
HETATM   65  O   HOH A  13      -2.657   2.115  10.792  1.00 39.15           O
HETATM   66  O   HOH A  14      -1.249  -0.406  12.085  1.00 43.49           O
END
"""

pdb_str_poor = """
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.459   3.510   5.663  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.653   3.104   4.720  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.668   2.542   3.775  1.00 16.16           C
ATOM      4  O   GLY A   1      -8.189   1.413   4.141  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.954   2.958   3.012  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.648   2.516   2.359  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.308   2.977   3.034  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.485   3.810   2.358  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.927   3.037   1.000  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.737   2.513   0.015  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -7.952   1.065  -0.629  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.904   3.178   0.091  1.00 11.72           N
ATOM     13  N   ASN A   3      -5.106   2.552   4.335  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.877   2.912   5.042  1.00 11.74           C
ATOM     15  C   ASN A   3      -2.738   1.912   4.250  1.00 11.10           C
ATOM     16  O   ASN A   3      -2.997   0.697   3.873  1.00 10.42           O
ATOM     17  CB  ASN A   3      -4.314   2.410   6.465  1.00 12.15           C
ATOM     18  CG  ASN A   3      -3.143   2.944   7.265  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -3.629   4.284   7.835  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -2.301   2.414   7.792  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.773   2.834   4.524  1.00 10.29           N
ATOM     22  CA  GLN A   4      -0.433   2.235   4.100  1.00 10.53           C
ATOM     23  C   GLN A   4       0.863   2.542   4.879  1.00 10.24           C
ATOM     24  O   GLN A   4       1.071   3.734   4.766  1.00  8.86           O
ATOM     25  CB  GLN A   4      -0.197   2.257   2.491  1.00  9.80           C
ATOM     26  CG  GLN A   4       0.875   1.525   1.853  1.00 10.25           C
ATOM     27  CD  GLN A   4       1.432   2.530   0.445  1.00 12.43           C
ATOM     28  OE1 GLN A   4       1.398   3.391   0.339  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.188   1.698  -0.122  1.00  9.05           N
ATOM     30  N   GLN A   5       1.429   2.044   5.626  1.00 10.38           N
ATOM     31  CA  GLN A   5       2.717   2.470   6.324  1.00 11.39           C
ATOM     32  C   GLN A   5       3.689   1.806   5.912  1.00 11.52           C
ATOM     33  O   GLN A   5       4.075   0.311   6.368  1.00 12.05           O
ATOM     34  CB  GLN A   5       2.487   2.139   8.044  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.286   2.887   8.602  1.00 10.81           C
ATOM     36  CD  GLN A   5       1.002   2.237   9.881  1.00 13.10           C
ATOM     37  OE1 GLN A   5      -0.061   1.468  10.062  1.00 10.65           O
ATOM     38  NE2 GLN A   5       1.201   3.083  11.024  1.00 12.30           N
ATOM     39  N   ASN A   6       4.778   2.294   5.275  1.00 11.99           N
ATOM     40  CA  ASN A   6       5.848   1.828   4.380  1.00 12.30           C
ATOM     41  C   ASN A   6       6.960   2.296   5.679  1.00 13.40           C
ATOM     42  O   ASN A   6       7.808   3.414   5.526  1.00 13.92           O
ATOM     43  CB  ASN A   6       6.158   2.620   3.220  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.088   1.892   1.935  1.00 12.77           C
ATOM     45  OD1 ASN A   6       4.679   1.048   2.113  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.040   2.792   1.222  1.00 10.07           N
ATOM     47  N   TYR A   7       7.375   1.303   6.416  1.00 14.70           N
ATOM     48  CA  TYR A   7       8.154   1.178   7.346  1.00 15.18           C
ATOM     49  C   TYR A   7       9.518   1.374   6.758  1.00 15.91           C
ATOM     50  O   TYR A   7       9.814   0.541   5.442  1.00 15.76           O
ATOM     51  CB  TYR A   7       8.159   0.426   8.488  1.00 15.35           C
ATOM     52  CG  TYR A   7       6.902   0.252   9.128  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.660   1.302  10.378  1.00 14.80           C
ATOM     54  CD2 TYR A   7       5.876  -0.539   8.736  1.00 15.68           C
ATOM     55  CE1 TYR A   7       5.585   1.138  10.877  1.00 14.33           C
ATOM     56  CE2 TYR A   7       4.628  -0.455   9.476  1.00 13.46           C
ATOM     57  CZ  TYR A   7       4.358   0.552  10.357  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.172   0.550  11.401  1.00 14.39           O
ATOM     59  OXT TYR A   7      10.717   1.652   7.258  1.00 17.49           O
TER
HETATM   60  O   HOH A   8      -6.602   4.759   6.855  1.00 22.62           O
HETATM   61  O   HOH A   9       9.693   2.396   2.540  1.00 19.71           O
HETATM   62  O   HOH A  10     -11.726   2.953  -1.969  1.00 17.08           O
HETATM   63  O   HOH A  11      12.407   4.738   8.725  1.00 23.99           O
HETATM   64  O   HOH A  12      14.188   2.403   8.259  1.00 26.17           O
HETATM   65  O   HOH A  13      -2.729   2.068  10.904  1.00 39.15           O
HETATM   66  O   HOH A  14      -1.123  -0.396  12.207  1.00 43.49           O
END
"""

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Refinement with CCTBX restraints, clustering is True, no H
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