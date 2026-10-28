import os
import iotbx.pdb
from scitbx.array_family import flex
from libtbx import easy_pickle
import time, sys
import libtbx.load_env
from libtbx.test_utils import approx_equal
from qrefine.command_line import granalyse
from qrefine.tests.unit import run_qrefine_in

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

pdb_str_in = """
CRYST1  109.965  109.703  107.335  90.00  90.00  90.00 P 1
SCALE1      0.009094  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009116  0.000000        0.00000
SCALE3      0.000000  0.000000  0.009317        0.00000
ATOM      1  N   CYS A   3      50.936  52.907  54.091  1.00 76.61           N
ATOM      2  CA  CYS A   3      51.734  52.240  55.098  1.00 78.13           C
ATOM      3  C   CYS A   3      51.581  52.965  56.433  1.00 83.42           C
ATOM      4  O   CYS A   3      50.648  53.748  56.632  1.00 87.64           O
ATOM      5  CB  CYS A   3      53.184  52.133  54.624  1.00 78.38           C
ATOM      6  SG  CYS A   3      54.247  53.578  54.761  1.00 63.83           S
ATOM      7  OXT CYS A   3      52.397  52.773  57.335  1.00 83.42           O
ATOM      8  H   CYS A   3      51.219  53.690  53.877  0.00 76.61           H
ATOM      9  H2  CYS A   3      50.929  52.360  53.254  0.00 76.61           H
ATOM     10  H3  CYS A   3      50.000  53.016  54.425  0.00 76.61           H
ATOM     11  HA  CYS A   3      51.403  51.215  55.263  0.00 78.13           H
ATOM     12  HB2 CYS A   3      53.661  51.342  55.203  0.00 78.38           H
ATOM     13  HB3 CYS A   3      53.165  51.862  53.568  0.00 78.38           H
ATOM     14  N   CYS A   8      52.848  58.796  54.595  1.00 96.46           N
ATOM     15  CA  CYS A   8      52.885  57.791  53.550  1.00 84.79           C
ATOM     16  C   CYS A   8      52.541  58.484  52.264  1.00 88.83           C
ATOM     17  O   CYS A   8      51.477  59.096  52.149  1.00 93.29           O
ATOM     18  CB  CYS A   8      51.916  56.638  53.808  1.00 85.00           C
ATOM     19  SG  CYS A   8      51.677  55.521  52.362  1.00 78.50           S
ATOM     20  OXT CYS A   8      53.324  58.444  51.315  1.00 88.83           O
ATOM     21  H   CYS A   8      52.199  58.721  55.154  0.00 96.46           H
ATOM     22  H2  CYS A   8      53.691  58.751  55.131  0.00 96.46           H
ATOM     23  H3  CYS A   8      52.764  59.703  54.183  0.00 96.46           H
ATOM     24  HA  CYS A   8      53.871  57.328  53.502  0.00 84.79           H
ATOM     25  HB2 CYS A   8      52.301  56.036  54.632  0.00 85.00           H
ATOM     26  HB3 CYS A   8      50.943  57.052  54.071  0.00 85.00           H
TER
ATOM     27  N   HIS A  30      56.097  51.132  53.603  1.00 66.84           N
ATOM     28  CA  HIS A  30      56.736  51.176  52.294  1.00 70.57           C
ATOM     29  C   HIS A  30      58.177  51.644  52.440  1.00 70.51           C
ATOM     30  O   HIS A  30      58.458  52.590  53.165  1.00 73.92           O
ATOM     31  CB  HIS A  30      56.046  52.104  51.302  1.00 70.79           C
ATOM     32  CG  HIS A  30      54.584  51.861  51.105  1.00 68.99           C
ATOM     33  ND1 HIS A  30      53.610  52.691  51.609  1.00 77.35           N
ATOM     34  CD2 HIS A  30      53.937  50.854  50.478  1.00 70.51           C
ATOM     35  CE1 HIS A  30      52.419  52.230  51.257  1.00 79.15           C
ATOM     36  NE2 HIS A  30      52.592  51.119  50.562  1.00 70.69           N
ATOM     37  OXT HIS A  30      59.076  51.070  51.825  1.00 70.51           O
ATOM     38  H   HIS A  30      55.423  51.655  53.708  0.00 66.84           H
ATOM     39  H2  HIS A  30      55.762  50.206  53.780  0.00 66.84           H
ATOM     40  H3  HIS A  30      56.760  51.387  54.307  0.00 66.84           H
ATOM     41  HA  HIS A  30      56.677  50.169  51.882  0.00 70.57           H
ATOM     42  HB2 HIS A  30      56.160  53.129  51.655  0.00 70.79           H
ATOM     43  HB3 HIS A  30      56.528  51.988  50.331  0.00 70.79           H
ATOM     44  HD2 HIS A  30      54.393  50.000  50.000  0.00 70.51           H
ATOM     45  HE1 HIS A  30      51.469  52.684  51.497  0.00 79.15           H
ATOM     46  HE2 HIS A  30      51.849  50.552  50.155  0.00 70.69           H
ATOM     47  N   CYS A  33      58.584  55.428  51.537  1.00 77.84           N
ATOM     48  CA  CYS A  33      58.095  56.294  52.613  1.00 77.49           C
ATOM     49  C   CYS A  33      59.049  56.395  53.796  1.00 82.56           C
ATOM     50  O   CYS A  33      58.904  57.293  54.632  1.00 86.85           O
ATOM     51  CB  CYS A  33      56.772  55.769  53.146  1.00 74.32           C
ATOM     52  SG  CYS A  33      55.522  55.599  51.926  1.00 77.49           S
ATOM     53  OXT CYS A  33      59.965  55.582  53.919  1.00 82.56           O
ATOM     54  H   CYS A  33      58.045  54.788  51.340  0.00 77.84           H
ATOM     55  H2  CYS A  33      58.720  55.970  50.708  0.00 77.84           H
ATOM     56  H3  CYS A  33      59.452  55.014  51.810  0.00 77.84           H
ATOM     57  HA  CYS A  33      57.981  57.289  52.183  0.00 77.49           H
ATOM     58  HB2 CYS A  33      56.939  54.786  53.587  0.00 74.32           H
ATOM     59  HB3 CYS A  33      56.403  56.460  53.904  0.00 74.32           H
TER
HETATM   60 ZN    ZN A 101      53.756  54.381  52.546  1.00 76.89          Zn
END
"""


def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Exercise gradients match: using clustering vs not using clustering.
  This is the case of several CYS coordinating ZN. EXTRA-HUGE box.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  pdb_in = "%s.pdb"%prefix
  open(pdb_in, "w").write(pdb_str_in)
  #
  #for two_buffers in [False, True]:
  for two_buffers in [False, ]:
    print("    two_buffers=", two_buffers)
    for clustering in ["False", "True"]:
      print("      clustering=", clustering)
      args = [
        pdb_in,
        "mode=opt",
        "stpmax=0.2",
        "gradient_only=true",
        "minimizer=lbfgs",
        "clustering=%s"%clustering,
        "dump_gradients=cluster_%s.pkl"%clustering,
        "restraints=qm",
        "quantum.engine_name=mopac",
        
        #"distance_buffer=True","distance_buffer_radius=6",
        
        #"restraints=cctbx",
         
        "debug=True",
        "number_of_micro_cycles=1",
        "max_iterations_refine=1",
        "two_buffers=%s"%str(two_buffers),
        ]
      r_qr = run_qrefine_in(prefix = prefix, args=args)
  g1 = easy_pickle.load(os.path.join(r_qr.folder,
    "debug_optimization_clusteringFalse.pkl")).gradients
  g2 = easy_pickle.load(os.path.join(r_qr.folder,
    "debug_optimization_clusteringTrue.pkl")).gradients
  assert len(g1) == len(g2)
  for g1i, g2i in zip(g1, g2):
    print()
    diffs = flex.sqrt((g1i-g2i).dot())
    norms = flex.sqrt(g1i.dot())
    print("|g1-g2| per atom min/max/mean:", diffs.min_max_mean().as_tuple())
    print("|g1|    per atom min/max/mean:", norms.min_max_mean().as_tuple())
    d = max(granalyse.get_grad_wdelta(ref=g1i, g=g2i))
    assert d < 1.e-6
    print(d, "get_grad_wdelta")
    print()

if(__name__ == '__main__'):
  run()
