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
from qrefine.command_line import granalyse

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

pdb_str_in = """
CRYST1   70.364   66.358   73.268  90.00  90.00  90.00 P 1
SCALE1      0.014212  0.000000  0.000000        0.00000
SCALE2      0.000000  0.015070  0.000000        0.00000
SCALE3      0.000000  0.000000  0.013649        0.00000
ATOM      1  N   CYS A   3      27.274  34.055  37.102  1.00 76.61           N
ATOM      2  CA  CYS A   3      28.072  33.388  38.109  1.00 78.13           C
ATOM      3  C   CYS A   3      27.919  34.113  39.444  1.00 83.42           C
ATOM      4  O   CYS A   3      26.986  34.896  39.643  1.00 87.64           O
ATOM      5  CB  CYS A   3      29.522  33.281  37.635  1.00 78.38           C
ATOM      6  SG  CYS A   3      30.585  34.726  37.772  1.00 63.83           S
ATOM      7  OXT CYS A   3      28.735  33.921  40.346  1.00 83.42           O
ATOM      8  H   CYS A   3      27.557  34.838  36.888  0.00 76.61           H
ATOM      9  H2  CYS A   3      27.267  33.508  36.265  0.00 76.61           H
ATOM     10  H3  CYS A   3      26.338  34.164  37.436  0.00 76.61           H
ATOM     11  HA  CYS A   3      27.741  32.363  38.274  0.00 78.13           H
ATOM     12  HB2 CYS A   3      29.999  32.490  38.214  0.00 78.38           H
ATOM     13  HB3 CYS A   3      29.503  33.010  36.579  0.00 78.38           H
ATOM     14  N   CYS A   8      29.186  39.944  37.606  1.00 96.46           N
ATOM     15  CA  CYS A   8      29.223  38.939  36.561  1.00 84.79           C
ATOM     16  C   CYS A   8      28.879  39.632  35.275  1.00 88.83           C
ATOM     17  O   CYS A   8      27.815  40.244  35.160  1.00 93.29           O
ATOM     18  CB  CYS A   8      28.254  37.786  36.819  1.00 85.00           C
ATOM     19  SG  CYS A   8      28.015  36.669  35.373  1.00 78.50           S
ATOM     20  OXT CYS A   8      29.662  39.592  34.326  1.00 88.83           O
ATOM     21  H   CYS A   8      28.537  39.869  38.165  0.00 96.46           H
ATOM     22  H2  CYS A   8      30.029  39.899  38.142  0.00 96.46           H
ATOM     23  H3  CYS A   8      29.102  40.851  37.194  0.00 96.46           H
ATOM     24  HA  CYS A   8      30.209  38.476  36.513  0.00 84.79           H
ATOM     25  HB2 CYS A   8      28.639  37.184  37.643  0.00 85.00           H
ATOM     26  HB3 CYS A   8      27.281  38.200  37.082  0.00 85.00           H
TER
ATOM     27  N   HIS A  30      32.435  32.280  36.614  1.00 66.84           N
ATOM     28  CA  HIS A  30      33.074  32.324  35.305  1.00 70.57           C
ATOM     29  C   HIS A  30      34.515  32.792  35.451  1.00 70.51           C
ATOM     30  O   HIS A  30      34.796  33.738  36.176  1.00 73.92           O
ATOM     31  CB  HIS A  30      32.384  33.252  34.313  1.00 70.79           C
ATOM     32  CG  HIS A  30      30.922  33.009  34.116  1.00 68.99           C
ATOM     33  ND1 HIS A  30      29.948  33.839  34.620  1.00 77.35           N
ATOM     34  CD2 HIS A  30      30.275  32.002  33.489  1.00 70.51           C
ATOM     35  CE1 HIS A  30      28.757  33.378  34.268  1.00 79.15           C
ATOM     36  NE2 HIS A  30      28.930  32.267  33.573  1.00 70.69           N
ATOM     37  OXT HIS A  30      35.414  32.218  34.836  1.00 70.51           O
ATOM     38  H   HIS A  30      31.761  32.803  36.719  0.00 66.84           H
ATOM     39  H2  HIS A  30      32.100  31.354  36.791  0.00 66.84           H
ATOM     40  H3  HIS A  30      33.098  32.535  37.318  0.00 66.84           H
ATOM     41  HA  HIS A  30      33.015  31.317  34.893  0.00 70.57           H
ATOM     42  HB2 HIS A  30      32.498  34.277  34.666  0.00 70.79           H
ATOM     43  HB3 HIS A  30      32.866  33.136  33.342  0.00 70.79           H
ATOM     44  HD2 HIS A  30      30.731  31.148  33.011  0.00 70.51           H
ATOM     45  HE1 HIS A  30      27.807  33.832  34.508  0.00 79.15           H
ATOM     46  HE2 HIS A  30      28.187  31.700  33.166  0.00 70.69           H
ATOM     47  N   CYS A  33      34.922  36.576  34.548  1.00 77.84           N
ATOM     48  CA  CYS A  33      34.433  37.442  35.624  1.00 77.49           C
ATOM     49  C   CYS A  33      35.387  37.543  36.807  1.00 82.56           C
ATOM     50  O   CYS A  33      35.242  38.441  37.643  1.00 86.85           O
ATOM     51  CB  CYS A  33      33.110  36.917  36.157  1.00 74.32           C
ATOM     52  SG  CYS A  33      31.860  36.747  34.937  1.00 77.49           S
ATOM     53  OXT CYS A  33      36.303  36.730  36.930  1.00 82.56           O
ATOM     54  H   CYS A  33      34.383  35.936  34.351  0.00 77.84           H
ATOM     55  H2  CYS A  33      35.058  37.118  33.719  0.00 77.84           H
ATOM     56  H3  CYS A  33      35.790  36.162  34.821  0.00 77.84           H
ATOM     57  HA  CYS A  33      34.319  38.437  35.194  0.00 77.49           H
ATOM     58  HB2 CYS A  33      33.277  35.934  36.598  0.00 74.32           H
ATOM     59  HB3 CYS A  33      32.741  37.608  36.915  0.00 74.32           H
TER
HETATM   60 ZN    ZN A 101      30.094  35.529  35.557  1.00 76.89          Zn
"""


def run(prefix):
  """
  Exercise gradients match: using clustering vs not using clustering.
  This is the case of several CYS coordinating ZN. Gradients match if two
  buffers are used, and they don't if one buffer is used.
  QM only.
  """
  pdb_in = "%s.pdb"%prefix
  open(pdb_in, "w").write(pdb_str_in)
  #
  for two_buffers in [False, True]:
    print("    two_buffers=", two_buffers)
    for clustering in ["true", "false"]:
      print("      clustering=", clustering)
      cmd = " ".join([
        "qr.refine",
        pdb_in,
        "mode=opt",
        "stpmax=0.2",
        "gradient_only=true",
        "clustering=%s"%clustering,
        "dump_gradients=cluster_%s.pkl"%clustering,
        "restraints=qm",
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
    print("        min/max/mean of (gradient1 - gradient2):", \
      diff.min_max_mean().as_tuple())
    d = max(granalyse.get_grad_wdelta(ref=g1, g=g2))
    print("get_grad_wdelta:", d)
    assert d < 1.e-4

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
