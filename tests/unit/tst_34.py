from __future__ import division
from __future__ import absolute_import

import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import Sorry, null_out
from qrefine import qr, refine
from qrefine.tests.unit import run_qrefine

# Three CYS and a HIS coordinating ZN in a tiny P1 cell: symmetry copies of the
# molecule are only 3.9 A away from it.
pdb_str_small_box = """
CRYST1   12.965   12.703   10.335  90.00  90.00  90.00 P 1
SCALE1      0.077131  0.000000  0.000000        0.00000
SCALE2      0.000000  0.078722  0.000000        0.00000
SCALE3      0.000000  0.000000  0.096759        0.00000
ATOM      1  N   CYS A   3       2.436   4.407   5.591  1.00 76.61           N
ATOM      2  CA  CYS A   3       3.234   3.740   6.598  1.00 78.13           C
ATOM      3  C   CYS A   3       3.081   4.465   7.933  1.00 83.42           C
ATOM      4  O   CYS A   3       2.148   5.248   8.132  1.00 87.64           O
ATOM      5  CB  CYS A   3       4.684   3.633   6.124  1.00 78.38           C
ATOM      6  SG  CYS A   3       5.747   5.078   6.261  1.00 63.83           S
ATOM      7  OXT CYS A   3       3.897   4.273   8.835  1.00 83.42           O
ATOM      8  H   CYS A   3       2.719   5.190   5.377  0.00 76.61           H
ATOM      9  H2  CYS A   3       2.429   3.860   4.754  0.00 76.61           H
ATOM     10  H3  CYS A   3       1.500   4.516   5.925  0.00 76.61           H
ATOM     11  HA  CYS A   3       2.903   2.715   6.763  0.00 78.13           H
ATOM     12  HB2 CYS A   3       5.161   2.842   6.703  0.00 78.38           H
ATOM     13  HB3 CYS A   3       4.665   3.362   5.068  0.00 78.38           H
ATOM     14  N   CYS A   8       4.348  10.296   6.095  1.00 96.46           N
ATOM     15  CA  CYS A   8       4.385   9.291   5.050  1.00 84.79           C
ATOM     16  C   CYS A   8       4.041   9.984   3.764  1.00 88.83           C
ATOM     17  O   CYS A   8       2.977  10.596   3.649  1.00 93.29           O
ATOM     18  CB  CYS A   8       3.416   8.138   5.308  1.00 85.00           C
ATOM     19  SG  CYS A   8       3.177   7.021   3.862  1.00 78.50           S
ATOM     20  OXT CYS A   8       4.824   9.944   2.815  1.00 88.83           O
ATOM     21  H   CYS A   8       3.699  10.221   6.654  0.00 96.46           H
ATOM     22  H2  CYS A   8       5.191  10.251   6.631  0.00 96.46           H
ATOM     23  H3  CYS A   8       4.264  11.203   5.683  0.00 96.46           H
ATOM     24  HA  CYS A   8       5.371   8.828   5.002  0.00 84.79           H
ATOM     25  HB2 CYS A   8       3.801   7.536   6.132  0.00 85.00           H
ATOM     26  HB3 CYS A   8       2.443   8.552   5.571  0.00 85.00           H
TER
ATOM     27  N   HIS A  30       7.597   2.632   5.103  1.00 66.84           N
ATOM     28  CA  HIS A  30       8.236   2.676   3.794  1.00 70.57           C
ATOM     29  C   HIS A  30       9.677   3.144   3.940  1.00 70.51           C
ATOM     30  O   HIS A  30       9.958   4.090   4.665  1.00 73.92           O
ATOM     31  CB  HIS A  30       7.546   3.604   2.802  1.00 70.79           C
ATOM     32  CG  HIS A  30       6.084   3.361   2.605  1.00 68.99           C
ATOM     33  ND1 HIS A  30       5.110   4.191   3.109  1.00 77.35           N
ATOM     34  CD2 HIS A  30       5.437   2.354   1.978  1.00 70.51           C
ATOM     35  CE1 HIS A  30       3.919   3.730   2.757  1.00 79.15           C
ATOM     36  NE2 HIS A  30       4.092   2.619   2.062  1.00 70.69           N
ATOM     37  OXT HIS A  30      10.576   2.570   3.325  1.00 70.51           O
ATOM     38  H   HIS A  30       6.923   3.155   5.208  0.00 66.84           H
ATOM     39  H2  HIS A  30       7.262   1.706   5.280  0.00 66.84           H
ATOM     40  H3  HIS A  30       8.260   2.887   5.807  0.00 66.84           H
ATOM     41  HA  HIS A  30       8.177   1.669   3.382  0.00 70.57           H
ATOM     42  HB2 HIS A  30       7.660   4.629   3.155  0.00 70.79           H
ATOM     43  HB3 HIS A  30       8.028   3.488   1.831  0.00 70.79           H
ATOM     44  HD2 HIS A  30       5.893   1.500   1.500  0.00 70.51           H
ATOM     45  HE1 HIS A  30       2.969   4.184   2.997  0.00 79.15           H
ATOM     46  HE2 HIS A  30       3.349   2.052   1.655  0.00 70.69           H
ATOM     47  N   CYS A  33      10.084   6.928   3.037  1.00 77.84           N
ATOM     48  CA  CYS A  33       9.595   7.794   4.113  1.00 77.49           C
ATOM     49  C   CYS A  33      10.549   7.895   5.296  1.00 82.56           C
ATOM     50  O   CYS A  33      10.404   8.793   6.132  1.00 86.85           O
ATOM     51  CB  CYS A  33       8.272   7.269   4.646  1.00 74.32           C
ATOM     52  SG  CYS A  33       7.022   7.099   3.426  1.00 77.49           S
ATOM     53  OXT CYS A  33      11.465   7.082   5.419  1.00 82.56           O
ATOM     54  H   CYS A  33       9.545   6.288   2.840  0.00 77.84           H
ATOM     55  H2  CYS A  33      10.220   7.470   2.208  0.00 77.84           H
ATOM     56  H3  CYS A  33      10.952   6.514   3.310  0.00 77.84           H
ATOM     57  HA  CYS A  33       9.481   8.789   3.683  0.00 77.49           H
ATOM     58  HB2 CYS A  33       8.439   6.286   5.087  0.00 74.32           H
ATOM     59  HB3 CYS A  33       7.903   7.960   5.404  0.00 74.32           H
TER
HETATM   60 ZN    ZN A 101       5.256   5.881   4.046  1.00 76.89          Zn
END
"""

def exercise_check(pdb_file):
  """
  The check itself. It raises for a crystal with neighbours in range and is
  silent when the neighbours are seen (clustering or expansion), when there
  are none in range, or for real-space refinement where the lattice is not
  physical.
  """
  pdb_inp = iotbx.pdb.input(file_name = pdb_file)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  def get_params(clustering, expansion, restraints="qm", radius=10.):
    params = qr.get_default_params()
    params.restraints = restraints
    params.cluster.clustering = clustering
    params.expansion = expansion
    params.cluster.select_within_radius = radius
    return params
  try:
    refine.check_qm_sees_crystal_neighbours(
      model    = model,
      params   = get_params(clustering=False, expansion=False),
      map_data = None)
  except Sorry as e:
    assert str(e).find("expansion=True") > -1, str(e)
  else: raise AssertionError("Expected Sorry")
  for params, map_data in [
      (get_params(clustering=True,  expansion=False), None),
      (get_params(clustering=False, expansion=True),  None),
      (get_params(clustering=False, expansion=False), 1), # real-space
      (get_params(clustering=False, expansion=False, radius=1.0), None),
      (get_params(clustering=False, expansion=False, restraints="cctbx"), None)]:
    refine.check_qm_sees_crystal_neighbours(
      model = model, params = params, map_data = map_data)

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  restraints=qm with clustering=False and expansion=False hands the bare
  asymmetric unit to the QM engine, so crystal neighbours are ignored. qr.refine
  must refuse this when neighbours are within cluster.select_within_radius.
  Dry runs: no QM calculations.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  pdb_small_box = "%s_small_box.pdb"%prefix
  with open(pdb_small_box, "w") as fo:
    fo.write(pdb_str_small_box)
  #
  pdb_big_box = "%s_big_box.pdb"%prefix
  with open(pdb_big_box, "w") as fo:
    for l in pdb_str_small_box.splitlines():
      if(l.startswith("SCALE")): continue
      if(l.startswith("CRYST1")):
        l = "CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1"
      fo.write(l+"\n")
  #
  exercise_check(pdb_file = pdb_small_box)
  #
  common = ["mode=opt", "dry_run=True", "restraints=qm",
            "quantum.engine_name=mopac"]
  # Neighbours in range, QM engine would only see the asymmetric unit: stops.
  r_qr = run_qrefine(prefix = prefix, args=[pdb_small_box,
    "clustering=False", "expansion=False"]+common, sorry_expected=True)
  assert r_qr.std_out_text.find("expansion=True") > -1, r_qr.std_out_text
  # Neighbours are part of the fragments: runs.
  r_qr = run_qrefine(prefix = prefix, args=[pdb_small_box,
    "clustering=True"]+common)
  # No neighbours in range: runs.
  r_qr = run_qrefine(prefix = prefix, args=[pdb_big_box,
    "clustering=False", "expansion=False"]+common)

if(__name__ == "__main__"):
  run()
