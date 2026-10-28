from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
from scitbx.array_family import flex
from qrefine import conformers
from qrefine.tests.unit import run_fmodel
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.tst_06 import r_work_of, refine_fast
from qrefine.tests.unit.tst_11 import get_restraints_manager

# data_files/altlocs.pdb: one residue, labels A and B only, no blank atoms
pdb_str = """
CRYST1   31.286   16.151   22.951  90.00  90.00  90.00 P 1
SCALE1      0.031963  0.000000  0.000000        0.00000
SCALE2      0.000000  0.061916  0.000000        0.00000
SCALE3      0.000000  0.000000  0.043571        0.00000
ATOM      1  N  AGLY A   1       5.732  10.151  12.637  1.00 16.77           N
ATOM      2  CA AGLY A   1       5.644   9.586  11.259  1.00 16.57           C
ATOM      3  C  AGLY A   1       6.751   8.596  10.960  1.00 16.16           C
ATOM      4  O  AGLY A   1       7.045   7.715  11.768  1.00 16.78           O
ATOM      5  H1 AGLY A   1       4.783  10.182  13.060  1.00 16.77           H
ATOM      6  H2 AGLY A   1       6.352   9.551  13.218  1.00 16.77           H
ATOM      7  H3 AGLY A   1       6.124  11.113  12.592  1.00 16.77           H
ATOM      8  HA2AGLY A   1       4.686   9.079  11.137  1.00 16.57           H
ATOM      9  HA3AGLY A   1       5.700  10.396  10.532  1.00 16.57           H
ATOM      1  N  BGLY A   1       5.700  10.100  12.600  1.00 16.77           N
ATOM      2  CA BGLY A   1       5.600   9.500  11.200  1.00 16.57           C
ATOM      3  C  BGLY A   1       6.700   8.500  10.900  1.00 16.16           C
ATOM      4  O  BGLY A   1       7.000   7.700  11.700  1.00 16.78           O
ATOM      5  H1 BGLY A   1       4.700  10.100  13.000  1.00 16.77           H
ATOM      6  H2 BGLY A   1       6.300   9.500  13.200  1.00 16.77           H
ATOM      7  H3 BGLY A   1       6.100  11.100  12.500  1.00 16.77           H
ATOM      8  HA2BGLY A   1       4.600   9.000  11.100  1.00 16.57           H
ATOM      9  HA3BGLY A   1       5.700  10.300  10.500  1.00 16.57           H
TER
END
"""

def exercise_routes(pdb_file):
  """
  Plain, expansion and clustering routes on a labels-only model give a
  from_conformers manager (multiplicity 1 for every atom) and the same
  gradients.
  """
  gradients = []
  for expansion, clustering in [(False, False), (True, False), (False, True)]:
    rm, model = get_restraints_manager(pdb_file, expansion=expansion,
      clustering=clustering, handling="conformers")
    assert isinstance(rm, conformers.from_conformers)
    assert rm.multiplicity.all_eq(1.0)
    sites_cart = model.get_sites_cart()
    g = rm.target_and_gradients(sites_cart=sites_cart)[1]
    assert g.size() == sites_cart.size()
    gradients.append(g)
  for g in gradients[1:]:
    diff = flex.max(flex.abs(g.as_double()-gradients[0].as_double()))
    assert diff < 1.e-6, diff

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine on a model whose atoms all carry altloc labels (no blank
  atoms): refinement works and the labels are kept.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb = os.path.abspath("%s.pdb"%prefix)
  with open(pdb, "w") as fo:
    fo.write(pdb_str)
  r_fm = run_fmodel(prefix=prefix, args=[pdb, "high_res=2"])
  mtz = os.path.abspath(r_fm.mtz)
  r_qr = run_qrefine_in(prefix=prefix+"_refine",
    args=[pdb, mtz, "restraints=cctbx", "clustering=False"] + refine_fast)
  r_work = r_work_of(r_qr.pdb, r_fm)
  print("r_work: %6.4f" % r_work)
  assert r_work < 0.05, r_work
  h = r_qr.model.get_hierarchy()
  assert sorted([k for k in h.altloc_indices().keys() if k.strip()]) == \
    ["A", "B"]
  exercise_routes(pdb_file=pdb)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
