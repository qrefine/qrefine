from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
from qrefine.tests.unit import run_fmodel
from qrefine.tests.unit.tst_06 import write_inputs, r_work_of
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.tst_06 import check_refined
from qrefine.tests.unit.tst_06 import refine_fast

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine with alternative conformations (altlocs.handling=conformers,
  cctbx restraints) with expansion, and with expansion plus refine.exclude on
  a residue without altlocs.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb_good, pdb_poor = write_inputs(prefix)
  r_fm = run_fmodel(prefix=prefix, args=[pdb_good, "high_res=2"])
  mtz = os.path.abspath(r_fm.mtz)
  r_start = r_work_of(pdb_poor, r_fm)
  routes = [
    ("expansion", ["expansion=True", "clustering=False"]),
    ("exclude",   ["expansion=True", "clustering=False",
                   "refine.exclude='resseq 99'"]),
  ]
  for name, args in routes:
    r_qr = run_qrefine_in(folder=name, prefix=prefix+"_"+name,
      args=[pdb_poor, mtz, "restraints=cctbx"] + args + refine_fast)
    check_refined(r_qr=r_qr, pdb_poor=pdb_poor, r_fm=r_fm, r_start=r_start)
    assert r_qr.std_out_text.find(
      "alternative conformations require minimizer=lbfgs") > -1

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
