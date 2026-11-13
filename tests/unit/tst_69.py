from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.tst_06 import refine_fast
from qrefine.tests.unit.tst_68 import rmsd
from qrefine.tests.unit.run_tests import qr_unit_tests_data

def r_work(cli, which):
  """R-work printed on the 'Best Rwork, Rfree (<which>)' line."""
  key = "Best Rwork, Rfree (%s):" % which
  lines = [l for l in cli.std_out_text.splitlines() if l.startswith(key)]
  assert len(lines) > 0, key
  return float(lines[-1].split(":")[1].split()[0])

def exercise_cctbx():
  """
  6uf7 (real data, ILE 4 A/B, ring closed through a symmetry mate), cctbx
  restraints: (a) native, (b) conformers, (c) conformers with clustering.
  """
  mtz = os.path.join(qr_unit_tests_data, "6uf7.mtz")
  pdb = os.path.join(qr_unit_tests_data, "6uf7_qmready.pdb")
  runs = {}
  for tag, extra in [
      ("a", ["clustering=False", "altlocs.handling=native"]),
      ("b", ["clustering=False", "altlocs.handling=conformers"]),
      ("c", ["clustering=True", "altlocs.handling=conformers",
             "maxnum_residues_in_cluster=2", "parallel.nproc=1"])]:
    t0 = time.time()
    r = run_qrefine_in(folder="cctbx_"+tag, prefix="cctbx_"+tag,
      args=[mtz, pdb, "restraints=cctbx"] + extra + refine_fast)
    r_start = r_work(r, "at refinement start")
    r_end = r_work(r, "after refinement")
    print("(%s) Rwork start %.4f end %.4f  wall %.1f s" % (
      tag, r_start, r_end, time.time()-t0))
    assert r_end <= r_start + 0.002, (tag, r_start, r_end)
    runs[tag] = r
  assert runs["c"].std_out_text.find("per-fragment variants") > -1
  assert runs["b"].std_out_text.find("per-fragment variants") == -1
  assert not os.path.exists(os.path.join(runs["c"].folder, "ase", "A"))
  d_bc = rmsd(runs["b"].model, runs["c"].model)
  d_ac = rmsd(runs["a"].model, runs["c"].model)
  print("rmsd(b,c): %8.5f rmsd(a,c): %8.5f" % (d_bc, d_ac))
  assert d_bc < 0.05, d_bc
  assert d_ac < 0.2, d_ac
  # non-vacuous: refinement moved the atoms in every run
  start = mmtbx.model.manager(model_input = iotbx.pdb.input(file_name=pdb),
    log = null_out())
  for tag in sorted(runs):
    d = rmsd(start, runs[tag].model)
    print("(%s) rmsd(start, end): %8.5f" % (tag, d))
    assert d > 1.e-4, (tag, d)

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  """
  6uf7, real data: the clustered altloc route (per-fragment conformer
  variants) against whole-model routes with cctbx restraints. MOPAC with
  altlocs and clustering is tst_67 (small model): on 6uf7 an optimization
  (10 PM7 jobs of 0.3-5 s per evaluation, many evaluations) is too long for
  the suite and is a manual acceptance run.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_cctbx()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
