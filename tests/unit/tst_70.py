from __future__ import absolute_import, division, print_function
import os
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.run_tests import qr_unit_tests_data

# Split from tst_63 (minimizer stop reasons): the 6uf7 weight-search series.

def parse_trials(text):
  """[(scale, r_free, n_fev, stop)] of the 'Optimal weight search' block."""
  import re
  trials = []
  inside = False
  for l in text.splitlines():
    if(l.startswith("Optimal weight search:")): inside = True
    if(l.startswith("Weight search:")): inside = False
    if(not inside): continue
    m = re.match(
      r"^\d+: Rw: [\d.]+ Rf: ([\d.]+).* rws: (\S+) n_fev: (\d+) stop: (\w+)",
      l)
    if(m is not None):
      trials.append((float(m.group(2)), float(m.group(1)), int(m.group(3)),
        m.group(4)))
  return trials

def assert_unimodal(rf, tolerance=0.0003):
  """Before the minimum no rise, after it no fall, beyond tolerance."""
  k = rf.index(min(rf))
  assert all(rf[i] >= rf[i+1] - tolerance for i in range(k)), rf
  assert all(rf[i] <= rf[i+1] + tolerance for i in range(k, len(rf)-1)), rf

def exercise_6uf7_series():
  """6uf7 with cctbx restraints (gradient-only path: A/B conformers): every
  trial reports its stop reason and the R-free series is unimodal."""
  r_qr = run_qrefine_in(folder="run_6uf7", prefix="w", args=[
    os.path.join(qr_unit_tests_data, "6uf7.mtz"),
    os.path.join(qr_unit_tests_data, "6uf7_qmready.pdb"),
    "restraints=cctbx", "r_tolerance=0.0001", "number_of_refine_cycles=1"])
  trials = parse_trials(r_qr.std_out_text)
  assert len(trials) >= 3, trials
  assert all(t[3] == "cap" for t in trials), trials
  assert_unimodal([t[1] for t in trials])

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_6uf7_series()
  print("OK")

if(__name__ == "__main__"):
  run()
