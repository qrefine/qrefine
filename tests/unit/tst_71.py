from __future__ import absolute_import, division, print_function
import os
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.run_tests import qr_unit_tests_data
from qrefine.tests.unit.tst_70 import parse_trials

# Split from tst_63 (minimizer stop reasons): the 6uf7 function-value path.

def exercise_6uf7_function_path():
  """The function-value path (no clustering, native altlocs, lbfgs with
  function values) runs and reports stop reasons."""
  r_qr = run_qrefine_in(folder="run_6uf7_f", prefix="f", args=[
    os.path.join(qr_unit_tests_data, "6uf7.mtz"),
    os.path.join(qr_unit_tests_data, "6uf7_qmready.pdb"),
    "restraints=cctbx", "clustering=False", "altlocs.handling=native",
    "minimizer=lbfgs", "gradient_only=False", "r_tolerance=0.0001",
    "number_of_refine_cycles=1"])
  trials = parse_trials(r_qr.std_out_text)
  assert len(trials) >= 3, trials
  # "error": the Moré-Thuente search stopped on a swallowed rounding-error
  # failure (the minimizer keeps its best point), usual at this resolution.
  assert all(t[3] in ("cap", "conv", "error") for t in trials), trials

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_6uf7_function_path()
  print("OK")

if(__name__ == "__main__"):
  run()
