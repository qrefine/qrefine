from __future__ import absolute_import, division, print_function
import os
from scitbx.array_family import flex
from mmtbx.hydrogens.tst_xh_mapper import pdb_str
from qrefine import qr, driver
from qrefine.tests.unit import run_fmodel
from qrefine.tests.unit.tst_54 import fixture

# The 6uf7 weight-search series and function-value path are tst_70, tst_71.

def minimize_recording(c, params, mode="weight"):
  """run_minimize with the total target and the largest atom shift recorded
  at every accepted iteration (callback_after_step)."""
  fs, shifts = [], []
  prev = [c.sites_cart_nuclear.deep_copy()]
  def rec(minimizer=None):
    fs.append(c.f)
    s = c.sites_cart_nuclear
    shifts.append(flex.max(flex.sqrt((s - prev[0]).dot())))
    prev[0] = s.deep_copy()
    return False
  c.callback_after_step = rec
  driver.run_minimize(calculator=c, params=params, mode=mode)
  return fs, shifts

def function_mode_params():
  params = qr.get_default_params()
  params.refine.minimizer = "lbfgs"
  params.refine.gradient_only = False
  params.refine.max_iterations_weight = 30
  return params

def exercise_function_mode_descends(r):
  """lbfgs with gradient_only=False: no stpmax crash; the total target is
  non-increasing over the accepted iterates and ends below its start."""
  c = fixture(r).calculator()
  fs, shifts = minimize_recording(c, function_mode_params())
  assert len(fs) > 2, fs
  assert all(b <= a + 1e-9 for a, b in zip(fs, fs[1:])), fs
  assert fs[-1] < fs[0], fs

def exercise_function_mode_cap(r):
  """No accepted step moves an atom by more than refine.stpmax."""
  c = fixture(r).calculator()
  params = function_mode_params()
  params.refine.stpmax = 0.05
  fs, shifts = minimize_recording(c, params)
  assert len(shifts) > 2
  assert max(shifts) <= 0.05 + 1e-9, max(shifts)

def exercise_no_early_stop(r):
  """The sites calculator does not stop the minimizer: a gradient-only
  minimization uses its whole iteration budget and reports stop 'cap'."""
  assert not hasattr(fixture(r).calculator(), "callback_after_step")
  c = fixture(r).calculator()
  params = qr.get_default_params()
  params.refine.minimizer = "lbfgs"
  params.refine.gradient_only = True
  params.refine.max_iterations_weight = 5
  m = driver.run_minimize(calculator=c, params=params, mode="weight")
  assert m.stop_reason == "cap", m.stop_reason
  assert m.minimized.minimizer.iter() == 5, m.minimized.minimizer.iter()

def exercise_trial_line_has_stop_reason():
  """A trial cut by max_iterations says 'stop: cap' on its line."""
  from qrefine.tests.unit.run_tests import setup_helix_example
  from qrefine.tests.unit.run_tests import qr_unit_tests_data
  from qrefine.tests.unit import run_qrefine_in
  setup_helix_example()
  r_qr = run_qrefine_in(folder="run_stop", prefix="stop", args=[
    "m00_good.mtz", os.path.join(qr_unit_tests_data, "m00_poor.pdb"),
    "restraints=cctbx", "clustering=False", "max_iterations_weight=3",
    "number_of_weight_search_cycles=1", "number_of_refine_cycles=1",
    "max_iterations_refine=3"])
  lines = [l for l in r_qr.std_out_text.splitlines() if l.startswith("0:")]
  assert len(lines) == 2, lines
  assert all(" stop: cap" in l for l in lines), lines

def exercise_stop_reason_values():
  """driver.stop_reason: 'error' when the minimizer swallowed a line-search
  error, else 'cap' at the iteration budget, else 'conv'."""
  from libtbx import group_args
  assert driver.stop_reason(group_args(error=None), 50, 50) == "cap"
  assert driver.stop_reason(group_args(error=None), 12, 50) == "conv"
  assert driver.stop_reason(group_args(error="Rounding errors"), 50, 50) \
    == "error"

def exercise_real_space_cap():
  """The real-space calculator bounds the largest atom shift of a
  minimizer step by refine.stpmax (gradient-only mode)."""
  import scitbx.lbfgs
  from libtbx.utils import null_out
  from qrefine import calculator
  from qrefine.tests.unit.run_tests import setup_helix_example
  from qrefine.tests.unit.run_tests import qr_unit_tests_data
  from qrefine.tests.unit.tst_54 import get_model, get_restraints_manager
  xrs_good = setup_helix_example()[0]
  fft_map = xrs_good.structure_factors(d_min=2.0).f_calc().fft_map(
    resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  model = get_model("n_gaussian", d_min=2.0,
    file_name=os.path.join(qr_unit_tests_data, "m00_poor.pdb"))
  o = calculator.sites_real_space(
    model=model, geometry_rmsd_manager=calculator.cctbx_geometry(model),
    max_bond_rmsd=0.03, map_data=fft_map.real_map_unpadded(),
    data_weight=1.0, refine_cycles=1, stpmax=0.02, gradient_only=True,
    line_search=True, number_of_weight_search_cycles=1,
    restraints_manager=get_restraints_manager(model, handling="native"),
    max_iterations=1, log=null_out())
  o.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations=1)
  start = model.get_sites_cart().deep_copy()
  refined = o.run_one(weight=1.0)
  shift = flex.max(flex.sqrt((refined.get_sites_cart() - start).dot()))
  assert 0 < shift <= 0.02 + 1e-9, shift

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  open("in.pdb", "w").write(pdb_str)
  r = run_fmodel(args=["in.pdb", "high_res=2.5"], prefix="data")
  exercise_function_mode_descends(r)
  exercise_function_mode_cap(r)
  exercise_no_early_stop(r)
  exercise_real_space_cap()
  exercise_stop_reason_values()
  exercise_trial_line_has_stop_reason()
  print("OK")

if(__name__ == "__main__"):
  run()
