from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from six.moves import cStringIO as StringIO
import mmtbx.f_model
from qrefine import qr, api, refine, conformers
from mmtbx.hydrogens.tst_fit_rotatable2 import (pdb_str, get_model,
  rotated_tyr_hh)

def get_inputs():
  """
  Model with the TYR HH rotated by 120 degrees off its true site and fmodel
  with f_obs = |Fcalc| of the true model at 1.15 A (below 1.2 A qr.refine
  fits rotatable H before refinement).
  """
  model = get_model(pdb_str)
  xrs = model.get_xray_structure().deep_copy_scatterers()
  f_obs = abs(xrs.structure_factors(d_min=1.15).f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.05)
  model.set_sites_cart(rotated_tyr_hh(model, 120))
  fmodel = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
    xray_structure=model.get_xray_structure())
  fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
  return model, fmodel

def get_params(model):
  p = qr.get_default_params()
  p.restraints = "cctbx"
  p.cluster.clustering = False
  p.expansion = False
  p.refine.number_of_weight_search_cycles = 1
  p.refine.max_iterations_weight = 2
  p.refine.number_of_refine_cycles = 1
  p.refine.max_iterations_refine = 2
  api.finalize_params(p, d_min=1.15, has_data=True,
    n_conformers=conformers.n_conformers(model.get_hierarchy()))
  return p

def h_fit_lines(max_increase):
  """
  The "Optimize H with d.o.f." block of the api.refine log, run with
  refine.H_FIT_MAX_R_WORK_INCREASE = max_increase: {label: r_work} for
  the start, fitted (only if reverted) and final lines, and whether the
  fit was reverted.
  """
  model, fmodel = get_inputs()
  log = StringIO()
  saved = refine.H_FIT_MAX_R_WORK_INCREASE
  refine.H_FIT_MAX_R_WORK_INCREASE = max_increase
  try:
    api.refine(model, get_params(model), fmodel=fmodel, working_folder=None,
      log=log)
  finally:
    refine.H_FIT_MAX_R_WORK_INCREASE = saved
  lines = log.getvalue().splitlines()
  i = lines.index("Optimize H with d.o.f.")
  r = {}
  for l in lines[i+1:]:
    if(not l.startswith("  ")): break
    for label in ["start", "fitted", "final"]:
      if(l.strip().startswith(label+":")):
        r[label] = float(l.split("r_work=")[1].split()[0])
  reverted = len([l for l in lines if "fit reverted" in l]) > 0
  return r, reverted

def test_fit_kept():
  """The fit puts HH back: R-work goes down, nothing is reverted."""
  assert refine.H_FIT_MAX_R_WORK_INCREASE == 1.e-4
  r, reverted = h_fit_lines(max_increase = refine.H_FIT_MAX_R_WORK_INCREASE)
  assert not reverted
  assert sorted(r.keys()) == ["final", "start"], r
  assert r["final"] < r["start"] - 1.e-3, r

def test_fit_reverted():
  """
  A fit that changes R-work by more than the limit (here: any fit, the
  limit is -1) is reverted: R-work is the one before the fit.
  """
  r, reverted = h_fit_lines(max_increase = -1)
  assert reverted
  assert sorted(r.keys()) == ["final", "fitted", "start"], r
  assert abs(r["fitted"] - r["start"]) > 1.e-3, r
  assert abs(r["final"] - r["start"]) < 1.5e-4, r

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine keeps the fit of rotatable H ("Optimize H with d.o.f.") only if
  it does not increase R-work by more than refine.H_FIT_MAX_R_WORK_INCREASE.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  test_fit_kept()
  test_fit_reverted()

if(__name__ == "__main__"):
  run()
  print("OK")
