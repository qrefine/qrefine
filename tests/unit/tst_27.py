from __future__ import division
from __future__ import print_function

import os
import iotbx.pdb
import libtbx.load_env

from scitbx.array_family import flex
import mmtbx.model

from qrefine import qr, refine, calculator
from libtbx.utils import null_out


qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests", "unit")


def get_model(file_name):
  file_name = os.path.join(qr_unit_tests, "data_files", file_name)
  pdb_inp = iotbx.pdb.input(file_name)

  model = mmtbx.model.manager(
    model_input=pdb_inp,
    log=null_out())

  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False

  model.process(
    make_restraints=True,
    grm_normalization=False,
    pdb_interpretation_params=params)

  return model


def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):

  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)

  model = get_model("helix.pdb")

  params = qr.get_default_params()
  params.restraints = "cctbx"
  params.expansion = False
  params.cluster.clustering = False

  #
  # Exclude GLY 87 from the restraints system.
  #
  exclude_selection = model.selection(
    string="chain A and resseq 87")

  keep_selection = ~exclude_selection

  freeze_selection = flex.bool(model.size(), False)
  freeze_selection[list(keep_selection).index(True)] = True

  assert exclude_selection.count(True) > 0
  assert keep_selection.count(True) > 0

  #
  # Restraints manager operates on the reduced model.
  #
  model_for_restraints = model.select(keep_selection)

  restraints_manager = refine.create_restraints_manager(
    params=params,
    model=model_for_restraints,
    altlocs_present=False)

  calc = calculator.sites_opt(
    model=model,
    max_shift=0.5,
    restraints_manager=restraints_manager,
    shift_eval="max",
    exclude_selection=exclude_selection,
    debug=False)

  #
  # The restraints manager must receive only the non-excluded coordinates.
  #
  f, g = calc.target_and_gradients()

  #
  # Gradient returned by sites_opt must be restored to full-model size.
  #
  assert g.size() == model.size()*3, (g.size(), model.size()*3)

  g = flex.vec3_double(g)

  assert g.size() == model.size()

  #
  # Excluded atoms are absent from restraints and therefore have zero gradient.
  #
  for gi in g.select(exclude_selection):
    assert gi == (0,0,0), gi

  #
  # At least one non-excluded atom must have a non-zero gradient.
  #
  non_zero = False
  for gi in g.select(keep_selection):
    if gi != (0,0,0):
      non_zero = True
      break

  assert non_zero

  # Combining exclude and freeze must preserve both semantics: excluded atoms
  # are absent from restraints, while frozen atoms remain in restraints and
  # are masked only after the reduced gradient is expanded.
  calc = calculator.sites_opt(
    model=model,
    max_shift=0.5,
    restraints_manager=restraints_manager,
    shift_eval="max",
    exclude_selection=exclude_selection,
    freeze_selection=freeze_selection,
    debug=False)
  _, combined_g = calc.target_and_gradients()
  combined_g = flex.vec3_double(combined_g)
  for gi in combined_g.select(exclude_selection):
    assert gi == (0,0,0), gi
  for gi in combined_g.select(freeze_selection):
    assert gi == (0,0,0), gi

  movable_selection = keep_selection & ~freeze_selection
  assert any(gi != (0,0,0) for gi in combined_g.select(movable_selection))

  print("OK")


if(__name__ == "__main__"):
  run()
