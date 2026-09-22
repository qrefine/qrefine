from __future__ import division
from __future__ import print_function

import os

import iotbx.pdb
import libtbx.load_env
import mmtbx.model
import scitbx.lbfgs

from libtbx.utils import null_out
from scitbx import minimizers
from scitbx.array_family import flex

from qrefine import calculator, qr, refine


qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests", "unit")


class recording_restraints_manager(object):
  def __init__(self, restraints_manager, expected_size, freeze_selection):
    self.restraints_manager = restraints_manager
    self.expected_size = expected_size
    self.freeze_selection = freeze_selection
    self.sizes_seen = []
    self.max_raw_frozen_gradient = 0

  def target_and_gradients(self, sites_cart):
    self.sizes_seen.append(sites_cart.size())
    assert sites_cart.size() == self.expected_size
    result = self.restraints_manager.target_and_gradients(
      sites_cart=sites_cart)
    frozen_gradients = result[1].select(self.freeze_selection)
    self.max_raw_frozen_gradient = max(
      self.max_raw_frozen_gradient,
      flex.max(flex.sqrt(frozen_gradients.dot())))
    return result


def get_model(file_name):
  file_name = os.path.join(qr_unit_tests, "data_files", file_name)
  pdb_inp = iotbx.pdb.input(file_name)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
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
  freeze_selection = model.selection(string="chain A and resseq 87")
  move_selection = ~freeze_selection
  assert freeze_selection.count(True) > 0
  assert move_selection.count(True) > 0

  # Distort both a frozen and a movable atom so the unmasked restraints
  # gradients are demonstrably non-zero in both parts of the model.
  sites_cart = model.get_sites_cart().deep_copy()
  i_frozen = list(freeze_selection).index(True)
  i_movable = list(move_selection).index(True)
  x, y, z = sites_cart[i_frozen]
  sites_cart[i_frozen] = (x+0.20, y-0.10, z+0.15)
  x, y, z = sites_cart[i_movable]
  sites_cart[i_movable] = (x-0.15, y+0.20, z-0.10)
  model.set_sites_cart(sites_cart=sites_cart)

  params = qr.get_default_params()
  assert params.refine.freeze is None
  params.restraints = "cctbx"
  params.expansion = False
  params.cluster.clustering = False
  restraints_manager = refine.create_restraints_manager(
    params=params,
    model=model,
    altlocs_present=False)
  recording_manager = recording_restraints_manager(
    restraints_manager=restraints_manager,
    expected_size=model.size(),
    freeze_selection=freeze_selection)

  sites_start = model.get_sites_cart().deep_copy()
  calc = calculator.sites_opt(
    model=model,
    max_shift=0.5,
    restraints_manager=recording_manager,
    shift_eval="max",
    freeze_selection=freeze_selection,
    debug=False)

  core_params = scitbx.lbfgs.core_parameters(
    stpmin=1.e-9,
    stpmax=0.5)
  minimizers.lbfgs(
    calculator=calc,
    mode="lbfgs",
    gradient_only=True,
    core_params=core_params,
    max_iterations=25)

  # The raw full-system restraints gradient really acts on the frozen atoms;
  # sites_opt, rather than the restraints manager, is what masks it.
  assert recording_manager.max_raw_frozen_gradient > 0

  calc.apply_x()
  sites_final = model.get_sites_cart()
  shifts = flex.sqrt((sites_final-sites_start).dot())

  assert recording_manager.sizes_seen
  assert set(recording_manager.sizes_seen) == set([model.size()])
  assert flex.max(shifts.select(freeze_selection)) < 1.e-12
  assert flex.max(shifts.select(move_selection)) > 1.e-6

  print("OK")


if __name__ == "__main__":
  run()
