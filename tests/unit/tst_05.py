from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from qrefine import conformers, restraints
from qrefine.tests.unit.tst_11 import get_restraints_manager

def cctbx_opt(model, restraints_manager, max_shift=0.2, max_iterations=25):
  from scitbx import minimizers
  from qrefine import calculator
  import scitbx.lbfgs
  s1 = model.get_sites_cart().deep_copy()
  C = calculator.sites_opt(model=model, max_shift=max_shift,
    restraints_manager=restraints_manager, shift_eval="mean")
  core_params = scitbx.lbfgs.core_parameters(stpmin=1.e-9, stpmax=max_shift)
  minimized = minimizers.lbfgs(calculator=C, mode="lbfgs", gradient_only=True,
    core_params=core_params, max_iterations=max_iterations)
  C.apply_x()
  print("Moved by:", flex.mean(flex.sqrt((s1 - C.model.get_sites_cart()).dot())))
  return C.model.get_sites_cart()

def gradients(file_name, handling, expansion, clustering):
  rm, model = get_restraints_manager(file_name, expansion=expansion,
    clustering=clustering, handling=handling)
  assert isinstance(rm, conformers.from_conformers) == \
    (handling == "conformers")
  g = rm.target_and_gradients(sites_cart=model.get_sites_cart())[1]
  return g, model

def exercise_gradients():
  """
  Conformers vs native handling, all routes: altloc atoms agree; blank atoms
  satisfy g_c * m == g_n + g_blank_only * (m - 1).
  """
  for fn in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
             "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:
    for expansion, clustering in [(False, False), (True, False), (False, True)]:
      print(fn, "expansion:", expansion, "clustering:", clustering)
      g_c, model = gradients(fn, "conformers", expansion, clustering)
      g_n, _     = gradients(fn, "native",     expansion, clustering)
      ph = model.get_hierarchy()
      cs = model.crystal_symmetry()
      assert g_c.size() == g_n.size() == ph.atoms_size()
      ci = ph.get_conformer_indices().conformer_indices
      altloc = ci != 0
      blank  = ci == 0
      assert approx_equal(g_c.select(altloc), g_n.select(altloc))
      if(not expansion and not clustering):
        assert approx_equal(g_n,
          restraints.get_cctbx_gradients(ph=ph, cs=cs).gradients)
      if(blank.count(True) == 0): continue
      m = conformers.multiplicity(ph).select(blank)
      g_blank_only = restraints.get_cctbx_gradients(
        ph=ph.select(blank), cs=cs).gradients
      assert approx_equal(g_c.select(blank) * m,
        g_n.select(blank) + g_blank_only * (m - 1))

def exercise_minimization():
  """
  h_altconf_2_complete_modified.pdb minimized once per handling (plain
  route): the two results are close.
  """
  results = {}
  for handling in ["conformers", "native"]:
    rm, model = get_restraints_manager("h_altconf_2_complete_modified.pdb",
      expansion=False, clustering=False, handling=handling)
    assert isinstance(rm, conformers.from_conformers) == \
      (handling == "conformers")
    results[handling] = cctbx_opt(model=model, restraints_manager=rm)
  d = flex.sqrt((results["conformers"] - results["native"]).dot())
  print("mean/max |conformers-native|:", flex.mean(d), flex.max(d))
  assert flex.mean(d) < 0.1 and flex.max(d) < 0.3, (flex.mean(d), flex.max(d))

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Altlocs: conformers vs native (cctbx) handling, gradients and minimization.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  exercise_gradients()
  exercise_minimization()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK time: %8.3f"%(time.time()-t0))
