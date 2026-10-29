from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
import random
from libtbx.utils import null_out
from scitbx.array_family import flex
from qrefine import refine, conformers
from qrefine.tests.unit.tst_11 import get_model, get_restraints_manager
from qrefine.tests.unit.tst_40 import get_params as get_mopac_params

def shaken(sites_cart, amplitude=0.05, seed=0):
  """
  sites_cart with every coordinate shifted by a uniform random number in
  [-amplitude, amplitude] (so no conformer coincides with another).
  """
  random.seed(seed)
  return flex.vec3_double([
    tuple([x + random.uniform(-amplitude, amplitude) for x in s])
    for s in sites_cart])

def finite_differences(energy, sites_cart, components, delta):
  """
  Central finite differences of energy(sites_cart) for every (i_atom, k) in
  components (k: 0, 1, 2 for x, y, z); returns a dict (i_atom, k) -> dE/dx.
  """
  result = {}
  for i, k in components:
    e = []
    for sign in (1, -1):
      s = sites_cart.deep_copy()
      x = list(s[i])
      x[k] += sign*delta
      s[i] = tuple(x)
      e.append(energy(s))
    result[(i, k)] = (e[0]-e[1])/(2*delta)
  return result

def check(fd, g, multiplicity, rel_tol, abs_tol, label):
  """
  fd[(i, k)] must agree with multiplicity[i]*g[i][k].
  """
  worst, worst_diff = 0, 0
  for (i, k), d in fd.items():
    a = multiplicity[i]*g[i][k]
    diff = abs(d - a)
    tol = max(abs_tol, rel_tol*abs(a))
    worst = max(worst, diff/tol)
    worst_diff = max(worst_diff, diff)
    assert diff <= tol, (label, i, k, d, a, diff, tol)
  print("%s: %d components, max |fd-analytical| = %.2e, max ratio to "
    "tolerance = %.3f" % (label, len(fd), worst_diff, worst))

def all_components(n_atoms):
  return [(i, k) for i in range(n_atoms) for k in range(3)]

def pick_components(model, multiplicity):
  """
  All three components of the first blank C, blank H, A N, A H, B C and B H
  atom of gly2_2 (blank atoms have multiplicity 2, the others 1).
  """
  atoms = list(model.get_hierarchy().atoms())
  picked = []
  for altloc, element in [("", "C"), ("", "H"), ("A", "N"), ("A", "H"),
                          ("B", "C"), ("B", "H")]:
    for i, a in enumerate(atoms):
      if(a.parent().altloc == altloc and a.element.strip() == element):
        picked.append(i)
        break
  assert len(picked) == 6, picked
  assert sorted(set([multiplicity[i] for i in picked])) == [1.0, 2.0]
  return [(i, k) for i in picked for k in range(3)]

def exercise_conformers(rm, model, components, delta, rel_tol, abs_tol,
                        label):
  """
  altlocs.handling=conformers: refinement uses g = rm.target_and_gradients()
  with g_i = (1/m_i) * sum over conformers c of g_c,i. The energies of the
  conformer managers are consistent with their gradients, so the finite
  differences of sum_c E_c(sites) must equal m_i*g_i. (Only this sum is
  differentiated; it is a test function, not a model energy.)
  """
  assert isinstance(rm, conformers.from_conformers)
  sites_cart = shaken(model.get_sites_cart())
  def energy(sites):
    return sum([c.restraints_manager.target_and_gradients(
      sites_cart=sites.select(c.iselection))[0] for c in rm.conformers])
  t, g = rm.target_and_gradients(sites_cart=sites_cart)
  assert t == 0 # energy undefined on this path
  fd = finite_differences(energy, sites_cart, components, delta)
  check(fd, g, rm.multiplicity, rel_tol, abs_tol, label)

def exercise_native(rm, model, delta, rel_tol, abs_tol, label):
  """
  altlocs.handling=native (cctbx): the target of rm.target_and_gradients() is
  the cctbx energy of the model with altlocs; finite differences of it must
  equal the gradients.
  """
  assert not isinstance(rm, conformers.from_conformers)
  sites_cart = shaken(model.get_sites_cart())
  def energy(sites):
    return rm.target_and_gradients(sites_cart=sites)[0]
  g = rm.target_and_gradients(sites_cart=sites_cart)[1]
  fd = finite_differences(energy, sites_cart,
    all_components(sites_cart.size()), delta)
  check(fd, g, flex.double(sites_cart.size(), 1), rel_tol, abs_tol, label)

def exercise_cctbx():
  """
  cctbx restraints, plain route, both handlings, every component of every
  atom. Models: blank atoms plus two conformers (gly2_2), conformers only
  (altlocs), a larger model with several residues in two conformers
  (h_altconf_complete).
  """
  for fn, n_conf in [("gly2_2.pdb", 2), ("altlocs.pdb", 2),
                     ("h_altconf_complete.pdb", 2)]:
    for handling in ["conformers", "native"]:
      rm, model = get_restraints_manager(fn, expansion=False,
        clustering=False, handling=handling)
      assert conformers.n_conformers(model.get_hierarchy()) == n_conf
      label = "cctbx %s %s" % (fn, handling)
      if(handling == "conformers"):
        m = rm.multiplicity
        assert set(m) == set([1.0, float(n_conf)]) or fn == "altlocs.pdb"
        exercise_conformers(rm=rm, model=model,
          components=all_components(model.get_number_of_atoms()),
          delta=1.e-4, rel_tol=1.e-5, abs_tol=1.e-4, label=label)
      else:
        exercise_native(rm=rm, model=model, delta=1.e-4, rel_tol=1.e-5,
          abs_tol=1.e-4, label=label)

def exercise_mopac():
  """
  MOPAC (PM7) restraints, plain route, altlocs.handling=conformers, on
  gly2_2 (blank residue plus a residue in two conformers): a few atoms of
  every kind (blank/A/B, heavy/hydrogen), all three components. The energies
  of different conformers are not additive as model energies, but their sum
  is a valid test function for the gradient combination.
  With the COSMO solvent used in refinement (EPS=78.4) MOPAC's energy is not
  smooth at the 0.01 kcal/mol level, so its own gradients agree with finite
  differences only to a few percent (checked loosely: a wrong combination
  would be off by a factor); in vacuum (EPS=None) they agree to 0.2 % (a
  systematic difference of MOPAC's own gradient, independent of the step).
  """
  model = get_model("gly2_2.pdb")
  params = get_mopac_params(expansion=False, clustering=False)
  rm = refine.create_restraints_manager(params=params, model=model)
  components = pick_components(model, rm.multiplicity)
  for eps, delta, rel_tol, abs_tol in [(78.4, 1.e-2, 1.e-1, 10.),
                                       (None, 2.e-3, 5.e-3, 1.e-1)]:
    for c in rm.conformers:
      c.restraints_manager.qm_engine.set(EPS=eps)
    exercise_conformers(rm=rm, model=model, components=components,
      delta=delta, rel_tol=rel_tol, abs_tol=abs_tol,
      label="mopac gly2_2.pdb EPS=%s" % eps)

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Finite-difference check of the gradients refinement uses with alternative
  conformations (refine.create_restraints_manager + target_and_gradients):
  cctbx (both handlings) and MOPAC (conformers), plain route. Expansion and
  clustering are not checked: their energies are undefined (0).
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_cctbx()
  exercise_mopac()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Total time (all tests): %6.2f"%(time.time()-t0))
  print("OK")
