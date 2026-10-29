from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from qrefine import qr, refine, conformers, restraints, cluster_restraints
from qrefine import fragment
import iotbx.pdb
import mmtbx.model
from cctbx import crystal
from libtbx.utils import null_out
from qrefine.tests.unit.tst_11 import qr_unit_tests
from qrefine.tests.unit.tst_40 import get_params as get_mopac_params
from qrefine.tests.unit.tst_45 import shaken, check, all_components
from qrefine.tests.unit.tst_45 import pick_components

def get_cctbx_params(expansion, clustering, handling, maxnum):
  params = qr.get_default_params()
  params.restraints = "cctbx"
  params.expansion = expansion
  params.cluster.clustering = clustering
  params.altlocs.handling = handling
  params.cluster.maxnum_residues_in_cluster = maxnum
  params.parallel.nproc = 1
  return params

def process(model):
  """
  As qr.py (and the super-sphere restraints of super_cell.py): the test
  helper of tst_11 differs (const_shrink_donor_acceptor), which matters for
  contacts between the model and its images.
  """
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params)
  return model

def closest_image_distance(xyz, box):
  """
  Smallest distance between an atom and an image of another atom in the
  P1 box (the 26 neighbouring cells).
  """
  result = None
  for i in (-1, 0, 1):
    for j in (-1, 0, 1):
      for l in (-1, 0, 1):
        if((i, j, l) == (0, 0, 0)): continue
        t = (i*box[0], j*box[1], l*box[2])
        images = xyz + flex.vec3_double(xyz.size(), t)
        for x in xyz:
          d = flex.min((images - flex.vec3_double(xyz.size(), x)).norms())
          if(result is None or d < result): result = d
  return result

def get_model(file_name, contact):
  """
  Processed model of data_files/file_name in a P1 box chosen so that the
  closest contact between the model and its images is at most contact A
  (the data files' boxes keep the images beyond the nonbonded cutoff): the
  largest of a series of boxes around the model (bounding box plus a gap)
  satisfying it. Written as <file_name>_contact<contact>.pdb in the current
  folder.
  """
  lines = open(os.path.join(qr_unit_tests, "data_files", file_name)).read()
  lines = [l for l in lines.splitlines()
           if not l.startswith("CRYST1") and not l.startswith("SCALE")]
  pdb_inp = iotbx.pdb.input(source_info=None, lines=lines)
  xyz = pdb_inp.atoms().extract_xyz()
  extent = [b-a for a, b in zip(xyz.min(), xyz.max())]
  for gap in [25.0] + [6.0 - 0.5*i for i in range(19)]:
    box = [e+gap for e in extent]
    d = closest_image_distance(xyz, box)
    if(d <= contact): break
  assert d <= contact and d > 1.0, (gap, d)
  print("%s: box gap %.1f A, closest model-image distance %.2f A" % (
    file_name, gap, d))
  cryst1 = "CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1" % tuple(box)
  file_name = "%s_contact%.1f.pdb" % (file_name[:-4], contact)
  with open(file_name, "w") as fo:
    fo.write("\n".join([cryst1] + lines) + "\n")
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  return process(model)

def supercell(model, sites_cart, n):
  """
  Processed n x n x n P1 supercell of model at sites_cart: the master copy
  first, in its atom order, then the translated copies as chains <id><k>.
  """
  ph = model.get_hierarchy().deep_copy()
  ph.atoms().set_xyz(sites_cart)
  a, b, c = model.crystal_symmetry().unit_cell().parameters()[:3]
  root = iotbx.pdb.hierarchy.root()
  mdl = iotbx.pdb.hierarchy.model()
  root.append_model(mdl)
  k = 0
  for i in range(n):
    for j in range(n):
      for l in range(n):
        for chain in ph.models()[0].chains():
          copy = chain.detached_copy()
          if(k > 0): copy.id = "%s%d" % (chain.id, k)
          copy.atoms().set_xyz(copy.atoms().extract_xyz() +
            flex.vec3_double(copy.atoms_size(), (i*a, j*b, l*c)))
          mdl.append_chain(copy)
        k += 1
  root.atoms_reset_serial()
  root.atoms().reset_i_seq()
  cs = crystal.symmetry(unit_cell=(n*a, n*b, n*c, 90, 90, 90),
    space_group_symbol="P1")
  return process(mmtbx.model.manager(model_input=None, pdb_hierarchy=root,
    crystal_symmetry=cs, log=null_out()))

def check_supercell(route, model, sites_cart, label):
  """
  Independent reference for a cctbx expansion or clustering manager built
  for model: the gradients of the master copy of a 2x2x2 supercell (plain
  cctbx restraints, every copy shaken alike) must equal the route's. The
  supercell is large enough for cctbx's own symmetry handling, and the
  frozen-environment gradients are exact for crystal contacts of distinct
  atoms (an atom never meets its own image here).
  """
  n_at = sites_cart.size()
  sc = supercell(model, sites_cart, 2)
  g_ref = sc.get_restraints_manager().energies_sites(
    sites_cart=sc.get_sites_cart(), compute_gradients=True).gradients[:n_at]
  g = route.target_and_gradients(sites_cart=sites_cart)[1]
  diff = flex.max((g-g_ref).norms())
  print("%s: max |route - supercell| = %.2e" % (label, diff))
  assert diff < 1.e-6, diff

def moved(sites_cart, j, k, delta):
  s = sites_cart.deep_copy()
  p = list(s[j])
  p[k] += delta
  s[j] = tuple(p)
  return s

class frozen_environment(object):
  """
  Energies of an expansion (restraints.from_expansion) or clustering
  (cluster_restraints.from_cluster) manager with one atom of its model moved
  and everything else (symmetry copies, buffer atoms, caps) held where the
  last target_and_gradients() call of the manager put them. The gradient
  these routes return for an atom is the partial derivative of that energy:
  for expansion the energy of the completed super sphere, for clustering the
  energy of the fragment the atom is a core atom of.
  The master copy is the first atoms of the super sphere (checked against
  sites_cart by the caller).
  """
  def __init__(self, route):
    self.route = route
    if(isinstance(route, restraints.from_expansion)):
      self.kind = "expansion"
      self.sites_super = \
        route.pdb_hierarchy_super_completed.atoms().extract_xyz()
      self.master = route.selection.iselection()
    else:
      assert isinstance(route, cluster_restraints.from_cluster)
      self.kind = "cluster"
      fm = route.fragment_manager
      self.sites_super = fm.pdb_hierarchy_super.atoms().extract_xyz()
      n = fm.pdb_hierarchy.atoms_size()
      self.master = flex.size_t(range(n))
      # QM fragments take their coordinates from this hierarchy
      assert route.restraints_manager.fragment_extracts.pdb_hierarchy_super \
        is fm.pdb_hierarchy_super
      # every model atom is a core atom of exactly one fragment
      self.fragment_index = {}
      for index, sel in enumerate(fm.fragment_selections):
        core = sel & ~fm.buffer_selections[index]
        for j in core.iselection():
          assert j not in self.fragment_index
          self.fragment_index[j] = index
      assert len(self.fragment_index) == n
      self.n_fragments = len(fm.fragment_selections)
      self.qm = isinstance(route.restraints_manager, restraints.from_qm)
      if(self.qm):
        self.cap_sites = dict(
          [(index, self._cap_sites(index)) for index in range(self.n_fragments)])

  def _cap_sites(self, index):
    fe = self.route.restraints_manager.fragment_extracts
    _, ph, fragment_selection = fragment.get_qm_file_name_and_pdb_hierarchy(
      fragment_extracts=fe, index=index)
    return ph.atoms().extract_xyz().select(~fragment_selection)

  def model_sites(self):
    return self.sites_super.select(self.master)

  def n_neighbours(self):
    """
    Atoms of the super sphere besides the master copy (and caps).
    """
    return self.sites_super.size() - self.master.size()

  def energy(self, j, k, delta):
    """
    Energy with atom j of the model moved by delta along k. None if the
    move re-places a capping hydrogen (QM clustering; the gradient does not
    account for that).
    """
    sites = moved(self.sites_super, self.master[j], k, delta)
    if(self.kind == "expansion"):
      return self.route.restraints_manager.target_and_gradients(
        sites_cart=sites)[0]
    fm = self.route.fragment_manager
    index = self.fragment_index[j]
    fm.pdb_hierarchy_super.atoms().set_xyz(sites)
    try:
      if(self.qm and self.cap_sites[index].size() > 0):
        d = flex.max((self._cap_sites(index) - self.cap_sites[index]).norms())
        if(d > 1.e-9): return None
      return self.route.restraints_manager(
        [fm.fragment_selections[index], sites, index])[0]
    finally:
      fm.pdb_hierarchy_super.atoms().set_xyz(self.sites_super)

def exercise(rm, model, components, delta, rel_tol, abs_tol, label,
             eps=78.4):
  """
  Finite differences of the frozen-environment energies versus the gradient
  rm.target_and_gradients() returns (the one refinement uses). For
  altlocs.handling=conformers rm is a from_conformers manager: the test
  function of atom i is the sum over the conformers containing i, and the
  finite differences must equal multiplicity times the combined gradient.
  """
  sites_cart = shaken(model.get_sites_cart())
  if(isinstance(rm, conformers.from_conformers)):
    routes = [(c.restraints_manager, c.iselection) for c in rm.conformers]
    multiplicity = rm.multiplicity
  else:
    routes = [(rm, flex.size_t(range(sites_cart.size())))]
    multiplicity = flex.double(sites_cart.size(), 1)
  for route, _ in routes:
    if(isinstance(route.restraints_manager, restraints.from_qm)):
      route.restraints_manager.qm_engine.set(EPS=eps)
  g = rm.target_and_gradients(sites_cart=sites_cart)[1]
  envs = []
  for route, isel in routes:
    env = frozen_environment(route)
    assert approx_equal(env.model_sites(), sites_cart.select(isel), 1.e-9)
    assert env.n_neighbours() > 0 # images within select_within_radius
    envs.append((env, dict([(i, pos) for pos, i in enumerate(isel)])))
  fd = {}
  skipped = []
  for i, k in components:
    d = 0
    for env, position in envs:
      if(i not in position): continue
      e = [env.energy(position[i], k, sign*delta) for sign in (1, -1)]
      if(None in e):
        d = None
        break
      d += (e[0]-e[1])/(2*delta)
    if(d is None): skipped.append((i, k))
    else:          fd[(i, k)] = d
  if(len(skipped) > 0):
    print("%s: %d components skipped (move a cap): %s" % (
      label, len(skipped), skipped))
  assert len(fd) >= len(components)//2, (len(fd), len(components))
  print("%s: %d neighbours" % (label, envs[0][0].n_neighbours()))
  check(fd, g, multiplicity, rel_tol, abs_tol, label)
  # the routes write sites_cart into the model's hierarchy: callers reuse
  # these sites rather than shaking the model again
  return envs, sites_cart

def exercise_cctbx():
  """
  cctbx restraints, expansion and clustering (several fragments), both
  handlings, every component of every atom, boxes small enough for an
  image to be within 2 A of the model (cctbx nonbonded energy is zero
  beyond the van der Waals distance); then the gradients against an
  explicit supercell.
  """
  for fn, maxnum in [("gly2_2.pdb", 1), ("h_altconf_complete.pdb", 3)]:
    for expansion, clustering in [(True, False), (False, True)]:
      for handling in ["conformers", "native"]:
        model = get_model(fn, contact=2.0)
        params = get_cctbx_params(expansion=expansion, clustering=clustering,
          handling=handling, maxnum=maxnum)
        rm = refine.create_restraints_manager(params=params, model=model)
        assert isinstance(rm, conformers.from_conformers) == \
          (handling == "conformers")
        label = "cctbx %s %s %s" % (fn,
          "expansion" if expansion else "clustering", handling)
        envs, sites_cart = exercise(rm=rm, model=model,
          components=all_components(model.get_number_of_atoms()),
          delta=1.e-4, rel_tol=1.e-5, abs_tol=1.e-4, label=label)
        if(clustering):
          for env, _ in envs: assert env.n_fragments > 1, env.n_fragments
        # the images do contribute: the same model with its images beyond
        # select_within_radius gives different gradients
        model_far = get_model(fn, contact=30.0) # images beyond the radius
        rm_far = refine.create_restraints_manager(params=params,
          model=model_far)
        diff = flex.max(
          (rm.target_and_gradients(sites_cart=sites_cart)[1] -
           rm_far.target_and_gradients(sites_cart=sites_cart)[1]).norms())
        assert diff > 1., diff
        # independent reference: explicit supercell, per conformer
        for env, _ in envs:
          route = env.route
          if(handling == "conformers"):
            c = [c for c in rm.conformers if c.restraints_manager is route][0]
            check_supercell(route=route, model=c.model,
              sites_cart=sites_cart.select(c.iselection),
              label=label + " conformer " + c.label)
          else:
            check_supercell(route=route, model=model,
              sites_cart=sites_cart, label=label)

def exercise_mopac():
  """
  MOPAC (PM7, vacuum: see tst_45), altlocs.handling=conformers, expansion
  and clustering (two fragments) on gly2_2 in a small box, a few atoms of
  every kind.
  """
  for expansion, clustering in [(True, False), (False, True)]:
    model = get_model("gly2_2.pdb", contact=4.0) # no clashes for the SCF
    params = get_mopac_params(expansion=expansion, clustering=clustering)
    params.cluster.maxnum_residues_in_cluster = 1
    # a few dozen image atoms within 5 A of the 17-atom conformers (hundreds
    # within 10 A: minutes per MOPAC energy)
    params.cluster.select_within_radius = 5
    rm = refine.create_restraints_manager(params=params, model=model)
    assert isinstance(rm, conformers.from_conformers)
    label = "mopac gly2_2.pdb %s" % (
      "expansion" if expansion else "clustering")
    envs, _ = exercise(rm=rm, model=model,
      components=pick_components(model, rm.multiplicity),
      delta=2.e-3, rel_tol=5.e-3, abs_tol=1.e-1, label=label, eps=None)
    if(clustering):
      for env, _ in envs: assert env.n_fragments == 2, env.n_fragments

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Finite-difference check of the gradients refinement uses with alternative
  conformations on the expansion and clustering routes. These routes have
  no model energy: expansion returns the gradient of the super-sphere energy
  on the master atoms with the symmetry copies and caps frozen, clustering
  the gradient of each atom's fragment energy with buffer and caps frozen.
  Both are partial derivatives of energies the managers hold, which are
  differentiated here with the environment frozen. tst_45 checks the plain
  route.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  t0 = time.time()
  exercise_cctbx()
  print("cctbx: %.1f s" % (time.time()-t0))
  t0 = time.time()
  exercise_mopac()
  print("mopac: %.1f s" % (time.time()-t0))

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Total time (all tests): %6.2f"%(time.time()-t0))
  print("OK")
