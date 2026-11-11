from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import iotbx.pdb
import mmtbx.f_model
import mmtbx.model
from cctbx import xray
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from scitbx.array_family import flex
from mmtbx.hydrogens import xh_mapper
from mmtbx.hydrogens.tst_xh_mapper import pdb_str
from qrefine import calculator, qr, refine
from qrefine.tests.unit import run_fmodel

def process(model):
  """Process the model as qr.py does."""
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params=params)
  return model

def get_model(scattering_table, d_min, file_name="in.pdb"):
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(file_name=file_name), log=null_out())
  process(model)
  model.setup_scattering_dictionaries(scattering_table=scattering_table,
    d_min=d_min)
  return model

def get_restraints_manager(model, handling="conformers"):
  """
  handling="native": cctbx restraints on the whole model, with a target
  value (the conformers route returns target 0, gradients only).
  """
  params = qr.get_default_params()
  params.restraints = "cctbx"
  params.expansion = False
  params.cluster.clustering = False
  params.altlocs.handling = handling
  params.parallel.nproc = 1
  return refine.create_restraints_manager(params=params, model=model,
    log=null_out())

def get_fmodel(model, r):
  """fmodel with exact (direct) structure factors and site gradient flags."""
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(f_obs=r.f_obs, r_free_flags=r.flags,
    xray_structure=model.get_xray_structure().deep_copy_scatterers(),
    target_name="ls_wunit_k1", sf_and_grads_accuracy_params=params)
  fmodel.xray_structure.scatterers().flags_set_grads(state=False)
  xray.set_scatterer_grad_flags(
    scatterers=fmodel.xray_structure.scatterers(), site=True)
  return fmodel

class fixture(object):
  """
  The tst_xh_mapper model with the n_gaussian table: active mapper, nuclear
  start sites (expanded), fmodel holding the shrunk sites.
  """
  def __init__(self, r, exclude_selection=None):
    self.model = get_model("n_gaussian", d_min=r.f_obs.d_min())
    self.mp = xh_mapper.mapper(self.model)
    assert not self.mp.is_identity
    self.nuclear = self.mp.expand(self.model.get_sites_cart())
    self.fmodel = get_fmodel(self.model, r)
    self.fmodel.xray_structure.set_sites_cart(self.mp.shrink(self.nuclear))
    self.fmodel.update_xray_structure(update_f_calc=True)
    self.grm = calculator.cctbx_geometry(self.model)
    self.exclude_selection = exclude_selection
    model_for_restraints = self.model
    if exclude_selection is not None:
      model_for_restraints = self.model.select(~exclude_selection)
    # The target of the native route makes the FD of the total possible.
    self.rm = get_restraints_manager(model_for_restraints, handling="native")
    assert self.rm.target_and_gradients(
      sites_cart=model_for_restraints.get_sites_cart())[0] > 0

  def calculator(self, **kwargs):
    c = calculator.sites(fmodel=self.fmodel, geometry_rmsd_manager=self.grm,
      max_bond_rmsd=0.03, max_shift=0.2, restraints_manager=self.rm,
      exclude_selection=self.exclude_selection, **kwargs)
    c.setw(data_weight=1.0, restraints_weight=1.0,
      restraints_weight_scale=1.0)
    return c

def value_at(c, sites_cart):
  """Total target of calculator c at the nuclear sites sites_cart."""
  c.update(x=(sites_cart - c.sites_cart).as_double())
  return c.f

def fd(func, sites, i, k, eps):
  values = []
  for sign in [1, -1]:
    s = sites.deep_copy()
    x = list(s[i])
    x[k] += sign * eps
    s[i] = tuple(x)
    values.append(func(s))
  return (values[0] - values[1]) / (2 * eps)

def check_fd(c, mp, nuclear, h_indices):
  """
  FD of the calculator target on the H in h_indices and their parents:
  the total target (data and restraints, weights 1) and the data term alone.
  The data gradient on H is small (~1e-5) next to the restraint gradient, so
  the data term is checked on its own with a relative tolerance, and the
  unmapped data gradient must fail that check.
  """
  atoms = set()
  for i_h, i_x in zip(mp.i_h, mp.i_x):
    if i_h in h_indices: atoms.update([i_h, i_x])
  assert len(atoms) == 2*len(h_indices), atoms
  x0 = flex.double(nuclear.size()*3, 0)
  # Total target. eps=1e-5: the restraint term has large third derivatives
  # at H (central differences with 1e-4 are off by ~1e-4 relative).
  c.setw(data_weight=1.0, restraints_weight=1.0)
  c.update(x=x0)
  g = c.g.deep_copy()
  for i in sorted(atoms):
    for k in range(3):
      gp = fd(lambda s: value_at(c, s), nuclear, i, k, 1.e-5)
      assert approx_equal(g[3*i+k], gp, max(1.e-5*abs(gp), 1.e-6)), (
        i, k, g[3*i+k], gp)
  # Data term alone.
  c.setw(data_weight=1.0, restraints_weight=0.0)
  c.update(x=x0)
  g_data = c.g.deep_copy()
  dg = flex.vec3_double(c.x_target_functor(compute_gradients=True).\
    gradients_wrt_atomic_parameters(site=True).packed()).as_double()
  n_unmapped_off = 0
  for i in sorted(atoms):
    for k in range(3):
      gp = fd(lambda s: value_at(c, s), nuclear, i, k, 1.e-4)
      tol = max(1.e-5*abs(gp), 1.e-12)
      assert approx_equal(g_data[3*i+k], gp, tol), (i, k, g_data[3*i+k], gp)
      if abs(dg[3*i+k] - gp) > 100*tol: n_unmapped_off += 1
  assert n_unmapped_off > 0
  c.setw(data_weight=1.0, restraints_weight=1.0)
  c.update(x=x0)
  return g

# H of different kinds, each with its own parent: N-terminal NH3 H (chain
# start), backbone HA, the capping H (occupancy 0) and side-chain CH2 H on
# the residue truncated at CB, the peptide D, the altloc hydroxyl H and a
# water H (separate molecule, end of the model).
diverse_h_labels = [("ARG", "H1"), ("ARG", "HA"), ("ARG", "HCG"),
  ("SER", "D"), ("SER", "HB2"), ("SER", "HG"), ("HOH", "H1")]

def diverse_h(f, labels=diverse_h_labels, excluded=()):
  """
  Indices of mapped H matching labels (resname, name), in that order, each
  with a parent not used by an earlier one (check_fd needs distinct
  parents); of an altloc pair only the first. All labels must be found.
  """
  atoms = f.model.get_hierarchy().atoms()
  parent = dict(zip(f.mp.i_h, f.mp.i_x))
  result, parents = [], set()
  for resname, name in labels:
    for i_h in f.mp.i_h:
      a = atoms[i_h]
      if((a.parent().resname, a.name.strip()) != (resname, name)): continue
      if(i_h in excluded or parent[i_h] in parents): continue
      result.append(i_h)
      parents.add(parent[i_h])
      break
    else:
      raise AssertionError("no mapped H %s %s" % (resname, name))
  return result

def check_invariant(c, mp):
  assert approx_equal(c.fmodel.xray_structure.sites_cart(),
    mp.shrink(c.sites_cart_nuclear), 1.e-9)

def test_fd_xray(r):
  """
  With an active mapper the calculator gradient is the derivative of its
  target with respect to the nuclear sites, and fmodel holds the shrunk
  sites.
  """
  f = fixture(r)
  c = f.calculator(sites_cart=f.nuclear, mapper=f.mp)
  assert c.mapper is f.mp
  # The given nuclear sites are kept as they are.
  assert c.sites_cart_nuclear.rms_difference(f.nuclear) == 0
  check_invariant(c, f.mp)
  check_fd(c, f.mp, f.nuclear, diverse_h(f))
  check_invariant(c, f.mp)
  # Required nuclear sites with an active mapper.
  try:
    f.calculator(mapper=f.mp)
  except AssertionError:
    pass
  else:
    raise AssertionError("sites_cart must be required")

def test_set_sites_round_trip(r):
  f = fixture(r)
  c = f.calculator(sites_cart=f.nuclear, mapper=f.mp)
  # Every C-H stretched by 0.1 A.
  nuclear2 = f.nuclear.deep_copy()
  elements = [e.strip() for e in f.model.get_hierarchy().atoms().
    extract_element()]
  n_ch = 0
  for i_h, i_x in zip(f.mp.i_h, f.mp.i_x):
    if elements[i_x] != "C": continue
    v = flex.vec3_double([f.nuclear[i_h]]) - flex.vec3_double(
      [f.nuclear[i_x]])
    u = v * (1./v.norms()[0])
    nuclear2[i_h] = (flex.vec3_double([f.nuclear[i_h]]) + u*0.1)[0]
    n_ch += 1
  assert n_ch > 0
  c.set_sites(nuclear2)
  assert c.sites_cart_nuclear.rms_difference(nuclear2) == 0
  assert c.sites_cart.rms_difference(nuclear2) == 0
  assert c.x.all_eq(0)
  # fmodel stores fractional coordinates: equal up to round-off.
  assert c.fmodel.xray_structure.sites_cart().rms_difference(
    f.mp.shrink(nuclear2)) < 1.e-12
  # A stretch along the bond is invisible to fmodel: the nuclear sites are
  # owned by the calculator and cannot be read back from fmodel.
  assert nuclear2.rms_difference(f.nuclear) > 0.01
  assert c.fmodel.xray_structure.sites_cart().rms_difference(
    f.mp.shrink(f.nuclear)) < 1.e-9
  # f_calc follows the sites.
  fc = c.fmodel.f_calc().deep_copy()
  check = c.fmodel.deep_copy()
  check.update_xray_structure(update_f_calc=True)
  assert approx_equal(fc.data(), check.f_calc().data(), 1.e-9)
  # reset_fmodel keeps the nuclear sites.
  c.reset_fmodel(c.fmodel.deep_copy())
  assert c.sites_cart_nuclear.rms_difference(nuclear2) == 0
  check_invariant(c, f.mp)
  # reset_fmodel with explicit nuclear sites.
  c.reset_fmodel(c.fmodel.deep_copy(), sites_cart=f.nuclear)
  assert c.sites_cart_nuclear.rms_difference(f.nuclear) == 0
  check_invariant(c, f.mp)
  # An fmodel with other sites and no nuclear sites is refused.
  other = c.fmodel.deep_copy()
  other.xray_structure.set_sites_cart(nuclear2)
  try:
    c.reset_fmodel(other)
  except AssertionError:
    pass
  else:
    raise AssertionError("reset_fmodel must refuse an fmodel with other sites")
  # After a step, update keeps the invariant.
  x = flex.double(f.nuclear.size()*3, 0.01)
  c.update(x=x)
  check_invariant(c, f.mp)
  assert c.sites_cart_nuclear.rms_difference(
    f.nuclear + flex.vec3_double(x)) == 0

def test_exclude(r):
  """
  Water excluded: its gradient is zero; FD on three non-excluded H of
  different kinds (chain start, peptide D, altloc hydroxyl) and their
  parents.
  """
  model = get_model("n_gaussian", d_min=r.f_obs.d_min())
  exclude_selection = model.selection("water")
  assert exclude_selection.count(True) == 3
  f = fixture(r, exclude_selection=exclude_selection)
  c = f.calculator(sites_cart=f.nuclear, mapper=f.mp)
  c.update(x=flex.double(f.nuclear.size()*3, 0))
  g = flex.vec3_double(c.g)
  for i in exclude_selection.iselection():
    assert g[i] == (0, 0, 0), (i, g[i])
  # The mapped water H and their parent O are excluded.
  excluded = set(exclude_selection.iselection())
  assert len([i for i in f.mp.i_h if i in excluded]) == 2
  h_indices = diverse_h(f, excluded=excluded,
    labels=[("ARG", "H1"), ("SER", "D"), ("SER", "HG")])
  check_fd(c, f.mp, f.nuclear, h_indices)
  check_invariant(c, f.mp)

def filtered_norm(g):
  d = flex.sqrt(g.dot())
  sel = d > flex.mean(d)*6
  return g.select(~sel).norm()

def test_compute_weight_maps(r):
  f = fixture(r)
  # Grad flags as calculator.sites sets them.
  f.calculator(sites_cart=f.nuclear, mapper=f.mp)
  fmodel_start = f.fmodel.deep_copy()
  w = calculator.compute_weight(fmodel=f.fmodel, restraints_manager=f.rm,
    sites_cart=f.nuclear, mapper=f.mp)
  # fmodel is not changed.
  assert f.fmodel.xray_structure.sites_cart().rms_difference(
    fmodel_start.xray_structure.sites_cart()) == 0
  # By hand: data gradient at the shrunk sites, mapped to the nuclear sites;
  # restraint gradient at the nuclear sites.
  fm = f.fmodel.deep_copy()
  fm.xray_structure.set_sites_cart(f.mp.shrink(f.nuclear))
  fm.update_xray_structure(update_f_calc=True)
  gx = flex.vec3_double(fm.target_functor()(compute_gradients=True).\
    gradients_wrt_atomic_parameters(site=True).packed())
  gx = f.mp.gradients(f.nuclear, gx)
  gc = f.rm.target_and_gradients(sites_cart=f.nuclear)[1]
  w_hand = filtered_norm(gc) / filtered_norm(gx)
  assert approx_equal(w, w_hand, 1.e-9), (w, w_hand)
  # The old-style call (no mapper, sites from fmodel) gives another weight.
  w_old = calculator.compute_weight(fmodel=f.fmodel, restraints_manager=f.rm)
  assert abs(w - w_old) > 1.e-3 * abs(w), (w, w_old)
  # The identity mapper and explicit fmodel sites give the old-style value.
  identity = xh_mapper.mapper(f.model, scattering_table="neutron")
  assert identity.is_identity
  w_id = calculator.compute_weight(fmodel=f.fmodel, restraints_manager=f.rm,
    sites_cart=f.fmodel.xray_structure.sites_cart(), mapper=identity)
  assert w_id == w_old, (w_id, w_old)
  # Shaking applies to the nuclear sites: same as shaking by hand.
  w_s = calculator.compute_weight(fmodel=f.fmodel, restraints_manager=f.rm,
    shake_sites=True, sites_cart=f.nuclear, mapper=f.mp)
  w_s2 = calculator.compute_weight(fmodel=f.fmodel, restraints_manager=f.rm,
    shake_sites=True, sites_cart=f.nuclear, mapper=f.mp)
  assert w_s == w_s2
  assert abs(w_s - w) > 1.e-6

def compare(c1, c2, x):
  """Bit-identical targets and gradients of two calculators after update."""
  c1.update(x=x)
  c2.update(x=x)
  assert c1.f == c2.f, (c1.f, c2.f)
  assert c1.g.all_eq(c2.g)
  assert c1.fmodel.r_work() == c2.fmodel.r_work()

def test_neutron_unchanged(r):
  """
  An identity mapper changes nothing: neutron table without H/D pairs, and
  the tst_25 exchangeable H/D model with hd_mapper.
  """
  model = get_model("neutron", d_min=r.f_obs.d_min())
  mp = xh_mapper.mapper(model)
  assert mp.is_identity
  grm = calculator.cctbx_geometry(model)
  rm = get_restraints_manager(model)
  fmodel = get_fmodel(model, r)
  def make(**kwargs):
    c = calculator.sites(fmodel=fmodel.deep_copy(), geometry_rmsd_manager=grm,
      max_bond_rmsd=0.03, max_shift=0.2, restraints_manager=rm, **kwargs)
    c.setw(data_weight=1.3, restraints_weight=0.7,
      restraints_weight_scale=1.1)
    return c
  x0 = flex.double(model.size()*3, 0)
  flex.set_random_seed(0)
  x1 = (flex.random_double(model.size()*3)-0.5)*0.1
  for kwargs in [dict(mapper=mp),
                 dict(mapper=mp, sites_cart=model.get_sites_cart())]:
    c_ref = make()
    c = make(**kwargs)
    assert c.sites_cart_nuclear.rms_difference(c_ref.sites_cart_nuclear) == 0
    for x in [x0, x1, x0]:
      compare(c_ref, c, x)
    c_ref.reset_fmodel(c_ref.fmodel.deep_copy())
    c.reset_fmodel(c.fmodel.deep_copy())
    compare(c_ref, c, x1)
  assert calculator.compute_weight(fmodel=fmodel, restraints_manager=rm) == \
    calculator.compute_weight(fmodel=fmodel, restraints_manager=rm,
      sites_cart=model.get_sites_cart(), mapper=mp)
  # Without an active mapper, given sites must agree with fmodel's.
  other = model.get_sites_cart() + (0.01, 0, 0)
  for call in [
      lambda: calculator.compute_weight(fmodel=fmodel, restraints_manager=rm,
        sites_cart=other, mapper=mp),
      lambda: make(sites_cart=other),
      lambda: make(sites_cart=other, mapper=mp),
      lambda: make().reset_fmodel(fmodel.deep_copy(), sites_cart=other)]:
    try:
      call()
    except AssertionError:
      pass
    else:
      raise AssertionError("sites different from fmodel's must be refused")
  # Exchangeable H/D sites (tst_25 fixture) with hd_mapper.
  from qrefine.tests.unit.tst_25 import read_hierarchy, split_hd
  from qrefine.tests.unit.tst_25 import get_model as get_model_25
  ph, cs = read_hierarchy()
  model_h = get_model_25(ph, cs)
  model_hd, original = split_hd(model_h, labels=("A","B"))
  with open("model_hd.pdb", "w") as fo:
    fo.write(model_hd.model_as_pdb())
  r_hd = run_fmodel(args=["model_hd.pdb", "high_resolution=2.0"],
    prefix="model_hd")
  model_hd.setup_scattering_dictionaries(scattering_table="neutron",
    d_min=r_hd.f_obs.d_min())
  mp = xh_mapper.mapper(model_hd)
  assert mp.is_identity
  hdm = refine.hd_mapper(model=model_hd)
  assert len(hdm.ehd) == 3
  rm = get_restraints_manager(hdm.get_single_model())
  grm = calculator.cctbx_geometry(model_hd)
  fmodel = get_fmodel(model_hd, r_hd)
  def make_hd(**kwargs):
    c = calculator.sites(fmodel=fmodel.deep_copy(), geometry_rmsd_manager=grm,
      max_bond_rmsd=0.03, max_shift=0.2, restraints_manager=rm, hdm=hdm,
      **kwargs)
    c.setw(data_weight=1.3, restraints_weight=0.7,
      restraints_weight_scale=1.1)
    return c
  x0 = flex.double(model_hd.size()*3, 0)
  x1 = (flex.random_double(model_hd.size()*3)-0.5)*0.1
  for kwargs in [dict(mapper=mp),
                 dict(mapper=mp, sites_cart=model_hd.get_sites_cart())]:
    c_ref = make_hd()
    c = make_hd(**kwargs)
    for x in [x0, x1]:
      compare(c_ref, c, x)
  assert calculator.compute_weight(fmodel=fmodel, restraints_manager=rm,
    hdm=hdm) == calculator.compute_weight(fmodel=fmodel,
      restraints_manager=rm, hdm=hdm, sites_cart=model_hd.get_sites_cart(),
      mapper=mp)

def test_special_position(r):
  """
  With an active mapper an atom off its special position is snapped on it
  at construction, as fmodel's sites are without one. The fixture in a P 2 2 2
  cell of 100 A (far from its symmetry elements) plus a water O 0.01 A off
  the two-fold along z at (1/2, 1/2, z).
  """
  lines = []
  for line in pdb_str.splitlines():
    if line.startswith("CRYST1"):
      line = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 2 2 2"
    if line.startswith("END"):
      line = ("HETATM   28  O   HOH S   2      50.010  50.010  30.000  1.00"
        "  8.00           O\nEND")
    lines.append(line)
  with open("in_special.pdb", "w") as fo:
    fo.write("\n".join(lines)+"\n")
  model = get_model("n_gaussian", d_min=2.0, file_name="in_special.pdb")
  assert str(model.crystal_symmetry().space_group_info()) == "P 2 2 2"
  i_o = model.size()-1
  assert model.get_hierarchy().atoms()[i_o].pdb_label_columns() == \
    " O   HOH S   2 "
  assert model.get_xray_structure().site_symmetry_table().is_special_position(
    i_o)
  mp = xh_mapper.mapper(model)
  assert not mp.is_identity
  nuclear = mp.expand(model.get_sites_cart())
  # Processing put the water O on the axis; move it off, as a refinement
  # step would.
  assert approx_equal(nuclear[i_o], (50.0, 50.0, 30.0), 1.e-6)
  nuclear[i_o] = (50.01, 50.01, 30.0)
  # apply_symmetry_sites moves the water O on the axis.
  xrs = model.get_xray_structure().deep_copy_scatterers()
  xrs.set_sites_cart(nuclear)
  xrs.apply_symmetry_sites()
  expected = xrs.sites_cart()
  assert approx_equal(expected[i_o], (50.0, 50.0, 30.0), 1.e-6)
  assert nuclear[i_o] != expected[i_o]
  f_obs = abs(model.get_xray_structure().structure_factors(d_min=2.0,
    algorithm="direct").f_calc())
  class data(object): pass
  d = data()
  d.f_obs = f_obs
  d.flags = f_obs.generate_r_free_flags(fraction=0.1)
  fmodel = get_fmodel(model, d)
  fmodel.xray_structure.set_sites_cart(mp.shrink(nuclear))
  fmodel.update_xray_structure(update_f_calc=True)
  rm = get_restraints_manager(model, handling="native")
  c = calculator.sites(fmodel=fmodel,
    geometry_rmsd_manager=calculator.cctbx_geometry(model),
    max_bond_rmsd=0.03, max_shift=0.2, restraints_manager=rm,
    sites_cart=nuclear, mapper=mp)
  assert approx_equal(c.sites_cart_nuclear, expected, 1.e-9)
  assert approx_equal(c.sites_cart_nuclear[i_o], (50.0, 50.0, 30.0), 1.e-9)
  check_invariant(c, mp)

def test_check_model_fmodel_order(r):
  """
  refine.check_model_fmodel_order: nothing on matching order, Sorry when
  two scatterers of fmodel are swapped with respect to the model atoms.
  """
  f = fixture(r)
  fmodel = get_fmodel(f.model, r)
  refine.check_model_fmodel_order(model=f.model, fmodel=fmodel)
  sites = f.model.get_sites_cart()
  swapped = sites.deep_copy()
  swapped[0], swapped[1] = sites[1], sites[0]
  fmodel.xray_structure.set_sites_cart(swapped)
  try:
    refine.check_model_fmodel_order(model=f.model, fmodel=fmodel)
  except Sorry as e:
    assert "atom order differ" in str(e), str(e)
  else:
    raise AssertionError("Sorry expected")

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  calculator.sites and compute_weight with the X-H mapper: nuclear sites in
  the calculator, shrunk sites in fmodel, mapped data gradients (FD), the
  set_sites/reset_fmodel round trip, the exclude selection, identity
  mappers bit-identical to no mapper, special positions with an active
  mapper, and the model/fmodel atom order check.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  with open("in.pdb", "w") as fo:
    fo.write(pdb_str)
  r = run_fmodel(args=["in.pdb", "high_res=2.5"], prefix="data")
  test_fd_xray(r)
  test_set_sites_round_trip(r)
  test_exclude(r)
  test_compute_weight_maps(r)
  test_neutron_unchanged(r)
  test_special_position(r)
  test_check_model_fmodel_order(r)

if(__name__ == "__main__"):
  run()
  print("OK")
