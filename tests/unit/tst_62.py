from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from cctbx import crystal
from cctbx.crystal import super_cell as cctbx_super_cell
from qrefine import conformers, fragment_capping, qr, refine, restraints
from qrefine.fragment import fragments
from qrefine.super_cell import expand

qr_unit_tests_data = os.path.join(libtbx.env.find_in_repositories("qrefine"),
  "tests", "unit", "data_files")

def n_bonds_in(hierarchy, geometry, chain, resseq, name):
  atoms = hierarchy.atoms()
  t = [a for a in atoms if a.name.strip() == name and
    a.parent().parent().resseq.strip() == resseq and
    a.parent().parent().parent().id == chain][0]
  table = geometry.bond_params_table
  n = len(table[t.i_seq])
  for i in range(t.i_seq):
    if(t.i_seq in table[i]): n += 1
  return n

def n_bonds_capper(c, chain, resseq, name):
  """Bonds of an atom of the capper's fragment, from its cached bond table."""
  t = [a for a in c.fragment_hierarchy.atoms() if a.name.strip() == name and
    a.parent().parent().resseq.strip() == resseq and
    a.parent().parent().parent().id == chain][0]
  return len(c.dangling.bond_table.get(t.i_seq, ()))

def get_model():
  pdb_inp = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,
    "6uf7_qmready.pdb"))
  return mmtbx.model.manager(model_input=pdb_inp, log=null_out())

def exercise_symmetry_closure():
  """
  6uf7 (half a cyclic peptide in the asymmetric unit, the ring closed through
  a symmetry mate): in the super sphere the closure is a C-N bond between
  the master chain and a copy with another chain id. The capper restores
  that bond, so PRO 1 N has three bonds and is neither dangling nor capped.
  """
  model = get_model()
  cs = model.crystal_symmetry()
  for conf in conformers.split(model):
    fm = fragments(working_folder=None, make_working_folder=False,
      clustering_method="bcc", maxnum_residues_in_cluster=15,
      pdb_hierarchy=conf.model.get_hierarchy(), qm_engine_name="mopac",
      crystal_symmetry=cs, select_within_radius=10, cap_fragments=True)
    assert len(fm.fragment_cappers) == 1
    c = fm.fragment_cappers[0]
    assert n_bonds_capper(c, "A", "1", "N") == 3
    assert c.n_closures >= 1, c.n_closures
    # no cap near PRO A 1 N
    capped = c.cap(sites_cart=fm.pdb_hierarchy_super.select(
      fm.fragment_super_selections[0]).atoms().extract_xyz()).pdb_hierarchy
    n = [a for a in capped.atoms() if a.name.strip() == "N" and
      a.parent().parent().resseq.strip() == "1" and
      a.parent().parent().parent().id == "A"][0]
    d = (capped.atoms().extract_xyz() - n.xyz).norms()
    close_h = [capped.atoms()[i] for i in (d < 1.2).iselection()
      if capped.atoms()[i].element_is_hydrogen()]
    assert len(close_h) == 0, [a.quote() for a in close_h]

def exercise_super_sphere_restraints():
  """
  The cctbx restraints of the super sphere (cluster route, restraints=cctbx)
  have the closure bond too: PRO 1 N of the master has three bonds.
  """
  model = get_model()
  for conf in conformers.split(model):
    e = expand(pdb_hierarchy=conf.model.get_hierarchy().deep_copy(),
      crystal_symmetry=model.crystal_symmetry(), select_within_radius=10)
    geometry = e.super_sphere_geometry_restraints_manager.geometry
    assert n_bonds_in(e.ph_super_sphere, geometry, "A", "1", "N") == 3

def exercise_expansion_restraints():
  """
  The expansion route with restraints=cctbx: the restraints of the completed
  super sphere have the closure bond.
  """
  model = get_model()
  params = qr.get_default_params()
  params.restraints = "cctbx"
  params.expansion = True
  params.cluster.clustering = False
  params.cluster.select_within_radius = 10
  params.parallel.nproc = 1
  rm = refine.create_restraints_manager(params=params, model=model,
    working_folder=None)
  routes = [c.restraints_manager for c in rm.conformers]
  assert len(routes) == 2, len(routes)
  for route in routes:
    assert isinstance(route, restraints.from_expansion)
    geometry = route.restraints_manager.geometry_restraints_manager.geometry
    assert n_bonds_in(route.pdb_hierarchy_super_completed, geometry,
      "A", "1", "N") == 3

def exercise_expansion_qm_capping():
  """
  The expansion route with restraints=qm (phenix.aquaref with data): the
  whole super sphere is capped. A copy at the sphere edge is closed through
  a copy outside the sphere; that copy is in the super cell (a cut bond, one
  cap), so no PRO N is an unsatisfied valence (was a Sorry at radius 7 and
  10). No QM job runs.
  """
  model = get_model()
  for radius in [7, 10]:
    params = qr.get_default_params()
    params.restraints = "qm"
    params.quantum.engine_name = "mopac"
    params.expansion = True
    params.cluster.clustering = False
    params.cluster.select_within_radius = radius
    params.parallel.nproc = 1
    rm = refine.create_restraints_manager(params=params,
      model=model.deep_copy(), working_folder=None)
    routes = [c.restraints_manager for c in rm.conformers]
    assert len(routes) == 2, len(routes)
    for route in routes:
      assert isinstance(route, restraints.from_expansion)
      assert route.capper is not None
      assert route.capper.n_caps > 0, route.capper.n_caps
      assert n_bonds_capper(route.capper, "A", "1", "N") == 3

def heavy_neighbour_counts(pdb_hierarchy, crystal_symmetry, cutoff):
  """
  Heavy atoms within cutoff of every heavy atom, symmetry mates included.
  """
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  heavy = [not a.element_is_hydrogen() for a in atoms]
  asu_mappings = crystal_symmetry.special_position_settings().asu_mappings(
    buffer_thickness = cutoff,
    sites_cart       = sites_cart)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(distance_cutoff=cutoff)
  pair_sym_table = pair_asu_table.extract_pair_sym_table(
    skip_j_seq_less_than_i_seq=False)
  counts = flex.int(sites_cart.size(), 0)
  for i_seq, row in enumerate(pair_sym_table):
    if(not heavy[i_seq]): continue
    for j_seq, rt_mx_list in row.items():
      if(heavy[j_seq]): counts[i_seq] += len(rt_mx_list)
  return counts

def exercise_bonded_copies(cutoff=2.5):
  """
  expand adds to the super cell the copies covalently bonded to a super
  sphere atom: every heavy sphere atom has in the super cell the heavy
  neighbours it has in the crystal. Only those copies are added, and the
  super sphere is that of the cctbx manager.
  """
  model = get_model()
  cs = model.crystal_symmetry()
  n_added = 0
  for conf in conformers.split(model):
    ph = conf.model.get_hierarchy()
    n = ph.atoms_size()
    expected = heavy_neighbour_counts(ph, cs, cutoff)
    for radius in [5, 7, 10]:
      plain = cctbx_super_cell.manager(pdb_hierarchy=ph.deep_copy(),
        crystal_symmetry=cs, select_within_radius=radius)
      n_plain = plain.super_cell_hierarchy.atoms_size()
      e = expand(pdb_hierarchy=ph.deep_copy(), crystal_symmetry=cs,
        select_within_radius=radius)
      atoms = e.SCM.super_cell_hierarchy.atoms()
      sites_cart = atoms.extract_xyz()
      heavy = flex.bool([not a.element_is_hydrogen() for a in atoms])
      touched = flex.bool(atoms.size(), False)
      for k in e.SCM.keep_selection.iselection():
        if(not heavy[k]): continue
        near = ((sites_cart - sites_cart[k]).dot() < cutoff**2) & heavy
        touched |= near
        # super cell atom k is a copy of master atom k % n
        assert near.count(True) - 1 == expected[k % n], (radius,
          atoms[k].id_str(), near.count(True) - 1, expected[k % n])
      # every added copy has an atom next to a sphere atom
      assert atoms.size() > n_plain
      for k in range(n_plain // n, atoms.size() // n):
        assert touched[k*n:(k+1)*n].count(True) > 0, (radius, k)
        n_added += 1
      # the sphere is not changed, and update() moves the added copies
      assert e.SCM.keep_selection[:n_plain].all_eq(plain.keep_selection)
      assert e.SCM.keep_selection[n_plain:].count(True) == 0
      assert approx_equal(e.ph_super_sphere.atoms().extract_xyz(),
        plain.super_sphere_hierarchy.atoms().extract_xyz())
      e.update(sites_cart=ph.atoms().extract_xyz())
      assert approx_equal(
        e.SCM.super_cell_hierarchy.atoms().extract_xyz(), sites_cart)
  assert n_added > 0

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_symmetry_closure()
  exercise_super_sphere_restraints()
  exercise_expansion_restraints()
  exercise_expansion_qm_capping()
  exercise_bonded_copies()

if(__name__ == "__main__"):
  run()
  print("OK")
