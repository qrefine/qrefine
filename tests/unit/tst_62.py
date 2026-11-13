from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
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

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_symmetry_closure()
  exercise_super_sphere_restraints()
  exercise_expansion_restraints()

if(__name__ == "__main__"):
  run()
  print("OK")
