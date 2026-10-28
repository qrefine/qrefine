from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import libtbx.load_env
import iotbx.pdb
import iotbx.pdb.hierarchy
import mmtbx.model
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex
from qrefine import restraints
from qrefine import conformers

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def build_cctbx(conformer_model, working_folder):
  return restraints.from_cctbx(restraints_manager=restraints.get_cctbx_gradients(
    ph=conformer_model.get_hierarchy(), cs=conformer_model.crystal_symmetry()
    ).model.get_restraints_manager())

def gradients_from_conformers(model):
  fc = conformers.from_conformers(model=model, build=build_cctbx)
  return fc.target_and_gradients(sites_cart=model.get_sites_cart())[1]

def cctbx_calculator(ph, crystal_symmetry):
  """
  Gradients from cctbx restraints built for the given hierarchy.
  """
  return restraints.get_cctbx_gradients(ph = ph, cs = crystal_symmetry).gradients

def duplicate(model, labels, residue_indices):
  """
  Copy of model where each residue group in residue_indices (each having one
  atom group) is replaced by identical conformers, one per label, with
  occupancies 1/len(labels). Returns the new model and, per atom of it, the
  i_seq of the original atom.
  """
  ph = model.get_hierarchy().deep_copy()
  for a in ph.atoms(): a.tmp = a.i_seq
  rgs = list(ph.residue_groups())
  for i in residue_indices:
    rg = rgs[i]
    ag = rg.only_atom_group()
    for label in labels:
      ag_c = ag.detached_copy()
      for a_c, a in zip(ag_c.atoms(), ag.atoms()): # detached_copy resets tmp
        a_c.tmp = a.tmp
      ag_c.altloc = label
      ag_c.atoms().set_occ(flex.double(ag_c.atoms_size(), 1.0/len(labels)))
      rg.append_atom_group(ag_c)
    rg.remove_atom_group(ag)
  ph.atoms_reset_serial()
  ph.atoms().reset_i_seq()
  original = flex.size_t([a.tmp for a in ph.atoms()])
  dup_model = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=model.crystal_symmetry(), log=null_out())
  return dup_model, original

def exercise_vs_cctbx():
  """
  cctbx counts each restraint once (blank-blank once, plus every blank-X
  term), which cannot be recovered from whole-conformer gradients. So:
    - altloc atoms: must equal cctbx exactly;
    - blank atoms: m*g = g_cctbx + (m-1)*g_blank_only, which also means
      equal to cctbx for blank atoms that do not interact with altlocs.
  g_blank_only is available from cctbx (not from QM) and is used only here.
  """
  for file_name in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
                    "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:
    pi = iotbx.pdb.input(file_name = os.path.join(qr_unit_tests_data,file_name))
    model = mmtbx.model.manager(model_input=pi, log=null_out())
    ph = model.get_hierarchy()
    cs = model.crystal_symmetry()
    g_cctbx = cctbx_calculator(ph, cs)
    g = gradients_from_conformers(model)
    # energies_sites (real-space refinement): energy 0, same gradients
    fc = conformers.from_conformers(model=model, build=build_cctbx)
    es = fc.energies_sites(sites_cart=model.get_sites_cart())
    assert es.target == 0
    assert approx_equal(es.gradients, fc.target_and_gradients(
      sites_cart=model.get_sites_cart())[1], 0)
    assert approx_equal(es.gradients, g, 0)
    ci = ph.get_conformer_indices().conformer_indices
    n_conformers = flex.max(ci)
    assert n_conformers == 2
    blank = ci == 0
    altloc = ~blank
    assert approx_equal(g.select(altloc), g_cctbx.select(altloc))
    if(blank.count(True) == 0): continue
    bsel = blank.iselection()
    g_blank_only = cctbx_calculator(ph.select(bsel), cs)
    m = n_conformers
    assert approx_equal(g.select(bsel)*m,
      g_cctbx.select(bsel) + g_blank_only*(m-1))
    # Blank atoms bonded to altloc atoms really differ from cctbx
    assert not approx_equal(g.select(bsel), g_cctbx.select(bsel), out=None)

def exercise_no_altlocs():
  """
  Without altlocs there are no conformer labels, so the caller
  (create_restraints_manager) never builds from_conformers.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  ph = pi.construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  assert conformers.n_conformers(ph) == 0

def exercise_duplicate(labels):
  """
  Duplicate a residue into identical conformers (one per label): every atom
  must get the gradient of the original single-conformer model.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  ph = pi.construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  model = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=pi.crystal_symmetry(), log=null_out())
  g_single = cctbx_calculator(ph, model.crystal_symmetry())
  dup_model, original = duplicate(model=model, labels=labels,
    residue_indices=[2])
  ph = dup_model.get_hierarchy()
  assert conformers.conformer_labels(ph) == labels
  g = gradients_from_conformers(dup_model)
  assert approx_equal(g, g_single.select(original))
  assert flex.max(conformers.multiplicity(ph)) == len(labels)

class _dummy_fragment_manager(object):
  def __init__(self):
    self.n_calls = 0
    self.interaction_list = [(1,2)]
  def set_up_cluster_qm(self):
    self.n_calls += 1

class _dummy_manager(object):
  def __init__(self, with_fragment_manager):
    if(with_fragment_manager):
      self.fragment_manager = _dummy_fragment_manager()
  def target_and_gradients(self, sites_cart):
    return 0, flex.vec3_double(sites_cart.size(), (0,0,0))

def exercise_proxy():
  """
  fragment_manager of from_conformers forwards set_up_cluster_qm to every
  conformer and concatenates their interaction lists.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  model = mmtbx.model.manager(model_input=pi, log=null_out())
  managers = []
  def build(conformer_model, working_folder):
    managers.append(_dummy_manager(with_fragment_manager=True))
    return managers[-1]
  fc = conformers.from_conformers(model=model, build=build)
  assert len(managers) == 2
  fc.fragment_manager.set_up_cluster_qm()
  assert [m.fragment_manager.n_calls for m in managers] == [1, 1]
  assert fc.fragment_manager.interaction_list == [(1,2),(1,2)]
  fc = conformers.from_conformers(model=model,
    build=lambda conformer_model, working_folder: _dummy_manager(False))
  assert fc.fragment_manager is None

def exercise_microheterogeneity():
  """
  Conformers with different residue names (microheterogeneity): in conformer
  X's model the residue takes the resname of the atom group with label X,
  blank atoms included. GLY 95 has blank, A and B atom groups; TYR 94 has
  A and B only (no merge needed).
  """
  file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb")
  for i_rg, original in [(8, "GLY"), (7, "TYR")]:
    pi = iotbx.pdb.input(file_name = file_name)
    ph = pi.construct_hierarchy()
    rg = list(ph.residue_groups())[i_rg]
    altlocs = [ag.altloc for ag in rg.atom_groups()]
    assert altlocs == {8: ["", "A", "B"], 7: ["A", "B"]}[i_rg]
    rg.atom_groups()[-1].resname = "ALA"
    ph.atoms().reset_i_seq()
    model = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
      crystal_symmetry=pi.crystal_symmetry(), log=null_out())
    rg_isel = set(rg.atoms().extract_i_seq())
    parts = conformers.split(model)
    assert [p.label for p in parts] == ["A", "B"]
    for p, resname in zip(parts, [original, "ALA"]):
      h = p.model.get_hierarchy()
      assert list(h.altloc_indices().keys()) == [""]
      assert max(r.atom_groups_size() for r in h.residue_groups()) == 1
      assert list(h.atoms().extract_name()) == list(
        ph.atoms().select(p.iselection).extract_name())
      assert h.atoms().extract_xyz().max_distance(
        ph.atoms().extract_xyz().select(p.iselection)) < 1.e-6
      ag = list(h.residue_groups())[i_rg].only_atom_group()
      assert ag.resname == resname
      # every atom of this residue in this conformer, blank atoms included
      rg_atoms = [k for k, i in enumerate(p.iselection) if i in rg_isel]
      assert rg_atoms == list(ag.atoms().extract_i_seq())
      for a in h.atoms().select(flex.size_t(rg_atoms)):
        assert a.parent().resname == resname
      if(i_rg == 8):
        n_blank = rg.atom_groups()[0].atoms_size()
        assert len(rg_atoms) == n_blank + 1

def exercise_incomplete_conformers():
  """
  Residue with A/B/C in a model where other residues have only A/B: conformer
  C would lack atoms, which QM cannot handle.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  model = mmtbx.model.manager(model_input=pi, log=null_out())
  ph = model.get_hierarchy()
  rgs = [rg for rg in ph.residue_groups() if rg.atom_groups_size() > 1]
  assert len(rgs) > 1
  ag_c = rgs[0].atom_groups()[-1].detached_copy()
  ag_c.altloc = "C"
  rgs[0].append_atom_group(ag_c)
  ph.atoms().reset_i_seq()
  try:
    conformers.split(model)
  except Sorry as e:
    assert str(e).find("conformers would be incomplete") > -1
  else: raise AssertionError("Expected Sorry for incomplete conformers")

def exercise_split():
  """
  Conformer labels, multiplicity and single-conformer models of
  conformers.split.
  """
  pi = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data, "h_altconf_complete.pdb"))
  model = mmtbx.model.manager(model_input=pi, log=null_out())
  ph = model.get_hierarchy()
  assert conformers.conformer_labels(ph) == ["A", "B"]
  assert conformers.n_conformers(ph) == 2
  m = conformers.multiplicity(ph)
  ci = ph.get_conformer_indices().conformer_indices
  assert m.select(ci == 0).all_eq(2.0) and m.select(ci != 0).all_eq(1.0)
  parts = conformers.split(model)
  assert [p.label for p in parts] == ["A", "B"]
  for p, c in zip(parts, [1, 2]):
    assert p.iselection.size() == ((ci == c) | (ci == 0)).count(True)
    assert list(p.iselection) == list(((ci == c) | (ci == 0)).iselection())
    h = p.model.get_hierarchy()
    assert list(h.altloc_indices().keys()) == [""]
    assert max(rg.atom_groups_size() for rg in h.residue_groups()) == 1
    assert h.atoms().extract_xyz().max_distance(
      ph.atoms().extract_xyz().select(p.iselection)) < 1.e-6
    assert list(h.atoms().extract_name()) == list(ph.atoms().select(p.iselection).extract_name())
  # no conformer labels: multiplicity 1, nothing to split
  h = ph.deep_copy()
  h.remove_alt_confs(always_keep_one_conformer=True)
  assert conformers.n_conformers(h) == 0
  assert conformers.multiplicity(h).all_eq(1.0)
  model_h = mmtbx.model.manager(model_input=None, pdb_hierarchy=h,
    crystal_symmetry=model.crystal_symmetry(), log=null_out())
  assert conformers.split(model_h) == []
  # hydrogen-only labels are not conformers
  rg = list(ph.residue_groups())[1] # GLY 88: first residue with atom H
  ag = rg.only_atom_group()
  h_atom = [a for a in ag.atoms() if a.name.strip() == "H"][0]
  ag_h = iotbx.pdb.hierarchy.atom_group(altloc="C", resname=ag.resname)
  ag_h.append_atom(h_atom.detached_copy())
  ag.remove_atom(h_atom)
  rg.append_atom_group(ag_h)
  ph.atoms().reset_i_seq()
  assert conformers.conformer_labels(ph) == ["A", "B"]
  assert conformers.multiplicity(ph)[ag_h.atoms()[0].i_seq] == 2.0
  for p in conformers.split(model):
    assert ag_h.atoms()[0].i_seq in p.iselection

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Combine gradients from single-conformer calculations (the QM route) for
  models with altlocs. cctbx is used as a stand-in for the QM calculator.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_split()
  exercise_vs_cctbx()
  exercise_no_altlocs()
  for labels in [["A","B"], ["A","B","C"], ["A","B","C","D"]]:
    exercise_duplicate(labels=labels)
  exercise_proxy()
  exercise_incomplete_conformers()
  exercise_microheterogeneity()

if(__name__ == '__main__'):
  run()
