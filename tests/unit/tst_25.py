from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import iotbx.pdb
import mmtbx.model
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx.array_family import flex
from qrefine import refine, conformers
from qrefine.tests.unit.tst_11 import (qr_unit_tests, process,
  get_restraints_manager)
from qrefine.tests.unit import run_fmodel
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.tst_06 import refine_fast

def read_hierarchy(file_name="h_altconf_complete.pdb"):
  pdb_inp = iotbx.pdb.input(os.path.join(qr_unit_tests, "data_files",
    file_name))
  return pdb_inp.construct_hierarchy(), pdb_inp.crystal_symmetry()

def get_model(ph, crystal_symmetry):
  model = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=crystal_symmetry, log=null_out())
  return process(model)

def split_hd(model_h, n_residues=3, labels=("C","D")):
  """
  Copy of model_h where atom H of the first n_residues blank-only residue
  groups (that have an H) is moved to atom group altloc labels[0] (occ 0.5)
  and duplicated as D in atom group altloc labels[1] (occ 0.5), same
  coordinates. labels=("A","B") is the standard neutron labelling.
  Returns the processed model and original: model_h i_seq of each atom of
  it (H and its D copy both map to the H).
  """
  ph = model_h.get_hierarchy().deep_copy()
  ph.reset_i_seq_if_necessary()
  for a in ph.atoms(): a.tmp = a.i_seq
  n_done = 0
  for rg in ph.residue_groups():
    if(n_done == n_residues): break
    ags = rg.atom_groups()
    if(len(ags) != 1 or ags[0].altloc != ""): continue
    ag = ags[0]
    hs = [a for a in ag.atoms() if a.name.strip() == "H"]
    if(len(hs) != 1): continue
    h = hs[0]
    i_seq = h.tmp
    ag.remove_atom(h)
    h_copy = h.detached_copy()
    h_copy.occ = 0.5
    h_copy.tmp = i_seq
    d_copy = h.detached_copy()
    d_copy.name = " D  "
    d_copy.element = " D"
    d_copy.occ = 0.5
    d_copy.tmp = i_seq
    for altloc, atom in [(labels[0], h_copy), (labels[1], d_copy)]:
      new_ag = iotbx.pdb.hierarchy.atom_group(altloc=altloc,
        resname=ag.resname)
      new_ag.append_atom(atom)
      rg.append_atom_group(new_ag)
    n_done += 1
  assert n_done == n_residues
  original = flex.size_t([a.tmp for a in ph.atoms()])
  ph.atoms_reset_serial()
  ph.atoms().reset_i_seq()
  names = [a.name for a in ph.atoms()]
  model_hd = get_model(ph, model_h.crystal_symmetry())
  # Processing keeps the atom order (sort_atoms=False), so original holds.
  assert [a.name for a in model_hd.get_hierarchy().atoms()] == names
  return model_hd, original

def check_gradients(model_h, model_hd, original):
  """
  hd_mapper composed with the restraints manager of the single model gives
  the gradients of model_h; H and its D copy both carry the H gradient.
  Returns the restraints manager of the single model.
  """
  hdm = refine.hd_mapper(model=model_hd)
  assert len(hdm.ehd) == 3
  single = hdm.get_single_model()
  elements = [e.strip() for e in single.get_hierarchy().atoms().extract_element()]
  assert elements.count("D") == 0
  assert single.size() == model_h.size()
  rm, _ = get_restraints_manager(None, expansion=False, clustering=False,
    handling="conformers", model=single)
  g = hdm.expand(rm.target_and_gradients(
    sites_cart=hdm.shrink(model_hd.get_sites_cart()))[1])
  rm_ref, _ = get_restraints_manager(None, expansion=False, clustering=False,
    handling="conformers", model=model_h)
  g_ref = rm_ref.target_and_gradients(sites_cart=model_h.get_sites_cart())[1]
  assert g.size() == model_hd.size()
  for i in range(g.size()):
    assert approx_equal(g[i], g_ref[original[i]]), i
  for pair in hdm.ehd:
    assert g[pair[0][0]] == g[pair[1][0]]
  return rm

def exercise_mixed(labels=("C","D")):
  """
  H/D exchange sites (altlocs labels, H-only atom groups) together with
  heavy-atom altlocs A/B. With labels A/B the H/D sites share the labels of
  the heavy-atom conformers; the exchanged H must still belong to both.
  """
  ph, cs = read_hierarchy()
  model_h = get_model(ph, cs)
  assert conformers.conformer_labels(model_h.get_hierarchy()) == ["A", "B"]
  model_hd, original = split_hd(model_h, labels=labels)
  h = model_hd.get_hierarchy()
  assert sorted(h.altloc_indices().keys()) == sorted(set(
    ["", "A", "B"] + list(labels)))
  assert conformers.conformer_labels(h) == ["A", "B"]
  rm = check_gradients(model_h, model_hd, original)
  assert isinstance(rm, conformers.from_conformers)
  return model_hd

def exercise_hd_only(labels=("C","D")):
  """
  H/D altlocs (labels) without heavy-atom altlocs: no conformers, plain
  route.
  """
  ph, cs = read_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  model_h = get_model(ph, cs)
  model_hd, original = split_hd(model_h, labels=labels)
  h = model_hd.get_hierarchy()
  assert sorted(h.altloc_indices().keys()) == ["", labels[0], labels[1]]
  assert conformers.conformer_labels(h) == []
  hdm = refine.hd_mapper(model=model_hd)
  assert conformers.n_conformers(hdm.get_single_model().get_hierarchy()) == 0
  rm = check_gradients(model_h, model_hd, original)
  assert not isinstance(rm, conformers.from_conformers)

def exercise_mixed_in_conformer_groups():
  """
  H/D exchange site inside the heavy-atom A/B atom groups of a disordered
  side chain: TYR 94 rebuilt with a blank backbone (from conformer A) and A/B
  side chains; its backbone H (blank parent N) is moved into the A atom group
  and a D copy (same position) is added to the B atom group. hd_mapper must
  give the H to both conformers: the composed gradient equals that of the
  reference model (H blank), H and D equal.
  """
  ph, cs = read_hierarchy()
  backbone = ["N", "CA", "C", "O", "H", "HA"]
  rg = [rg for rg in ph.residue_groups() if rg.resseq.strip() == "94"][0]
  ag_a, ag_b = rg.atom_groups()
  assert [ag_a.altloc, ag_b.altloc] == ["A", "B"]
  blank = iotbx.pdb.hierarchy.atom_group(altloc="", resname=ag_a.resname)
  for ag in [ag_a, ag_b]:
    for a in ag.atoms():
      if(a.name.strip() not in backbone): continue
      if(ag.altloc == "A"):
        a_copy = a.detached_copy()
        a_copy.occ = 1.0
        blank.append_atom(a_copy)
      ag.remove_atom(a)
  rg.insert_atom_group(0, blank)
  ph.atoms_reset_serial()
  ph.atoms().reset_i_seq()
  model_h = get_model(ph, cs)
  assert conformers.conformer_labels(model_h.get_hierarchy()) == ["A", "B"]
  # Mixed model: blank H -> A atom group, D copy -> B atom group.
  ph = model_h.get_hierarchy().deep_copy()
  ph.reset_i_seq_if_necessary()
  for a in ph.atoms(): a.tmp = a.i_seq
  rg = [rg for rg in ph.residue_groups() if rg.resseq.strip() == "94"][0]
  blank, ag_a, ag_b = rg.atom_groups()
  h = [a for a in blank.atoms() if a.name.strip() == "H"][0]
  i_seq = h.tmp
  blank.remove_atom(h)
  h_copy = h.detached_copy()
  h_copy.occ = 0.5
  h_copy.tmp = i_seq
  d_copy = h.detached_copy()
  d_copy.name = " D  "
  d_copy.element = " D"
  d_copy.occ = 0.5
  d_copy.tmp = i_seq
  ag_a.append_atom(h_copy)
  ag_b.append_atom(d_copy)
  original = flex.size_t([a.tmp for a in ph.atoms()])
  ph.atoms_reset_serial()
  ph.atoms().reset_i_seq()
  names = [a.name for a in ph.atoms()]
  model_hd = get_model(ph, cs)
  assert [a.name for a in model_hd.get_hierarchy().atoms()] == names
  assert len(model_hd.get_hierarchy().exchangeable_hd_selections()) == 0
  hdm = refine.hd_mapper(model=model_hd)
  assert len(hdm.ehd) == 1
  single = hdm.get_single_model()
  assert single.size() == model_h.size()
  rg = [rg for rg in single.get_hierarchy().residue_groups()
    if rg.resseq.strip() == "94"][0]
  assert [ag.altloc for ag in rg.atom_groups()] == ["", "A", "B"]
  assert "H" in [a.name.strip() for a in rg.atom_groups()[0].atoms()]
  rm, _ = get_restraints_manager(None, expansion=False, clustering=False,
    handling="conformers", model=single)
  assert isinstance(rm, conformers.from_conformers)
  g = hdm.expand(rm.target_and_gradients(
    sites_cart=hdm.shrink(model_hd.get_sites_cart()))[1])
  rm_ref, _ = get_restraints_manager(None, expansion=False, clustering=False,
    handling="conformers", model=model_h)
  g_ref = rm_ref.target_and_gradients(sites_cart=model_h.get_sites_cart())[1]
  assert g.size() == model_hd.size()
  for i in range(g.size()):
    assert approx_equal(g[i], g_ref[original[i]]), i
  for pair in hdm.ehd:
    assert g[pair[0][0]] == g[pair[1][0]]

def exercise_qrefine(model_hd, labels):
  """
  qr.refine on a model with H/D exchange sites (altlocs labels) and heavy-atom
  altlocs A/B: the hd_mapper is used and the H/D atom groups and A/B labels
  are kept.
  mode=refine (with data): mode=opt does not pass hdm to calculator.sites_opt
  (see create_calculator), so H/D models fail there regardless of altlocs.
  """
  tag = "".join(labels)
  pdb = os.path.abspath("model_hd_%s.pdb" % tag)
  with open(pdb, "w") as fo:
    fo.write(model_hd.model_as_pdb())
  r_fm = run_fmodel(args=[pdb, "high_resolution=2.0"], prefix="model_hd_%s" % tag)
  r_qr = run_qrefine_in(prefix="hd_%s" % tag, args=[pdb,
    os.path.abspath(r_fm.mtz), "restraints=cctbx",
    "number_of_micro_cycles=1", "max_iterations_refine=2"] + refine_fast)
  assert r_qr.std_out_text.find("Using hd_mapper") > -1
  assert r_qr.model is not None
  h = r_qr.model.get_hierarchy()
  assert sorted(h.altloc_indices().keys()) == sorted(set(
    ["", "A", "B"] + list(labels)))
  assert conformers.conformer_labels(h) == ["A", "B"]
  n_hd = 0
  for rg in h.residue_groups():
    altlocs = [ag.altloc for ag in rg.atom_groups()]
    if(altlocs != ["", labels[0], labels[1]]): continue
    # GLY 95 of the input has its own A/B H atoms (no D)
    if([a.element.strip() for a in rg.atom_groups()[2].atoms()] != ["D"]):
      continue
    ag_c, ag_d = rg.atom_groups()[1:]
    assert [a.name.strip() for a in ag_c.atoms()] == ["H"]
    assert [a.name.strip() for a in ag_d.atoms()] == ["D"]
    assert [a.element.strip() for a in ag_d.atoms()] == ["D"]
    assert approx_equal(ag_c.atoms()[0].xyz, ag_d.atoms()[0].xyz)
    n_hd += 1
  assert n_hd == 3

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  H/D exchange (hd_mapper) together with heavy-atom altlocs.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  for labels in [("C","D"), ("A","B")]:
    model_hd = exercise_mixed(labels=labels)
    exercise_hd_only(labels=labels)
  exercise_mixed_in_conformer_groups()
  exercise_qrefine(model_hd, labels=("A","B"))

if(__name__ == "__main__"):
  run()
