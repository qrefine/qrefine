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

def _model_with_tyr94_c(file_name="h_altconf_complete.pdb"):
  """
  file_name with a C atom group on TYR 94 only (copy of its B group);
  occupancies of TYR 94: A 0.4, B 0.35, C 0.25.
  """
  pi = iotbx.pdb.input(file_name = os.path.join(qr_unit_tests_data,file_name))
  ph = pi.construct_hierarchy()
  rg = [rg for rg in ph.residue_groups() if rg.resseq_as_int() == 94][0]
  assert [ag.altloc for ag in rg.atom_groups()] == ["A", "B"]
  ag_c = rg.atom_groups()[-1].detached_copy()
  ag_c.altloc = "C"
  rg.append_atom_group(ag_c)
  for ag, occ in zip(rg.atom_groups(), [0.4, 0.35, 0.25]):
    ag.atoms().set_occ(flex.double(ag.atoms_size(), occ))
  ph.atoms().reset_i_seq()
  ph.atoms_reset_serial()
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=pi.crystal_symmetry(), log=null_out())

def exercise_fallback():
  """
  Residue group lacking a conformer label X takes the atoms of its first own
  label into conformer X: TYR 94 has A/B/C, GLY 95 only A/B (H atoms), so
  conformer C takes GLY 95's A atoms. Multiplicity is membership.
  """
  model = _model_with_tyr94_c()
  ph = model.get_hierarchy()
  sel = conformers.conformer_selections(ph)
  assert [s.label for s in sel] == ["A", "B", "C"]
  assert [s.n_fallback for s in sel] == [0, 0, 1]
  def i_seqs(resseq, altloc):
    return flex.size_t([a.i_seq for a in ph.atoms()
      if a.parent().parent().resseq_as_int() == resseq and
         a.parent().altloc == altloc])
  gly95_a, gly95_b = i_seqs(95, "A"), i_seqs(95, "B")
  assert gly95_a.size() == gly95_b.size() == 1
  tyr94 = flex.size_t([a.i_seq for a in ph.atoms()
    if a.parent().parent().resseq_as_int() == 94])
  blank = flex.size_t([a.i_seq for a in ph.atoms() if a.parent().altloc == ""])
  isel_c = set(sel[2].iselection)
  assert all(i in isel_c for i in gly95_a)
  assert not any(i in isel_c for i in gly95_b)
  for s in sel:
    assert list(s.iselection) == sorted(s.iselection)
  m = conformers.multiplicity(ph)
  assert m.size() == ph.atoms_size()
  assert m.select(gly95_a).all_eq(2.0) and m.select(gly95_b).all_eq(1.0)
  assert m.select(tyr94).all_eq(1.0) and m.select(blank).all_eq(3.0)
  # the old Sorry ("conformers would be incomplete") is gone
  parts = conformers.split(model)
  assert len(parts) == 3
  assert [p.label for p in parts] == ["A", "B", "C"]
  for p, s in zip(parts, sel):
    assert list(p.iselection) == list(s.iselection)
    assert p.model.get_hierarchy().atoms_size() == s.iselection.size()
  # from_conformers reports the fallback
  from six.moves import cStringIO as StringIO
  log = StringIO()
  conformers.from_conformers(model=model, log=log,
    build=lambda conformer_model, working_folder: None)
  text = log.getvalue()
  assert text.find("1 residue groups without label C take their first own "
    "label") > -1, text
  assert text.find("without label A") == -1 and \
    text.find("without label B") == -1

def _water_rg(resseq, labels, o_xyz):
  """
  HOH residue group: one atom group (O, H1, H2) per label, all at the same
  place (a blank label "" gives an ordinary water).
  """
  rg = iotbx.pdb.hierarchy.residue_group(resseq="%4d" % resseq, icode=" ")
  x, y, z = o_xyz
  for label in labels:
    ag = iotbx.pdb.hierarchy.atom_group(altloc=label, resname="HOH")
    for name, element, xyz in [
        (" O  ", " O", (x, y, z)),
        (" H1 ", " H", (x-0.96, y, z)),
        (" H2 ", " H", (x+0.24, y-0.93, z))]:
      a = iotbx.pdb.hierarchy.atom()
      a.name = name
      a.element = element
      a.xyz = xyz
      a.occ = 1.0 if label == "" else 1.0/len(labels)
      a.b = 30.0
      a.hetero = True
      ag.append_atom(a)
    rg.append_atom_group(ag)
  return rg

def model_with_b_only_water(labels=["B"], bonded=False, tyr94_c=False):
  """
  h_altconf_complete (with tyr94_c: plus a C conformer of TYR 94 as in
  _model_with_tyr94_c) and a chain W of two waters: HOH 100 with only the
  given labels, overlapping TYR 94 A (O 0.9 A from CD1 A, 3.4 A from CD1
  B; H2 0.6 A from CE1 A), as the B-only waters of 1PWM; with bonded its
  O is 1.3 A from CA of GLY 99 (blank) instead. HOH 101 blank, more than
  4 A from any atom. Unprocessed model, crystal symmetry of the file.
  """
  if(tyr94_c):
    model = _model_with_tyr94_c()
    ph = model.get_hierarchy()
    cs = model.crystal_symmetry()
  else:
    pi = iotbx.pdb.input(file_name = os.path.join(qr_unit_tests_data,
      "h_altconf_complete.pdb"))
    ph = pi.construct_hierarchy()
    cs = pi.crystal_symmetry()
  o_xyz = (7.697, 5.715, 13.492)
  if(bonded): o_xyz = (12.354+1.3, 8.727, 24.053)
  chain = iotbx.pdb.hierarchy.chain(id="W")
  chain.append_residue_group(_water_rg(100, labels, o_xyz))
  chain.append_residue_group(_water_rg(101, [""], (15.5, 12.0, 20.0)))
  ph.only_model().append_chain(chain)
  ph.atoms().reset_i_seq()
  ph.atoms_reset_serial()
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=cs, log=null_out())

def _residue_i_seqs(ph, resseq, altloc=None):
  return flex.size_t([a.i_seq for a in ph.atoms()
    if a.parent().parent().resseq_as_int() == resseq and
       (altloc is None or a.parent().altloc == altloc)])

def _hand_conformer(ph, label):
  """
  Independent single-conformer hierarchy of conformer label: atom groups
  with another conformer label removed (empty residue groups too), then
  iotbx remove_alt_confs. No fallback needed in the models used here.
  """
  h = ph.deep_copy()
  for rg in h.residue_groups():
    for ag in rg.atom_groups():
      if(ag.altloc not in ("", label)): rg.remove_atom_group(ag)
    if(rg.atom_groups_size() == 0): rg.parent().remove_residue_group(rg)
  h.remove_alt_confs(always_keep_one_conformer=True)
  h.atoms().reset_i_seq()
  return h

def _key(atom):
  return (atom.parent().parent().parent().id,
    atom.parent().parent().resseq_as_int(), atom.name.strip())

def exercise_single_conformer_only_water():
  """
  A residue group without blank atoms, with a proper subset of the model's
  conformer labels and no bonded partner (a water labelled B only, 1PWM)
  exists only in the conformers of its own labels: no fallback.
  """
  model = model_with_b_only_water()
  ph = model.get_hierarchy()
  cs = model.crystal_symmetry()
  w_b, w_0 = _residue_i_seqs(ph, 100), _residue_i_seqs(ph, 101)
  assert w_b.size() == w_0.size() == 3
  sel = conformers.conformer_selections(ph)
  assert [s.label for s in sel] == ["A", "B"]
  assert [s.n_fallback for s in sel] == [0, 0]
  assert [s.n_absent for s in sel] == [1, 0]
  isel_a, isel_b = set(sel[0].iselection), set(sel[1].iselection)
  assert not any(i in isel_a for i in w_b)
  assert all(i in isel_b for i in w_b)
  assert all(i in isel_a and i in isel_b for i in w_0)
  m = conformers.multiplicity(ph)
  assert m.select(w_b).all_eq(1.0) and m.select(w_0).all_eq(2.0)
  # split: A without HOH 100, B with it, HOH 101 in both
  parts = conformers.split(model)
  assert [p.n_absent for p in parts] == [1, 0]
  for p, s in zip(parts, sel):
    h = p.model.get_hierarchy()
    assert list(p.iselection) == list(s.iselection)
    resseqs = set([rg.resseq_as_int() for rg in h.residue_groups()])
    assert (100 in resseqs) == (p.label == "B")
    assert 101 in resseqs
    assert list(h.altloc_indices().keys()) == [""]
    assert h.atoms_size() == s.iselection.size()
    assert h.is_similar_hierarchy(_hand_conformer(ph, p.label))
  # gradients against cctbx on independent single-conformer hierarchies:
  # A atoms get conformer A's gradient (no water term), the B-only water
  # and B atoms conformer B's, blank atoms the mean
  g = gradients_from_conformers(model)
  refs = {}
  for label in ["A", "B"]:
    h = _hand_conformer(ph, label)
    g_h = cctbx_calculator(h, cs)
    refs[label] = dict([(_key(a), g_h[a.i_seq]) for a in h.atoms()])
  n_checked = 0
  for a in ph.atoms():
    k, altloc = _key(a), a.parent().altloc
    if(altloc == ""):
      expected = [(x+y)/2 for x, y in zip(refs["A"][k], refs["B"][k])]
    else:
      expected = refs[altloc][k]
    assert approx_equal(g[a.i_seq], expected, eps=1.e-6), (a.id_str())
    n_checked += 1
  assert n_checked == ph.atoms_size()
  assert ("W", 100, "O") in refs["B"] and ("W", 100, "O") not in refs["A"]
  # the water does touch TYR 94 A: with it in conformer A (the old
  # fallback) the A atoms would differ from the reference
  tyr_a = _residue_i_seqs(ph, 94, "A")
  h_old = _hand_conformer(ph, "A")
  rg_w = [rg for rg in ph.residue_groups() if rg.resseq_as_int() == 100][0]
  chain = [c for c in h_old.chains() if c.id == "W"][0]
  rg_copy = rg_w.detached_copy()
  rg_copy.only_atom_group().altloc = ""
  chain.insert_residue_group(0, rg_copy)
  h_old.atoms().reset_i_seq()
  g_old = cctbx_calculator(h_old, cs)
  old = dict([(_key(a), g_old[a.i_seq]) for a in h_old.atoms()])
  assert not approx_equal([old[_key(ph.atoms()[i])] for i in tyr_a],
    [g[i] for i in tyr_a], out=None)
  # from_conformers reports it
  from six.moves import cStringIO as StringIO
  log = StringIO()
  conformers.from_conformers(model=model, log=log,
    build=lambda conformer_model, working_folder: None)
  text = log.getvalue()
  assert text.find("1 residue groups absent from conformer A") > -1, text
  assert text.find("absent from conformer B") == -1, text
  # bonded to a blank atom (O 1.3 A from CA of GLY 99): the fallback stays
  model = model_with_b_only_water(bonded=True)
  ph = model.get_hierarchy()
  w_b = _residue_i_seqs(ph, 100)
  sel = conformers.conformer_selections(ph)
  assert [s.n_fallback for s in sel] == [1, 0]
  assert [s.n_absent for s in sel] == [0, 0]
  assert all(i in set(sel[0].iselection) for i in w_b)
  assert conformers.multiplicity(ph).select(w_b).all_eq(2.0)
  # three labels (TYR 94 A/B/C), water with B and C: absent from A only;
  # GLY 95 (A/B, bonded) keeps its fallback in C
  model = model_with_b_only_water(labels=["B", "C"], tyr94_c=True)
  ph = model.get_hierarchy()
  sel = conformers.conformer_selections(ph)
  assert [s.label for s in sel] == ["A", "B", "C"]
  assert [s.n_fallback for s in sel] == [0, 0, 1]
  assert [s.n_absent for s in sel] == [1, 0, 0]
  w_b, w_c = _residue_i_seqs(ph, 100, "B"), _residue_i_seqs(ph, 100, "C")
  for s, has_b, has_c in zip(sel, [False, True, False], [False, False, True]):
    isel = set(s.iselection)
    assert all((i in isel) == has_b for i in w_b)
    assert all((i in isel) == has_c for i in w_c)
  m = conformers.multiplicity(ph)
  assert m.select(w_b).all_eq(1.0) and m.select(w_c).all_eq(1.0)
  # a water labelled B 1.0 A from a chloride labelled A (1PWM HOH 1322,
  # CL 1321), away from the protein: each absent from the other conformer
  model = model_with_b_only_water()
  ph = model.get_hierarchy()
  rg_w = [rg for rg in ph.residue_groups() if rg.resseq_as_int() == 100][0]
  for a in rg_w.atoms():
    a.xyz = (a.xyz[0] + 15.5 - 7.697, a.xyz[1] + 12.0 - 5.715,
      a.xyz[2] + 4.0 - 13.492)
  rg_w.parent().remove_residue_group(
    [rg for rg in ph.residue_groups() if rg.resseq_as_int() == 101][0])
  rg_cl = iotbx.pdb.hierarchy.residue_group(resseq=" 102", icode=" ")
  ag = iotbx.pdb.hierarchy.atom_group(altloc="A", resname=" CL")
  a = iotbx.pdb.hierarchy.atom()
  a.name, a.element, a.xyz, a.occ, a.b = "CL  ", "CL", (15.5, 12.0, 5.0), \
    0.8, 30.
  ag.append_atom(a)
  rg_cl.append_atom_group(ag)
  rg_w.parent().append_residue_group(rg_cl)
  # a blank water H 2.1 A from the chloride (1PWM: HIS 110 HD1): the
  # distance test calls it a bond, but bonds of hydrogens do not count
  rg_w.parent().append_residue_group(_water_rg(103, [""], (15.5, 15.05, 5.0)))
  ph.atoms().reset_i_seq()
  h2 = [a for a in ph.atoms() if a.parent().parent().resseq_as_int() == 103
    and a.name.strip() == "H2"][0]
  assert abs(h2.distance(ag.atoms()[0]) - 2.13) < 0.01
  sel = conformers.conformer_selections(ph)
  assert [s.n_absent for s in sel] == [1, 1]
  assert [s.n_fallback for s in sel] == [0, 0]
  cl, w_b = _residue_i_seqs(ph, 102), _residue_i_seqs(ph, 100)
  assert ph.atoms()[cl[0]].distance(ph.atoms()[w_b[0]]) < 1.01
  assert set(cl).issubset(set(sel[0].iselection))
  assert not set(cl) & set(sel[1].iselection)
  assert set(w_b).issubset(set(sel[1].iselection))
  assert not set(w_b) & set(sel[0].iselection)
  # the same chloride blank: the water is bonded to it, the fallback stays
  ag.altloc = ""
  sel = conformers.conformer_selections(ph)
  assert [s.n_absent for s in sel] == [0, 0]
  assert [s.n_fallback for s in sel] == [1, 0]
  # a water with every model label (A and B, no blank atoms) is in both
  model = model_with_b_only_water(labels=["A", "B"])
  sel = conformers.conformer_selections(model.get_hierarchy())
  assert [s.n_absent for s in sel] == [0, 0]
  assert [s.n_fallback for s in sel] == [0, 0]

def _one_atom_rg(resseq, label, xyz):
  """LIG residue group of one C atom in one atom group (altloc label)."""
  rg = iotbx.pdb.hierarchy.residue_group(resseq="%4d" % resseq, icode=" ")
  ag = iotbx.pdb.hierarchy.atom_group(altloc=label, resname="LIG")
  a = iotbx.pdb.hierarchy.atom()
  a.name, a.element, a.xyz, a.occ, a.b = " C1 ", " C", xyz, 0.5, 30.
  a.hetero = True
  ag.append_atom(a)
  rg.append_atom_group(ag)
  return rg

def exercise_single_conformer_only_pair():
  """
  Two residue groups labelled B only, bonded to each other (1.5 A) and to
  nothing else (a B-only two-residue ligand or glycan chain): the pair is
  absent from A as a whole. Control: one member also bonded to a blank
  residue group, the whole pair falls back.
  """
  model = model_with_b_only_water()
  ph = model.get_hierarchy()
  chain = [c for c in ph.chains() if c.id == "W"][0]
  chain.remove_residue_group(
    [rg for rg in chain.residue_groups() if rg.resseq_as_int() == 100][0])
  chain.append_residue_group(_one_atom_rg(104, "B", (15.5, 12.0, 5.0)))
  chain.append_residue_group(_one_atom_rg(105, "B", (17.0, 12.0, 5.0)))
  ph.atoms().reset_i_seq()
  sel = conformers.conformer_selections(ph)
  assert [s.n_absent for s in sel] == [2, 0], [s.n_absent for s in sel]
  assert [s.n_fallback for s in sel] == [0, 0]
  pair = _residue_i_seqs(ph, 104).concatenate(_residue_i_seqs(ph, 105))
  assert not set(pair) & set(sel[0].iselection)
  assert set(pair).issubset(set(sel[1].iselection))
  assert conformers.multiplicity(ph).select(pair).all_eq(1.0)
  # control: residue 105 bonded (1.5 A) to a blank residue group
  chain.append_residue_group(_one_atom_rg(106, "", (18.5, 12.0, 5.0)))
  ph.atoms().reset_i_seq()
  sel = conformers.conformer_selections(ph)
  assert [s.n_absent for s in sel] == [0, 0], [s.n_absent for s in sel]
  assert [s.n_fallback for s in sel] == [2, 0]
  assert set(pair).issubset(set(sel[0].iselection))

pdb_str_link_guard = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N  AALA A   1      -1.521   3.898   0.713  0.40 20.00           N
ATOM      2  CA AALA A   1      -1.879   2.678   1.370  0.40 20.00           C
ATOM      3  CB AALA A   1      -3.346   2.715   1.747  0.40 20.00           C
ATOM      4  C  AALA A   1      -1.551   1.486   0.458  0.40 20.00           C
ATOM      5  O  AALA A   1      -2.007   1.439  -0.662  0.40 20.00           O
ATOM      6  N  BALA A   1      -1.521   3.898   4.713  0.40 20.00           N
ATOM      7  CA BALA A   1      -1.879   2.678   5.370  0.40 20.00           C
ATOM      8  CB BALA A   1      -3.346   2.715   5.747  0.40 20.00           C
ATOM      9  C  BALA A   1      -1.551   1.486   4.458  0.40 20.00           C
ATOM     10  O  BALA A   1      -2.007   1.439   3.338  0.40 20.00           O
ATOM     11  N  CALA A   1      -1.521   3.898   4.713  0.20 20.00           N
ATOM     12  CA CALA A   1      -1.879   2.678   5.370  0.20 20.00           C
ATOM     13  CB CALA A   1      -3.346   2.715   5.747  0.20 20.00           C
ATOM     14  C  CALA A   1      -1.551   1.486   4.458  0.20 20.00           C
ATOM     15  O  CALA A   1      -2.007   1.439   3.338  0.20 20.00           O
ATOM     16  N  %sALA A   2      -0.744   0.404   0.952  0.50 20.00           N
ATOM     17  CA %sALA A   2      -0.493  -0.749   0.119  0.50 20.00           C
ATOM     18  CB %sALA A   2      -1.476  -1.885   0.331  0.50 20.00           C
ATOM     19  C  %sALA A   2       0.972  -1.193   0.100  0.50 20.00           C
ATOM     20  O  %sALA A   2       1.603  -1.231   1.129  0.50 20.00           O
ATOM     21  N  %sALA A   2      -0.744   0.404   4.952  0.50 20.00           N
ATOM     22  CA %sALA A   2      -0.493  -0.749   4.119  0.50 20.00           C
ATOM     23  CB %sALA A   2      -1.476  -1.885   4.331  0.50 20.00           C
ATOM     24  C  %sALA A   2       0.972  -1.193   4.100  0.50 20.00           C
ATOM     25  O  %sALA A   2       1.603  -1.231   5.129  0.50 20.00           O
END
"""

def exercise_link_guard():
  """
  ALA-ALA dipeptide (as tst_41.pdb_str_swapped): residue 1 under A/B/C with
  B and C the A conformer shifted 4 A along z; residue 2 under A/B only.
  With residue 2's B shifted, conformer C of residue 1 bonds to residue 2's
  B, but residue 2 takes its first label A into conformer C: the peptide
  bond would be broken (Sorry). With residue 2's A shifted: no Sorry.
  check_all_links (the clustering route's check) raises the same Sorry.
  """
  for unshifted, shifted, fails in [("A", "B", True), ("B", "A", False)]:
    pdb_str = pdb_str_link_guard % tuple([unshifted]*5 + [shifted]*5)
    pi = iotbx.pdb.input(source_info=None, lines=pdb_str)
    model = mmtbx.model.manager(model_input=pi, log=null_out())
    sel = conformers.conformer_selections(model.get_hierarchy())
    assert [s.n_fallback for s in sel] == [0, 0, 1]
    message = None
    try:
      conformers.check_all_links(model.get_hierarchy())
    except Sorry as e:
      message = str(e)
    assert (message is not None) == fails
    try:
      parts = conformers.split(model)
    except Sorry as e:
      assert fails
      e = str(e)
      assert e == message, (e, message)
      assert e.find("conformer C") > -1, e
      assert e.find("residue") > -1 and e.find("not bonded") > -1, e
      assert e.find(" ALA A   1 ") > -1 and e.find(" ALA A   2 ") > -1, e
      assert e.find(" C ") > -1 and e.find(" N ") > -1, e
    else:
      assert not fails
      assert len(parts) == 3

def exercise_check_completed():
  """
  check_completed: Sorry naming the atoms reduce2 qmready would complete.
  TYR 94 B made identical to A, then its OH removed: OH of A is bonded to
  CZ of B, so it would be shared.
  """
  file_name = os.path.join(qr_unit_tests_data, "h_altconf_complete.pdb")
  ph = iotbx.pdb.input(file_name = file_name).construct_hierarchy()
  conformers.check_completed(ph) # complete model: no exception
  rg = [rg for rg in ph.residue_groups() if rg.resseq_as_int() == 94][0]
  ag_a, ag_b = rg.atom_groups()
  xyz_a = dict([(a.name, a.xyz) for a in ag_a.atoms()])
  for a in ag_b.atoms():
    a.xyz = xyz_a[a.name]
  ag_b.remove_atom([a for a in ag_b.atoms() if a.name.strip() == "OH"][0])
  ph.atoms().reset_i_seq()
  try:
    conformers.check_completed(ph)
  except Sorry as e:
    e = str(e)
    assert e.find("mmtbx.reduce2 qmready=true") > -1, e
    assert e.find("OH") > -1, e
    assert e.find("Altloc completion needed: 1 atoms to share and 0 to "
      "copy ( TYR A  94  OH)") > -1, e
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
  # A hydrogen-only label is a conformer like any other (a methyl in two
  # orientations, 6uf7 ILE 4), unless its atoms are H/D exchange sites: an
  # H and a D of equal name.strip()[1:] in atom groups with different
  # altlocs (iotbx exchangeable_hd_selections).
  rg = list(ph.residue_groups())[1] # GLY 88: first residue with atom H
  ag = rg.only_atom_group()
  h_atom = [a for a in ag.atoms() if a.name.strip() == "H"][0]
  ag_h = iotbx.pdb.hierarchy.atom_group(altloc="C", resname=ag.resname)
  ag_h.append_atom(h_atom.detached_copy())
  ag.remove_atom(h_atom)
  rg.append_atom_group(ag_h)
  ph.atoms().reset_i_seq()
  assert conformers.conformer_labels(ph) == ["A", "B", "C"]
  # the H/D partner under label D makes C and D exchange sites, no conformers
  ag_d = iotbx.pdb.hierarchy.atom_group(altloc="D", resname=ag.resname)
  d_atom = h_atom.detached_copy()
  d_atom.name = " D  "
  d_atom.element = " D"
  ag_d.append_atom(d_atom)
  rg.append_atom_group(ag_d)
  ph.atoms().reset_i_seq()
  assert conformers.conformer_labels(ph) == ["A", "B"]
  assert conformers.multiplicity(ph)[ag_h.atoms()[0].i_seq] == 2.0
  for p in conformers.split(model):
    assert ag_h.atoms()[0].i_seq in p.iselection
  # two H orientations of one methyl under A and B (both H): conformers
  rg.remove_atom_group(ag_d)
  ag_h.altloc = "A"
  h2 = h_atom.detached_copy()
  h2.xyz = (h2.xyz[0] + 0.5, h2.xyz[1], h2.xyz[2])
  ag_b = iotbx.pdb.hierarchy.atom_group(altloc="B", resname=ag.resname)
  ag_b.append_atom(h2)
  rg.append_atom_group(ag_b)
  ph.atoms().reset_i_seq()
  assert conformers.conformer_labels(ph) == ["A", "B"]
  assert conformers.multiplicity(ph)[ag_h.atoms()[0].i_seq] == 1.0
  assert conformers.multiplicity(ph)[ag_b.atoms()[0].i_seq] == 1.0

def exercise_super_selections():
  # h_altconf_complete.pdb in the small box of tst_46.get_model(fn, 2.0)
  # (several images in the sphere); for every label the tiled sphere
  # selection equals the per-residue rule applied directly to the sphere
  # hierarchy (same label set, conformer_selections on the sphere)
  from cctbx.crystal import super_cell
  # tst_46 imports tst_11, which imports tst_31: import it here
  from qrefine.tests.unit import tst_46
  model = tst_46.get_model("h_altconf_complete.pdb", contact=2.0)
  ph = model.get_hierarchy()
  scm = super_cell.manager(pdb_hierarchy=ph,
    crystal_symmetry=model.crystal_symmetry(), select_within_radius=10)
  sels = conformers.conformer_selections(ph)
  ss = conformers.super_selections(sels, ph.atoms_size(), scm)
  direct = dict([(s.label, s.iselection) for s in
    conformers.conformer_selections(scm.super_sphere_hierarchy)])
  assert sorted(ss.keys()) == ["A", "B"]
  for label in ss:
    assert ss[label].sphere.size() == scm.super_sphere_hierarchy.atoms_size()
    assert list(ss[label].sphere.iselection()) == list(direct[label])
    assert ss[label].cell.size() == scm.super_cell_hierarchy.atoms_size()
    # master block first
    assert list(ss[label].cell[:ph.atoms_size()].iselection()) == \
      list([s for s in sels if s.label == label][0].iselection)

def exercise_single_conformer_hierarchy():
  # public helper: same result as split() for every conformer
  pi = iotbx.pdb.input(os.path.join(qr_unit_tests_data,
    "h_altconf_complete.pdb"))
  model = mmtbx.model.manager(model_input=pi, log=null_out())
  ph = model.get_hierarchy()
  labels = set(conformers.conformer_labels(ph))
  for c, s in zip(conformers.split(model), conformers.conformer_selections(ph)):
    h = conformers.single_conformer_hierarchy(ph, s.iselection, s.label, labels)
    assert h.is_similar_hierarchy(c.model.get_hierarchy())
    assert list(h.altloc_indices().keys()) == [""]

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
  exercise_fallback()
  exercise_single_conformer_only_water()
  exercise_single_conformer_only_pair()
  exercise_link_guard()
  exercise_check_completed()
  exercise_microheterogeneity()
  exercise_single_conformer_hierarchy()
  exercise_super_selections()

if(__name__ == '__main__'):
  run()
