from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import re
import libtbx.load_env
import iotbx.pdb
import mmtbx.model
from cctbx import crystal
from libtbx import group_args
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex
from qrefine import conformers, fragment, fragment_capping
from qrefine.fragment import fragments
from qrefine.super_cell import expand
from qrefine.tests.unit import altloc_scenarios

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine, "tests", "unit", "data_files")

def boxed(pdb_hierarchy, gap=30.0):
  """
  Unprocessed model of pdb_hierarchy in a P1 box of its extent plus gap
  (tst_46.get_model style): a gap of 30 A puts no image within the 10 A
  super sphere.
  """
  xyz = pdb_hierarchy.atoms().extract_xyz()
  box = [b-a+gap for a, b in zip(xyz.min(), xyz.max())]
  cs = crystal.symmetry(unit_cell=box+[90, 90, 90], space_group_symbol="P1")
  pdb_hierarchy.atoms().reset_i_seq()
  pdb_hierarchy.atoms_reset_serial()
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=pdb_hierarchy,
    crystal_symmetry=cs, log=null_out())

def peptide():
  """
  Unprocessed model of altloc_scenarios.base_model() without the waters.
  """
  ph = altloc_scenarios.base_model().get_hierarchy()
  ph = ph.select(flex.bool([a.parent().resname != "HOH" for a in ph.atoms()]),
    copy_atoms=True)
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=altloc_scenarios.base_model().crystal_symmetry(),
    log=null_out())

def three_label_model():
  """
  The peptide in a P1 box (no images): ASN 2 with side-chain conformers A,
  B (chi1 +0, +120), GLN 5 with B, C (chi1 +0, -60: +120 puts OE1 0.65 A
  from TYR 7 OH), the other residues blank. Global labels A, B, C.
  """
  m = peptide()
  altloc_scenarios.side_chain_conformers(altloc_scenarios.residue(m, 2),
    ["A", "B"], [0.5, 0.5], [0, 120])
  altloc_scenarios.side_chain_conformers(altloc_scenarios.residue(m, 5),
    ["B", "C"], [0.5, 0.5], [0, -60])
  return boxed(m.get_hierarchy())

def hetero_model():
  """
  The peptide plus water 8 with H1, H2, moved 15 A away from the peptide
  (its fragment is the water alone), in a P1 box; the water atoms are
  label A and a copy shifted 0.5 A along x is label B.
  """
  ph = altloc_scenarios.base_model().get_hierarchy()
  far = ph.select(flex.bool([a.parent().resname != "HOH"
    for a in ph.atoms()])).atoms().extract_xyz().max()
  keep = flex.bool([a.parent().resname != "HOH" or
    a.parent().parent().resseq_as_int() == 8 for a in ph.atoms()])
  ph = ph.select(keep, copy_atoms=True)
  rg = [rg for rg in ph.residue_groups() if rg.resseq_as_int() == 8][0]
  ag = rg.only_atom_group()
  o = ag.only_atom()
  o.xyz = (far[0]+15, far[1], far[2])
  for name, d in [("H1", (0.96, 0, 0)), ("H2", (-0.24, 0.93, 0))]:
    h = o.detached_copy()
    h.name = " %s " % name
    h.element = " H"
    h.xyz = tuple([x+dx for x, dx in zip(o.xyz, d)])
    ag.append_atom(h)
  ph.atoms().reset_i_seq()
  altloc_scenarios.label_atoms(rg, ["O", "H1", "H2"], "A", 0.5)
  for ag in rg.atom_groups():
    if(ag.atoms_size() == 0): rg.remove_atom_group(ag)
  ag_b = rg.only_atom_group().detached_copy()
  ag_b.altloc = "B"
  for a in ag_b.atoms():
    a.xyz = (a.xyz[0]+0.5, a.xyz[1], a.xyz[2])
  rg.append_atom_group(ag_b)
  return boxed(ph)

def make_fragments(model, with_conformers=True, cap_fragments=True):
  ph = model.get_hierarchy()
  selections = None
  if(with_conformers):
    selections = conformers.conformer_selections(ph)
  return fragments(working_folder=None, make_working_folder=False,
    clustering_method="bcc", maxnum_residues_in_cluster=1,
    pdb_hierarchy=ph, crystal_symmetry=model.crystal_symmetry(),
    qm_engine_name="mopac", cap_fragments=cap_fragments,
    conformer_selections=selections)

def jobs_of(f, i):
  return [j for j, k in enumerate(f.variant_fragment_index) if k == i]

def residues_of(h):
  return set([rg.resseq_as_int() for rg in h.residue_groups()])

def labelled(h):
  """(resseq, altloc) of the atoms of h with a non-blank altloc"""
  return set([(a.parent().parent().resseq_as_int(), a.parent().altloc)
    for a in h.atoms() if a.parent().altloc.strip() != ""])

def check_common(f, n_conf):
  """
  Brief assertions: weights per fragment sum to n_conf, at most labels+1
  variants (and at most n_conf), relabeled variant hierarchies of the
  variant's atoms, cappers and charges present.
  """
  n = f.pdb_hierarchy.atoms_size()
  lists = [f.variant_fragment_index, f.variant_weights, f.variant_labels,
    f.fragment_selections, f.fragment_super_selections, f.cluster_selections,
    f.buffer_selections, f.fragment_cappers, f.fragment_charges,
    f.variant_hierarchies]
  assert len(set([len(l) for l in lists])) == 1
  assert f.n_fragments == len(f.clusters)
  assert len(f.fragment_cluster_selections) == f.n_fragments
  assert len(f.fragment_full_super_selections) == f.n_fragments
  assert sorted(set(f.variant_fragment_index)) == list(range(f.n_fragments))
  for i in range(f.n_fragments):
    jobs = jobs_of(f, i)
    assert sum(f.variant_weights[j] for j in jobs) == n_conf
    full = f.fragment_full_super_selections[i]
    labels = set(conformers.conformer_labels(
      f.pdb_hierarchy_super.select(full)))
    assert len(jobs) <= len(labels) + 1 and len(jobs) <= n_conf
    if(len(labels) == 0):
      assert len(jobs) == 1 and f.variant_labels[jobs[0]] == ""
    union = flex.bool(full.size(), False)
    for j in jobs:
      v = f.fragment_super_selections[j]
      assert (v & ~full).count(True) == 0
      union = union | v
      assert list(f.fragment_selections[j]) == list(v[:n])
      assert list(f.cluster_selections[j]) == list(
        f.fragment_cluster_selections[i])
      h = f.variant_hierarchies[j]
      assert list(h.altloc_indices().keys()) == [""]
      assert max([rg.atom_groups_size() for rg in h.residue_groups()]) == 1
      assert h.atoms_size() == v.count(True)
      assert list(h.atoms().extract_name()) == list(
        f.pdb_hierarchy_super.select(v).atoms().extract_name())
      assert h.atoms().extract_xyz().max_distance(
        f.pdb_hierarchy_super.select(v).atoms().extract_xyz()) < 1.e-6
      if(f.cap_fragments):
        assert f.fragment_cappers[j] is not None
        assert f.fragment_charges[j] is not None
      assert (len(jobs) == 1) == (f.variant_labels[j] == "")
    assert list(union) == list(full)

def expected_variants(residues):
  """
  Hand-derived variants of a fragment of three_label_model with these
  residues: [(label, weight, labelled atoms)]. Conformer X takes the X atoms
  of a residue with X, else (fallback) those of its first own label; equal
  atom sets collapse, the smallest producing label represents them.
  """
  has_2, has_5 = 2 in residues, 5 in residues
  if(has_2 and has_5):
    return [("A", 1, set([(2, "A"), (5, "B")])),
            ("B", 1, set([(2, "B"), (5, "B")])),
            ("C", 1, set([(2, "A"), (5, "C")]))]
  if(has_2):
    return [("A", 2, set([(2, "A")])), ("B", 1, set([(2, "B")]))]
  if(has_5):
    return [("A", 2, set([(5, "B")])), ("C", 1, set([(5, "C")]))]
  return [("", 3, set())]

def exercise_three_labels():
  model = three_label_model()
  f = make_fragments(model)
  check_common(f, n_conf=3)
  assert f.pdb_hierarchy_super.atoms_size() == f.pdb_hierarchy.atoms_size()
  # conformer spheres are built for charge embedding only
  assert f.conformer_spheres == {}
  n_jobs = 0
  cases = set()
  for i in range(f.n_fragments):
    residues = residues_of(f.pdb_hierarchy_super.select(
      f.fragment_full_super_selections[i]))
    expected = expected_variants(residues)
    cases.add((2 in residues, 5 in residues))
    n_jobs += len(expected)
    got = []
    for j in jobs_of(f, i):
      got.append((f.variant_labels[j], f.variant_weights[j], labelled(
        f.pdb_hierarchy_super.select(f.fragment_super_selections[j]))))
    assert sorted(got) == sorted(expected), (i, residues, got)
  # fragments with residues 2 and 5, 2 only, 5 only (every fragment of this
  # peptide reaches 2 or 5; hetero_model has fragments without altlocs)
  assert cases == set([(True, True), (True, False), (False, True)]), cases
  assert len(f.variant_weights) == n_jobs

def residue_keys(h):
  return set([(rg.parent().id, rg.resid()) for rg in h.residue_groups()])

def select_residues(h, keys):
  sel = flex.bool()
  for rg in h.residue_groups():
    sel.extend(flex.bool(rg.atoms_size(), (rg.parent().id, rg.resid()) in keys))
  return sel

def exercise_variant_atoms_match_split():
  """
  Every variant's master atoms are exactly the atoms of the fragment's
  residues in the conformer model conformers.split gives for its label
  (any label for a single-variant fragment: all give the same atoms).
  """
  model = three_label_model()
  f = make_fragments(model)
  parts = dict([(p.label, p) for p in conformers.split(model)])
  for j, label in enumerate(f.variant_labels):
    label = label or "A"
    p = parts[label]
    h_c = p.model.get_hierarchy()
    i = f.variant_fragment_index[j]
    keys = residue_keys(f.pdb_hierarchy.select(f.fragment_selections[j]))
    sel = select_residues(h_c, keys)
    expected = p.iselection.select(sel.iselection())
    assert list(f.fragment_selections[j].iselection()) == list(expected)
    h = f.variant_hierarchies[j]
    assert h.is_similar_hierarchy(h_c.select(sel))
    assert [a.parent().resname for a in h.atoms()] == \
      [a.parent().resname for a in h_c.select(sel).atoms()]

def backbone_model():
  """
  The peptide in a P1 box with backbone conformers A (the original) and B
  (rotated 30 degrees about C 2 ... N 4) of ASN 3: the C of ASN 3 and the N
  of GLN 4 (blank) are bonded under both labels.
  """
  m = peptide()
  altloc_scenarios.backbone_conformers(m, 3, "B", 0.5, 30)
  return boxed(m.get_hierarchy())

def check_capper_parent(model):
  """
  For every variant: as many caps (and the same charge) as
  fragment_capping.capper gives for the same residue groups selected in the
  super cell of the conformer model of its label (A, the smallest, for a
  single-variant fragment: every label gives its atoms). Returns the number
  of variants that capping against the altloc super cell (the wrong parent)
  caps more: a cut bond to an atom of the other conformer.
  """
  f = make_fragments(model)
  cs = model.crystal_symmetry()
  expansions = {}
  for p in conformers.split(model):
    expansions[p.label] = expand(pdb_hierarchy=p.model.get_hierarchy(),
      crystal_symmetry=cs, select_within_radius=10)
  SCM = f.expansion.SCM
  keep = SCM.keep_selection.iselection()
  n_more_with_altloc_parent = 0
  for j, label in enumerate(f.variant_labels):
    h = f.variant_hierarchies[j]
    e = expansions[label or "A"]
    cell = e.SCM.super_cell_hierarchy
    sel = select_residues(cell, residue_keys(h))
    c = fragment_capping.capper(parent_hierarchy=cell, selection=sel,
      crystal_symmetry=e.cs_box)
    assert c.fragment_hierarchy.is_similar_hierarchy(h)
    assert c.n_caps == f.fragment_cappers[j].n_caps, (j, label)
    assert c.charge == f.fragment_charges[j]
    # the capper of the variant against the altloc super cell
    sel_alt = flex.bool(SCM.super_cell_hierarchy.atoms_size(),
      keep.select(f.fragment_super_selections[j].iselection()))
    try:
      c_alt = fragment_capping.capper(
        parent_hierarchy=SCM.super_cell_hierarchy, selection=sel_alt,
        crystal_symmetry=f.expansion.cs_box)
    except Sorry:
      continue # two labels in one fragment, or a cap clash
    assert c_alt.n_caps >= c.n_caps
    if(c_alt.n_caps > c.n_caps): n_more_with_altloc_parent += 1
  return n_more_with_altloc_parent

def exercise_capper_parent():
  """
  Side-chain conformers (three_label_model) have no bonds to other residues:
  either parent gives the same caps. Backbone conformers of ASN 3: a
  fragment ending at GLN 4 (or starting after ASN 3) cut from the altloc
  super cell would cap its N twice, for C of A and C of B.
  """
  assert check_capper_parent(three_label_model()) == 0
  assert check_capper_parent(backbone_model()) > 0

def exercise_single_conformer_unchanged():
  """
  A model without altlocs (P 1 21 1, with images): one variant per
  fragment, weight 1, label "", the selection lists of the rule before
  variants (fragment_super_atoms), and the same lists with
  conformer_selections=[] as with None.
  """
  model = altloc_scenarios.base_model()
  ph = model.get_hierarchy()
  assert conformers.conformer_selections(ph) == []
  f0 = make_fragments(model, with_conformers=False)
  f1 = make_fragments(model, with_conformers=True)
  n, n_super = ph.atoms_size(), f0.pdb_hierarchy_super.atoms_size()
  assert n_super > n
  assert list(f0.conformer_spheres.keys()) == [""]
  assert f0.conformer_spheres[""] is f0.pdb_hierarchy_super
  for f in [f0, f1]:
    check_common(f, n_conf=1)
    assert f.n_fragments == len(f.fragment_super_atoms)
    assert f.variant_fragment_index == list(range(f.n_fragments))
    assert f.variant_weights == [1]*f.n_fragments
    assert f.variant_labels == [""]*f.n_fragments
    for i in range(f.n_fragments):
      atoms = f.fragment_super_atoms[i]
      sel = flex.bool(n_super, False)
      sel.set_selected(flex.size_t([k-1 for k in atoms]), True)
      assert list(f.fragment_super_selections[i]) == list(sel)
      assert list(f.fragment_full_super_selections[i]) == list(sel)
      assert list(f.fragment_selections[i]) == list(sel[:n])
      assert f.variant_hierarchies[i].is_similar_hierarchy(
        f.pdb_hierarchy_super.select(sel))
  fe0, fe1 = f0.get_fragment_extracts(), f1.get_fragment_extracts()
  for name in ["fragment_selections", "fragment_super_selections",
               "cluster_selections", "buffer_selections"]:
    l0, l1 = getattr(fe0, name), getattr(fe1, name)
    assert [list(s) for s in l0] == [list(s) for s in l1], name
  assert fe0.fragment_charges == fe1.fragment_charges
  assert [c.n_caps for c in fe0.fragment_cappers] == \
    [c.n_caps for c in fe1.fragment_cappers]
  for i in range(f0.n_fragments):
    assert fragment.job_folder_name(fe0, i) == str(i)

def exercise_copy_in_buffer():
  """
  h_altconf_complete in a box with images 2 A away: a fragment whose master
  residues carry no altloc but whose buffer holds an image of an altloc
  residue has two variants (labels A, B, weight 1 each).
  """
  from qrefine.tests.unit.tst_46 import get_model
  model = get_model("h_altconf_complete.pdb", contact=2.0)
  f = make_fragments(model, cap_fragments=False)
  check_common(f, n_conf=2)
  n = f.pdb_hierarchy.atoms_size()
  n_found = 0
  for i in range(f.n_fragments):
    full = f.fragment_full_super_selections[i]
    master = f.pdb_hierarchy_super.select(full.iselection().select(
      full.iselection() < n))
    images = f.pdb_hierarchy_super.select(full.iselection().select(
      full.iselection() >= n))
    if(conformers.conformer_labels(master) != []): continue
    if(conformers.conformer_labels(images) == []):
      assert len(jobs_of(f, i)) == 1
      continue
    jobs = jobs_of(f, i)
    assert [f.variant_labels[j] for j in jobs] == ["A", "B"]
    assert [f.variant_weights[j] for j in jobs] == [1, 1]
    for j in jobs:
      assert list(f.fragment_selections[j]) == list(full[:n])
    n_found += 1
  assert n_found > 0

def exercise_hetero():
  """
  Water with conformers A and B: its fragment has two variants of three
  atoms each, capped with zero caps and charge 0.
  """
  model = hetero_model()
  f = make_fragments(model)
  check_common(f, n_conf=2)
  n_found = 0
  for i in range(f.n_fragments):
    h = f.pdb_hierarchy_super.select(f.fragment_full_super_selections[i])
    if(residues_of(h) != set([8])): continue
    jobs = jobs_of(f, i)
    assert [f.variant_labels[j] for j in jobs] == ["A", "B"]
    for j, label in zip(jobs, ["A", "B"]):
      assert labelled(f.pdb_hierarchy_super.select(
        f.fragment_super_selections[j])) == set([(8, label)])
      assert f.variant_hierarchies[j].atoms_size() == 3
      assert f.fragment_cappers[j].n_caps == 0
      assert f.fragment_charges[j] == 0
    n_found += 1
  assert n_found == 1

def exercise_micro():
  """
  Microheterogeneity (tst_31.exercise_microheterogeneity): GLY 95 of
  h_altconf_complete has blank, A and B atom groups, the B group renamed
  ALA. In every variant containing residue 95, it is ALA in label B
  variants and GLY in label A variants.
  """
  pi = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,
    "h_altconf_complete.pdb"))
  ph = pi.construct_hierarchy()
  rg = list(ph.residue_groups())[8]
  assert [ag.altloc for ag in rg.atom_groups()] == ["", "A", "B"]
  rg.atom_groups()[-1].resname = "ALA"
  ph.atoms().reset_i_seq()
  model = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=pi.crystal_symmetry(), log=null_out())
  f = make_fragments(model, cap_fragments=False)
  check_common(f, n_conf=2)
  n_checked = 0
  for j, label in enumerate(f.variant_labels):
    rgs = [rg for rg in f.variant_hierarchies[j].residue_groups()
      if rg.resseq_as_int() == 95]
    if(len(rgs) == 0): continue
    assert label in ["A", "B"], label
    expected = {"A": "GLY", "B": "ALA"}[label]
    for rg in rgs:
      assert rg.only_atom_group().resname == expected
      assert rg.atoms_size() == 7
    n_checked += 1
  assert n_checked > 0

def exercise_job_folder_name():
  fe = group_args(variant_fragment_index=[0, 1, 1, 2],
    variant_labels=["", "A", "C", ""])
  names = [fragment.job_folder_name(fe, j) for j in range(4)]
  assert names == ["0", "1_A", "1_C", "2"], names
  f = make_fragments(three_label_model())
  fe = f.get_fragment_extracts()
  for name in ["variant_fragment_index", "variant_weights", "variant_labels",
               "variant_hierarchies", "conformer_spheres", "n_fragments",
               "fragment_cluster_selections",
               "fragment_full_super_selections"]:
    assert getattr(fe, name) is getattr(f, name), name
  for j in range(len(f.variant_labels)):
    i = f.variant_fragment_index[j]
    expected = str(i)
    if(len(jobs_of(f, i)) > 1):
      expected = "%d_%s" % (i, f.variant_labels[j])
    assert fragment.job_folder_name(fe, j) == expected

def exercise_charge_embedding():
  """
  Point charges of a variant (xtb, charge embedding; no engine is run) come
  from its conformer's relabeled super sphere minus the variant: no altloc,
  no residue of the variant, every atom one of the conformer model.
  """
  model = three_label_model()
  ph = model.get_hierarchy()
  f = fragments(working_folder="ase", clustering_method="bcc",
    maxnum_residues_in_cluster=1, pdb_hierarchy=ph,
    crystal_symmetry=model.crystal_symmetry(), qm_engine_name="xtb",
    cap_fragments=True, charge_embedding=True,
    conformer_selections=conformers.conformer_selections(ph))
  assert sorted(f.conformer_spheres.keys()) == ["A", "B", "C"]
  parts = dict([(p.label, p) for p in conformers.split(model)])
  fe = f.get_fragment_extracts()
  for j, label in enumerate(f.variant_labels):
    name = fragment.job_folder_name(fe, j)
    file_name = fragment.write_mm_charge_file(fe, j)
    assert os.path.basename(file_name) == name + "_qxyz_cctbx.dat"
    env = iotbx.pdb.input(file_name=os.path.join("ase", name,
      name + "_mm.pdb")).construct_hierarchy()
    assert list(env.altloc_indices().keys()) == [""]
    assert env.atoms_size() > 0
    assert residues_of(env) & residues_of(f.variant_hierarchies[j]) == set()
    h_c = parts[label or "A"].model.get_hierarchy()
    xyz = dict([((a.parent().parent().resseq_as_int(), a.name), a.xyz)
      for a in h_c.atoms()])
    for a in env.atoms():
      key = (a.parent().parent().resseq_as_int(), a.name)
      assert a.distance(xyz[key]) < 1.e-3, (name, key)

def exercise_sorry_prefix():
  """
  Residue 5 (labels B, C) renamed XXX (no dictionary): the Sorry of the
  first variant capped names its fragment and conformer (A: the first
  variant of a fragment with residue 5 is that of label A, which takes B
  there by the fallback).
  """
  model = three_label_model()
  for ag in altloc_scenarios.residue(model, 5).atom_groups():
    ag.resname = "XXX"
  try:
    make_fragments(model)
  except Sorry as e:
    msg = str(e)
    # the capper interprets the variant first: pdb_interpretation refuses
    # the unknown residue before the capper's own no-dictionary check
    assert re.match(r"fragment \d+, conformer [ABC]: (Capping: no "
      r"dictionary for XXX|Fatal problems interpreting model file)", msg), msg
  else:
    raise AssertionError("no Sorry")

def sphere_keys(f, original, selection):
  """
  (copy index, original atom index) of every atom of selection (flex.bool
  over the super sphere of fragments f): the super cell is the master
  followed by copies of the whole model, so cell position k*n+i is copy k
  of model atom i; original maps model atoms to the single-conformer
  model's (None: identity).
  """
  n = f.pdb_hierarchy.atoms_size()
  cell = f.expansion.SCM.keep_selection.iselection()
  result = set()
  for i in selection.iselection():
    k, j = divmod(cell[i], n)
    if(original is not None): j = original[j]
    result.add((k, j))
  return frozenset(result)

def exercise_interactions_match_single_conformer():
  """
  h_altconf_complete without altlocs, images within 2 A, and its copies
  with residue 2 duplicated into identical A/B, residues 2 and 3 into
  A/B/C: fragments with conformer selections on a duplicate gives the
  interaction list, bonded pairs, clusters and full fragment selections
  (through the original atoms, images included) of the single-conformer
  model. (pair_interaction on the altloc hierarchy would lose the
  interactions of the duplicated residues.)
  """
  from qrefine.tests.unit.tst_31 import duplicate
  from qrefine.tests.unit.tst_65 import contact_model
  ph = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,
    "h_altconf_complete.pdb")).construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  single = contact_model(ph, contact=2.0, name="single_conformer")
  f_s = make_fragments(single, with_conformers=False, cap_fragments=False)
  def fragment_keys(f, original):
    return sorted([sphere_keys(f, original, s)
      for s in f.fragment_full_super_selections], key=sorted)
  keys_s = fragment_keys(f_s, None)
  for labels, residue_indices in [(["A","B"], [2]), (["A","B","C"], [2,3])]:
    dup, original = duplicate(model=single, labels=labels,
      residue_indices=residue_indices)
    f_d = make_fragments(dup, cap_fragments=False)
    assert f_d.labels == labels
    assert sorted(f_d.interaction_list) == sorted(f_s.interaction_list)
    assert f_d.bonded_pairs == f_s.bonded_pairs
    assert sorted([sorted(c) for c in f_d.clusters]) == \
      sorted([sorted(c) for c in f_s.clusters])
    assert sphere_keys(f_d, original, flex.bool(
      f_d.pdb_hierarchy_super.atoms_size(), True)) == sphere_keys(f_s, None,
      flex.bool(f_s.pdb_hierarchy_super.atoms_size(), True))
    assert fragment_keys(f_d, original) == keys_s
    assert f_d.conformer_spheres == {} # released after setup
    print("duplicate %s of residues %s: %d interactions, %d clusters, "
      "%d fragments as single conformer" % ("".join(labels),
      residue_indices, len(f_d.interaction_list), len(f_d.clusters),
      len(keys_s)))

def exercise_single_conformer_only_water():
  """
  tst_31.model_with_b_only_water (HOH 100 labelled B only, overlapping
  TYR 94 A; no bonded partner) clustered one residue per cluster and capped:
  the water exists in conformer B only. Its own fragment has one job, B,
  weight 1 (no A variant: it would hold no cluster atom); a fragment with
  the water in its buffer has variant A without the water atoms and
  variant B with them, weights 1 and 1. The capper refuses nothing.
  """
  from qrefine.tests.unit import tst_31
  model = boxed(tst_31.model_with_b_only_water().get_hierarchy())
  f = make_fragments(model)
  ph = f.pdb_hierarchy
  n = ph.atoms_size()
  assert f.pdb_hierarchy_super.atoms_size() == n # no images
  water = flex.bool(n, False)
  water.set_selected(flex.size_t([a.i_seq for a in ph.atoms()
    if a.parent().parent().resseq_as_int() == 100]), True)
  assert water.count(True) == 3
  n_own, n_buffer = 0, 0
  for i in range(f.n_fragments):
    jobs = jobs_of(f, i)
    full = f.fragment_full_super_selections[i]
    cluster = f.fragment_cluster_selections[i]
    got = sorted([(f.variant_labels[j], f.variant_weights[j],
      (f.fragment_selections[j] & water).count(True)) for j in jobs])
    for j in jobs:
      assert f.fragment_cappers[j] is not None
      assert f.fragment_charges[j] is not None
      h = f.variant_hierarchies[j]
      assert list(h.altloc_indices().keys()) == [""]
      assert h.atoms_size() == f.fragment_super_selections[j].count(True)
    if((cluster & water).count(True) > 0):
      n_own += 1
      assert (cluster & ~water).count(True) == 0 # maxnum 1
      assert got == [("B", 1, 3)], got
    elif((full & water).count(True) > 0):
      n_buffer += 1
      assert got == [("A", 1, 0), ("B", 1, 3)], (i, got)
    else:
      assert sum([w for l, w, nw in got]) == 2, got
  assert n_own == 1 and n_buffer > 0, (n_own, n_buffer)

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  """
  Conformer variants of the fragments (fragment.fragments with
  conformer_selections): variant atom sets, weights and labels, capper
  parents, single-conformer models unchanged, job folder names, single-label
  residue groups, interactions of altloc residues found by the single-
  conformer views.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_three_labels()
  exercise_variant_atoms_match_split()
  exercise_capper_parent()
  exercise_single_conformer_unchanged()
  exercise_copy_in_buffer()
  exercise_hetero()
  exercise_micro()
  exercise_job_folder_name()
  exercise_charge_embedding()
  exercise_sorry_prefix()
  exercise_interactions_match_single_conformer()
  exercise_single_conformer_only_water()

if(__name__ == "__main__"):
  run()
  print("OK")
