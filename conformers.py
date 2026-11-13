"""
Conformer models of a model with alternative conformations (altlocs).

Conformer labels are the altloc labels that have at least one non-hydrogen
atom or a hydrogen that is not an H/D exchange site; the labels of neutron
H/D exchange sites (an H and a D of the same name) are not conformers, a
methyl in two H orientations is.
The own labels of a residue group are the conformer labels present in it.
Conformer X takes, from every residue group: all its atoms if it has no own
labels; its blank atoms, X atoms and hydrogen-only-label atoms if X is one
of its own labels; otherwise (fallback) its blank atoms, hydrogen-only-label
atoms and the atoms of its parent label, the first of its sorted own labels.
Exception to the fallback (absent residue groups): a residue group with no
blank or hydrogen-only-label atoms, whose own labels are a proper subset of
the model's conformer labels and which is covalently bonded to no other
residue group (a partial-occupancy water or free ligand labelled B only)
exists only in the conformers of its own labels; every other conformer
takes none of its atoms. Bonded is by distance
(cctbx.crystal.distance_based_connectivity, as fragment.bonded_residue_pairs)
in one of its own conformers: its non-hydrogen atom with label X and a
non-hydrogen atom of another residue group that conformer X contains (by
the fallback rule, or by its own label for an absent one); atoms of
different conformers never meet, so a water labelled B overlapping an A
side chain or an A chloride is unbonded.
A peptide link to a chain neighbour under any labels (the check_links test)
counts too. Bonds to a symmetry mate are not seen (the test is on the model
alone). Computed once per model.
Candidates bonded to each other by atoms of the same label (a B-only
two-residue ligand or glycan chain) are one entity: bonds inside it do not
count, it is absent or keeps the fallback as a whole. Residue groups with
blank atoms or a bonded partner outside their entity (protein residues,
bonded ligands, glycans) keep the fallback.
Atoms are in model order. split() returns each conformer as a genuine
single-conformer model: every atom group relabeled to blank and merged into
one atom group per residue group, atom order preserved, so no downstream
code (fragments, capper, charges, QM engines) sees an altloc.
single_conformer_hierarchy() is that construction for any selection;
super_selections() tiles the conformer selections of a model onto the super
cell and super sphere of a cctbx.crystal.super_cell.manager.

The fallback could break a peptide bond: a residue group lacking X whose
parent link atom (C, or N) is labelled, next to a residue group with X whose
link atom (N, or C) under X is labelled too. split() then requires the two
link atoms to be bonded (iotbx.pdb.altloc_completion.bonded), otherwise
Sorry (check_all_links does the same check without splitting). Pairs of
residue groups without any bonded C-N pair (chain breaks, ligands) are not
checked.

The multiplicity of an atom is the number of conformer models it belongs to
(membership): blank atoms belong to all, an atom with a conformer label to
its own conformer plus every conformer its residue group takes it into by
the fallback (none for an absent residue group: the atoms of a water
labelled B only have multiplicity 1). Occupancies are not used.

from_conformers combines the gradients of one restraints manager per
conformer: g_i = (1/m_i) * sum over conformers c containing atom i of g_c,i.
It serves the plain and expansion routes; the clustering route
(refine.create_restraints_manager) clusters the altloc model itself, with
per-fragment conformer variants (fragment.fragments conformer_selections,
cluster_restraints.from_cluster).

check_completed() requires the model to be completed already
(iotbx.pdb.altloc_completion, done by mmtbx.reduce2 qmready=true).
"""
from __future__ import absolute_import, division, print_function
import os
import mmtbx.model
from iotbx.pdb import altloc_completion
from libtbx import group_args
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex

def _is_hydrogen(atom):
  return atom.element.strip().upper() in ("H", "D")

def _hd_names_by_label(residue_group):
  """
  {altloc: {name.strip()[1:]: set of elements}} of the H and D atoms of
  residue_group (the name rule of iotbx exchangeable_hd_selections).
  """
  result = {}
  for ag in residue_group.atom_groups():
    d = result.setdefault(ag.altloc, {})
    for a in ag.atoms():
      e = a.element.strip()
      if(e in ("H", "D")):
        d.setdefault(a.name.strip()[1:], set()).add(e)
  return result

def _is_exchange_site(atom, altloc, hd_names_by_label):
  """
  True if atom (H or D, in the atom group with altloc) has its other
  isotope of equal name.strip()[1:] in an atom group with another altloc.
  """
  other = "D" if atom.element.strip() == "H" else "H"
  key = atom.name.strip()[1:]
  for label, names in hd_names_by_label.items():
    if(label == altloc): continue
    if(other in names.get(key, ())): return True
  return False

def conformer_labels(pdb_hierarchy):
  """
  Sorted altloc labels that are conformers: labels with at least one
  non-hydrogen atom, or with an H or D atom that is not an H/D exchange
  site (_is_exchange_site). A methyl in two H orientations under A and B
  makes A and B conformers; the H and D of a neutron exchange site do not.
  """
  labels = set()
  for rg in pdb_hierarchy.residue_groups():
    hd_names = None
    for ag in rg.atom_groups():
      if(ag.altloc in labels or ag.altloc.strip() == ""): continue
      for a in ag.atoms():
        if(not _is_hydrogen(a)):
          labels.add(ag.altloc)
          break
        if(hd_names is None): hd_names = _hd_names_by_label(rg)
        if(not _is_exchange_site(a, ag.altloc, hd_names)):
          labels.add(ag.altloc)
          break
  return sorted(labels)

def n_conformers(pdb_hierarchy):
  return len(conformer_labels(pdb_hierarchy))

def _own_labels(residue_group, labels):
  """
  Sorted conformer labels (labels: set of the model's) of a residue group:
  those of its atom groups with a non-hydrogen atom. An H-only atom group on
  a blank parent (1pq7 CYS A 57: HG B only, the disulfide to CYS A 41 SG A
  takes HG A) is that conformer's H, not a conformer of the residue, as in
  iotbx.pdb.altloc_completion.
  """
  return sorted(set([ag.altloc for ag in residue_group.atom_groups()
    if any([not a.element_is_hydrogen() for a in ag.atoms()])]) & labels)

def _residue_key(residue_group):
  """
  (model id, chain id, resseq, icode) of a residue group: the residue.
  """
  chain = residue_group.parent()
  return (chain.parent().id, chain.id, residue_group.resseq,
    residue_group.icode)

def _membership(rgs, labels):
  """
  Per atom (model order): the set of conformer labels whose conformer
  contains it by the fallback rule (no absent residue groups); None for all
  (blank and hydrogen-only-label atoms, residue groups without own labels).
  rgs: (residue group, own labels) of every residue group, labels: sorted
  conformer labels.
  """
  result = []
  for rg, own in rgs:
    lacking = [l for l in labels if l not in own]
    for ag in rg.atom_groups():
      m = None
      if(ag.altloc in own):
        m = set([ag.altloc])
        if(ag.altloc == own[0]): m.update(lacking)
      result.extend([m]*ag.atoms_size())
  return result

def _candidates(rgs, labels):
  """
  Indices (of rgs: (residue group, own labels) of every residue group) of
  the residue groups without blank or hydrogen-only-label atoms whose own
  labels are a proper subset of labels (the sorted conformer labels).
  """
  label_set = set(labels)
  result = set()
  for i, (rg, own) in enumerate(rgs):
    if(len(own) == 0 or len(own) == len(labels)): continue
    if(all([ag.altloc in label_set for ag in rg.atom_groups()])):
      result.add(i)
  return result

def _absent_residue_groups(pdb_hierarchy, rgs, labels):
  """
  Indices of the absent residue groups (module docstring) of pdb_hierarchy;
  rgs: (residue group, own labels) of each of its residue groups, labels
  its sorted conformer labels. A candidate (_candidates) is bonded if
  peptide-linked to a chain neighbour under any labels (_peptide_linked,
  the check_links test), or if one of its non-hydrogen atoms, label X, is
  bonded by distance (cctbx.crystal.distance_based_connectivity, as
  fragment.bonded_residue_pairs) to a non-hydrogen atom of another residue
  group that conformer X contains: by the fallback rule (_membership), for
  an absent residue group its own label only. Candidates bonded to each
  other by atoms of the same label form one entity (connected components,
  a B-only two-residue ligand or glycan chain); bonds inside an entity do
  not count and the entity is absent or falls back as a whole. Absent starts
  as all unlinked entities and loses the bonded ones until none is left
  (memberships only grow, so this ends): a water labelled B next to a
  chloride labelled A (1PWM HOH 1322, CL 1321) or overlapping an A side
  chain is unbonded.
  """
  candidates = _candidates(rgs, labels)
  if(len(candidates) == 0): return candidates
  from cctbx.crystal import distance_based_connectivity
  atoms = pdb_hierarchy.atoms()
  residue_of = flex.size_t()
  for i, (rg, own) in enumerate(rgs):
    residue_of.extend(flex.size_t(rg.atoms_size(), i))
  assert residue_of.size() == atoms.size()
  absent = set(candidates)
  index = dict([(rg.memory_id(), i) for i, (rg, own) in enumerate(rgs)])
  for chain in pdb_hierarchy.chains():
    chain_rgs = chain.residue_groups()
    for rg_1, rg_2 in zip(chain_rgs[:-1], chain_rgs[1:]):
      i_1, i_2 = index[rg_1.memory_id()], index[rg_2.memory_id()]
      if(not (i_1 in absent or i_2 in absent)): continue
      if(_peptide_linked(rg_1, rg_2)):
        absent.difference_update([i_1, i_2])
  bond_sets = distance_based_connectivity.build_simple_two_way_bond_sets(
    sites_cart = atoms.extract_xyz(),
    elements   = atoms.extract_element())
  # (residue, label of its atom, partner atom) of every bond of a candidate
  # atom to another residue group; hydrogens excluded (the distance test
  # takes a Cl...H-N hydrogen bond of 2.1 A for a bond, 1PWM CL 1321)
  bonds = []
  for i_seq in range(atoms.size()):
    r = residue_of[i_seq]
    if(r not in candidates or _is_hydrogen(atoms[i_seq])): continue
    label = atoms[i_seq].parent().altloc
    for j_seq in bond_sets[i_seq]:
      if(residue_of[j_seq] != r and not _is_hydrogen(atoms[j_seq])):
        bonds.append((r, label, j_seq))
  # entities: connected components of candidates over same-label bonds
  entity = dict([(i, i) for i in candidates])
  def root(i):
    while entity[i] != i: i = entity[i]
    return i
  for r, label, j_seq in bonds:
    r_j = residue_of[j_seq]
    if(r_j in candidates and atoms[j_seq].parent().altloc == label):
      entity[root(r)] = root(r_j)
  members = {}
  for i in candidates: members.setdefault(root(i), set()).add(i)
  def whole(residues):
    result = set()
    for i in residues: result.update(members[root(i)])
    return result
  absent -= whole(candidates - absent)
  membership = _membership(rgs, labels)
  while True:
    bonded = set()
    for r, label, j_seq in bonds:
      if(r in bonded or r not in absent): continue
      r_j = residue_of[j_seq]
      if(r_j in candidates and root(r_j) == root(r)): continue
      if(r_j in absent):
        meets = atoms[j_seq].parent().altloc == label
      else:
        m = membership[j_seq]
        meets = m is None or label in m
      if(meets): bonded.add(r)
    if(len(bonded) == 0): return absent
    absent -= whole(bonded)

def _selections(pdb_hierarchy):
  """
  Per conformer label: (group_args(label, iselection, n_fallback, n_absent),
  parents) with parents mapping residue group memory_id to the parent label
  of every residue group lacking this label and taking the fallback;
  n_absent is the number of absent residue groups lacking this label. The
  residue groups of one residue share their labels: a residue with two
  residue names (1ejg PRO A / SER B, C 22) is read as two residue groups,
  and the second must not fall back to its own first label (SER B next to
  PRO A in conformer A).
  """
  labels = conformer_labels(pdb_hierarchy)
  label_set = set(labels)
  own_of_residue = {}
  for rg in pdb_hierarchy.residue_groups():
    own_of_residue.setdefault(_residue_key(rg), set()).update(
      _own_labels(rg, label_set))
  rgs = [(rg, sorted(own_of_residue[_residue_key(rg)]))
    for rg in pdb_hierarchy.residue_groups()]
  absent = set()
  if(len(labels) > 1): absent = _absent_residue_groups(pdb_hierarchy, rgs,
    labels)
  result = []
  for label in labels:
    parents = {}
    n_absent = 0
    keep = flex.bool()
    for i, (rg, own) in enumerate(rgs):
      if(i in absent and label not in own):
        n_absent += 1
        keep.extend(flex.bool(rg.atoms_size(), False))
        continue
      take = label
      if(len(own) > 0 and label not in own):
        take = own[0]
        parents[rg.memory_id()] = take
      for ag in rg.atom_groups():
        keep.extend(flex.bool(ag.atoms_size(),
          ag.altloc not in label_set or ag.altloc == take))
    assert keep.size() == pdb_hierarchy.atoms_size()
    result.append((group_args(
      label      = label,
      iselection = keep.iselection(),
      n_fallback = len(parents),
      n_absent   = n_absent), parents))
  return result

def conformer_selections(pdb_hierarchy):
  """
  One group_args(label, iselection, n_fallback, n_absent) per conformer
  label, in conformer_labels order: iselection (ascending) selects the atoms
  of conformer label (module docstring), n_fallback is the number of residue
  groups lacking label that take their parent label, n_absent the number of
  absent residue groups lacking label (none of their atoms selected).
  """
  return [s for s, parents in _selections(pdb_hierarchy)]

def multiplicity(pdb_hierarchy):
  """
  Per atom number of conformer models it belongs to (conformer_selections).
  Never 0: a model without conformer labels gives 1 for every atom.
  """
  selections = conformer_selections(pdb_hierarchy)
  if(len(selections) == 0):
    return flex.double(pdb_hierarchy.atoms_size(), 1.0)
  result = flex.double(pdb_hierarchy.atoms_size(), 0.0)
  for s in selections:
    result.set_selected(s.iselection, result.select(s.iselection) + 1.0)
  return result

def _atom(residue_group, altloc, name):
  for ag in residue_group.atom_groups():
    if(ag.altloc != altloc): continue
    for a in ag.atoms():
      if(a.name.strip() == name): return a
  return None

def _peptide_linked(rg_1, rg_2):
  """
  True if some atom C of rg_1 (any label) is bonded to some N of rg_2.
  """
  cs = [a for a in rg_1.atoms() if a.name.strip() == "C"]
  ns = [a for a in rg_2.atoms() if a.name.strip() == "N"]
  for c in cs:
    for n in ns:
      if(altloc_completion.bonded(c, n)): return True
  return False

def _residue_id(atom):
  """Residue id with a blank altloc position, e.g. ' TYR A  94 '"""
  return " " + atom.parent().id_str()[1:]

def check_links(pdb_hierarchy, label, parents):
  """
  Sorry if conformer label breaks a peptide bond by the fallback. parents
  maps residue group memory_id to the parent label of every residue group
  lacking label. For consecutive residue groups of a chain, one lacking
  label (parent P) and the other having it, the link atom of the lacking one
  under P (C if it comes first, N otherwise) and the link atom of the other
  under label must be bonded when both exist. Not checked when no C-N pair
  of the two residue groups is bonded (chain break, ligand).
  """
  for chain in pdb_hierarchy.chains():
    rgs = chain.residue_groups()
    for rg_1, rg_2 in zip(rgs[:-1], rgs[1:]):
      lacks_1 = rg_1.memory_id() in parents
      lacks_2 = rg_2.memory_id() in parents
      if(lacks_1 == lacks_2): continue
      if(lacks_1): lacking, other, name_l, name_o = rg_1, rg_2, "C", "N"
      else:        lacking, other, name_l, name_o = rg_2, rg_1, "N", "C"
      parent = parents[lacking.memory_id()]
      a = _atom(lacking, parent, name_l)
      b = _atom(other, label, name_o)
      if(a is None or b is None): continue
      if(altloc_completion.bonded(a, b)): continue
      if(not _peptide_linked(rg_1, rg_2)): continue
      raise Sorry("conformer %s: residue %s takes its first label %s, but "
        "its %s is not bonded to %s of residue %s (%.2f A)" % (label,
        _residue_id(a), parent, name_l, name_o, _residue_id(b),
        a.distance(b)))

def _check_all_links(pdb_hierarchy, selections):
  for s, parents in selections:
    check_links(pdb_hierarchy=pdb_hierarchy, label=s.label, parents=parents)

def check_all_links(pdb_hierarchy):
  """
  check_links for every conformer label (the check split runs first).
  """
  _check_all_links(pdb_hierarchy, _selections(pdb_hierarchy))

def check_completed(pdb_hierarchy):
  """
  Sorry if iotbx.pdb.altloc_completion would share or copy atoms, i.e. the
  model was not completed by mmtbx.reduce2 qmready=true.
  """
  r = altloc_completion.run(pdb_hierarchy, check_only=True)
  if(len(r.shared) + len(r.copied) == 0): return
  atoms = ["%s %s" % (e[0], e[1]) for e in list(r.shared) + list(r.copied)]
  listing = ", ".join(atoms[:5])
  if(len(atoms) > 5): listing += ", ..."
  raise Sorry("Altloc completion needed: %d atoms to share and %d to copy "
    "(%s): run mmtbx.reduce2 qmready=true on the model first" % (
    len(r.shared), len(r.copied), listing))

def single_conformer_hierarchy(pdb_hierarchy, iselection, label, labels):
  """
  pdb_hierarchy.select(iselection) with all atom groups of a residue group
  relabeled to blank and merged into the first one (atom order preserved).
  The merged atom group takes the resname of the atom group with this label,
  if any, else that of the atom group with another conformer label (the
  parent label, the only one selected), else the first (microheterogeneity).
  """
  ph = pdb_hierarchy.deep_copy()
  for i, a in enumerate(ph.atoms()):
    a.tmp = i
  # Without copy_atoms the selection shares (and keeps the tmp of) the atoms
  # of ph, so atom order can be checked; deep_copy below detaches it from ph.
  h = ph.select(iselection)
  for rg in h.residue_groups():
    ags = rg.atom_groups()
    resnames = [ag.resname for ag in ags if ag.altloc == label] + \
      [ag.resname for ag in ags if ag.altloc in labels] + [ags[0].resname]
    resname = resnames[0]
    for ag in ags:
      ag.altloc = ""
      # merge_atom_groups asserts equal resname
      ag.resname = resname
    for ag in ags[1:]:
      rg.merge_atom_groups(ags[0], ag)
      rg.remove_atom_group(ag)
  assert list(h.atoms().extract_tmp_as_size_t()) == list(iselection)
  h = h.deep_copy()
  h.atoms().reset_i_seq()
  h.atoms_reset_serial()
  return h

def super_selections(conformer_selections, n_atoms, scm):
  """
  {label: group_args(cell, sphere)} for every conformer selection (as returned
  by conformer_selections of a hierarchy of n_atoms atoms): the selection of
  that conformer as flex.bool over scm.super_cell_hierarchy and over
  scm.super_sphere_hierarchy, scm being a cctbx.crystal.super_cell.manager of
  that hierarchy. The super cell is the master model followed by chain-wise
  copies of the whole model in model order, so atom k*n_atoms+i of it is
  master atom i and the master selection is tiled; the sphere selection is
  the cell one restricted to scm.keep_selection.
  """
  n_cell = scm.super_cell_hierarchy.atoms_size()
  assert n_cell % n_atoms == 0
  result = {}
  for s in conformer_selections:
    master = flex.bool(n_atoms, False)
    master.set_selected(s.iselection, True)
    cell = flex.bool()
    for i in range(n_cell // n_atoms): cell.extend(master)
    result[s.label] = group_args(cell = cell,
      sphere = cell.select(scm.keep_selection))
  return result

def split(model):
  """
  One group_args(label, iselection, model, n_fallback, n_absent) per
  conformer label, in conformer_labels order (the model of a conformer has
  no residue group absent from it). iselection is that of conformer_selections;
  model is the unprocessed single-conformer mmtbx.model.manager of these
  atoms, in iselection order. Runs check_links for every label first.
  Microheterogeneity (atom groups of a residue group with different resnames)
  is supported: in conformer X's model the residue takes the resname of the
  atom group with label X (or of its parent label); blank and
  hydrogen-only-label atoms adopt it.
  """
  ph = model.get_hierarchy()
  labels = set(conformer_labels(ph))
  selections = _selections(ph)
  _check_all_links(ph, selections)
  result = []
  for s, parents in selections:
    h = single_conformer_hierarchy(ph, s.iselection, s.label, labels)
    m = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = h,
      crystal_symmetry  = model.crystal_symmetry(),
      restraint_objects = model.get_restraint_objects(),
      log               = null_out())
    result.append(group_args(label = s.label, iselection = s.iselection,
      model = m, n_fallback = s.n_fallback, n_absent = s.n_absent))
  return result

class from_conformers(object):
  """
  Restraints manager for a model with conformer labels: one manager per
  conformer, built by build(conformer_model, working_folder) in
  os.path.join(working_folder, label) (None when working_folder is None:
  no files). Energy is undefined (0). The plain and expansion routes only
  (docs/superpowers/specs/2026-11-07-altloc-clustering-variants-design.md,
  section 3.1); clustering uses per-fragment variants instead.
  """
  def __init__(self, model, build, working_folder="ase", log=None):
    if(log is None): log = null_out()
    ph = model.get_hierarchy()
    self.size = ph.atoms_size()
    self.multiplicity = multiplicity(ph)
    self.conformers = split(model)
    for c in self.conformers:
      folder = None
      if(working_folder is not None):
        folder = os.path.join(working_folder, c.label)
      print("Conformer %s: %d atoms, folder %s" % (
        c.label, c.iselection.size(), folder), file=log)
      if(c.n_fallback > 0):
        print("  %d residue groups without label %s take their first own "
          "label" % (c.n_fallback, c.label), file=log)
      if(c.n_absent > 0):
        print("  %d residue groups absent from conformer %s (single-label, "
          "unbonded)" % (c.n_absent, c.label), file=log)
      c.restraints_manager = build(c.model, folder)

  def target_and_gradients(self, sites_cart):
    assert sites_cart.size() == self.size
    g = flex.vec3_double(self.size, (0,0,0))
    for c in self.conformers:
      isel = c.iselection
      g_c = c.restraints_manager.target_and_gradients(
        sites_cart = sites_cart.select(isel))[1]
      assert g_c.size() == isel.size()
      g = g.set_selected(isel, g.select(isel) + g_c)
    return 0, g * (1/self.multiplicity)

  def energies_sites(self, sites_cart, compute_gradients=True):
    tg = self.target_and_gradients(sites_cart=sites_cart)
    return group_args(
      target    = tg[0],
      gradients = tg[1])
