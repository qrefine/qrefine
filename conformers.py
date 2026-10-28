"""
Conformer models of a model with alternative conformations (altlocs).

Conformer labels are the altloc labels that have at least one non-hydrogen
atom; hydrogen-only labels (neutron H/D exchange sites) are not conformers.
Conformer X is the set of atoms with blank altloc, with label X, or with a
hydrogen-only label, in model order. split() returns each conformer as a
genuine single-conformer model: every atom group relabeled to blank and
merged into one atom group per residue group, atom order preserved, so no
downstream code (fragments, capper, charges, QM engines) sees an altloc.

Every residue group that has conformer labels must have all of the model's
conformer labels, otherwise conformer models would lack atoms (Sorry).

The multiplicity of an atom is the number of conformers it belongs to: 1 for
atoms with a conformer label, the number of conformers (at least 1) for all
others.
Occupancies are not used.

from_conformers combines the gradients of one restraints manager per
conformer: g_i = (1/m_i) * sum over conformers c containing atom i of g_c,i.
"""
from __future__ import absolute_import, division, print_function
import os
import mmtbx.model
from libtbx import group_args
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex

def _is_hydrogen(atom):
  return atom.element.strip().upper() in ("H", "D")

def conformer_labels(pdb_hierarchy):
  """
  Sorted altloc labels that have at least one non-hydrogen atom.
  """
  labels = set()
  for ag in pdb_hierarchy.atom_groups():
    if(ag.altloc in labels or ag.altloc.strip() == ""): continue
    for a in ag.atoms():
      if(not _is_hydrogen(a)):
        labels.add(ag.altloc)
        break
  return sorted(labels)

def n_conformers(pdb_hierarchy):
  return len(conformer_labels(pdb_hierarchy))

def multiplicity(pdb_hierarchy):
  """
  Per atom number of conformers it belongs to: 1 for atoms with a conformer
  label, max(n_conformers, 1) for blank and hydrogen-only-label atoms. Never
  0: a model without conformer labels gives 1 for every atom.
  """
  labels = set(conformer_labels(pdb_hierarchy))
  n = float(max(len(labels), 1))
  return flex.double(
    [1.0 if a.parent().altloc in labels else n for a in pdb_hierarchy.atoms()])

def _check_complete(pdb_hierarchy, all_labels):
  all_labels = set(all_labels)
  for rg in pdb_hierarchy.residue_groups():
    labels = set([ag.altloc for ag in rg.atom_groups()]) & all_labels
    if(len(labels) > 0 and labels != all_labels):
      raise Sorry("Residue %s has altlocs %s but model has %s: conformers "
        "would be incomplete." % (rg.id_str(), sorted(labels),
        sorted(all_labels)))

def _single_conformer_hierarchy(pdb_hierarchy, iselection, label):
  """
  pdb_hierarchy.select(iselection) with all atom groups of a residue group
  relabeled to blank and merged into the first one (atom order preserved).
  The merged atom group takes the resname of the atom group with this label,
  if any (microheterogeneity).
  """
  ph = pdb_hierarchy.deep_copy()
  for i, a in enumerate(ph.atoms()):
    a.tmp = i
  # Without copy_atoms the selection shares (and keeps the tmp of) the atoms
  # of ph, so atom order can be checked; deep_copy below detaches it from ph.
  h = ph.select(iselection)
  for rg in h.residue_groups():
    ags = rg.atom_groups()
    resnames = [ag.resname for ag in ags if ag.altloc == label]
    resname = resnames[0] if len(resnames) > 0 else ags[0].resname
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

def split(model):
  """
  One group_args(label, iselection, model) per conformer label, in
  conformer_labels order. iselection (ascending) selects blank atoms, atoms
  with this label and atoms with hydrogen-only labels; model is the
  unprocessed single-conformer mmtbx.model.manager of these atoms, in
  iselection order.
  Microheterogeneity (atom groups of a residue group with different resnames)
  is supported: in conformer X's model the residue takes the resname of the
  atom group with label X; blank and hydrogen-only-label atoms adopt it.
  """
  ph = model.get_hierarchy()
  labels = conformer_labels(ph)
  _check_complete(ph, labels)
  altlocs = [a.parent().altloc for a in ph.atoms()]
  result = []
  for label in labels:
    iselection = flex.bool(
      [not (al in labels and al != label) for al in altlocs]).iselection()
    h = _single_conformer_hierarchy(ph, iselection, label)
    m = mmtbx.model.manager(
      model_input       = None,
      pdb_hierarchy     = h,
      crystal_symmetry  = model.crystal_symmetry(),
      restraint_objects = model.get_restraint_objects(),
      log               = null_out())
    result.append(group_args(label = label, iselection = iselection, model = m))
  return result

class fragment_manager_proxy(object):
  """
  Fragment manager facade over the fragment managers of all conformers.
  """
  def __init__(self, fragment_managers):
    self.fragment_managers = fragment_managers

  def set_up_cluster_qm(self):
    for fm in self.fragment_managers:
      fm.set_up_cluster_qm()

  @property
  def interaction_list(self):
    result = []
    for fm in self.fragment_managers:
      result.extend(fm.interaction_list)
    return result

class from_conformers(object):
  """
  Restraints manager for a model with conformer labels: one manager per
  conformer, built by build(conformer_model, working_folder) in
  os.path.join(working_folder, label). Energy is undefined (0).
  """
  def __init__(self, model, build, working_folder="ase", log=None):
    if(log is None): log = null_out()
    ph = model.get_hierarchy()
    self.size = ph.atoms_size()
    self.multiplicity = multiplicity(ph)
    self.conformers = split(model)
    fragment_managers = []
    for c in self.conformers:
      folder = os.path.join(working_folder, c.label)
      print("Conformer %s: %d atoms, folder %s" % (
        c.label, c.iselection.size(), folder), file=log)
      c.restraints_manager = build(c.model, folder)
      fm = getattr(c.restraints_manager, "fragment_manager", None)
      if(fm is not None): fragment_managers.append(fm)
    self.fragment_manager = None
    if(len(fragment_managers) > 0):
      self.fragment_manager = fragment_manager_proxy(fragment_managers)

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
