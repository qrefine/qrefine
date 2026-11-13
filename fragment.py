from __future__ import print_function
from __future__ import absolute_import
import os
from libtbx.utils import Sorry
from scitbx.array_family import flex
from cctbx import uctbx
from cctbx import crystal
from .utils import fragment_utils
from libtbx import group_args
from qrefine.super_cell import expand
from qrefine.utils.hierarchy_utils import write_hierarchy
from . import fragment_capping
from . import conformers
from .charges import charges_class
from mmtbx.pair_interaction import pair_interaction

def check_residue_integrity(pdb_hierarchy, selection):
  """
  Assert that selection (flex.bool over the atoms of pdb_hierarchy) takes
  every residue group of pdb_hierarchy either whole or not at all: clusters
  and fragments are made of whole residues.
  """
  assert selection.size() == pdb_hierarchy.atoms_size()
  offset = 0
  for rg in pdb_hierarchy.residue_groups():
    n = rg.atoms_size()
    n_selected = selection[offset:offset+n].count(True)
    assert n_selected in (0, n), \
      'error in cluster: residue %s has %d of %d atoms selected' % (
        rg.id_str(), n_selected, n)
    offset += n

def bonded_residue_pairs(pdb_hierarchy):
  """
  Pairs (i, j), i < j, of 1-based residue group indices of pdb_hierarchy
  with a covalent bond between the two residues, found by distance
  (cctbx.crystal.distance_based_connectivity: any residue type).
  """
  from cctbx.crystal import distance_based_connectivity
  atoms = pdb_hierarchy.atoms()
  bond_sets = distance_based_connectivity.build_simple_two_way_bond_sets(
    sites_cart = atoms.extract_xyz(),
    elements   = atoms.extract_element())
  residue_of = []
  for i, rg in enumerate(pdb_hierarchy.residue_groups()):
    residue_of.extend([i+1]*rg.atoms_size())
  result = set()
  for i, bonded in enumerate(bond_sets):
    for j in bonded:
      if(residue_of[i] != residue_of[j]):
        result.add((min(residue_of[i], residue_of[j]),
                    max(residue_of[i], residue_of[j])))
  return sorted(result)

def _view_residues(pdb_hierarchy, iselection):
  """
  1-based indices of the residue groups of pdb_hierarchy with an atom in
  iselection, in order: residue k+1 of pdb_hierarchy.select(iselection) is
  residue result[k] of pdb_hierarchy.
  """
  selected = flex.bool(pdb_hierarchy.atoms_size(), False)
  selected.set_selected(iselection, True)
  result = []
  offset = 0
  for i, rg in enumerate(pdb_hierarchy.residue_groups()):
    n = rg.atoms_size()
    if(selected[offset:offset+n].count(True) > 0): result.append(i+1)
    offset += n
  return result

def _map_pairs(pairs, residues):
  """
  Residue index pairs of a view mapped to the model (residues:
  _view_residues of the view; None: unchanged).
  """
  if(residues is None): return pairs
  return [(residues[i-1], residues[j-1]) for i, j in pairs]

def residues_within_radius(pdb_hierarchy, selection, radius):
  """
  Return positional indices (into pdb_hierarchy.atoms()) of all atoms of the
  residue groups having at least one atom within radius of the selection.
  Non-periodic: the hierarchy is placed in a P1 box with a buffer layer
  larger than radius so that no images are picked up.
  """
  atoms = pdb_hierarchy.atoms()
  xyz = atoms.extract_xyz()
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart = xyz, buffer_layer = radius+1.0)
  cs = crystal.symmetry(unit_cell=box.unit_cell, space_group_symbol="P1")
  tmp = pdb_hierarchy.deep_copy()
  tmp.atoms().set_xyz(box.sites_cart)
  xrs = tmp.extract_xray_structure(crystal_symmetry=cs)
  within = xrs.selection_within(radius=radius, selection=selection)
  result = flex.size_t()
  offset = 0
  for rg in pdb_hierarchy.residue_groups():
    n = rg.atoms().size()
    if(within[offset:offset+n].count(True)>0):
      result.extend(flex.size_t(range(offset, offset+n)))
    offset += n
  return result

class fragments(object):

  def __init__(self,
      working_folder             = "ase",
      make_working_folder        = True, # A way to not make it if not needed
      clustering_method          = None,
      maxnum_residues_in_cluster = 20,
      bcc_threshold              = 9,
      charge_embedding           = False,
      two_buffers                = False,
      distance_buffer            = False,
      distance_buffer_radius     = 5.0,
      pdb_hierarchy              = None,
      qm_engine_name             = None,
      crystal_symmetry           = None,
      cif_objects                = None,
      debug                      = False,
      charge_cutoff              = 8.0,
      save_clusters              = False,
      select_within_radius       = 10.0,
      clusters_only              = False,
      cap_fragments              = True,
      conformer_selections       = None):
    # cap_fragments: cap the fragments and compute their charges (needed
    # for QM restraints only; cctbx restraints use the fragment selections)
    self.cap_fragments = cap_fragments
    # conformer_selections: conformers.conformer_selections(pdb_hierarchy);
    # None or empty: single-conformer behaviour (one variant per fragment)
    self.conformer_selections = conformer_selections or None
    self.labels = []
    if(self.conformer_selections is not None):
      self.labels = [s.label for s in self.conformer_selections]
    self.select_within_radius = select_within_radius
    self.charge_embedding = charge_embedding
    self.two_buffers = two_buffers
    self.distance_buffer = distance_buffer
    self.distance_buffer_radius = distance_buffer_radius
    self.crystal_symmetry = crystal_symmetry
    self.cif_objects = cif_objects
    # working_folder None: no file is written (expansion, fragments)
    self.working_folder = None
    if(working_folder is not None):
      self.working_folder = os.path.abspath(working_folder)
    self.pdb_hierarchy = pdb_hierarchy
    self.charge_cutoff = charge_cutoff
    self.qm_engine_name = qm_engine_name
    self.clustering_method = clustering_method
    self.debug = debug
    self.maxnum_residues_in_cluster =  maxnum_residues_in_cluster
    self.bcc_threshold = bcc_threshold
    self.save_clusters = save_clusters
    self.expansion = None
    self.expansion_file = None
    self.pdb_hierarchy_super = None
    self.clusters = None
    self.clusters_only = clusters_only
    self.charge_service = None
    if(self.charge_embedding):
      if(self.qm_engine_name in ["mopac"]):
        raise Sorry(
          "cluster.charge_embedding=True is not supported by the %s engine: "
          "it ignores external point charges." % self.qm_engine_name)
      # Point charges for charge embedding. Only the monomer library server
      # and ion charge table are reused; the hierarchy is replaced per
      # fragment in write_mm_charge_file.
      self.charge_service = charges_class(
        pdb_hierarchy    = pdb_hierarchy,
        crystal_symmetry = crystal_symmetry,
        cif_objects      = cif_objects)
    #
    if(make_working_folder and self.working_folder is not None):
      os.makedirs(self.working_folder, exist_ok=True)
    self._expand()
    self.set_up_cluster_qm()

  def _expand(self):
    self.expansion = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.select_within_radius)
    self.pdb_hierarchy_super = self.expansion.ph_super_sphere
    if(self.working_folder is None): return
    # .cif instead of .pdb when the super sphere does not fit PDB format
    self.expansion_file = self.expansion.write_super_cell_selected_in_sphere(
      file_name=os.path.join(self.working_folder, "expansion.pdb"))

  def update_xyz(self, sites_cart):
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    self.expansion = self.expansion.update_xyz(sites_cart=sites_cart)
    self.pdb_hierarchy_super = self.expansion.ph_super_sphere

  def set_up_cluster_qm(self):
    self.get_clusters()
    if not self.clusters_only:
      self.set_up_conformers()
      self.get_fragments()
      self._release_sphere_views()
      self.get_fragment_hierarchies_and_charges()

  def _model_views(self):
    """
    (hierarchy, residues) pairs the interaction search runs on: the model
    without conformer selections (residues None), else its single-conformer
    view per label, residues mapping its 1-based residue group indices to
    those of the model (_view_residues: the relabel keeps residue group
    order; residue groups absent from a conformer are missing from its
    view). pair_interaction on an altloc hierarchy assigns a grid point
    near an altloc residue to its own A and B atoms (same residue: no
    pair), so interactions of altloc residues would be lost.
    """
    if(self.conformer_selections is None): return [(self.pdb_hierarchy, None)]
    result = []
    for s in self.conformer_selections:
      h = conformers.single_conformer_hierarchy(self.pdb_hierarchy,
        s.iselection, s.label, self.labels)
      residues = _view_residues(self.pdb_hierarchy, s.iselection)
      assert len(list(h.residue_groups())) == len(residues)
      result.append((h, residues))
    return result

  def _release_sphere_views(self):
    """
    The single-conformer super-sphere views are setup-time objects (buffer
    search); they are kept after setup as conformer_spheres only with
    charge embedding.
    """
    if(self.conformer_selections is None): return
    if(not self.charge_embedding): self.conformer_spheres = {}

  def set_up_conformers(self):
    """
    Per conformer label: its selection over the super cell and the super
    sphere (self._super_selections) and the relabeled single-conformer
    super sphere (conformer_spheres: the buffer search of get_fragments;
    kept after setup for charge embedding only, _release_sphere_views).
    Built once per (re-)clustering, shared by all variants of a label.
    """
    SCM = self.expansion.SCM
    self._super_selections = {}
    self.conformer_spheres = {}
    if(self.conformer_selections is None):
      self.conformer_spheres = {"": self.pdb_hierarchy_super}
      return
    self._super_selections = conformers.super_selections(
      conformer_selections = self.conformer_selections,
      n_atoms              = self.pdb_hierarchy.atoms_size(),
      scm                  = SCM)
    for label in self.labels:
      s = self._super_selections[label]
      self.conformer_spheres[label] = \
        conformers.single_conformer_hierarchy(self.pdb_hierarchy_super,
          s.sphere.iselection(), label, self.labels)

  def _set_up_conformer_cells(self):
    """
    Per conformer label (capping with conformers only): the relabeled
    single-conformer super cell, capper parent of the variants of that
    label (conformer_cells), and the rank array mapping super cell
    positions to its atoms (_cell_ranks). Each is a copy of most of the
    super cell: built just before the cappers and released after them
    (_release_conformer_cells).
    """
    self.conformer_cells = {}
    self._cell_ranks = {}
    if(self.conformer_selections is None or not self.cap_fragments): return
    SCM = self.expansion.SCM
    for label in self.labels:
      s = self._super_selections[label]
      isel = s.cell.iselection()
      self.conformer_cells[label] = conformers.single_conformer_hierarchy(
        SCM.super_cell_hierarchy, isel, label, self.labels)
      rank = flex.size_t(s.cell.size(), 0)
      rank.set_selected(isel, flex.size_t(range(isel.size())))
      self._cell_ranks[label] = rank

  def _release_conformer_cells(self):
    self.conformer_cells = {}
    self._cell_ranks = {}

  def _variants(self, fragment_super_selection, cluster_selection=None):
    """
    Conformer variants of a fragment (fragment_super_selection: flex.bool
    over the super sphere): group_args(selection, weight, label, parent) per
    distinct atom set fragment & conformer, in label order. weight: number
    of conformers giving that atom set; parent: the smallest of them (the
    label the variant is relabeled and capped with); label: parent, or ""
    when the fragment has one variant from all conformers. A conformer
    without any atom of cluster_selection (flex.bool over the master atoms;
    its cluster residues are absent from that conformer, conformers module
    docstring) gives no variant: it would contribute no gradient. Without
    conformers: the fragment itself, weight 1, label and parent "".
    """
    if(self.conformer_selections is None):
      return [group_args(selection = fragment_super_selection, weight = 1,
        label = "", parent = "")]
    result = []
    by_atoms = {}
    n_master = self.pdb_hierarchy.atoms_size()
    for label in self.labels:
      v = fragment_super_selection & self._super_selections[label].sphere
      if(cluster_selection is not None and
         (v[:n_master] & cluster_selection).count(True) == 0):
        continue
      key = tuple(v.iselection())
      if(key in by_atoms):
        by_atoms[key].weight += 1
        continue
      by_atoms[key] = group_args(selection = v, weight = 1, label = label,
        parent = label)
      result.append(by_atoms[key])
    if(len(result) == 1 and result[0].weight == len(self.labels)):
      result[0].label = ""
    return result

  def _variant_hierarchy(self, selection, parent):
    """
    The variant selection (over the super sphere) as a single-conformer
    hierarchy at the current coordinates: plain selection without
    conformers, else relabeled with the parent label (only the variant's
    atoms are copied, not the super sphere).
    """
    if(self.conformer_selections is None):
      return self.pdb_hierarchy_super.select(selection)
    h = self.pdb_hierarchy_super.select(selection, copy_atoms=True)
    result = conformers.single_conformer_hierarchy(h,
      flex.size_t(range(h.atoms_size())), parent, self.labels)
    assert list(result.altloc_indices().keys()) == [""]
    return result

  def get_clusters(self):
    """
    Interacting and bonded residue pairs (of the model, or the union over
    its single-conformer views: _model_views), then the clusters.
    """
    n_residues=len(list(self.pdb_hierarchy.residue_groups()))
    if(self.conformer_selections is None):
      self.interaction_list = pair_interaction.run(self.pdb_hierarchy)
      bonded_pairs = bonded_residue_pairs(self.pdb_hierarchy)
    else:
      interactions = set()
      bonded = set()
      for h, residues in self._model_views():
        interactions.update(_map_pairs(pair_interaction.run(h), residues))
        bonded.update(_map_pairs(bonded_residue_pairs(h), residues))
      self.interaction_list = sorted(interactions)
      bonded_pairs = sorted(bonded)
    self.bonded_pairs = bonded_pairs
    from . import clustering
    clusters = clustering.betweenness_centrality_clustering(
      self.interaction_list,
      size = n_residues,
      maxnum_residues_in_cluster = self.maxnum_residues_in_cluster,
      bcc_threshold = self.bcc_threshold,
      bonded_pairs = bonded_pairs).get_clusters()
    self.clusters = sorted(clusters, key=len, reverse=True)

  def get_fragments(self):
    """
    Cluster atoms, fragment (cluster plus buffer) atoms and fragment
    residues of every cluster, from pair_interaction.run on the super sphere.
    Atoms are 1-based positions in pdb_hierarchy_super (the master copy is
    its first block, so the cluster atoms are positions in pdb_hierarchy too);
    residues are 1-based residue group indices.
    """
    self.pdb_hierarchy_super.atoms_reset_serial() # for the written PDB files
    ph = self.pdb_hierarchy_super
    clusters = self.clusters
    cluster_atoms_in_ph = []
    fragment_super_atoms_in_ph = []
    molecules_in_fragments = []
    # The super sphere (each of its single-conformer views) is the same for
    # every cluster: prepared once.
    prepared = self._prepared_sphere_views()
    for cluster in clusters:
      atoms_in_one_cluster, atoms_in_one_fragment, molecules_in_one_fragment = \
        self._run_prepared(prepared, cluster)
      cluster_atoms_in_ph.append(atoms_in_one_cluster)
      fragment_super_atoms_in_ph.append(atoms_in_one_fragment)
      molecules_in_fragments.append(molecules_in_one_fragment)
    if(self.two_buffers): # a second buffer layer around the first
      fragment_super_atoms_in_ph = []
      for molecules in molecules_in_fragments:
        atoms_in_one_fragment = self._run_prepared(prepared, molecules)[1]
        fragment_super_atoms_in_ph.append(atoms_in_one_fragment)
    if(self.distance_buffer):
      fragment_super_atoms_in_ph = self.add_distance_buffer(
        ph                         = ph,
        cluster_atoms_in_ph        = cluster_atoms_in_ph,
        fragment_super_atoms_in_ph = fragment_super_atoms_in_ph)
    #
    self.cluster_atoms = []
    self.fragment_super_atoms = []
    for i_cluster in range(len(clusters)):
      self.collect_cluster_and_fragment(cluster_atoms_in_ph,
                                    fragment_super_atoms_in_ph, i_cluster)

  def _prepared_sphere_views(self):
    """
    pair_interaction.prepare of the super sphere, or of each of its
    single-conformer views (conformer_spheres, set_up_conformers), with
    the view's iselection over the super sphere and its residues
    (_view_residues: 1-based residue group indices of the super sphere).
    """
    if(self.conformer_selections is None):
      return [(pair_interaction.prepare(self.pdb_hierarchy_super), None,
        None)]
    result = []
    for label in self.labels:
      h = self.conformer_spheres[label]
      isel = self._super_selections[label].sphere.iselection()
      residues = _view_residues(self.pdb_hierarchy_super, isel)
      assert len(list(h.residue_groups())) == len(residues)
      result.append((pair_interaction.prepare(h), isel, residues))
    return result

  def _run_prepared(self, prepared, residues):
    """
    pair_interaction.run_prepared(core=residues) on the super sphere: core
    atoms, fragment atoms (1-based positions in the super sphere) and
    fragment residues (1-based residue group indices). With views: the
    union over the views, core atoms mapped back to the super sphere,
    fragment atoms all atoms of the union of the fragment residues (whole
    residues of the super sphere; variants cut them per conformer). A view
    lacking every residue of the core (absent from that conformer) is
    skipped.
    """
    if(len(prepared) == 1 and prepared[0][1] is None):
      return pair_interaction.run_prepared(prepared[0][0], residues)
    core, molecules = set(), set()
    for p, isel, view_residues in prepared:
      in_view = dict([(r, k+1) for k, r in enumerate(view_residues)])
      core_v = [in_view[r] for r in residues if r in in_view]
      if(len(core_v) == 0): continue
      c, f, m = pair_interaction.run_prepared(p, core_v)
      core.update([isel[i-1]+1 for i in c])
      molecules.update([view_residues[i-1] for i in m])
    residue_atoms = self._sphere_residue_atoms()
    atoms = []
    for m in molecules:
      atoms.extend(residue_atoms[m-1])
    return sorted(core), sorted(atoms), sorted(molecules)

  def _sphere_residue_atoms(self):
    """1-based atom positions of every residue group of the super sphere."""
    result = []
    offset = 1
    for rg in self.pdb_hierarchy_super.residue_groups():
      n = rg.atoms_size()
      result.append(list(range(offset, offset+n)))
      offset += n
    return result

  def add_distance_buffer(self, ph, cluster_atoms_in_ph,
                          fragment_super_atoms_in_ph):
    """
    Extend each fragment with all residues of ph (the super-sphere hierarchy)
    having an atom within self.distance_buffer_radius of the cluster. Atoms
    in the input and output lists are 1-based positions in ph.
    """
    result = []
    for cluster, fragment in zip(cluster_atoms_in_ph, fragment_super_atoms_in_ph):
      selection = flex.bool(ph.atoms_size(), False)
      selection = selection.set_selected(
        flex.size_t([i-1 for i in cluster]), True)
      positions = residues_within_radius(
        pdb_hierarchy = ph,
        selection     = selection,
        radius        = self.distance_buffer_radius)
      extra = [i+1 for i in positions]
      result.append(sorted(set(fragment) | set(extra)))
    return result

  def collect_cluster_and_fragment(self, cluster_atoms_in_ph,
                                fragment_super_atoms_in_ph, i_cluster):
    self.cluster_atoms.append(cluster_atoms_in_ph[i_cluster])
    self.fragment_super_atoms.append(fragment_super_atoms_in_ph[i_cluster])

  def make_capper(self, fragment_super_selection, parent="",
                  variant_hierarchy=None):
    """
    The capper of the fragment (or variant) fragment_super_selection (a
    selection on pdb_hierarchy_super, the super sphere). Cut bonds are found
    against the super cell, not the sphere: a symmetry-copy residue at the
    edge of the sphere has its peptide partner outside the sphere, a cut
    bond (one H), not an unsatisfied valence. A variant of conformer parent
    is cut out of that conformer's relabeled super cell (conformer_cells),
    so atoms of other conformers are no bond partners; variant_hierarchy is
    its relabeled hierarchy.
    """
    SCM = self.expansion.SCM
    cell_positions = SCM.keep_selection.iselection().select(
      fragment_super_selection.iselection())
    parent_hierarchy = SCM.super_cell_hierarchy
    if(parent != ""):
      # the conformer cells exist only while
      # get_fragment_hierarchies_and_charges builds the cappers
      # (_set_up_conformer_cells / _release_conformer_cells)
      assert parent in getattr(self, "conformer_cells", {}), (parent,
        "conformer cells are built only during capper construction")
      assert self._super_selections[parent].cell.select(
        cell_positions).all_eq(True)
      parent_hierarchy = self.conformer_cells[parent]
      cell_positions = self._cell_ranks[parent].select(cell_positions)
    if(variant_hierarchy is None):
      variant_hierarchy = self.pdb_hierarchy_super.select(
        fragment_super_selection)
    selection = flex.bool(parent_hierarchy.atoms_size(), cell_positions)
    c = fragment_capping.capper(
      parent_hierarchy = parent_hierarchy,
      selection        = selection,
      crystal_symmetry = self.expansion.cs_box,
      cif_objects      = self.cif_objects)
    assert c.fragment_hierarchy.is_similar_hierarchy(variant_hierarchy)
    return c

  def get_fragment_hierarchies_and_charges(self):
    """
    Per fragment: its cluster selection and full super-sphere selection;
    per variant (job): the variant's selections, relabeled hierarchy,
    capper and charge, with the fragment's cluster and buffer selections
    repeated (module docstring of conformers for the conformer rule).
    """

    def pdb_hierarchy_select(atoms_size, selection):
      selection_array = flex.bool(atoms_size, False)
      for item in selection:
        if(item<=atoms_size):
          selection_array[item-1] = True
      return selection_array

    n_master = self.pdb_hierarchy.atoms_size()
    self.fragment_cluster_selections = []
    self.fragment_full_super_selections = []
    self.variant_fragment_index = []
    self.variant_weights = []
    self.variant_labels = []
    self.variant_hierarchies = []
    self.fragment_selections = []
    self.fragment_super_selections = []
    self.fragment_charges = []
    self.cluster_selections = []
    self.buffer_selections = []
    self.fragment_cappers = []
    self.n_fragments = len(self.fragment_super_atoms)
    self._set_up_conformer_cells()
    try:
      for i in range(self.n_fragments):
        fragment_selection = pdb_hierarchy_select(n_master,
          self.fragment_super_atoms[i])
        ## QM part is fragment_super
        fragment_super_selection = pdb_hierarchy_select(
          self.pdb_hierarchy_super.atoms_size(),
          self.fragment_super_atoms[i])
        cluster_selection = pdb_hierarchy_select(n_master,
          self.cluster_atoms[i])
        s = fragment_selection==cluster_selection
        buffer_selection = fragment_selection.deep_copy().set_selected(s, False)
        check_residue_integrity(self.pdb_hierarchy_super,
          fragment_super_selection)
        check_residue_integrity(self.pdb_hierarchy, cluster_selection)
        self.fragment_cluster_selections.append(cluster_selection)
        self.fragment_full_super_selections.append(fragment_super_selection)
        for v in self._variants(fragment_super_selection, cluster_selection):
          variant_hierarchy = self._variant_hierarchy(v.selection, v.parent)
          # caps (and the charge of the capped variant) as reduce2
          # capping=True places them; the capper re-places them at every
          # gradient evaluation
          if self.cap_fragments:
            try:
              c = self.make_capper(v.selection, parent=v.parent,
                variant_hierarchy=variant_hierarchy)
            except Sorry as e:
              raise Sorry("fragment %d, conformer %s: %s" % (i,
                v.label or "-", str(e))) from None
            self.fragment_cappers.append(c)
            fragment_charge = c.charge
          else:
            self.fragment_cappers.append(None)
            fragment_charge = None
          self.variant_fragment_index.append(i)
          self.variant_weights.append(v.weight)
          self.variant_labels.append(v.label)
          self.variant_hierarchies.append(variant_hierarchy)
          self.fragment_super_selections.append(v.selection)
          # master atoms of the variant: the master copy is the first block
          # of the super sphere
          self.fragment_selections.append(v.selection[:n_master])
          self.fragment_charges.append(fragment_charge)
          self.cluster_selections.append(cluster_selection)
          self.buffer_selections.append(buffer_selection)
    finally:
      self._release_conformer_cells()

  def get_fragment_extracts(self, with_super_sphere_grm=False):
    # with_super_sphere_grm: include the cctbx restraints manager of the
    # super sphere (built on first use; cctbx restraints of fragments only)
    super_sphere_grm = None
    if(with_super_sphere_grm):
      super_sphere_grm = self.expansion.super_sphere_geometry_restraints_manager
    return group_args(
      cluster_selections        = self.cluster_selections,
      fragment_charges          = self.fragment_charges,
      fragment_selections       = self.fragment_selections,
      fragment_super_selections = self.fragment_super_selections,
      fragment_cappers          = self.fragment_cappers,
      working_folder            = self.working_folder,
      fragment_super_atoms      = self.fragment_super_atoms,
      cluster_atoms             = self.cluster_atoms,
      qm_engine_name            = self.qm_engine_name,
      charge_embedding          = self.charge_embedding,
      crystal_symmetry          = self.crystal_symmetry,
      pdb_hierarchy             = self.pdb_hierarchy,
      pdb_hierarchy_super       = self.pdb_hierarchy_super,
      expansion_cs              = self.expansion.cs_box,
      buffer_selections         = self.buffer_selections,
      debug                     = self.debug,
      charge_service            = self.charge_service,
      charge_cutoff             = self.charge_cutoff,
      expansion_file            = self.expansion_file,
      save_clusters             = self.save_clusters,
      n_fragments               = self.n_fragments,
      fragment_cluster_selections    = self.fragment_cluster_selections,
      fragment_full_super_selections = self.fragment_full_super_selections,
      variant_fragment_index    = self.variant_fragment_index,
      variant_weights           = self.variant_weights,
      variant_labels            = self.variant_labels,
      variant_hierarchies       = self.variant_hierarchies,
      conformer_spheres         = self.conformer_spheres,
      conformer_sphere_selections = dict([(label, s.sphere)
        for label, s in self._super_selections.items()]),
      super_sphere_geometry_restraints_manager = super_sphere_grm)

def job_folder_name(fragment_extracts, index):
  """
  Folder (under the working folder) and file stem of job index: the
  fragment index for a fragment with one variant of full weight (all
  conformers, so always for a single-conformer model), else
  <fragment index>_<label>.
  """
  i = fragment_extracts.variant_fragment_index[index]
  label = fragment_extracts.variant_labels[index]
  if(label == ""): return str(i)
  return "%d_%s" % (i, label)

def variant_hierarchy(fragment_extracts, index):
  """
  The relabeled hierarchy of job index (a copy) at the current coordinates
  of the super sphere. Without conformers: the selection of the current
  super sphere, as before variants.
  """
  fe = fragment_extracts
  selection = fe.fragment_super_selections[index]
  if(len(fe.conformer_sphere_selections) == 0):
    return fe.pdb_hierarchy_super.select(selection)
  xyz = fe.pdb_hierarchy_super.atoms().extract_xyz().select(selection)
  result = fe.variant_hierarchies[index].deep_copy()
  assert result.atoms_size() == xyz.size()
  result.atoms().set_xyz(xyz)
  return result

def get_qm_file_name_and_pdb_hierarchy(fragment_extracts, index):
  fragment_hierarchy = variant_hierarchy(fragment_extracts, index)
  if(fragment_extracts.working_folder is None):
    # no files: the capped fragment only (file name None)
    r = fragment_extracts.fragment_cappers[index].cap(
      sites_cart = fragment_hierarchy.atoms().extract_xyz())
    return None, r.pdb_hierarchy, r.fragment_selection
  name = job_folder_name(fragment_extracts, index)
  sub_working_folder = fragment_extracts.working_folder + "/" + name + "/"
  if (not os.path.isdir(sub_working_folder)):
    os.mkdir(sub_working_folder)
  qm_pdb_file = sub_working_folder + name + ".pdb"
  complete_qm_pdb_file = sub_working_folder + name + "_capping.pdb"
  # The fragment (user information) and the capped fragment (user
  # information; its name is the QM engine's label) are written at every
  # call, as .cif when they do not fit PDB format; the QM engine gets the
  # capped hierarchy, not the file.
  write_hierarchy(fragment_hierarchy, qm_pdb_file,
    crystal_symmetry=fragment_extracts.expansion_cs)
  # re-capping because geometry of the fragment has changed.
  r = fragment_extracts.fragment_cappers[index].cap(
    sites_cart = fragment_hierarchy.atoms().extract_xyz())
  complete_qm_pdb_file = write_hierarchy(r.pdb_hierarchy,
    complete_qm_pdb_file, crystal_symmetry=fragment_extracts.expansion_cs)
  return os.path.abspath(complete_qm_pdb_file), r.pdb_hierarchy, \
    r.fragment_selection

def charge(fragment_extracts, index):
  return fragment_extracts.fragment_charges[index]

def charge_environment(fragment_extracts, index):
  """
  (sphere, selection): the hierarchy the point charges of job index come
  from and the job's selection over it. Without conformers the super
  sphere; else the relabeled super sphere of the variant's conformer (the
  smallest label for a single-variant fragment: all give its atoms), its
  coordinates set to the current ones in place. Known limitation: a merged
  job (label "") thus has conformer A's field even where the environment
  within charge_cutoff differs between conformers.
  """
  fe = fragment_extracts
  selection = fe.fragment_super_selections[index]
  if(len(fe.conformer_sphere_selections) == 0):
    return fe.pdb_hierarchy_super, selection
  label = fe.variant_labels[index]
  if(label == ""): label = sorted(fe.conformer_sphere_selections.keys())[0]
  sphere_selection = fe.conformer_sphere_selections[label]
  assert (selection & ~sphere_selection).count(True) == 0
  sphere = fe.conformer_spheres[label]
  isel = sphere_selection.iselection()
  sphere.atoms().set_xyz(fe.pdb_hierarchy_super.atoms().extract_xyz().select(
    isel))
  return sphere, selection.select(isel)

def write_mm_charge_file(fragment_extracts, index):
  file_name = None
  if (fragment_extracts.charge_embedding is True):
    if(fragment_extracts.working_folder is None):
      raise Sorry("cluster.charge_embedding needs a working folder (the "
        "point charges are passed to the QM engine in a file).")
    ph, fragment_selection = charge_environment(fragment_extracts, index)
    if(fragment_extracts.charge_cutoff is not None):
      if(fragment_extracts.debug):
        print("charge_cutoff: ",fragment_extracts.charge_cutoff)
      xrs_super = ph.extract_xray_structure()
      non_fragment_selection_super = xrs_super.selection_within(
        radius=fragment_extracts.charge_cutoff,
        selection=fragment_selection)
      non_fragment_selection_super = non_fragment_selection_super&~fragment_selection
      non_fragment_hierarchy_super = ph.select(non_fragment_selection_super)
    else:
      non_fragment_hierarchy_super = ph.select(~fragment_selection)
    non_fragment_hierarchy = non_fragment_hierarchy_super
    name = job_folder_name(fragment_extracts, index)
    sub_working_folder = fragment_extracts.working_folder + "/" + name + "/"
    if (not os.path.isdir(sub_working_folder)):
      os.mkdir(sub_working_folder)
    if(fragment_extracts.debug): print("write mm pdb file:", index)
    non_fragment_pdb_file = sub_working_folder + name + "_mm.pdb"
    write_hierarchy(non_fragment_hierarchy, non_fragment_pdb_file,
      crystal_symmetry=fragment_extracts.expansion_cs)
    charge_scaling_positions = fragment_utils.get_edge_atom_positions(
      ph, non_fragment_hierarchy, charge_embed=True)
    fragment_extracts.charge_service.update_pdb_hierarchy(
      non_fragment_hierarchy,
      fragment_extracts.expansion_cs,
    )
    if(fragment_extracts.qm_engine_name == "turbomole"):
      file_name = sub_working_folder + name + "_xyzq_cctbx.dat"
      fragment_extracts.charge_service.write_pdb_hierarchy_xyzq_file(
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    if(fragment_extracts.qm_engine_name in ["terachem","xtb",'mopac']):
      file_name = sub_working_folder + name + "_qxyz_cctbx.dat"
      fragment_extracts.charge_service.write_pdb_hierarchy_qxyz_file(
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    if(file_name is None):
      raise Sorry("There is no point charge file")
    file_name = os.path.abspath(file_name)
  return file_name

def write_cluster_and_fragments_pdbs(fragments,directory):
  # write current fragment and cluster PDBs into ./<directory>
  # makes a fresh(!) <directory>
  # fragments: fragment extracts; per fragment <i>_cluster.pdb and
  # <i>_frag.pdb, per variant <job_folder_name>_capped0.pdb
  F=fragments
  if not F.save_clusters:
    return
  from shutil import rmtree
  cwd = os.getcwd()
  frag_dir = os.path.join(cwd,directory)
  if os.path.exists(frag_dir):
    rmtree(frag_dir)
  os.mkdir(frag_dir)
  os.chdir(frag_dir)
  for index in range(F.n_fragments):
    index_cluster = F.pdb_hierarchy.select(
      F.fragment_cluster_selections[index])
    index_frag = F.pdb_hierarchy_super.select(
      F.fragment_full_super_selections[index])
    write_hierarchy(index_cluster, "%s_cluster.pdb" %(index),
      crystal_symmetry = F.expansion_cs)
    write_hierarchy(index_frag, "%s_frag.pdb" %(index),
      crystal_symmetry = F.expansion_cs)
  for index, capper in enumerate(F.fragment_cappers):
    if capper is None: continue # cctbx restraints: no caps
    h = variant_hierarchy(F, index)
    capped_hierarchy = capper.cap(
      sites_cart = h.atoms().extract_xyz()).pdb_hierarchy
    write_hierarchy(capped_hierarchy,
      "%s_capped0.pdb" % job_folder_name(F, index),
      crystal_symmetry=F.expansion_cs)

  log=open('fragment_info.txt','w')
  print('~  # clusters  : ',len(F.cluster_atoms), file=log)
  print('~  list of atoms per cluster:', file=log)
  print('~   ',[len(x) for x in F.cluster_atoms], file=log)
  print('~  list of atoms per fragment:', file=log)
  print('~   ',[len(x) for x in F.fragment_super_atoms], file=log)
  if(len(F.conformer_sphere_selections) > 0):
    print('~  variants per fragment (label: weight):', file=log)
    for index in range(F.n_fragments):
      jobs = [j for j, i in enumerate(F.variant_fragment_index) if i == index]
      print('~    %d: %s' % (index, ", ".join(["%s: %d" % (
        F.variant_labels[j] or "-", F.variant_weights[j]) for j in jobs])),
        file=log)
  log.close()
  os.chdir(cwd)
