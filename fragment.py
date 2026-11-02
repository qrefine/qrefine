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
from . import fragment_capping
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
      cap_fragments              = True):
    # cap_fragments: cap the fragments and compute their charges (needed
    # for QM restraints only; cctbx restraints use the fragment selections)
    self.cap_fragments = cap_fragments
    self.select_within_radius = select_within_radius
    self.charge_embedding = charge_embedding
    self.two_buffers = two_buffers
    self.distance_buffer = distance_buffer
    self.distance_buffer_radius = distance_buffer_radius
    self.crystal_symmetry = crystal_symmetry
    self.cif_objects = cif_objects
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
    if(make_working_folder):
      os.makedirs(self.working_folder, exist_ok=True)
    self._expand()
    self.set_up_cluster_qm()

  def _expand(self):
    self.expansion = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.select_within_radius)
    self.pdb_hierarchy_super = self.expansion.ph_super_sphere
    self.expansion_file = "expansion.pdb"
    self.expansion.write_super_cell_selected_in_sphere(
      file_name=self.expansion_file)

  def update_xyz(self, sites_cart):
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    self.expansion = self.expansion.update_xyz(sites_cart=sites_cart)
    self.pdb_hierarchy_super = self.expansion.ph_super_sphere

  def set_up_cluster_qm(self):
    self.get_clusters()
    if not self.clusters_only:
      self.get_fragments()
      self.get_fragment_hierarchies_and_charges()

  def get_clusters(self):
    n_residues=len(list(self.pdb_hierarchy.residue_groups()))
    # pair_interaction.run changes the coordinates of the hierarchy it gets
    # (Angstrom to Bohr): give it a copy.
    self.interaction_list = pair_interaction.run(self.pdb_hierarchy.deep_copy())
    from . import clustering
    clusters = clustering.betweenness_centrality_clustering(
      self.interaction_list,
      size = n_residues,
      maxnum_residues_in_cluster = self.maxnum_residues_in_cluster,
      bcc_threshold = self.bcc_threshold,
      bonded_pairs = bonded_residue_pairs(self.pdb_hierarchy)).get_clusters()
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
    for cluster in clusters:
      # pair_interaction.run changes the coordinates of the hierarchy it gets
      atoms_in_one_cluster, atoms_in_one_fragment, molecules_in_one_fragment = \
        pair_interaction.run(ph.deep_copy(), cluster)
      cluster_atoms_in_ph.append(atoms_in_one_cluster)
      fragment_super_atoms_in_ph.append(atoms_in_one_fragment)
      molecules_in_fragments.append(molecules_in_one_fragment)
    if(self.two_buffers): # a second buffer layer around the first
      fragment_super_atoms_in_ph = []
      for molecules in molecules_in_fragments:
        atoms_in_one_fragment = pair_interaction.run(ph.deep_copy(), molecules)[1]
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

  def make_capper(self, fragment_super_selection):
    """
    The capper of the fragment fragment_super_selection (a selection on
    pdb_hierarchy_super, the super sphere). Cut bonds are found against the
    super cell, not the sphere: a symmetry-copy residue at the edge of the
    sphere has its peptide partner outside the sphere, a cut bond (one H),
    not an unsatisfied valence.
    """
    SCM = self.expansion.SCM
    keep = SCM.keep_selection.iselection()
    selection = flex.bool(SCM.super_cell_hierarchy.atoms().size(),
      keep.select(fragment_super_selection.iselection()))
    c = fragment_capping.capper(
      parent_hierarchy = SCM.super_cell_hierarchy,
      selection        = selection,
      crystal_symmetry = self.expansion.cs_box,
      cif_objects      = self.cif_objects)
    assert c.fragment_hierarchy.is_similar_hierarchy(
      self.pdb_hierarchy_super.select(fragment_super_selection))
    return c

  def get_fragment_hierarchies_and_charges(self):

    def pdb_hierarchy_select(atoms_size, selection):
      selection_array = flex.bool(atoms_size, False)
      for item in selection:
        if(item<=atoms_size):
          selection_array[item-1] = True
      return selection_array

    self.fragment_selections = []
    self.fragment_super_selections = []
    self.fragment_charges = []
    self.cluster_selections = []
    self.buffer_selections = []
    self.fragment_capped_initial = []
    self.fragment_cappers = []
    for i in range(len(self.fragment_super_atoms)):
      fragment_selection = pdb_hierarchy_select(
          self.pdb_hierarchy.atoms_size(),
          self.fragment_super_atoms[i])
      ## QM part is fragment_super
      fragment_super_selection = pdb_hierarchy_select(
        self.pdb_hierarchy_super.atoms_size(),
        self.fragment_super_atoms[i])
      fragment_super_hierarchy = self.pdb_hierarchy_super.select(
        fragment_super_selection)
      # caps (and the charge of the capped fragment) as reduce2 capping=True
      # places them; the capper re-places them at every gradient evaluation
      if self.cap_fragments:
        c = self.make_capper(fragment_super_selection)
        self.fragment_cappers.append(c)
        self.fragment_capped_initial.append(c.cap().pdb_hierarchy)
        fragment_charge = c.charge
      else:
        self.fragment_cappers.append(None)
        self.fragment_capped_initial.append(fragment_super_hierarchy)
        fragment_charge = None
      self.fragment_super_selections.append(fragment_super_selection)
      #
      self.fragment_selections.append(fragment_selection)
      self.fragment_charges.append(fragment_charge)
      cluster_selection = pdb_hierarchy_select(
        self.pdb_hierarchy.atoms_size(), self.cluster_atoms[i])
      self.cluster_selections.append(cluster_selection)
      s = fragment_selection==cluster_selection
      buffer_selection = fragment_selection.deep_copy().set_selected(s, False)
      self.buffer_selections.append(buffer_selection)
      check_residue_integrity(self.pdb_hierarchy_super, fragment_super_selection)
      check_residue_integrity(self.pdb_hierarchy, cluster_selection)

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
      fragment_capped_initial   = self.fragment_capped_initial,
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
      super_sphere_geometry_restraints_manager = super_sphere_grm)

def get_qm_file_name_and_pdb_hierarchy(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  fragment_hierarchy = fragment_extracts.pdb_hierarchy_super.select(
    fragment_selection)
  sub_working_folder = fragment_extracts.working_folder + "/"+ str(index) + "/"
  if (not os.path.isdir(sub_working_folder)):
    os.mkdir(sub_working_folder)
  qm_pdb_file = sub_working_folder + str(index) + ".pdb"
  complete_qm_pdb_file = qm_pdb_file[:-4] + "_capping.pdb"
  # The fragment (user information) and the capped fragment (user
  # information; its name is the QM engine's label) are written at every
  # call; the QM engine gets the capped hierarchy, not the file.
  fragment_hierarchy.write_pdb_file(
    file_name=qm_pdb_file,
    crystal_symmetry=fragment_extracts.expansion_cs)
  # re-capping because geometry of the fragment has changed.
  r = fragment_extracts.fragment_cappers[index].cap(
    sites_cart = fragment_hierarchy.atoms().extract_xyz())
  r.pdb_hierarchy.write_pdb_file(file_name=complete_qm_pdb_file,
                    crystal_symmetry=fragment_extracts.expansion_cs)
  return os.path.abspath(complete_qm_pdb_file), r.pdb_hierarchy, \
    r.fragment_selection

def charge(fragment_extracts, index):
  return fragment_extracts.fragment_charges[index]

def write_mm_charge_file(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  file_name = None
  if (fragment_extracts.charge_embedding is True):
    if(fragment_extracts.charge_cutoff is not None):
      if(fragment_extracts.debug):
        print("charge_cutoff: ",fragment_extracts.charge_cutoff)
      xrs_super = fragment_extracts.pdb_hierarchy_super.extract_xray_structure()
      non_fragment_selection_super = xrs_super.selection_within(
        radius=fragment_extracts.charge_cutoff,
        selection=fragment_selection)
      non_fragment_selection_super = non_fragment_selection_super&~fragment_selection
      non_fragment_hierarchy_super = fragment_extracts.pdb_hierarchy_super.\
                       select(non_fragment_selection_super)
    else:
      non_fragment_hierarchy_super = fragment_extracts.pdb_hierarchy_super.\
                       select(~fragment_selection)
    non_fragment_hierarchy = non_fragment_hierarchy_super
    ph = fragment_extracts.pdb_hierarchy_super
    sub_working_folder = fragment_extracts.working_folder + "/" + str(index) + "/"
    if (not os.path.isdir(sub_working_folder)):
      os.mkdir(sub_working_folder)
    if(fragment_extracts.debug): print("write mm pdb file:", index)
    non_fragment_pdb_file = sub_working_folder + str(index) + "_mm.pdb"
    non_fragment_hierarchy.write_pdb_file(
      file_name=non_fragment_pdb_file,
      crystal_symmetry=fragment_extracts.expansion_cs)
    charge_scaling_positions = fragment_utils.get_edge_atom_positions(
      ph, non_fragment_hierarchy, charge_embed=True)
    fragment_extracts.charge_service.update_pdb_hierarchy(
      non_fragment_hierarchy,
      fragment_extracts.expansion_cs,
    )
    if(fragment_extracts.qm_engine_name == "turbomole"):
      file_name = sub_working_folder + str(index) + "_xyzq_cctbx.dat"
      fragment_extracts.charge_service.write_pdb_hierarchy_xyzq_file(
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    if(fragment_extracts.qm_engine_name in ["terachem","xtb",'mopac']):
      file_name = sub_working_folder + str(index) + "_qxyz_cctbx.dat"
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
  for index, selection_fragment in enumerate(F.fragment_selections):
    cluster_selection = F.cluster_selections[index]
    frag_selection = F.fragment_super_selections[index]
    index_cluster = F.pdb_hierarchy.select(cluster_selection)
    index_frag = F.pdb_hierarchy_super.select(frag_selection)
    filename_cluster = "%s_cluster.pdb" %(index)
    filename_frag = "%s_frag.pdb" %(index)
    filename_capped = "%s_capped0.pdb" %(index)
    index_cluster.write_pdb_file(
    file_name        = filename_cluster,
    crystal_symmetry = F.expansion_cs)
    index_frag.write_pdb_file(
    file_name        = filename_frag,
    crystal_symmetry = F.expansion_cs)
    capper = F.fragment_cappers[index]
    if capper is not None: # cctbx restraints: no caps (cap_fragments=False)
      capped_hierarchy = capper.cap(
        sites_cart = index_frag.atoms().extract_xyz()).pdb_hierarchy
      capped_hierarchy.write_pdb_file(file_name=filename_capped,
                crystal_symmetry=F.expansion_cs)

  log=open('fragment_info.txt','w')
  print('~  # clusters  : ',len(F.cluster_atoms), file=log)
  print('~  list of atoms per cluster:', file=log)
  print('~   ',[len(x) for x in F.cluster_atoms], file=log)
  print('~  list of atoms per fragment:', file=log)
  print('~   ',[len(x) for x in F.fragment_super_atoms], file=log)
  os.chdir(cwd)
