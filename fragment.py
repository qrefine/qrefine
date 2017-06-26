import os
import time
import libtbx.load_env
from scitbx.array_family import flex
from charges import get_total_charge_from_pdb
from charges import write_pdb_hierarchy_qxyz_file
from charges import write_pdb_hierarchy_xyzq_file
from utils import fragment_utils
from libtbx import group_args
from qrefine.super_cell import expand
import qrefine.completion as model_completion
from qrefine.plugin.yoink.pyoink import PYoink
from qrefine.utils.yoink_utils import write_yoink_infiles
import completion

qrefine = libtbx.env.find_in_repositories("qrefine")

class fragments(object):

  def __init__(self,
      working_folder             = "ase",
      clustering_method          = None,
      maxnum_residues_in_cluster = 20,
      charge_embedding           = False,
      two_buffers                = False,
      pdb_hierarchy              = None,
      qm_engine_name             = None,
      crystal_symmetry           = None,
      clustering                 = True):
    self.charge_embedding = charge_embedding
    self.two_buffers = two_buffers
    self.crystal_symmetry = crystal_symmetry
    self.working_folder = os.path.abspath(working_folder)
    self.pdb_hierarchy = pdb_hierarchy
    self.qm_engine_name = qm_engine_name
    self.clustering_method = clustering_method
    self.maxnum_residues_in_cluster =  maxnum_residues_in_cluster
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    self.backbone_connections = fragment_utils.get_backbone_connections(
      self.pdb_hierarchy)
    self.super_cell = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = 10.0)
    self.pdb_hierarchy_super = self.super_cell.ph_super_sphere
    if(clustering):
      self.yoink_dat_path = os.path.join(qrefine,"plugin","yoink","dat")
      self.pyoink = PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"))
      #t0 = time.time()
      self.set_up_cluster_qm()
      #print "time taken for interaction graph",(time.time() - t0)


  def set_up_cluster_qm(self, sites_cart=None):
    if(sites_cart is not None):
      ## update the selection of super_cell_sphere, and its
      ## geometry_restraints_manager and pdb_hierarchy
      self.pdb_hierarchy_super = self.super_cell.update(
        sites_cart=sites_cart).ph_super_sphere
    ###get clusters and their buffer regions using yoink and graph clustering.
    try:
      pre_clusters = self.clusters
    except:
      pre_clusters = None
    self.get_clusters()
    if pre_clusters!=self.clusters:
      self.get_fragments()
      self.get_fragment_hierarchies_and_charges()

  def get_clusters(self):
    self.cluster_file_name = self.working_folder + "/cluster.xml"
    self.qmmm_file_name = self.working_folder + "/qmmm.xml"
    ##  write yoink input file to get interactions
    write_yoink_infiles(self.cluster_file_name, self.qmmm_file_name,
                        self.pdb_hierarchy, self.yoink_dat_path)
    self.pyoink.input_file = self.cluster_file_name
    self.pyoink.update()
    self.interaction_list, weight = self.pyoink.get_interactions_list()
    self.interacting_pairs = len(self.interaction_list)
    self.interaction_list += self.backbone_connections

    import clustering
    #t0 = time.time()
    self.clustering = clustering.betweenness_centrality_clustering(
      self.interaction_list,
      size=len(self.pyoink.molecules),
      maxnum_residues_in_cluster=self.maxnum_residues_in_cluster)
    clusters = self.clustering.get_clusters()
    self.clusters = sorted(clusters,
      lambda x, y: 1 if len(x) < len(y) else -1 if len(x) > len(y) else 0)
    #print "time taken for clustering", (time.time() - t0)

  def get_fragments(self):

    def selected_atom_indices_in_entire_ph(selected_atom_indices_in_sub_ph, sub_ph):
      selected_atom_indices_in_entire_ph = []
      for index, number in enumerate(sub_ph.atoms().extract_serial()):
        if(index+1 in selected_atom_indices_in_sub_ph):
          selected_atom_indices_in_entire_ph.append(int(number))
      return selected_atom_indices_in_entire_ph
    self.pdb_hierarchy_super.atoms_reset_serial()
    phs = [self.pdb_hierarchy_super]
    altloc_size = self.pdb_hierarchy_super.altloc_indices().size()
    if(altloc_size>1):
      ## make pdb_hierarchy for each altloc case
      phs = []
      asc = self.pdb_hierarchy_super.atom_selection_cache()
      ## the first one altloc is " ", A B.. altlocs start from 1
      for altloc in self.pdb_hierarchy_super.altloc_indices().keys()[1:]:
        sel = asc.selection("altloc %s or altloc ' '"%altloc)
        ph_altloc = self.pdb_hierarchy_super.select(sel)
        phs.append(ph_altloc)
    cluster_atoms_in_phs = []
    fragment_super_atoms_in_phs = []
    clusters = self.clusters##from graph clustring, molecular indices
    ##loop over each cluster in every pdb_hierarchy to define buffer region
    ##fragment consists of cluster and buffer
    ##all pdb_hierarchies have the same clusters at molecular level
    for ph in phs:
      pyoink = self.pyoink
      cluster_atoms_in_ph = []
      fragment_super_atoms_in_ph = []
      ## write yoink input file to get fragment
      write_yoink_infiles(self.cluster_file_name, self.qmmm_file_name,
                          ph, self.yoink_dat_path)
      molecules_in_fragments = []
      for i in range(len(clusters)):
        pyoink.input_file = self.qmmm_file_name
        pyoink.update(clusters[i])
        atoms_in_one_cluster = pyoink.qm_core_fixed_indices
        atoms_in_one_cluster = selected_atom_indices_in_entire_ph(
                                                    atoms_in_one_cluster, ph)
        cluster_atoms_in_ph.append(atoms_in_one_cluster)
        atoms_in_one_fragment, molecules_in_one_fragment = pyoink.get_qm_indices()
        atoms_in_one_fragment = selected_atom_indices_in_entire_ph(
                                                     atoms_in_one_fragment, ph)
        fragment_super_atoms_in_ph.append(atoms_in_one_fragment)
        molecules_in_fragments.append(molecules_in_one_fragment)
        if(0):
          print i, "atoms in cluster: ", atoms_in_one_cluster
      if(self.two_buffers):## define a second buffer layer
        fragment_super_atoms_in_ph = []
        for molecules in molecules_in_fragments:
          pyoink.input_file = self.qmmm_file_name
          pyoink.update(molecules)
          atoms_in_one_fragment, junk = pyoink.get_qm_indices()
          atoms_in_one_fragment = selected_atom_indices_in_entire_ph(
                                                     atoms_in_one_fragment, ph)
          fragment_super_atoms_in_ph.append(atoms_in_one_fragment)
      cluster_atoms_in_phs.append(cluster_atoms_in_ph)
      fragment_super_atoms_in_phs.append(fragment_super_atoms_in_ph)
    ##check conformers and get all clusters and fragments
    self.cluster_atoms = []
    self.fragment_super_atoms = []
    self.fragment_scales = []
    for i_cluster in range(len(clusters)):
      ##always collect the clustering result from phs[0]
      self.collect_cluster_and_fragment(cluster_atoms_in_phs,
                                    fragment_super_atoms_in_phs, i_cluster, 0)
      if(len(phs)>1):
        for j_ph in xrange(1, len(phs)):
          fragment_same = (set(fragment_super_atoms_in_phs[0][i_cluster])==
               set(fragment_super_atoms_in_phs[j_ph][i_cluster]))
          if(fragment_same):continue
          self.collect_cluster_and_fragment(cluster_atoms_in_phs,
                              fragment_super_atoms_in_phs, i_cluster, j_ph)
          overlap_atoms_in_one_cluster = self.atoms_overlap(
                                        cluster_atoms_in_phs, i_cluster, j_ph)
          overlap_atoms_in_one_fragment = self.atoms_overlap(
                                  fragment_super_atoms_in_phs, i_cluster, j_ph)
          self.cluster_atoms.append(list(overlap_atoms_in_one_cluster))
          self.fragment_super_atoms.append(list(overlap_atoms_in_one_fragment))
          self.fragment_scales.append(-1.0)#substract the contribution from overlap

  def atoms_overlap(self, cluster_atoms_in_phs, i_cluster, j_ph):
    overlap_atoms_in_one_cluster = set(cluster_atoms_in_phs[0][i_cluster]) & \
                                   set(cluster_atoms_in_phs[j_ph][i_cluster])
    return overlap_atoms_in_one_cluster

  def collect_cluster_and_fragment(self, cluster_atoms_in_phs,
                                fragment_super_atoms_in_phs, i_cluster, j_ph):
    self.cluster_atoms.append(cluster_atoms_in_phs[j_ph][i_cluster])
    self.fragment_super_atoms.append(fragment_super_atoms_in_phs[j_ph][i_cluster])
    self.fragment_scales.append(1.0)

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
    for i in range(len(self.fragment_super_atoms)):
      fragment_selection = pdb_hierarchy_select(
          self.pdb_hierarchy.atoms_size(), self.fragment_super_atoms[i])
      ## QM part is fragment_super
      fragment_super_selection = pdb_hierarchy_select(
        self.pdb_hierarchy_super.atoms_size(), self.fragment_super_atoms[i])
      qm_pdb_hierarchy = self.pdb_hierarchy_super.select(fragment_super_selection)
      raw_records = qm_pdb_hierarchy.as_pdb_string(
        crystal_symmetry=self.super_cell.cs_box)
      ##TODO: remove if-esle statement
      ## finalise could not complete altloc yet
      ## charge calculation needs a completed structure
      ## temprary solution: skip charge calculation if conformers in PDB
      if(self.pdb_hierarchy_super.altloc_indices().size()>1):
        charge = 0
      else:
        charge = get_total_charge_from_pdb(raw_records=raw_records)
      self.fragment_super_selections.append(fragment_super_selection)
      #
      self.fragment_selections.append(fragment_selection)
      self.fragment_charges.append(charge)
      cluster_selection = pdb_hierarchy_select(
        self.pdb_hierarchy.atoms_size(), self.cluster_atoms[i])
      self.cluster_selections.append(cluster_selection)
      s = fragment_selection==cluster_selection
      buffer_selection = fragment_selection.deep_copy().set_selected(s, False)
      self.buffer_selections.append(buffer_selection)
      if(0):
        print "write pdb files for cluster and fragment"
        qm_pdb_hierarchy.write_pdb_file(file_name=str(i)+"_frag.pdb",
          crystal_symmetry=self.super_cell.cs_box)
        cluster_pdb_hierarchy = self.pdb_hierarchy_super.select(cluster_selection)
        cluster_pdb_hierarchy.write_pdb_file(file_name=str(i)+"_cluster.pdb",
          crystal_symmetry=self.super_cell.cs_box)

def get_qm_file_name_and_pdb_hierarchy(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  fragment_hierarchy = fragment_extracts.pdb_hierarchy_super.select(
    fragment_selection)
  sub_working_folder = fragment_extracts.working_folder + "/"+ str(index) + "/"
  if (not os.path.isdir(sub_working_folder)):
    os.mkdir(sub_working_folder)
  qm_pdb_file = sub_working_folder + str(index) + ".pdb"
  complete_qm_pdb_file = qm_pdb_file[:-4] + "_capping.pdb"
  ph = completion.run(pdb_hierarchy=fragment_hierarchy,
                      crystal_symmetry=fragment_extracts.super_cell.cs_box,
                      model_completion=False)
  ##for debugging
  if(0):  ## for degugging
    fragment_hierarchy.write_pdb_file(
      file_name=qm_pdb_file,
      crystal_symmetry=fragment_extracts.super_cell.cs_box)
    ph.write_pdb_file(file_name=complete_qm_pdb_file)
  return os.path.abspath(complete_qm_pdb_file), ph

def charge(fragment_extracts, index):
  return fragment_extracts.fragment_charges[index]

def write_mm_charge_file(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  file_name = None
  if (fragment_extracts.charge_embedding is True):
    non_fragment_hierarchy = \
      fragment_extracts.pdb_hierarchy_super.select(~fragment_selection)
    sub_working_folder = fragment_extracts.working_folder + str(index) + "/"
    if (not os.path.isdir(sub_working_folder)):
      os.mkdir(sub_working_folder)
    non_fragment_pdb_file = sub_working_folder + str(index) + "_mm.pdb"
    non_fragment_hierarchy.write_pdb_file(
      file_name=non_fragment_pdb_file,
      crystal_symmetry=fragment_extracts.crystal_symmetry)
    pdb_hierarchy = fragment_extracts.pdb_hierarchy_super
    non_qm_edge_positions = fragment_utils.get_edge_atom_positions(
      pdb_hierarchy, non_fragment_hierarchy, charge_embed=True)
    charge_scaling_positions = non_qm_edge_positions
    if (fragment_extracts.qm_engine_name == "turbomole"):
      file_name = sub_working_folder + str(index) + "_xyzq_cctbx.dat"
      write_pdb_hierarchy_xyzq_file(
        non_fragment_hierarchy,
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    if (fragment_extracts.qm_engine_name == "terachem"):
      file_name = sub_working_folder + str(index) + "_qxyz_cctbx.dat"
      write_pdb_hierarchy_qxyz_file(
        non_fragment_hierarchy,
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    file_name = os.path.abspath(file_name)
  return file_name

def fragment_extracts(fragments):
  return group_args(
    fragment_charges     = fragments.fragment_charges,
    fragment_selections  = fragments.fragment_selections,
    fragment_super_selections=fragments.fragment_super_selections,
    super_sphere_geometry_restraints_manager= \
          fragments.super_cell.super_sphere_geometry_restraints_manager,
    working_folder       = fragments.working_folder,
    fragment_super_atoms      = fragments.fragment_super_atoms,
    cluster_atoms        = fragments.cluster_atoms,
    qm_engine_name       = fragments.qm_engine_name,
    charge_embedding     = fragments.charge_embedding,
    crystal_symmetry     = fragments.crystal_symmetry,
    pdb_hierarchy        = fragments.pdb_hierarchy,
    pdb_hierarchy_super  = fragments.pdb_hierarchy_super,
    super_cell           = fragments.super_cell,
    buffer_selections    = fragments.buffer_selections,
    fragment_scales      = fragments.fragment_scales)
