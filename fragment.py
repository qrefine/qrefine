import os
import time
import itertools 
import libtbx.load_env
from libtbx.utils import Sorry
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
      altloc_method              = None,
      maxnum_residues_in_cluster = 20,
      charge_embedding           = False,
      two_buffers                = False,
      pdb_hierarchy              = None,
      qm_engine_name             = None,
      crystal_symmetry           = None,
      clustering                 = True,
      qm_run                     = True,
      debug                      = False):
    self.charge_embedding = charge_embedding
    self.two_buffers = two_buffers
    self.crystal_symmetry = crystal_symmetry
    self.working_folder = os.path.abspath(working_folder)
    self.pdb_hierarchy = pdb_hierarchy
    self.system_size = pdb_hierarchy.atoms_size()
    self.qm_engine_name = qm_engine_name
    self.clustering_method = clustering_method
    self.altloc_method =  altloc_method
    self.debug = debug
    self.maxnum_residues_in_cluster =  maxnum_residues_in_cluster
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    self.backbone_connections = fragment_utils.get_backbone_connections(
      self.pdb_hierarchy)
    self.get_altloc_molecular_indices()
    self.super_cell = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = 10.0)
    self.pdb_hierarchy_super = self.super_cell.ph_super_sphere
    ## write super_cell.pdb as the reference for capping 
    self.super_cell_file = "super_cell.pdb"
    self.super_cell.ph_super_cell.write_pdb_file(file_name=self.super_cell_file) 
    if(1):
      self.altloc_atoms = [atom for atom in list(pdb_hierarchy.atoms())
                           if atom.pdb_label_columns()[4]!=" "]
    if(clustering):
      self.yoink_dat_path = os.path.join(qrefine,"plugin","yoink","dat")
      self.pyoink = PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"))
      self.qm_run = qm_run
      #t0 = time.time()
      self.set_up_cluster_qm()
      #print "time taken for interaction graph",(time.time() - t0)
  
  def get_altloc_molecular_indices(self):
    self.altloc_molecular_indices=[]
    index = 0
    for chain in self.pdb_hierarchy.chains():
      for residue_group in chain.residue_groups():
        index +=1
        if(residue_group.have_conformers()):
          self.altloc_molecular_indices.append(index)
    
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
    if(self.qm_run is True and pre_clusters!=self.clusters):
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
    ## isolate altloc molecules
    new_interaction_list = []
    if(0):
      for item in self.interaction_list:
        contain_altlocs = set(self.altloc_molecular_indices)&set(item)
        if(len(contain_altlocs)==0):
          new_interaction_list.append(item)
      self.interaction_list = new_interaction_list
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
      ## generate pdb_hierarchy for each altloc case
      phs = []
      asc = self.pdb_hierarchy_super.atom_selection_cache()
      if(self.debug):self.pdb_hierarchy_super.write_pdb_file(file_name="super.pdb")
      ## the first one altloc is " ", A B.. altlocs start from 1
      for altloc in self.pdb_hierarchy_super.altloc_indices().keys()[1:]:
        sel = asc.selection("altloc %s or altloc ' '"%altloc)
        ph_altloc = self.pdb_hierarchy_super.select(sel)
        phs.append(ph_altloc)
        if(self.debug):ph_altloc.write_pdb_file(file_name="super-"+str(altloc)+".pdb")
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
    #
    ##always collect the clustering result from phs[0]
    self.cluster_atoms = []
    self.fragment_super_atoms = []
    self.fragment_scales = []
    for i_cluster in range(len(clusters)):
      self.collect_cluster_and_fragment(cluster_atoms_in_phs,
                                    fragment_super_atoms_in_phs, i_cluster, 0)
    ##check alternative locations and get all clusters and fragments
    overlap_clusters={}
    overlap_fragments_super = {}
    if(len(phs)>1):
      for i_cluster in range(len(clusters)):
        for j_ph in xrange(1, len(phs)):
          fragment_same = (set(fragment_super_atoms_in_phs[0][i_cluster])==
               set(fragment_super_atoms_in_phs[j_ph][i_cluster]))
          # two same fragments for same non-altloc clusters
          if(fragment_same):continue
          # check the overlap
          overlap_atoms_in_one_cluster = self.atoms_overlap(
                                        cluster_atoms_in_phs, i_cluster, j_ph)
          empty_overlap_cluster = (len(overlap_atoms_in_one_cluster)==0)
          #substract the contribution from overlap
          if(self.altloc_method=="subtract"):
            self.collect_cluster_and_fragment(cluster_atoms_in_phs,
                              fragment_super_atoms_in_phs, i_cluster, j_ph)
            # different fragments for different altloc clusters
            if(empty_overlap_cluster):continue
            else:
            # two same non-altloc clusters, the overlap is a cluster
            # two different altloc clusters, the overlap is part of a residue, even an atom
            # the cluster overlap will cause troubles for QM calculation, expecially when it is an atom
              overlap_atoms_in_one_fragment = self.atoms_overlap(
                                  fragment_super_atoms_in_phs, i_cluster, j_ph)
              self.cluster_atoms.append(list(overlap_atoms_in_one_cluster))
              self.fragment_super_atoms.append(list(overlap_atoms_in_one_fragment))
              scale_list = [-1.0]*sum(i <= self.system_size
                                      for i in overlap_atoms_in_one_fragment)
              self.fragment_scales.append(scale_list)
          ##average the contributions from overlap
          if(self.altloc_method=="average"):
            # different fragments for different altloc clusters
            if(empty_overlap_cluster):
              self.collect_cluster_and_fragment(cluster_atoms_in_phs,
                              fragment_super_atoms_in_phs, i_cluster, j_ph)
            else:
            # two same non-altloc clusters, the overlap is a cluster
            # two different altloc clusters, the overlap is part of a residue, even an atom
            # collect all overlap clusters and fragments
              try:
                overlap_clusters[i_cluster] = overlap_clusters[i_cluster].append(
                                   cluster_atoms_in_phs[j_ph][i_cluster])
                overlap_fragments_super[i_cluster] = \
                  overlap_fragments_super[i_cluster].append(
                                   fragment_super_atoms_in_phs[j_ph][i_cluster])
              except:
                overlap_clusters[i_cluster] = \
                  [cluster_atoms_in_phs[j_ph][i_cluster]]
                overlap_fragments_super[i_cluster] = \
                  [fragment_super_atoms_in_phs[j_ph][i_cluster]]
    overlap_atoms = []
    for i_cluster, overlap_cluster in overlap_clusters.items():
      overlap_atoms = overlap_atoms+list(itertools.chain.from_iterable(
        overlap_cluster +[self.cluster_atoms[i_cluster]]))#[atom_index, atom_index]
    frequency_overlap_atoms =  {x:overlap_atoms.count(x) for x in overlap_atoms}#{atom_index,frequency}
    for i_cluster, clusters in overlap_clusters.items():
      ## reset the fragment scale for the ith fragment in ph[0]
      for index, atom in enumerate([i for i in self.fragment_super_atoms[i_cluster] 
                                       if i <= self.system_size]):
          if(atom in self.cluster_atoms[i_cluster] and
             atom in frequency_overlap_atoms.keys() and
               not self.bond_with_altloc(atom)):
              self.fragment_scales[i_cluster][index] = \
                1.0/frequency_overlap_atoms[atom]
      ## add overlap clusters and fragments
      for index, fragment_super in  enumerate(overlap_fragments_super[i_cluster]):
         scale_list = []
         for atom in [i for i in fragment_super if i <= self.system_size]:
           if(atom in clusters[index] and atom in frequency_overlap_atoms.keys()
                and not self.bond_with_altloc(atom)):
             scale_list.append(1.0/frequency_overlap_atoms[atom])
           else: scale_list.append(1.0)
         self.cluster_atoms.append(clusters[index]) 
         self.fragment_super_atoms.append(fragment_super)
         self.fragment_scales.append(scale_list)
             
  def bond_with_altloc(self, atom_index):
    ph_atoms = list(self.pdb_hierarchy.atoms())
    ph_atom = ph_atoms[atom_index-1]
    bond = False
    for altloc_atom in self.altloc_atoms:
      distance = ph_atom.distance(altloc_atom)
      if(distance<1.7):
        bond =True
        break
      ##TODO
      ##check bond, better from bond topology
    return bond
    
  def atoms_overlap(self, cluster_atoms_in_phs, i_cluster, j_ph):
    overlap_atoms_in_one_cluster = set(cluster_atoms_in_phs[0][i_cluster]) & \
                                   set(cluster_atoms_in_phs[j_ph][i_cluster])
    return overlap_atoms_in_one_cluster

  def collect_cluster_and_fragment(self, cluster_atoms_in_phs,
                                fragment_super_atoms_in_phs, i_cluster, j_ph):
    self.cluster_atoms.append(cluster_atoms_in_phs[j_ph][i_cluster])
    self.fragment_super_atoms.append(fragment_super_atoms_in_phs[j_ph][i_cluster])
    scale_list = [1.0]*sum(i <= self.system_size
                          for i in fragment_super_atoms_in_phs[j_ph][i_cluster])
    self.fragment_scales.append(scale_list)

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
      fragment_super_hierarchy = self.pdb_hierarchy_super.select(fragment_super_selection)
      if(self.debug):
        fragment_super_hierarchy.write_pdb_file(file_name=str(i)+"-origin-cs.pdb")
        fragment_super_hierarchy.write_pdb_file(file_name=str(i)+".pdb",
          crystal_symmetry=self.super_cell.cs_box)
      charge_hierarchy = completion.run(pdb_hierarchy=fragment_super_hierarchy,
                      crystal_symmetry=self.super_cell.cs_box,
                      model_completion=False,
                      original_pdb_filename=self.super_cell_file)
      raw_records = charge_hierarchy.as_pdb_string(
        crystal_symmetry=self.super_cell.cs_box)
      if(self.debug):charge_hierarchy.write_pdb_file(file_name=str(i)+"_capping.pdb",
          crystal_symmetry=self.super_cell.cs_box)
      ##TODO: remove if-esle statement
      ## temprary solution: skip charge calculation for an altloc pdb
      if(self.pdb_hierarchy_super.altloc_indices().size()>1):
        charge = 0
      if(1):
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
      if(self.debug):
        fragment_super_hierarchy.write_pdb_file(file_name=str(i)+"_frag.pdb",
          crystal_symmetry=self.super_cell.cs_box)
        cluster_pdb_hierarchy = self.pdb_hierarchy.select(cluster_selection)
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
  if(fragment_extracts.debug):  ## for degugging
    fragment_hierarchy.write_pdb_file(
      file_name=qm_pdb_file,
      crystal_symmetry=fragment_extracts.super_cell_cs)
  ph = completion.run(pdb_hierarchy=fragment_hierarchy,
                      crystal_symmetry=fragment_extracts.super_cell_cs,
                      model_completion=False,
                      original_pdb_filename=self.super_cell_file) 
  ##for debugging
  if(fragment_extracts.debug):  
    fragment_hierarchy.write_pdb_file(
      file_name=qm_pdb_file,
      crystal_symmetry=fragment_extracts.super_cell_cs)
    ph.write_pdb_file(file_name=complete_qm_pdb_file)
  return os.path.abspath(complete_qm_pdb_file), ph

def charge(fragment_extracts, index):
  return fragment_extracts.fragment_charges[index]

def write_mm_charge_file(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  file_name = None
  if (fragment_extracts.charge_embedding is True):
    altlocs = fragment_extracts.pdb_hierarchy_super.altloc_indices().keys()
    non_fragment_hierarchy_super = fragment_extracts.pdb_hierarchy_super.\
                       select(~fragment_selection)
    # the pdb has no altlocs
    if (len(altlocs)==1):
      non_fragment_hierarchy = non_fragment_hierarchy_super
      ph = fragment_extracts.pdb_hierarchy_super
    # the pdb has  altlocs
    else:
      fragment_altlocs = fragment_extracts.pdb_hierarchy_super.\
                        select(fragment_selection).altloc_indices().keys()
      asc_ph = fragment_extracts.pdb_hierarchy_super.atom_selection_cache()
      asc_non_fragment = non_fragment_hierarchy_super.atom_selection_cache()
      # the fragment has  altlocs
      if (len(fragment_altlocs)==2):
        sel_non_fragment = asc_non_fragment.\
          selection("altloc %s or altloc ' '"%fragment_altlocs[1])
        sel_ph = asc_ph.selection("altloc %s or altloc ' '"%fragment_altlocs[1])
      # the fragment has no altlocs
      else:
        sel_non_fragment = asc_non_fragment.\
          selection("altloc %s or altloc ' '"%fragment_altlocs[0])
        sel_ph = asc_ph.selection("altloc %s or altloc ' '"%fragment_altlocs[0])
      non_fragment_hierarchy = non_fragment_hierarchy_super.select(sel_non_fragment)
      ph = fragment_extracts.pdb_hierarchy_super.select(sel_ph)
    sub_working_folder = fragment_extracts.working_folder + "/" + str(index) + "/"
    if (not os.path.isdir(sub_working_folder)):
      os.mkdir(sub_working_folder)
    if(fragment_extracts.debug): print "write mm pdb file:", index
    non_fragment_pdb_file = sub_working_folder + str(index) + "_mm.pdb"
    non_fragment_hierarchy.write_pdb_file(
      file_name=non_fragment_pdb_file,
      crystal_symmetry=fragment_extracts.super_cell_cs)
    non_qm_edge_positions = fragment_utils.get_edge_atom_positions(
      ph, non_fragment_hierarchy, charge_embed=True)
    charge_scaling_positions = non_qm_edge_positions
    if(fragment_extracts.qm_engine_name == "turbomole"):
      file_name = sub_working_folder + str(index) + "_xyzq_cctbx.dat"
      write_pdb_hierarchy_xyzq_file(
        non_fragment_hierarchy,
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    if(fragment_extracts.qm_engine_name == "terachem"):
      file_name = sub_working_folder + str(index) + "_qxyz_cctbx.dat"
      write_pdb_hierarchy_qxyz_file(
        non_fragment_hierarchy,
        file_name=file_name,
        exclude_water=False,
        charge_scaling_positions=charge_scaling_positions)
    if(file_name is None):
      raise Sorry("There is no point charge file") 
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
    fragment_super_atoms = fragments.fragment_super_atoms,
    cluster_atoms        = fragments.cluster_atoms,
    qm_engine_name       = fragments.qm_engine_name,
    charge_embedding     = fragments.charge_embedding,
    crystal_symmetry     = fragments.crystal_symmetry,
    pdb_hierarchy        = fragments.pdb_hierarchy,
    pdb_hierarchy_super  = fragments.pdb_hierarchy_super,
    super_cell_cs        = fragments.super_cell.cs_box,
    buffer_selections    = fragments.buffer_selections,
    fragment_scales      = fragments.fragment_scales,
    debug                = fragments.debug,
    super_cell_file      = fragments.super_cell_file)
