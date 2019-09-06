import os
import time
import copy
import itertools
import libtbx.load_env
from libtbx.utils import Sorry
from scitbx.array_family import flex
from utils import fragment_utils
from libtbx import group_args
from qrefine.super_cell import expand
import qrefine.completion as model_completion
import completion
from charges import charges_class
from mmtbx.pair_interaction import pair_interaction

qrefine = libtbx.env.find_in_repositories("qrefine")

def check_atoms_integrity(atoms, verbose=False):
  rc = {}
  for atom in atoms:
    resid = atom.parent().parent().id_str()
    rc.setdefault(resid, [])
    if atom.name in [' CA ', ' CB ']:
      if atom.quote().find('GLY')==-1:
        rc[resid].append(atom.quote())
  for key, item in rc.items():
    if verbose: print key, item
    assert len(item) in [0,2], 'error in cluster %s %s' % (key, item)

def check_selection_integrity(atoms, indices, verbose=False):
  selection = []
  for j in indices:
    atom = atoms[j-1]
    if verbose: print atom.quote()
    selection.append(atom)
  check_atoms_integrity(selection, verbose=verbose)

def check_hierarchy(hierarchy, verbose=False):
  check_atoms_integrity(hierarchy.atoms(), verbose=verbose)

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
      cif_objects                = None,
      debug                      = False,
      charge_cutoff              = 8.0,
      save_clusters              = False,
      fast_interaction           = False):
    self.clustering = clustering
    self.fast_interaction = fast_interaction
    self.charge_embedding = charge_embedding
    self.two_buffers = two_buffers
    self.crystal_symmetry = crystal_symmetry
    self.working_folder = os.path.abspath(working_folder)
    self.pdb_hierarchy = pdb_hierarchy
    self.charge_cutoff = charge_cutoff
    self.system_size = pdb_hierarchy.atoms_size()
    self.qm_engine_name = qm_engine_name
    self.clustering_method = clustering_method
    self.altloc_method =  altloc_method
    self.debug = debug
    self.maxnum_residues_in_cluster =  maxnum_residues_in_cluster
    self.save_clusters = save_clusters
    raw_records = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
    self.charge_service = charges_class(
        raw_records           = raw_records,
        ligand_cif_file_names = cif_objects)
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    self.backbone_connections = fragment_utils.get_backbone_connections(
      self.pdb_hierarchy)
    self.get_altloc_molecular_indices()
    if(1):
      self.altloc_atoms = [atom for atom in list(pdb_hierarchy.atoms())
                           if atom.pdb_label_columns()[4]!=" "]
    self.expansion = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = 10.0)
    self.pdb_hierarchy_super = self.expansion.ph_super_sphere
    ## write expansion.pdb as the reference for capping
    self.expansion_file = "expansion.pdb"
    self.expansion.write_super_cell_selected_in_sphere(file_name=self.expansion_file)
    if(not self.fast_interaction):
      from qrefine.plugin.yoink.pyoink import PYoink
      from qrefine.utils.yoink_utils import write_yoink_infiles
      self.yoink_dat_path = os.path.join(qrefine,"plugin","yoink","dat")
      self.pyoink = PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"))
    self.qm_run = qm_run
    #t0 = time.time()
    self.set_up_cluster_qm()
    #print "time taken for interaction graph",(time.time() - t0)
  
  def update_xyz(self,sites_cart):
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    pre_atoms=[atom.pdb_label_columns() for atom in self.pdb_hierarchy_super.atoms()]
    self.expansion = self.expansion.update_xyz(
                                 sites_cart=sites_cart)
    self.pdb_hierarchy_super = self.expansion.ph_super_sphere
    new_atoms=[ atom.pdb_label_columns() for atom in self.pdb_hierarchy_super.atoms()]
    #if(self.expansion.ph_super_sphere.atoms_size()!=pre_size): 
    if(pre_atoms!=new_atoms):
      if(self.debug):
        print("the content of the super sphere has been changed,reset up fragments")
      #Note: the atom size of self.expansion.ph_super_sphere gets changeed,
      #while the atom size in super_sphere_geometry_restraints_manager
      #does not get changed. Re-generate the object of expand
      if(1):
        self.expansion = expand(
          pdb_hierarchy        = self.pdb_hierarchy,
          crystal_symmetry     = self.crystal_symmetry,
          select_within_radius = 10.0)
        self.pdb_hierarchy_super = self.expansion.ph_super_sphere
      self.get_fragments()
      self.get_fragment_hierarchies_and_charges()

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
      ## update the selection of expansion_sphere, and its
      ## geometry_restraints_manager and pdb_hierarchy
      self.pdb_hierarchy_super = self.expansion.update(
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
    print(self.clustering)
    n_residues=len(list(self.pdb_hierarchy.residue_groups()))
    if(not self.clustering):
      return(range(1,n_residues+1,1) )
    if(self.fast_interaction):
      self.interaction_list = pair_interaction.run(copy.deepcopy(self.pdb_hierarchy))  ##deepcopy
    else: # to be deprecated.
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
    # t0 = time.time()
    self.clustering = clustering.betweenness_centrality_clustering(
      self.interaction_list,
      size = n_residues,
      maxnum_residues_in_cluster = self.maxnum_residues_in_cluster)
    clusters = self.clustering.get_clusters()
    self.clusters = sorted(clusters,
      lambda x, y: 1 if len(x) < len(y) else -1 if len(x) > len(y) else 0)
    # print "time taken for clustering", (time.time() - t0)

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
      altlocs = self.pdb_hierarchy_super.altloc_indices().keys()
      altlocs.sort()
      for altloc in altlocs[1:]:
        sel = asc.selection("altloc '%s' or altloc '' or altloc ' '"%altloc)
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
      cluster_atoms_in_ph = []
      fragment_super_atoms_in_ph = []
      molecules_in_fragments = []
      ## write yoink input file to get fragment
      if(not self.fast_interaction):
        pyoink = self.pyoink
        write_yoink_infiles(self.cluster_file_name,
                          self.qmmm_file_name,
                          ph,
                          self.yoink_dat_path)
      for i in range(len(clusters)):
        # print 'processing cluster', i
        if(self.fast_interaction):
          atoms_in_one_cluster, atoms_in_one_fragment, molecules_in_one_fragment = \
            pair_interaction.run(copy.deepcopy(ph), clusters[i])  ##deepcopy
          # print("clusters[i]",clusters[i])
          # print("molecules_in_one_fragment:", molecules_in_one_fragment)
          # print("atoms_in_one_fragment",atoms_in_one_fragment)
        else:
          pyoink.input_file = self.qmmm_file_name
          pyoink.update(clusters[i])
          atoms_in_one_cluster = pyoink.qm_core_fixed_indices
          atoms_in_one_fragment, molecules_in_one_fragment = pyoink.get_qm_indices()
      
        atoms_in_one_cluster = selected_atom_indices_in_entire_ph(
                                                    atoms_in_one_cluster, ph)
        cluster_atoms_in_ph.append(atoms_in_one_cluster)

        atoms_in_one_fragment = selected_atom_indices_in_entire_ph(
                                                     atoms_in_one_fragment, ph)
        fragment_super_atoms_in_ph.append(atoms_in_one_fragment)
        molecules_in_fragments.append(molecules_in_one_fragment)
        if(0):
          print i, "atoms in cluster: ", atoms_in_one_cluster
        if True:
          atoms = self.pdb_hierarchy_super.atoms()
          check_selection_integrity(atoms, atoms_in_one_cluster)
      # print "cluster->fragments done"
      if(self.two_buffers):## define a second buffer layer
        # print "adding second layer"
        fragment_super_atoms_in_ph = []
        for molecules in molecules_in_fragments:
          if(self.fast_interaction):
            junk1,atoms_in_one_fragment,junk2 = pair_interaction.run(copy.deepcopy(ph),molecules)
          else:
            pyoink.input_file = self.qmmm_file_name
            pyoink.update(list(molecules))
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
            # two different altloc clusters, the overlap is part of a residue,
            # even an atom
            # the cluster overlap will cause troubles for QM calculation,
            # expecially when it is an atom
              atoms = self.pdb_hierarchy_super.atoms()
              overlap_atoms_in_one_fragment = self.atoms_overlap(
                                  fragment_super_atoms_in_phs,
                                  i_cluster,
                                  j_ph)
              check_selection_integrity(atoms, overlap_atoms_in_one_fragment)
              self.cluster_atoms.append(list(overlap_atoms_in_one_cluster))
              self.fragment_super_atoms.append(list(overlap_atoms_in_one_fragment))
              scale_list = [-1.0]*sum(i <= self.system_size
                                      for i in overlap_atoms_in_one_fragment)
              self.fragment_scales.append(scale_list)
          ##average the contributions from overlap
          elif(self.altloc_method=="average"):
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
    self.cluster_selections = []
    self.fragment_capped_initial = []
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
      if(self.debug):
        fragment_super_hierarchy.write_pdb_file(file_name=str(i)+"-origin-cs.pdb")
        fragment_super_hierarchy.write_pdb_file(file_name=str(i)+".pdb",
          crystal_symmetry=self.expansion.cs_box)
      charge_hierarchy = completion.run(pdb_hierarchy=fragment_super_hierarchy,
                      crystal_symmetry=self.expansion.cs_box,
                      model_completion=False,
                      original_pdb_filename=self.expansion_file)
      self.fragment_capped_initial.append(charge_hierarchy)
      raw_records = charge_hierarchy.as_pdb_string(
        crystal_symmetry=self.expansion.cs_box)
      if(1):charge_hierarchy.write_pdb_file(file_name=str(i)+"_capping_tmp.pdb",
          crystal_symmetry=self.expansion.cs_box)

      self.charge_service.update_pdb_hierarchy(
        charge_hierarchy,
        self.expansion.cs_box,
      )
      #TODO: do not why self.charge_service could not right charge
      #charge = self.charge_service.get_total_charge()
      charge = charges_class(pdb_filename=str(i)+"_capping_tmp.pdb").get_total_charge()
      #the capping pdb file causes problem for tests, remove it
      os.remove(str(i)+"_capping_tmp.pdb")
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
      # if(self.debug):
      #   fragment_super_hierarchy.write_pdb_file(file_name=str(i)+"_frag.pdb",
      #     crystal_symmetry=self.expansion.cs_box)
      #   cluster_pdb_hierarchy = self.pdb_hierarchy.select(cluster_selection)
      #   cluster_pdb_hierarchy.write_pdb_file(file_name=str(i)+"_cluster.pdb",
      #     crystal_symmetry=self.expansion.cs_box)
      check_hierarchy(fragment_super_hierarchy)

def get_qm_file_name_and_pdb_hierarchy(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  fragment_hierarchy = fragment_extracts.pdb_hierarchy_super.select(
    fragment_selection)
  sub_working_folder = fragment_extracts.working_folder + "/"+ str(index) + "/"
  if (not os.path.isdir(sub_working_folder)):
    os.mkdir(sub_working_folder)
  qm_pdb_file = sub_working_folder + str(index) + ".pdb"
  complete_qm_pdb_file = qm_pdb_file[:-4] + "_capping.pdb"
  # if(fragment_extracts.debug):
  if(1): # at least for now
    fragment_hierarchy.write_pdb_file(
      file_name=qm_pdb_file,
      crystal_symmetry=fragment_extracts.expansion_cs)
  # re-capping because geometry of the fragment has changed.
  ph = completion.run(pdb_hierarchy=fragment_hierarchy,
                      crystal_symmetry=fragment_extracts.expansion_cs,
                      model_completion=False,
                      original_pdb_filename=fragment_extracts.expansion_file)
  # we now want this file by default  
  ph.write_pdb_file(file_name=complete_qm_pdb_file,
                    crystal_symmetry=fragment_extracts.expansion_cs)
  return os.path.abspath(complete_qm_pdb_file), ph

def charge(fragment_extracts, index):
  return fragment_extracts.fragment_charges[index]

def write_mm_charge_file(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  file_name = None
  if (fragment_extracts.charge_embedding is True):
    altlocs = fragment_extracts.pdb_hierarchy_super.altloc_indices().keys()
    altlocs.sort()
    if(fragment_extracts.charge_cutoff is not None):
      if(fragment_extracts.debug):
        print "charge_cutoff: ",fragment_extracts.charge_cutoff
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
    # the pdb has no altlocs
    if (len(altlocs)==1):
      non_fragment_hierarchy = non_fragment_hierarchy_super
      ph = fragment_extracts.pdb_hierarchy_super
    # the pdb has  altlocs
    else:
      fragment_altlocs = fragment_extracts.pdb_hierarchy_super.\
                        select(fragment_selection).altloc_indices().keys()
      fragment_altlocs.sort()
      asc_ph = fragment_extracts.pdb_hierarchy_super.atom_selection_cache()
      asc_non_fragment = non_fragment_hierarchy_super.atom_selection_cache()
      # the fragment has  altlocs
      if (len(fragment_altlocs)==2):
        sel_non_fragment = asc_non_fragment.\
          selection("altloc '%s' or altloc '' or altloc ' '"%fragment_altlocs[1])
        sel_ph = asc_ph.selection("altloc '%s' or altloc '' or altloc ' '"%fragment_altlocs[1])
      # the fragment has no altlocs
      else:
        sel_non_fragment = asc_non_fragment.\
          selection("altloc '%s' or altloc '' or altloc ' '"%fragment_altlocs[0])
        sel_ph = asc_ph.selection("altloc '%s' or altloc '' or altloc ' '"%fragment_altlocs[0])
      non_fragment_hierarchy = non_fragment_hierarchy_super.select(sel_non_fragment)
      ph = fragment_extracts.pdb_hierarchy_super.select(sel_ph)
    sub_working_folder = fragment_extracts.working_folder + "/" + str(index) + "/"
    if (not os.path.isdir(sub_working_folder)):
      os.mkdir(sub_working_folder)
    if(fragment_extracts.debug): print "write mm pdb file:", index
    non_fragment_pdb_file = sub_working_folder + str(index) + "_mm.pdb"
    non_fragment_hierarchy.write_pdb_file(
      file_name=non_fragment_pdb_file,
      crystal_symmetry=fragment_extracts.expansion_cs)
    non_qm_edge_positions = fragment_utils.get_edge_atom_positions(
      ph, non_fragment_hierarchy, charge_embed=True)
    charge_scaling_positions = non_qm_edge_positions
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

def fragment_extracts(fragments):
  return group_args(
    cluster_selections   = fragments.cluster_selections,
    fragment_charges     = fragments.fragment_charges,
    fragment_selections  = fragments.fragment_selections,
    fragment_super_selections=fragments.fragment_super_selections,
    fragment_capped_initial=fragments.fragment_capped_initial,
    super_sphere_geometry_restraints_manager= \
          fragments.expansion.super_sphere_geometry_restraints_manager,
    working_folder       = fragments.working_folder,
    fragment_super_atoms = fragments.fragment_super_atoms,
    cluster_atoms        = fragments.cluster_atoms,
    qm_engine_name       = fragments.qm_engine_name,
    charge_embedding     = fragments.charge_embedding,
    crystal_symmetry     = fragments.crystal_symmetry,
    pdb_hierarchy        = fragments.pdb_hierarchy,
    pdb_hierarchy_super  = fragments.pdb_hierarchy_super,
    expansion_cs         = fragments.expansion.cs_box,
    buffer_selections    = fragments.buffer_selections,
    fragment_scales      = fragments.fragment_scales,
    debug                = fragments.debug,
    charge_service       = fragments.charge_service,
    charge_cutoff        = fragments.charge_cutoff,
    expansion_file       = fragments.expansion_file,
    save_clusters        = fragments.save_clusters)

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
    # capped_hierarchy = completion.run(pdb_hierarchy=index_frag,
    #           crystal_symmetry=F.expansion_cs,
    #           model_completion=False,
    #           original_pdb_filename=filename_frag)
    capped_hierarchy = F.fragment_capped_initial[index]
    capped_hierarchy.write_pdb_file(file_name=filename_capped,
              crystal_symmetry=F.expansion_cs)

  log=open('fragment_info.txt','w')
  print >> log, '~  # clusters  : ',len(F.cluster_atoms)
  print >> log, '~  list of atoms per cluster:'
  print >> log, '~   ',[len(x) for x in F.cluster_atoms]
  print >> log, '~  list of atoms per fragment:'
  print >> log, '~   ',[len(x) for x in F.fragment_super_atoms]
  os.chdir(cwd)
