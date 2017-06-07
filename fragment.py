import os
import libtbx.load_env
from scitbx.array_family import flex
from charges import get_total_charge_from_pdb
from charges import write_pdb_hierarchy_qxyz_file
from charges import write_pdb_hierarchy_xyzq_file
from utils import fragment_utils
from libtbx import group_args
import completion
try:
  from plugin.yoink.pyoink import PYoink
  from utils.yoink_utils import write_yoink_infiles
except ImportError, e:
  print str(e)
  pass

if (os.getenv("QR_REPO_PARENT") is not None):
  qr_yoink_path = os.getenv("QR_REPO_PARENT")+"/./qr-core/plugin/yoink/"
else:
  qrefine_path = libtbx.env.find_in_repositories("qrefine")
  qr_path = os.path.join(qrefine_path, "core")
  qr_yoink_path =os.path.join(qr_path, "plugin/yoink/")

class fragments(object):

  def __init__(self,
      working_folder             = "./ase/",
      clustering_method          = None,
      maxnum_residues_in_cluster = 20,
      charge_embedding           = False,
      pdb_hierarchy              = None,
      qm_engine_name             = None,
      crystal_symmetry           = None,
      clustering                 = True):
    self.charge_embedding = charge_embedding
    self.crystal_symmetry = crystal_symmetry
    self.working_folder = working_folder
    self.pdb_hierarchy = pdb_hierarchy
    self.qm_engine_name = qm_engine_name
    self.clustering_method = clustering_method
    self.maxnum_residues_in_cluster =  maxnum_residues_in_cluster
    if(os.path.exists(self.working_folder) is not True):
      os.mkdir(self.working_folder)
    self.backbone_connections = fragment_utils.get_backbone_connections(
      self.pdb_hierarchy)
    from utils.super_cell import expand
    self.expand = expand(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = 10.0)
    if(clustering):
      self.yoink_dat_path = os.path.join(qr_path,"plugin/yoink/dat")
      self.pyoink = PYoink(os.path.join(qr_path,"plugin/yoink/Yoink-0.0.1.jar"))
      self.set_up_cluster_qm()

  def set_up_cluster_qm(self, sites_cart=None):
    if(sites_cart is not None):
      self.pdb_hierarchy_super = self.expand.update(sites_cart=sites_cart)
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
    self.cluster_file_name = self.working_folder + "./cluster.xml"
    self.qmmm_file_name = self.working_folder + "./qmmm.xml"
    ##  write yoink input file to get interactions
    write_yoink_infiles(self.cluster_file_name, self.qmmm_file_name,
                        self.pdb_hierarchy, self.yoink_dat_path)
    self.pyoink.input_file = self.cluster_file_name
    self.pyoink.update()
    self.interaction_list, weight = self.pyoink.get_interactions_list()
    self.interacting_pairs = len(self.interaction_list)
    self.interaction_list += self.backbone_connections
    self.clustering = self.clustering_method(
      self.interaction_list,
      size=len(self.pyoink.molecules),
      maxnum_residues_in_cluster=self.maxnum_residues_in_cluster)
    clusters = self.clustering.get_clusters()
    self.clusters = sorted(clusters,
      lambda x, y: 1 if len(x) < len(y) else -1 if len(x) > len(y) else 0)


  def get_fragments(self):
    self.cluster_atoms = []
    self.fragment_super_atoms = []
    pyoink = self.pyoink
    clusters = self.clusters##from graph clustring, molecular indices
    sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.pdb_hierarchy_super = self.expand.update(sites_cart=sites_cart)
    ## write yoink input file to get fragment
    write_yoink_infiles(self.cluster_file_name, self.qmmm_file_name,
                        self.pdb_hierarchy_super, self.yoink_dat_path)
    ##TODO conformer check
    for i in range(len(clusters)):
      pyoink.input_file = self.qmmm_file_name
      pyoink.update(clusters[i])
      atoms_in_one_cluster = pyoink.qm_core_fixed_indices
      self.cluster_atoms.append(atoms_in_one_cluster)
      atoms_in_one_fragment, qm_molecules = pyoink.get_qm_indices()
      self.fragment_super_atoms.append(atoms_in_one_fragment)

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
      qm_pdb_hierarchy.write_pdb_file(file_name=str(i)+".pdb",
        crystal_symmetry=self.expand.cs_box)
      raw_records = qm_pdb_hierarchy.as_pdb_string(
        crystal_symmetry=self.expand.cs_box)
      ## cell size self.expand.cs_p1 from expand, seems not right
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

def get_qm_file_name_and_pdb_hierarchy(fragment_extracts, index):
  fragment_selection = fragment_extracts.fragment_super_selections[index]
  fragment_hierarchy = fragment_extracts.pdb_hierarchy_super.select(
    fragment_selection)
  sub_working_folder = fragment_extracts.working_folder + str(index) + "/"
  if (not os.path.isdir(sub_working_folder)):
    os.mkdir(sub_working_folder)
  qm_pdb_file = sub_working_folder + str(index) + ".pdb"
  ##for debugging
  if(0):
    fragment_hierarchy.write_pdb_file(
      file_name=qm_pdb_file,
      crystal_symmetry=fragment_extracts.expand.cs_box)
  complete_qm_pdb_file = qm_pdb_file[:-4] + "_capping.pdb"
  # _capping_pdb_filename(qm_pdb_file)
  ph = completion.run(pdb_hierarchy=fragment_hierarchy,
                      crystal_symmetry=fragment_extracts.expand.cs_box,
                      model_completion=False)
  return os.path.abspath(complete_qm_pdb_file),ph

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
    fragment_charges    = fragments.fragment_charges,
    fragment_selections = fragments.fragment_selections,
    fragment_super_selections=fragments.fragment_super_selections,
    working_folder      = fragments.working_folder,
    fragment_super_atoms      = fragments.fragment_super_atoms,
    cluster_atoms       = fragments.cluster_atoms,
    qm_engine_name      = fragments.qm_engine_name,
    charge_embedding    = fragments.charge_embedding,
    crystal_symmetry    = fragments.crystal_symmetry,
    pdb_hierarchy       = fragments.pdb_hierarchy,
    pdb_hierarchy_super = fragments.pdb_hierarchy_super,
    expand              = fragments.expand,
    buffer_selections   = fragments.buffer_selections)
