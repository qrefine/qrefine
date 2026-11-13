from libtbx.utils import null_out
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from cctbx.crystal import super_cell as cctbx_super_cell
from qrefine.utils.hierarchy_utils import write_hierarchy, \
  peptide_closures, bond_edits_params

class expand(object):
  """
  Super sphere of pdb_hierarchy: the model and all its symmetry copies
  within select_within_radius (cctbx.crystal.super_cell.manager). The
  sphere starts with the master copy; update(sites_cart) moves the master
  and the copies. The cctbx restraints manager of the sphere is built on
  first use of super_sphere_geometry_restraints_manager.
  """
  def __init__(self, pdb_hierarchy, crystal_symmetry, select_within_radius=15):
    self.crystal_symmetry = crystal_symmetry
    self.select_within_radius = select_within_radius
    self.pdb_hierarchy = pdb_hierarchy
    self._super_sphere_geometry_restraints_manager = None
    self.SCM = cctbx_super_cell.manager(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.select_within_radius)
    self.ph_super_sphere = self.SCM.super_sphere_hierarchy
    self.cs_box = self.SCM.cs_super_sphere

  def update(self, sites_cart=None):
    if(sites_cart is not None):
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      self.ph_super_sphere = self.SCM.update(sites_cart=sites_cart)
    return self

  def update_xyz(self, sites_cart=None):
    # keep the selection of super_sphere, update its coordinates
    return self.update(sites_cart=sites_cart)

  @property
  def super_sphere_geometry_restraints_manager(self):
    if(self._super_sphere_geometry_restraints_manager is None):
      self._super_sphere_geometry_restraints_manager = \
        self._build_super_sphere_geometry_restraints_manager()
    return self._super_sphere_geometry_restraints_manager

  def _build_super_sphere_geometry_restraints_manager(self):
    # Same interpretation settings as qr.py uses for the model.
    mon_lib_srv = mmtbx.monomer_library.server.server()
    ener_lib    = mmtbx.monomer_library.server.ener_lib(
      use_neutron_distances = True)
    params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
    params.use_neutron_distances = True
    params.restraints_library.cdl = False
    params.sort_atoms = False
    processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv              = mon_lib_srv,
      ener_lib                 = ener_lib,
      params                   = params,
      pdb_hierarchy            = self.ph_super_sphere,
      strict_conflict_handling = False,
      crystal_symmetry         = self.cs_box,
      force_symmetry           = True,
      log                      = null_out())
    xrs = processed_pdb_file.xray_structure()
    sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
    # peptide bonds between the master and a copy or between copies (a
    # macrocycle closed through a symmetry mate): pdb_interpretation does
    # not link two chains
    closures = peptide_closures(self.ph_super_sphere)
    params_edits = None
    if(len(closures) > 0):
      params_edits = bond_edits_params(closures).geometry_restraints.edits
    # no plain_pairs_radius: the plain pair table serves ADP restraints only,
    # and was the largest part of every fragment's selected restraints
    geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      assume_hydrogens_all_missing = not has_hd,
      params_edits                 = params_edits)
    return mmtbx.restraints.manager(geometry = geometry, normalization = False)

  def write_super_cell_selected_in_sphere(self, file_name):
    """
    Write the super sphere into file_name (.cif instead of .pdb when it does
    not fit PDB format); nothing if it is None. Returns the name written.
    """
    if(file_name is None): return None
    # box around super-sphere, not original cs or cs in p1.
    return write_hierarchy(self.ph_super_sphere, file_name,
      crystal_symmetry = self.cs_box)
