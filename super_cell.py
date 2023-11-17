import iotbx.pdb
from scitbx.array_family import flex
from cctbx import uctbx
from cctbx import crystal
from libtbx.utils import null_out
import mmtbx
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from mmtbx import monomer_library
from libtbx import group_args
import math
import iotbx.pdb.utils
from cctbx.crystal import super_cell as cctbx_super_cell

mon_lib_srv = mmtbx.monomer_library.server.server()
ener_lib    = mmtbx.monomer_library.server.ener_lib(
  use_neutron_distances = True)

class expand(object):
  def __init__(self, pdb_hierarchy, crystal_symmetry, select_within_radius=15,
               create_restraints_manager=True):
    # fixed members
    self.create_restraints_manager = create_restraints_manager
    self.crystal_symmetry = crystal_symmetry
    self.select_within_radius = select_within_radius
    self.cs_p1 = crystal.symmetry(self.crystal_symmetry.unit_cell(), "P1")
    self.sel_str_focus = " or ".join([
      "chain %s"%c.id for c in pdb_hierarchy.chains()])
    # variable members (can change upon self.update() call) or by end of
    # constructor execution
    self.pdb_hierarchy = pdb_hierarchy
    self.super_sphere_geometry_restraints_manager = None
    self.ph_p1           = None
    self.selection_keep  = None
    self.ph_super_cell   = None
    self.cs_box          = None
    self.ph_super_sphere = None
    self.SCM = cctbx_super_cell.manager(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.select_within_radius)
    self.ph_super_sphere = self.SCM.super_sphere_hierarchy
    self.cs_box = self.SCM.cs_super_sphere
    #
    self.update()

  def update(self, sites_cart=None, redo_super_cell=False):
    if(sites_cart is not None):
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      if(not redo_super_cell):
        self.ph_super_sphere = self.SCM.update(sites_cart=sites_cart)
    if(redo_super_cell):
      self.SCM = cctbx_super_cell.manager(
        pdb_hierarchy        = self.pdb_hierarchy,
        crystal_symmetry     = self.crystal_symmetry,
        select_within_radius = self.select_within_radius)
      self.ph_super_sphere = self.SCM.super_sphere_hierarchy.deep_copy()
      self.cs_box = self.SCM.cs_super_sphere
    if(self.create_restraints_manager):
      if(self.super_sphere_geometry_restraints_manager is None or redo_super_cell):
        self.update_super_sphere_geometry_restraints_manager()
    return self

  def update_xyz(self, sites_cart=None):
    # keep the selection of super_sphere, update its coordinates
    self.update(sites_cart=sites_cart)
    return self

  def update_super_sphere_geometry_restraints_manager(self):
    # XXX Unify with process_model_file of qr.py
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
    geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      assume_hydrogens_all_missing = not has_hd,
      plain_pairs_radius           = 5.0)
    self.super_sphere_geometry_restraints_manager = mmtbx.restraints.manager(
       geometry = geometry, normalization = False)

  def write_super_cell_selected_in_sphere(self, file_name="super_sphere.pdb"):
    self.ph_super_sphere.write_pdb_file(file_name = file_name,
      crystal_symmetry = self.cs_box) # box around super-sphere, not original cs
                                      # or cs in p1.
