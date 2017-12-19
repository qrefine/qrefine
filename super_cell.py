import iotbx.pdb
from scitbx.array_family import flex
import string
from cctbx import uctbx
from cctbx import crystal
from libtbx.utils import null_out
import mmtbx
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from mmtbx import monomer_library
import iotbx.pdb.utils
from libtbx import group_args
import mmtbx.model

mon_lib_srv = mmtbx.monomer_library.server.server()
ener_lib    = mmtbx.monomer_library.server.ener_lib()

def create_super_sphere(pdb_hierarchy, crystal_symmetry, select_within_radius):
  pdb_str = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
  pdb_inp = iotbx.pdb.input(
    source_info = "pdb_hierarchy",
    lines       = flex.split_lines(pdb_str))
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.nonbonded_distance_cutoff = select_within_radius
  m = mmtbx.model.manager(
    model_input               = pdb_inp,
    pdb_interpretation_params = p,
    log                       = null_out())
  grm = m.get_restraints_manager()
  c = iotbx.pdb.hierarchy.chain(id = "SS")
  r = grm.geometry.pair_proxies().nonbonded_proxies.\
    get_symmetry_interacting_indices_unique(
      sites_cart = pdb_hierarchy.atoms().extract_xyz())
  fm = crystal_symmetry.unit_cell().fractionalization_matrix()
  om = crystal_symmetry.unit_cell().orthogonalization_matrix()
  selection = r.keys()
  for rg in pdb_hierarchy.residue_groups():
    keep=False
    for i in rg.atoms().extract_i_seq():
      if(i in selection):
        keep=True
        break
    if(keep):
      ops = r[i]
      for op in ops:
        rg_ = rg.detached_copy()
        xyz = rg_.atoms().extract_xyz()
        new_xyz = flex.vec3_double()
        for xyz_ in xyz:
          t1 = fm*flex.vec3_double([xyz_])
          t2 = op*t1[0]
          t3 = om*flex.vec3_double([t2])
          new_xyz.append(t3[0])
        rg_.atoms().set_xyz(new_xyz)
        rg_.link_to_previous=True
        c.append_residue_group(rg_)
  resseq = 0
  for rg in c.residue_groups():
    rg.resseq = "%4d" % resseq
    resseq+=1
  super_sphere_hierarchy = pdb_hierarchy.deep_copy()
  assert len(super_sphere_hierarchy.models())==1
  super_sphere_hierarchy.models()[0].append_chain(c)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = super_sphere_hierarchy.atoms().extract_xyz(),
    buffer_layer = 10)
  cs_super_sphere = box.crystal_symmetry()
  return group_args(
    hierarchy         = super_sphere_hierarchy,
    crystal_symmetry = cs_super_sphere)

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
    #
    self.update()

  def update(self, sites_cart=None):
    if(sites_cart is not None):
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    o = create_super_sphere(
      pdb_hierarchy        = self.pdb_hierarchy,
      crystal_symmetry     = self.crystal_symmetry,
      select_within_radius = self.select_within_radius)
    self.ph_super_sphere = o.hierarchy
    self.ph_super_cell = self.ph_super_sphere.deep_copy() # nonsense, we don't have super-cell anymore
    self.cs_box = o.crystal_symmetry
    if(self.create_restraints_manager):
      self.update_super_sphere_geometry_restraints_manager()
    return self

  def update_xyz(self, sites_cart=None):
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
      crystal_symmetry = self.cs_p1)
