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
from itertools import product
import string

mon_lib_srv = mmtbx.monomer_library.server.server()
ener_lib    = mmtbx.monomer_library.server.ener_lib()

def all_chain_ids(): # XXX This may go to cctbx/iotbx/pdb.
  chars = string.ascii_letters+string.digits
  def permutations(iterable, r=None):
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    for indices in product(range(n), repeat=r):
      if(len(set(indices)) == r):
        yield tuple(pool[i] for i in indices)
  result = permutations(iterable = chars, r = 2)
  return [c for c in chars]+["".join(p) for p in result]

class expand(object):
  def __init__(self, pdb_hierarchy, crystal_symmetry, select_within_radius=15):
    # fixed members
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
    if(self.ph_p1 is None or sites_cart is not None):
      self.ph_p1 = self.pdb_hierarchy.expand_to_p1(
        crystal_symmetry = self.crystal_symmetry)
    self._build_super_cell()
    # Box is needed for selection criteria to not select periodically
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart   = self.ph_super_cell.atoms().extract_xyz(),
      buffer_layer = 10)
    self.cs_box = box.crystal_symmetry()
    self.xrs_super = self.ph_super_cell.extract_xray_structure(
      crystal_symmetry = self.cs_box)
    #
    self._select_within()
    self.ph_super_sphere = self.ph_super_cell.select(self.selection_keep)
    self.update_super_sphere_geometry_restraints_manager()
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
      pdb_inp                  = self.ph_super_sphere.as_pdb_input(), # XXX Does this loose precision?
      strict_conflict_handling = False,
      crystal_symmetry         = self.cs_box, #XXX I hope this is correct
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

  def write_super_cell(self, file_name="super.pdb"):
    self.ph_super_cell.write_pdb_file(file_name=file_name)

  def write_p1(self, file_name="p1.pdb"):
    self.ph_p1.write_pdb_file(
      file_name        = file_name,
      crystal_symmetry = self.cs_p1)

  def _build_super_cell(self):
    om = self.cs_p1.unit_cell().orthogonalization_matrix()
    fm = self.cs_p1.unit_cell().fractionalization_matrix()
    #
    r = iotbx.pdb.hierarchy.root()
    m = iotbx.pdb.hierarchy.model()
    chains = list(self.ph_p1.chains())
    #
    original_chain_ids = set([c.id for c in self.pdb_hierarchy.chains()])
    all_ci = set(all_chain_ids())
    available_ci = list(original_chain_ids.symmetric_difference(all_ci))
    taken_ci = []
    # make sure original chains come first
    for chain in chains:
      m.append_chain(chain.detached_copy())
    #
    for x in [-1, 0, 1]:
      for y in [-1, 0, 1]:
        for z in [-1, 0, 1]:
          if([x,y,z]==[0,0,0]): continue # skip original chains
          for chain in chains:
            c_ = chain.detached_copy()
            sites_frac = fm*(c_.atoms().extract_xyz())
            sites_cart = om*(sites_frac+[x,y,z])
            c_.atoms().set_xyz(sites_cart)
            new_ci = None
            for ci in available_ci:
              if(not ci in taken_ci):
                new_ci = ci
                taken_ci.append(new_ci)
                break
            assert new_ci is not None
            c_.id = new_ci
            m.append_chain(c_)
    r.append_model(m)
    self.ph_super_cell = r
    self.ph_super_cell.atoms().reset_i_seq()

  def _select_within(self):
    #
    sel_model_orig = self.ph_super_cell.atom_selection_cache().selection(
      string = self.sel_str_focus)
    selection = self.xrs_super.selection_within(radius=self.select_within_radius,
      selection=sel_model_orig)
    #
    self.selection_keep = flex.size_t()
    for rg in self.ph_super_cell.residue_groups():
      keep=True
      sel = rg.atoms().extract_i_seq()
      for i in sel:
        if(not selection[i]):
          keep=False
          break
      if(keep): self.selection_keep.extend(sel)
