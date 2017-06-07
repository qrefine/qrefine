import iotbx.pdb
from scitbx.array_family import flex
import string
from cctbx import uctbx
from cctbx import crystal

# This must be a method of iotbx.pdb.hierarchy and used in phenix.models_as_chains
def models_as_chains(h, model_id_orig, original_chain_ids):
  r = iotbx.pdb.hierarchy.root()
  m = iotbx.pdb.hierarchy.model()
  idl = [i for i in string.ascii_lowercase]
  idu = [i for i in string.ascii_uppercase]
  taken = original_chain_ids[:]
  c1 = None
  c2 = None
  n_atoms = []
  original_chains = []
  neibor_chains = []
  for model_id, m_ in enumerate(list(h.models())):
    if(model_id !=model_id_orig):
      for c_ in m_.chains():
        n_at = len(c_.atoms())
        if(not n_at in n_atoms): n_atoms.append(n_at)
        c_ = c_.detached_copy()
        found = False
        for idu_ in idu:
          for idl_ in idl:
            id_ = idu_+idl_
            if(not id_ in taken):
              taken.append(id_)
              found = id_
              break
          if(found): break
        c_.id = found
        neibor_chains.append(c_)
        #m.append_chain(c_)
    else:
      for c_ in m_.chains():
        original_chains.append(c_.detached_copy())
        #m.append_chain(c_.detached_copy())
  ## put original pdb in the beginning and neibours at the end
  for chain in original_chains:
    m.append_chain(chain)
  for chain in neibor_chains:
    m.append_chain(chain)
  r.append_model(m)
  return r

class expand(object):
  def __init__(self, pdb_hierarchy, crystal_symmetry, select_within_radius=15):
    self.pdb_hierarchy = pdb_hierarchy
    self.select_within_radius = select_within_radius
    self.crystal_symmetry = crystal_symmetry
    self.cs_p1 = crystal.symmetry(self.crystal_symmetry.unit_cell(), "P1")
    self.ph_p1 = self.pdb_hierarchy.expand_to_p1(
      crystal_symmetry=self.crystal_symmetry)
    # original chain IDs
    self.original_chain_ids = []
    for c in self.pdb_hierarchy.chains():
      self.original_chain_ids.append(c.id)
    self.sel_str_focus = " or ".join([
      "chain %s"%c for c in self.original_chain_ids])
    #
    self.ph_super  = None
    self.ph_chains = None
    self._build_super_cell()
    # Box is needed for selection criteria to not select periodically
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart   = self.ph_chains.atoms().extract_xyz(),
      buffer_layer = 10)
    self.cs_box = box.crystal_symmetry()
    self.xrs_super = self.ph_chains.extract_xray_structure(
                              crystal_symmetry = self.cs_box)
    #
    self._select_within()

  def write_super_cell_selected_in_sphere(self, file_name="super_sphere.pdb"):
    self.ph_chains.select(self.selection_keep).write_pdb_file(
      file_name=file_name, crystal_symmetry = self.cs_p1)

  def write_super_cell(self, file_name="super.pdb"):
    self.ph_chains.write_pdb_file(
      file_name        = file_name,
      crystal_symmetry = self.cs_p1)

  def write_p1(self, file_name="p1.pdb"):
    self.ph_p1.write_pdb_file(file_name=file_name,
      crystal_symmetry=self.crystal_symmetry)

  def update(self, sites_cart):
    self.pdb_hierarchy.atoms().set_xyz(sites_cart)
    self.ph_p1 = self.pdb_hierarchy.expand_to_p1(
      crystal_symmetry=self.crystal_symmetry)
    self._build_super_cell()
    return self.ph_chains.select(self.selection_keep)

  def _build_super_cell(self):
    # super-cell (models)
    self.ph_super = iotbx.pdb.hierarchy.root()
    # INEFFICIENT XXXX
    xrs_p1 = self.ph_p1.extract_xray_structure(crystal_symmetry = self.cs_p1)
    sites_frac = xrs_p1.sites_frac()
    #
    cntr=0
    model_id_orig = None
    for x in [-1, 0, 1]:
      for y in [-1, 0, 1]:
        for z in [-1, 0, 1]:
          if([x,y,z]==[0,0,0]): model_id_orig = cntr
          ph_ = self.ph_p1.deep_copy()
          #
          tmp = sites_frac+[x,y,z]
          xrs_p1.set_sites_frac(tmp)
          tmp = xrs_p1.sites_cart()
          #
          ph_.atoms().set_xyz(tmp)
          models = ph_.models()
          md = models[0].detached_copy()
          md.id = str(cntr)
          self.ph_super.append_model(md)
          cntr+=1
    # super-cell (chains)
    self.ph_chains = models_as_chains(
      h                  = self.ph_super,
      model_id_orig      = model_id_orig,
      original_chain_ids = self.original_chain_ids)
    self.ph_chains.atoms().reset_i_seq()

  def _select_within(self):
    #
    sel_model_orig = self.ph_chains.atom_selection_cache().selection(
      string = self.sel_str_focus)
    selection = self.xrs_super.selection_within(radius=self.select_within_radius,
      selection=sel_model_orig)
    #
    self.selection_keep = flex.size_t()
    for rg in self.ph_chains.residue_groups():
      keep=True
      sel = rg.atoms().extract_i_seq()
      for i in sel:
        if(not selection[i]):
          keep=False
          break
      if(keep): self.selection_keep.extend(sel)

def exercise():
  pdb_inp = iotbx.pdb.input(file_name="2ona_complete_occH0_cctbx_refined.pdb")
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()

  o = expand(pdb_hierarchy = ph, crystal_symmetry = cs)
  o.write_super_cell()
  o.write_p1()
  o.write_super_cell_selected_in_sphere()
  sites_cart = ph.atoms().extract_xyz()
  o.update(sites_cart = sites_cart)
  o.write_super_cell()

if __name__ == "__main__":
  exercise()
