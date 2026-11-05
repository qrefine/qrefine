"""This module stores and handles all of the data
   needed for logging the results of the refinement.
"""
from __future__ import division
from __future__ import print_function
import os
from scitbx.array_family import flex

class manager(object):
  def __init__(self, model, geometry_rmsd_manager, log, output_mapper=None):
    """
    model holds the nuclear sites. output_mapper (None or an
    mmtbx.hydrogens.xh_mapper.mapper) is applied (shrink) to them when
    writing a PDB file (an identity mapper is not used); the model itself
    stays nuclear.
    """
    self.log = log
    self.output_mapper = output_mapper
    if(output_mapper is not None and output_mapper.is_identity):
      self.output_mapper = None
    self.geometry_rmsd_manager = geometry_rmsd_manager
    self.model = model#.deep_copy()
    self.sites_cart_start = self.model.get_sites_cart().deep_copy()
    self.sites_cart_prev  = self.sites_cart_start.deep_copy()
    self.sites_cart       = self.sites_cart_start.deep_copy()
    # Things we keep track of
    self.r_work = None
    self.r_free = None
    self.b_rmsd = self._get_b_rmsd()
    self.a_rmsd = self._get_a_rmsd()
    self.shift_start = 0
    self.shift_prev  = 0

  def _get_b_rmsd(self):
    return self.geometry_rmsd_manager.bond_rmsd(sites_cart = self.sites_cart)

  def _get_a_rmsd(self):
    return self.geometry_rmsd_manager.angle_rmsd(sites_cart = self.sites_cart)

  def update(self, fmodel=None, model=None, sites_cart=None):
    """
    Sites (nuclear) from sites_cart when given, else from model, else from
    fmodel; r_work and r_free from fmodel when given.
    """
    def _(x,y): return flex.mean(flex.sqrt((x - y).dot()))
    assert [fmodel, model, sites_cart].count(None) < 3
    if(sites_cart is not None):
      self.sites_cart = sites_cart.deep_copy()
    elif(model is not None):
      self.sites_cart = model.get_sites_cart()
    else:
      assert self.output_mapper is None, ("nuclear sites must be passed "
        "explicitly when an X-H mapper is active")
      self.sites_cart = fmodel.xray_structure.sites_cart()
    if(fmodel is not None):
      self.r_work = fmodel.r_work()
      self.r_free = fmodel.r_free()
    self.model.set_sites_cart(sites_cart = self.sites_cart)
    self.b_rmsd = self._get_b_rmsd()
    self.a_rmsd = self._get_a_rmsd()
    self.shift_start = _(self.sites_cart_start, self.sites_cart)
    self.shift_prev  = _(self.sites_cart_prev,  self.sites_cart)
    self.sites_cart_prev = self.sites_cart.deep_copy()

  def get_stats(self):
    return self.model.geometry_statistics(use_hydrogens=False)

  def g_info(self):
    stats = self.get_stats()
    s = stats.show_short()
    s = s.split()
    s = " ".join(s)
    s += " shift start/prev: %5.3f %5.3f"%(self.shift_start, self.shift_prev)
    return s

  def r_info(self):
    if self.r_work is None: return ""
    fmt="Rw: %6.4f Rf: %6.4f Rf-Rw: %6.4f"
    return fmt%(self.r_work, self.r_free, self.r_free-self.r_work)

  def show(self, prefix="", suffix=""):
    s = prefix + " " + self.r_info() + " " + self.g_info() + " " + suffix
    print(s, file = self.log)
    self.log.flush()

  def write_pdb_file(self, output_file_name, output_folder_name):
    """
    Write the model into output_folder_name (created if missing) as
    <stem>.pdb, or <stem>.cif when it does not fit PDB format, where stem is
    output_file_name without its extension. Returns the file name written;
    nothing written and None returned when output_folder_name is None.
    """
    if(output_folder_name is None): return None
    os.makedirs(output_folder_name, exist_ok=True)
    ext = ".pdb" if self.model.can_be_output_as_pdb() else ".cif"
    file_name = os.path.join(output_folder_name,
      os.path.splitext(output_file_name)[0] + ext)
    sites_cart = self.model.get_sites_cart()
    if(self.output_mapper is not None):
      self.model.set_sites_cart(self.output_mapper.shrink(sites_cart))
    try:
      with open(file_name,"w") as fo:
        fo.write(self.model.as_pdb_or_mmcif_string(target_format="pdb"))
    finally:
      if(self.output_mapper is not None):
        self.model.set_sites_cart(sites_cart)
    return file_name
