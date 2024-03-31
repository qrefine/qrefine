"""This module stores and handles all of the data
   needed for logging the results of the refinement.
"""
from __future__ import division
from __future__ import print_function
import os,sys
import mmtbx.utils
from scitbx.array_family import flex

class manager(object):
  def __init__(self, model, geometry_rmsd_manager, log):
    self.log = log
    self.geometry_rmsd_manager = geometry_rmsd_manager
    self.model = model.deep_copy()
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

  def update(self, fmodel=None, model=None):
    def _(x,y): return flex.mean(flex.sqrt((x - y).dot()))
    if fmodel is not None:
      self.sites_cart = fmodel.xray_structure.sites_cart()
      self.r_work = fmodel.r_work()
      self.r_free = fmodel.r_free()
    elif model is not None:
      self.sites_cart = model.get_sites_cart()
    else:
      assert [fmodel, model].count(None) == 1
    self.model.set_sites_cart(sites_cart = self.sites_cart)
    self.b_rmsd = self._get_b_rmsd()
    self.a_rmsd = self._get_a_rmsd()
    self.shift_start = _(self.sites_cart_start, self.sites_cart)
    self.shift_prev  = _(self.sites_cart_prev,  self.sites_cart)
    self.sites_cart_prev = self.sites_cart.deep_copy()

  def g_info(self):
    stats = self.model.geometry_statistics(use_hydrogens=False)
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

  def write_final_pdb_files(self, output_file_name, output_folder_name):
    if(os.path.exists(output_folder_name) is False):
      os.mkdir(output_folder_name)
    output_file_name = os.path.basename(output_file_name)
    self.write_pdb_file(output_file_name = output_file_name,
      output_folder_name = output_folder_name)

  def write_pdb_file(self, output_file_name, output_folder_name):
    file_name = "%s" % os.path.join(output_folder_name, output_file_name)
    with open(file_name,"w") as fo:
      fo.write(self.model.model_as_pdb())

  def finalize(self,
               input_file_name_prefix,
               output_file_name_prefix,
               output_folder_name):
    if(output_file_name_prefix is not None):
      file_name = "%s_refined.pdb"%output_file_name_prefix
    else:
      file_name = "%s_refined.pdb"%input_file_name_prefix
    self.write_final_pdb_files(
      output_file_name   = file_name,
      output_folder_name = output_folder_name)
    print("See %s in %s folder."%(file_name, output_folder_name), file=self.log)
