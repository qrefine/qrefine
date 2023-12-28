"""This module stores and handles all of the data
   needed for logging the results of the refinement.
"""
from __future__ import division
from __future__ import print_function
import os,sys
import mmtbx.utils
from scitbx.array_family import flex

log = sys.stdout

def selxrs(xrss, s):
  tmp = []
  for i,s_ in enumerate(s):
    if(s_):
      tmp.append(xrss[i].deep_copy_scatterers())
  return tmp[:]

class manager(object):
  def __init__(self, model, max_bond_rmsd,
                     restraints_weight_scale, max_r_work_r_free_gap,
                     mode, r_work=None, r_free=None):
    assert [r_work, r_free].count(None) in [0,2]
    xrs = model.get_xray_structure()
    self.mode = mode
    self.pdb_hierarchy = model.get_hierarchy()
    self.states = mmtbx.utils.states(pdb_hierarchy  = self.pdb_hierarchy)
    self.crystal_symmetry = xrs.crystal_symmetry()
    self.states.add(sites_cart = xrs.sites_cart())
    self.max_bond_rmsd = max_bond_rmsd
    self.max_r_work_r_free_gap = max_r_work_r_free_gap
    # Collectables
    self.r_works = flex.double()
    self.r_frees = flex.double()
    self.bs   = flex.double()
    self.xrss = []
    if([r_work, r_free].count(None) == 0):
      self.r_works.append(r_work)
      self.r_frees.append(r_free)
    self.bs.append(model.get_bonds_rmsd())
    self.xrss.append(xrs.deep_copy_scatterers())
    self.restraints_weight_scales = flex.double([restraints_weight_scale])
    self.n_fev = 0

  def reset_custom(self):
    # the first model is the initial input model
    self.r_works = self.r_works[0:1]
    self.r_frees = self.r_frees[0:1]
    self.bs   = self.bs[0:1]
    self.xrss = self.xrss[0:1]
    self.restraints_weight_scales = self.restraints_weight_scales[0:1]

  def update(self, r_work=None, r_free=None, b=None, xrs=None,
             restraints_weight_scale=None, n_fev=None):
    if(r_work is not None): self.r_works.append(r_work)
    if(r_free is not None): self.r_frees.append(r_free)
    if(b is not None): self.bs.append(b)
    if(xrs is not None):
      self.xrss.append(xrs.deep_copy_scatterers())
      self.states.add(sites_cart = xrs.sites_cart())
      self.pdb_hierarchy.adopt_xray_structure(xrs)
    if(restraints_weight_scale is not None):
      self.restraints_weight_scales.append(restraints_weight_scale)
    if(n_fev is not None):
      self.n_fev += n_fev

  def choose_last(self):
    xrs_best = self.xrss[len(self.xrss)-1].deep_copy_scatterers()
    self.pdb_hierarchy.adopt_xray_structure(xrs_best)

  def show(self, prefix):
    if(self.r_works.size()>0):
      fmt="%s %3d Rw: %6.4f Rf: %6.4f Rf-Rw: %6.4f rmsd(b): %7.4f rws: %6.3f n_fev: %d"
      i = self.r_works.size()-1
      print(fmt%(prefix, i, self.r_works[-1], self.r_frees[-1],
        self.r_frees[-1]-self.r_works[-1], self.bs[-1],
        self.restraints_weight_scales[-1], self.n_fev), file=log)
    else:
      fmt="%s %3d rmsd(b): %7.4f rws: %6.3f n_fev: %d"
      i = self.bs.size()-1
      print(fmt%(prefix, i, self.bs[-1],
        self.restraints_weight_scales[-1], self.n_fev), file=log)
    log.flush()

  def write_final_pdb_files(self, output_file_name, output_folder_name):
    if(os.path.exists(output_folder_name) is False):
      os.mkdir(output_folder_name)
    output_file_name = os.path.basename(output_file_name)
    self.pdb_hierarchy.write_pdb_file(
      file_name        = "%s" % os.path.join(output_folder_name, output_file_name),
      crystal_symmetry = self.crystal_symmetry)
    self.states.write(
      file_name = "%s"%os.path.join(output_folder_name,
                                    output_file_name+"_all_states.pdb"))

  def write_pdb_file(self, output_file_name, output_folder_name):
    self.pdb_hierarchy.write_pdb_file(
      file_name = "%s" % os.path.join(output_folder_name, output_file_name),
      crystal_symmetry = self.crystal_symmetry,
      )

  def finalize(self,
               input_file_name_prefix,
               output_file_name_prefix,
               output_folder_name,
               use_r_work):
    self.choose_last()
    if(output_file_name_prefix is not None):
      file_name = "%s_refined.pdb"%output_file_name_prefix
    else:
      file_name = "%s_refined.pdb"%input_file_name_prefix
    self.write_final_pdb_files(
      output_file_name   = file_name,
      output_folder_name = output_folder_name)
    print("See %s in %s folder."%(file_name, output_folder_name), file=log)
