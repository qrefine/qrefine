"""This module stores and handles all of the data
   needed for logging the results of the refinement.
"""
from __future__ import division
import os
import mmtbx.utils
from scitbx.array_family import flex

def selxrs(xrss, s):
  tmp = []
  for i,s_ in enumerate(s):
    if(s_):
      tmp.append(xrss[i].deep_copy_scatterers())
  return tmp[:]

class manager(object):
  def __init__(self, r_work, r_free, b, xrs, max_bond_rmsd,
                     restraints_weight_scale, max_r_work_r_free_gap,
                     pdb_hierarchy, mode, log):
    self.mode = mode
    self.log = log
    self.pdb_hierarchy = pdb_hierarchy
    self.states = mmtbx.utils.states(
      pdb_hierarchy  = self.pdb_hierarchy,
      xray_structure = xrs)
    self.crystal_symmetry = xrs.crystal_symmetry()
    self.states.add(sites_cart = xrs.sites_cart())
    self.max_bond_rmsd = max_bond_rmsd
    self.max_r_work_r_free_gap = max_r_work_r_free_gap
    self.r_works = flex.double()
    self.r_frees = flex.double()
    self.bs   = flex.double()
    self.xrss = []
    self.r_works.append(r_work)
    self.r_frees.append(r_free)
    self.bs.append(b)
    self.xrss.append(xrs.deep_copy_scatterers())
    self.restraints_weight_scales = flex.double([restraints_weight_scale])
    self.n_fev = 0

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

  def choose_best(self):
    # Do not inlude initial model in decision-making.
    rfs  = self.r_frees[1:]
    rws  = self.r_works[1:]
    bs   = self.bs    [1:]
    gaps = rfs - rws
    xrss = self.xrss  [1:]
    rewescas = self.restraints_weight_scales[1:]
    # Select all that satisfy bonds and Rfree-Rwork gap criteria
    s  = bs < self.max_bond_rmsd
    rfs  = rfs .select(s)
    rws  = rws .select(s)
    bs   = bs  .select(s)
    gaps = gaps.select(s)
    rewescas = rewescas.select(s)
    xrss = selxrs(xrss=xrss, s=s)
    # Rfree-Rwork gap
    filtered_by_gap=False
    s  = gaps>0
    s &= flex.abs(gaps)*100.<self.max_r_work_r_free_gap
    if(s.count(True)>0):
      rfs      = rfs .select(s)
      rws      = rws .select(s)
      bs       = bs  .select(s)
      gaps     = gaps.select(s)
      rewescas = rewescas.select(s)
      xrss     = selxrs(xrss=xrss, s=s)
      filtered_by_gap=True
    if(rfs.size()==0):
      return None, None, None, None
    else:
      # Choose the one that has lowest Rfree
      min_r = flex.min(rfs)
      min_gap = flex.min(gaps)
      index_best = None
      if(filtered_by_gap):
        for i in xrange(rfs.size()):
          if(abs(rfs[i]-min_r)<1.e-5):
            index_best = i
            break
      else:
        for i in xrange(gaps.size()):
          if(abs(gaps[i]-min_gap)<1.e-5):
            index_best = i
            break
      # This is the result
      self.pdb_hierarchy.adopt_xray_structure(xrss[index_best])
      return xrss[index_best], rws[index_best], rfs[index_best],\
        rewescas[index_best]

  def choose_last(self):
    xrs_best = self.xrss[len(self.xrss)-1].deep_copy_scatterers()
    self.pdb_hierarchy.adopt_xray_structure(xrs_best)

  def show(self, prefix):
    fmt="%s %3d Rw: %6.4f Rf: %6.4f Rf-Rw: %6.4f rmsd(b): %7.4f rws: %4.1f n_fev: %d"
    i = self.r_works.size()-1
    print >> self.log, fmt%(prefix, i, self.r_works[-1], self.r_frees[-1],
      self.r_frees[-1]-self.r_works[-1], self.bs[-1],
      self.restraints_weight_scales[-1], self.n_fev)
    self.log.flush()

  def write_final_pdb_files(self, output_file_name, output_folder_name):
    if(os.path.exists(output_folder_name) is False):
      os.mkdir(output_folder_name)
    self.pdb_hierarchy.write_pdb_file(
      file_name        = "%s/%s"%(output_folder_name, output_file_name),
      crystal_symmetry = self.crystal_symmetry)
    self.states.write(
      file_name = "%s/%s"%(
        output_folder_name, output_file_name+"_all_states.pdb"))

  def write_pdb_file(self, output_file_name, output_folder_name):
    self.pdb_hierarchy.write_pdb_file(
      file_name = "%s/%s"%(output_folder_name, output_file_name))

  def finalize(self,
               input_file_name_prefix,
               output_file_name_prefix,
               output_folder_name):
    xrs_best = None
    if(self.mode == "refine"):
      xrs_best, r_work, r_free, dummy = self.choose_best()
      if(xrs_best is not None):
        print >> self.log, "Best r_work: %6.4f r_free: %6.4f"%(r_work, r_free)
      else:
        print >> self.log, " r_factor (best): None"
        print >> self.log, " take the last structure"
        self.show(prefix="")
        xrs_best = self.xrss[0]
    if(self.mode == "opt"):
      self.choose_last()
    if(output_file_name_prefix is not None):
      file_name = "%s_refined.pdb"%output_file_name_prefix
    else:
      file_name = "%s_refined.pdb"%input_file_name_prefix
    self.write_final_pdb_files(
      output_file_name   = file_name,
      output_folder_name = output_folder_name)
    print >> self.log, "See %s in %s folder."%(file_name, output_folder_name)
    return xrs_best
