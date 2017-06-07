"""This module handles the weight factors and scaling.
   An adaptive restraints weight factor calculator is implemented, whereby the
   weight factor is doubled if a sufficiently large bond-RMSD is observed.
   Conversely, if a sufficiently small bond-RMSD is observed, then the weight
   factor is halved.
"""
from __future__ import division
import random
from cctbx import xray
from libtbx import adopt_init_args
from scitbx.array_family import flex

def get_bonds_angles_rmsd(restraints_manager, xrs):
  hd_sel = xrs.hd_selection()
  energies_sites = \
    restraints_manager.select(~hd_sel).energies_sites(
      sites_cart        = xrs.sites_cart().select(~hd_sel),
      compute_gradients = False)
  b_mean = energies_sites.bond_deviations()[2]
  return b_mean

class weights(object):
  def __init__(self,
               shake_sites             = True,
               restraints_weight       = None,
               data_weight             = None,
               restraints_weight_scale = 1.0):
    adopt_init_args(self, locals())
    if(self.data_weight is not None):
      self.weight_was_provided = True
    else:
      self.weight_was_provided = False
    self.restraints_weight_scales = flex.double([self.restraints_weight_scale])

  def scale_restraints_weight(self):
    if(self.weight_was_provided): return
    self.restraints_weight_scale *= 4.0

  def adjust_restraints_weight_scale(
        self,
        fmodel,
        geometry_rmsd_manager,
        max_bond_rmsd):
    adjusted = None
    if(self.weight_was_provided): return adjusted
    rw = fmodel.r_work()
    rf = fmodel.r_free()
    cctbx_rm_bonds_rmsd = get_bonds_angles_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = fmodel.xray_structure)
    r_work_and_r_free_are_ok = (rf>rw and abs(rf-rw)*100.<5.)
    if(cctbx_rm_bonds_rmsd<0.01 and r_work_and_r_free_are_ok):
      self.restraints_weight_scale /= 2.
      adjusted=True
    if(not r_work_and_r_free_are_ok):
      self.restraints_weight_scale *= 2.
      adjusted=True
    if(not adjusted and cctbx_rm_bonds_rmsd>max_bond_rmsd):
      self.restraints_weight_scale *= 2.
      adjusted=True
    return adjusted

  def add_restraints_weight_scale_to_restraints_weight_scales(self):
    if(self.weight_was_provided): return
    self.restraints_weight_scales.append(self.restraints_weight_scale)

  def compute_weight(self, fmodel, rm):
    if(self.weight_was_provided): return
    random.seed(1)
    flex.set_random_seed(1)
    #
    fmodel_dc = fmodel.deep_copy()
    xrs = fmodel_dc.xray_structure.deep_copy_scatterers()
    if(self.shake_sites):
      xrs.shake_sites_in_place(mean_distance=0.2)
    fmodel_dc.update_xray_structure(xray_structure=xrs, update_f_calc=True)
    x_target_functor = fmodel_dc.target_functor()
    tgx = x_target_functor(compute_gradients=True)
    gx = flex.vec3_double(tgx.\
            gradients_wrt_atomic_parameters(site=True).packed())
    tc, gc = rm.target_and_gradients(sites_cart=xrs.sites_cart())
    x = gc.norm()
    y = gx.norm()
    # filter out large contributions
    gx_d = flex.sqrt(gx.dot())
    sel = gx_d>flex.mean(gx_d)*6
    y = gx.select(~sel).norm()
    #
    gc_d = flex.sqrt(gc.dot())
    sel = gc_d>flex.mean(gc_d)*6
    x = gc.select(~sel).norm()
    ################
    if(y != 0.0): self.data_weight = x/y
    else:         self.data_weight = 1.0 # ad hoc default fallback

class calculator(object):
  def __init__(self,
               fmodel,
               restraints_weight_scale = 1.0):
    self.fmodel = fmodel
    self.restraints_weight_scale = restraints_weight_scale

  #def update_restraints_weight_scale(self, restraints_weight_scale):
  #  self.restraints_weight_scale = restraints_weight_scale

  def update_fmodel(self):
    self.fmodel.xray_structure.tidy_us()
    self.fmodel.xray_structure.apply_symmetry_sites()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True,
      update_f_mask  = True)
    self.fmodel.update_all_scales(remove_outliers=False)

  def target_and_gradients(self):
    f, g = self.target_and_gradients(x = self.x)
    return f, g.as_double()

class sites_opt(calculator):
  def __init__(self, restraints_manager, fmodel):
    self.restraints_manager = restraints_manager
    self.x = None
    self.fmodel = fmodel
    self.initialize(fmodel = self.fmodel)

  def initialize(self, fmodel=None):
    self.x = self.fmodel.xray_structure.sites_cart().as_double()

  def update(self, x):
    self.x = x

  def target_and_gradients(self, x):
    self.update(x = x)
    f, g = self.restraints_manager.target_and_gradients(
      sites_cart = flex.vec3_double(self.x))
    return f, g.as_double()

  def update_fmodel_opt(self):
    self.fmodel.xray_structure.set_sites_cart(
      sites_cart = flex.vec3_double(self.x))
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

class sites(calculator):
  def __init__(self,
               fmodel=None,
               restraints_manager=None,
               weights=None):
    adopt_init_args(self, locals())
    self.x = None
    self.x_target_functor = None
    self.initialize(fmodel = self.fmodel)

  def initialize(self, fmodel=None):
    assert fmodel is not None
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    xray.set_scatterer_grad_flags(
      scatterers = self.fmodel.xray_structure.scatterers(),
      site       = True)
    self.x = self.fmodel.xray_structure.sites_cart().as_double()
    self.x_target_functor = self.fmodel.target_functor()

  def calculate_weight(self):
    self.weights.compute_weight(
      fmodel = self.fmodel,
      rm     = self.restraints_manager)

  def reset_fmodel(self, fmodel=None):
    if(fmodel is not None):
      self.initialize(fmodel=fmodel)
      self.fmodel = fmodel
      self.update_fmodel()

  def update_restraints_weight_scale(self, restraints_weight_scale):
    self.weights.restraints_weight_scale = restraints_weight_scale

  def update(self, x):
    self.x = flex.vec3_double(x)
    self.fmodel.xray_structure.set_sites_cart(sites_cart = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def target_and_gradients(self, x):
    self.update(x = x)
    rt, rg = self.restraints_manager.target_and_gradients(sites_cart = self.x)
    tgx = self.x_target_functor(compute_gradients=True)
    dt = tgx.target_work()
    dg = flex.vec3_double(tgx.\
      gradients_wrt_atomic_parameters(site=True).packed())
    t = dt*self.weights.data_weight + \
      self.weights.restraints_weight*rt*self.weights.restraints_weight_scale
    g = dg*self.weights.data_weight + \
      self.weights.restraints_weight*rg*self.weights.restraints_weight_scale
    return t, g.as_double()

class adp(calculator):
  def __init__(self,
               fmodel=None,
               restraints_manager=None,
               restraints_weight=None,
               data_weight=None,
               restraints_weight_scale=None):
    adopt_init_args(self, locals())
    self.x = None
    self.x_target_functor = None
    self.initialize(fmodel = self.fmodel)

  def initialize(self, fmodel=None):
    assert fmodel is not None
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    assert self.fmodel.xray_structure.scatterers().size() == \
      self.fmodel.xray_structure.use_u_iso().count(True)
    sel = flex.bool(
      self.fmodel.xray_structure.scatterers().size(), True).iselection()
    self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(iselection=sel)
    self.x = fmodel.xray_structure.extract_u_iso_or_u_equiv()
    self.x_target_functor = self.fmodel.target_functor()

  def calculate_weight(self):
    raise RuntimeError("Not implemented.")
    self.data_weight = compute_weight(
      fmodel = self.fmodel,
      rm     = self.restraints_manager)

  def reset_fmodel(self, fmodel=None):
    if(fmodel is not None):
      self.initialize(fmodel=fmodel)
      self.fmodel = fmodel

  def update(self, x):
    self.x = x
    self.fmodel.xray_structure.set_u_iso(values = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def target_and_gradients(self, x):
    self.update(x = x)
    tgx = self.x_target_functor(compute_gradients=True)
    f = tgx.target_work()
    g = tgx.gradients_wrt_atomic_parameters(u_iso=True)
    return f, g
