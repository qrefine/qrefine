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
import cctbx.maptbx.real_space_refinement_simple
import scitbx.lbfgs
from libtbx import group_args
import qr
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.clashscore import clashscore
from libtbx.utils import null_out

def get_bonds_rmsd(restraints_manager, xrs):
  hd_sel = xrs.hd_selection()
  energies_sites = \
    restraints_manager.select(~hd_sel).energies_sites(
      sites_cart        = xrs.sites_cart().select(~hd_sel),
      compute_gradients = False)
  return energies_sites.bond_deviations()[2]

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
    self.r_frees = []
    self.r_works = []

  def scale_restraints_weight(self):
    if(self.weight_was_provided): return
    self.restraints_weight_scale *= 4.0

  def adjust_restraints_weight_scale(
        self,
        fmodel,
        geometry_rmsd_manager,
        max_bond_rmsd,
        scale):
    adjusted = None
    if(self.weight_was_provided): return adjusted
    rw = fmodel.r_work()
    rf = fmodel.r_free()
    cctbx_rm_bonds_rmsd = get_bonds_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = fmodel.xray_structure)
    ####
    adjusted = False
    if(cctbx_rm_bonds_rmsd>max_bond_rmsd):
      self.restraints_weight_scale *= scale
      adjusted = True
    if(not adjusted and rf<rw):
      self.restraints_weight_scale /= scale
      adjusted = True
    if(not adjusted and cctbx_rm_bonds_rmsd<max_bond_rmsd and rf>rw and
       abs(rf-rw)*100.<5.):
      self.restraints_weight_scale /= scale
      adjusted = True
    if(not adjusted and cctbx_rm_bonds_rmsd<max_bond_rmsd and rf>rw and
       abs(rf-rw)*100.>5.):
      self.restraints_weight_scale *= scale
      adjusted = True
    ####
    self.r_frees.append(round(rf,4))
    self.r_works.append(round(rw,4))
    return adjusted

  def add_restraints_weight_scale_to_restraints_weight_scales(self):
    if(self.weight_was_provided): return
    self.restraints_weight_scales.append(self.restraints_weight_scale)

  def compute_weight(self, fmodel, rm, verbose=False):
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
    if verbose: print '>>> gradient norms c,x %0.2f %0.2f' % (x, y)
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
    if verbose: print '>>> data_weight %0.2f' % self.data_weight

class calculator(object):
  def __init__(self,
               fmodel=None,
               xray_structure=None,
               restraints_weight_scale = 1.0):
    assert [fmodel, xray_structure].count(None)==1
    self.fmodel=None
    self.xray_structure=None
    if(fmodel is not None):
      self.fmodel = fmodel
    if(xray_structure is not None):
      self.xray_structure = xray_structure
    self.restraints_weight_scale = restraints_weight_scale

  def update_fmodel(self):
    if(self.fmodel is not None):
      self.fmodel.xray_structure.tidy_us()
      self.fmodel.xray_structure.apply_symmetry_sites()
      self.fmodel.update_xray_structure(
        xray_structure = self.fmodel.xray_structure,
        update_f_calc  = True,
        update_f_mask  = True)
      self.fmodel.update_all_scales(remove_outliers=False)
    else:
      self.xray_structure.tidy_us()
      self.xray_structure.apply_symmetry_sites()

class sites_opt(calculator):
  def __init__(self, restraints_manager, xray_structure, dump_gradients,
                     ase_atoms):
    self.dump_gradients = dump_gradients
    self.restraints_manager = restraints_manager
    self.x = None
    self.xray_structure = xray_structure
    self.not_hd_selection = None # XXX UGLY
    self.ase_atoms = ase_atoms
    self.initialize(xray_structure = self.xray_structure)

  def initialize(self, xray_structure=None):
    self.not_hd_selection = ~self.xray_structure.hd_selection() # XXX UGLY
    self.x = self.xray_structure.sites_cart().as_double()

  def update(self, x):
    self.x = x
    self.update_xray_structure()

  def target_and_gradients(self, x):
    self.update(x = x)
    f, g = self.restraints_manager.target_and_gradients(
      sites_cart = flex.vec3_double(self.x))
    if(self.dump_gradients is not None):
      from libtbx import easy_pickle
      easy_pickle.dump(self.dump_gradients, g)
      STOP()
    return f, g.as_double()

  def update_xray_structure(self):
    self.xray_structure.set_sites_cart(
      sites_cart = flex.vec3_double(self.x))

class sites(calculator):
  def __init__(self,
               fmodel=None,
               restraints_manager=None,
               weights=None,
               dump_gradients=None,
               ase_atoms=None):
    adopt_init_args(self, locals())
    self.x = None
    self.x_target_functor = None
    self.not_hd_selection = None # XXX UGLY
    self.initialize(fmodel = self.fmodel)

  def initialize(self, fmodel=None):
    self.not_hd_selection = ~self.fmodel.xray_structure.hd_selection() # XXX UGLY
    assert fmodel is not None
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    xray.set_scatterer_grad_flags(
      scatterers = self.fmodel.xray_structure.scatterers(),
      site       = True)
    self.x = self.fmodel.xray_structure.sites_cart().as_double()
    self.x_target_functor = self.fmodel.target_functor()

  def calculate_weight(self, verbose=False):
    self.weights.compute_weight(
      fmodel  = self.fmodel,
      rm      = self.restraints_manager,
      verbose = verbose)

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
    if(self.dump_gradients is not None):
      from libtbx import easy_pickle
      easy_pickle.dump(self.dump_gradients+"_dg", dg.as_double())
      easy_pickle.dump(self.dump_gradients+"_rg", rg.as_double())
      easy_pickle.dump(self.dump_gradients+"_g", g.as_double())
      STOP()
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

class sites_real_space(object):
  def __init__(self,
               model,
               geometry_rmsd_manager,
               max_bond_rmsd,
               stpmax,
               gradient_only,
               line_search,
               data_weight,
               refine_cycles,
               skip_weight_search,
               log,
               map_data=None,
               restraints_manager=None,
               max_iterations=None):
    adopt_init_args(self, locals())
    self.gradient_only = True
    self.max_iterations = 100
    self.weight = data_weight
    self.sites_cart_start = self.model.get_xray_structure().sites_cart()
    self.show(model=self.model)
    #
    self.rama_fav_best = None
    self.cbeta_best    = None
    self.rota_best     = None
    self.clash_best    = None
    #
    if(self.weight is None):
      self.weight = 1.
    self.refine_cycles = refine_cycles
    self.skip_weight_search = skip_weight_search
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = self.max_iterations)
    self.lbfgs_core_params = scitbx.lbfgs.core_parameters(
      stpmin = 1.e-9,
      stpmax = stpmax)
    self.lbfgs_exception_handling_params = scitbx.lbfgs.\
      exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = False,
        ignore_line_search_failed_step_at_upper_bound = False,
        ignore_line_search_failed_maxfev              = False)
    self.sites_cart_refined = None
    self.cctbx_rm_bonds_rmsd = get_bonds_rmsd(
      restraints_manager = self.geometry_rmsd_manager.geometry,
      xrs                = self.model.get_xray_structure())

  def get_shift(self, other):
    s1 = self.sites_cart_start
    s2 = other.sites_cart()
    return flex.mean(flex.sqrt((s1 - s2).dot()))

  def get_scores(self, model):
    rama_fav = ramalyze(
      pdb_hierarchy = model.get_hierarchy(),
      outliers_only = False).percent_favored
    cbeta = cbetadev(
      pdb_hierarchy = model.get_hierarchy(),
      outliers_only = True,
      out           = null_out()).get_outlier_percent()
    rota = rotalyze(
      pdb_hierarchy = model.get_hierarchy(),
      outliers_only = False).percent_outliers
    b_rmsd = get_bonds_rmsd(
      restraints_manager = self.geometry_rmsd_manager.geometry,
      xrs                = model.get_xray_structure())
    clash = clashscore(
    pdb_hierarchy = model.get_hierarchy(),
    keep_hydrogens = True).get_clashscore
    return group_args(
      rama_fav = rama_fav, cbeta = cbeta, rota = rota, b_rmsd = b_rmsd, clash = clash)

  def ready_to_stop(self, sc):
    return (sc.rama_fav < self.rama_fav_best and
            abs(sc.rama_fav-self.rama_fav_best)>1.) or \
           sc.cbeta > self.cbeta_best or \
           sc.rota > self.rota_best   or \
           sc.clash > self.clash_best

  def macro_cycle(self, weights):
    print "RSR: weights to try:", weights
    weight_best = None
    i_best = None
    model_best = None
    models = []
    for i, w in enumerate(weights):
      self.weight = w
      m = self.run_one()
      models.append(m.deep_copy())
      sc = self.get_scores(model = m)
      if(i==0 and self.rama_fav_best is None): # we assume best Rama favored with smallest weight
        self.rama_fav_best = sc.rama_fav
        self.cbeta_best    = sc.cbeta
        self.rota_best     = sc.rota
        self.clash_best    = sc.clash
      elif(i==0): # 2nd round: fine-tuning
        if(self.ready_to_stop(sc)):
          break
      if(sc.b_rmsd<self.max_bond_rmsd):
        weight_best = w
        i_best = i
        model_best = models[i_best]
      else:
        break
      #
      if(i>0):
        if(self.ready_to_stop(sc)):
          i_best = i-1
          weight_best = weights[i_best]
          model_best = models[i_best]
          break
      #
    print "RSR: weight_best:", weight_best
    return model_best, weight_best, i_best

  def run(self):
    weights = [0.1, 1.0, 10, 20, 30, 40, 50,  200]
    model, weight, i = self.macro_cycle(weights = weights)
    #
    if(weight==50.):
      new_weights = [50,60,70,80,90,100,110,120,130,140,150,160,170,180,190]
    elif(weight>1 and i!=len(weights)-1):
      new_weights = []
      w=weights[i]
      while w<weights[i+1]:
        w+=1
        new_weights.append(w)
    elif(weight == 1.0):
      new_weights = [1,2,3,4,5,6,7,8,9]
    elif(weight == 0.1):
      new_weights = [0.1,0.2,0.5,0.7]
    elif(weight == 0.01):
      new_weights = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]
    else:
      print "RSR: FALED TO FIND BEST WEIGHT"
      STOP()
    print "RSR: new_weights:", new_weights
    #
    model_, weight_, i_ = self.macro_cycle(weights = new_weights)
    self.model  = model
    self.weight = weight
    if(weight_ is not None):
      self.model  = model_
      self.weight = weight_
    #
    rmsd = get_bonds_rmsd(
      restraints_manager = self.geometry_rmsd_manager.geometry,
      xrs                = self.model.get_xray_structure())
    self.show(model=self.model, prefix="(start macro-cycles)")
    #
    for mc in [1,2,3,4,5]:
      self.model = self.run_one()
    #
    return self.model

  def show(self, model, prefix=""):
    s = model.geometry_statistics(use_hydrogens=False).show_short()
    s = s.split()
    s = " ".join(s)
    dist = self.get_shift(other=model.get_xray_structure())
    if(self.weight is not None): w = "%5.2f"%self.weight
    else:                        w = "%5s"%str(None)
    cc_mask = qr.show_cc(
      map_data=self.map_data, xray_structure=model.get_xray_structure())
    print "RSR", prefix, "weight=%s"%w, s, "shift=%6.4f"%dist, \
      "cc_mask=%6.4f"%cc_mask
    with open("weight_%s.pdb"%w.strip(), "w") as of:
      of.write(model.model_as_pdb())

  def run_one(self):
    model = self.model.deep_copy()
    xrs = model.get_xray_structure()
    uc = xrs.crystal_symmetry().unit_cell()
    refined = cctbx.maptbx.real_space_refinement_simple.lbfgs(
      unit_cell                       = uc,
      gradients_method                = "tricubic",
      sites_cart                      = xrs.sites_cart(),
      density_map                     = self.map_data,
      geometry_restraints_manager     = self.restraints_manager,
      real_space_target_weight        = self.weight,
      real_space_gradients_delta      = 0.25,
      gradient_only                   = self.gradient_only,
      line_search                     = self.line_search,
      lbfgs_core_params               = self.lbfgs_core_params,
      lbfgs_termination_params        = self.lbfgs_termination_params,
      lbfgs_exception_handling_params = self.lbfgs_exception_handling_params)
    model.set_sites_cart(sites_cart=refined.sites_cart)
    ####
    #rmsd = get_bonds_rmsd(
    #  restraints_manager = self.geometry_rmsd_manager.geometry,
    #  xrs                = model.get_xray_structure())
    self.show(model = model)
    return model
