"""This module handles the weight factors and scaling.
   An adaptive restraints weight factor calculator is implemented, whereby the
   weight factor is doubled if a sufficiently large bond-RMSD is observed.
   Conversely, if a sufficiently small bond-RMSD is observed, then the weight
   factor is halved.
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import random, time
from cctbx import xray
from libtbx import adopt_init_args
from scitbx.array_family import flex
import cctbx.maptbx.real_space_refinement_simple
import scitbx.lbfgs
from libtbx import group_args
from . import refine
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.clashscore import clashscore
from libtbx.utils import null_out
from cctbx import maptbx

def compute_weight(fmodel, restraints_manager, shake_sites=False):
  if(shake_sites):
    random.seed(1)
    flex.set_random_seed(1)
  fmodel_dc = fmodel.deep_copy()
  xrs = fmodel_dc.xray_structure.deep_copy_scatterers()
  if(shake_sites):
    xrs.shake_sites_in_place(mean_distance=0.2)
  fmodel_dc.update_xray_structure(xray_structure=xrs, update_f_calc=True)
  x_target_functor = fmodel_dc.target_functor()
  tgx = x_target_functor(compute_gradients=True)
  gx = flex.vec3_double(tgx.\
    gradients_wrt_atomic_parameters(site=True).packed())
  tc, gc = restraints_manager.target_and_gradients(sites_cart=xrs.sites_cart())
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
  #
  if(y != 0.0): data_weight = x/y
  else:         data_weight = 1.0 # ad hoc default fallback
  return data_weight

class cctbx_geometry(object):
  def __init__(self, model):
    restraints_manager = model.get_restraints_manager()
    self.hd_sel = model.get_xray_structure().hd_selection()
    self.restraints_manager = restraints_manager.select(~self.hd_sel)

  def bond_rmsd(self, sites_cart):
    assert self.hd_sel.size() == sites_cart.size()
    energies_sites = self.restraints_manager.geometry.energies_sites(
      sites_cart        = sites_cart.select(~self.hd_sel),
      compute_gradients = False)
    return energies_sites.bond_deviations()[2]

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
        #update_f_mask  = True
        )
      #self.fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
    else:
      self.xray_structure.tidy_us()
      self.xray_structure.apply_symmetry_sites()

class sites_opt(object):
  """
  General calculator for model geometry optimization. For native CCTBX
  restraints, restraints_manager and model.restraints_manager are the same
  things.
  However, restraints_manager can be an external entity, such as coming from
  external packeges (eg., QM).
  Ideally, and this is probably a TODO item for the future, any
  restraints_manager should always be in the model.
  dump_gradients is used for debugging.
  """

  def __init__(self, model, max_shift, restraints_manager, shift_eval,
               dump_gradients=False, convergence_threshold=1.e-3,
               convergence_reached_times=3):
    self.model = model
    self.restraints_manager = restraints_manager
    self.dump_gradients = dump_gradients
    self.convergence_threshold = convergence_threshold
    self.convergence_reached_times = convergence_reached_times
    self.meat_convergence_criteria = 0
    self.x = flex.double(self.model.size()*3, 0)
    self.n = self.x.size()
    self.f = None
    self.g = None
    self.f_start = None
    self.max_shift_between_resets = 0
    self.sites_cart = self.model.get_sites_cart()
    self.sites_cart_start = self.sites_cart.deep_copy()
    self.bound_flags = flex.int(self.n, 2)
    self.lower_bound = flex.double([-1*max_shift]*self.n)
    self.upper_bound = flex.double([   max_shift]*self.n)
    self.shift_eval = shift_eval
    assert self.shift_eval in ["max", "mean"]
    if(self.shift_eval == "mean"): self.shift_eval_func = flex.mean
    else:                          self.shift_eval_func = flex.max
    self.total_time = 0
    self.number_of_target_and_gradients_calls = 0

  def target_and_gradients(self):
    self.number_of_target_and_gradients_calls+=1
    t0=time.time()
    sites_plus_x = self.sites_cart+flex.vec3_double(self.x)
    self.f, self.g = self.restraints_manager.target_and_gradients(
      sites_cart = sites_plus_x)
    self.g = self.g.as_double()
    # For tests
    if(self.dump_gradients):
      from libtbx import easy_pickle
      easy_pickle.dump(self.dump_gradients, self.g)
      STOP()
    #
    if(self.f_start is None):
      self.f_start = self.f

    self.max_shift_between_resets = self.shift_eval_func(flex.sqrt((
      self.sites_cart - sites_plus_x).dot()))
    self.total_time += (time.time()-t0)
    return self.f, self.g

  def compute_functional_and_gradients(self):
    return self.target_and_gradients()

  def apply_x(self):
    self.f_start = self.f
    self.model.set_sites_cart(
      sites_cart = self.sites_cart+flex.vec3_double(self.x))
    self.x = flex.double(self.model.size()*3, 0)
    self.sites_cart = self.model.get_sites_cart()
    if(self.max_shift_between_resets < self.convergence_threshold):
      self.meat_convergence_criteria += 1

  def converged(self):
    if(self.meat_convergence_criteria >= self.convergence_reached_times):
      return True
    return False

  def mean_shift_from_start(self):
    # Assumes apply_x has been called before so that self.sites_cart are updated
    return flex.mean(flex.sqrt((
      self.sites_cart - self.sites_cart_start).dot()))

  def __call__(self):
    f, g = self.target_and_gradients()
    return self.x, f, g

class sites(calculator):
  def __init__(self,
               fmodel=None,
               restraints_manager=None,
               dump_gradients=None):
    adopt_init_args(self, locals())
    self.x = None
    self.x_target_functor = None
    self.not_hd_selection = None # XXX UGLY

    self.sites_cart_start = None

    self._initialize(fmodel = self.fmodel)
    self.total_time = 0
    self.number_of_target_and_gradients_calls = 0
    self.data_weight = 1.
    self.restraints_weight_scale = 1.
    self.restraints_weight = 1.
    self.r_works = flex.double()
    self.n_converged = 0

  def _initialize(self, fmodel):
    self.sites_cart_start = fmodel.xray_structure.sites_cart().deep_copy()
    self.r_works = flex.double()
    self.n_converged = 0
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

  def reset_fmodel(self, fmodel):
    self._initialize(fmodel=fmodel)
    self.update_fmodel()

  def setw(self,
           data_weight = None,
           restraints_weight_scale = None,
           restraints_weight = None):
    if(data_weight is not None):
      self.data_weight = data_weight
    if(restraints_weight_scale is not None):
      self.restraints_weight_scale = restraints_weight_scale
    if(restraints_weight is not None):
      self.restraints_weight = restraints_weight

  def scale_restraints_weight_down(self, scale=2):
    self.setw(restraints_weight_scale = self.restraints_weight_scale/scale)

  def scale_restraints_weight_up(self, scale=2):
    self.setw(restraints_weight_scale = self.restraints_weight_scale*scale)

  def update(self, x):
    self.x = flex.vec3_double(x)
    self.fmodel.xray_structure.set_sites_cart(sites_cart = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    # This is to monitor convergence
    r_work = self.fmodel.r_work()
    self.r_works.append(r_work)

    diff = flex.sqrt((self.x - self.sites_cart_start).dot())
    shift_me = flex.mean(diff)
    shift_ma = flex.max(diff)
    print("%7.5f %7.3f %7.3f"%(r_work, shift_me, shift_ma), self.r_works.size())
    self.sites_cart_start = self.x.deep_copy()

    # Early termination
    #
    #if self.r_works.size()>20:
    #  a = self.r_works[-3:]
    #  aa = list(set([abs(i-j) for i in a for j in a if i != j]))
    #  if len(aa)>0:
    #    d = max(aa)
    #    d = True if d<0.0003 else False
    #    if d: self.n_converged += 1
    #    else:
    #      if self.n_converged>0: self.n_converged -= 1

  def converged(self):
    if self.n_converged >= 3: return True
    else:                   return False

  def target_and_gradients(self, x):
    self.number_of_target_and_gradients_calls+=1
    t0=time.time()
    self.update(x = x)
    rt, rg = self.restraints_manager.target_and_gradients(sites_cart = self.x)
    tgx = self.x_target_functor(compute_gradients=True)
    dt = tgx.target_work()
    dg = flex.vec3_double(tgx.\
      gradients_wrt_atomic_parameters(site=True).packed())
    t = dt*self.data_weight + \
      self.restraints_weight*rt*self.restraints_weight_scale
    g = dg*self.data_weight + \
      self.restraints_weight*rg*self.restraints_weight_scale
    if(self.dump_gradients is not None):
      from libtbx import easy_pickle
      easy_pickle.dump(self.dump_gradients+"_dg", dg.as_double())
      easy_pickle.dump(self.dump_gradients+"_rg", rg.as_double())
      easy_pickle.dump(self.dump_gradients+"_g", g.as_double())
      STOP()
    self.total_time += (time.time()-t0)
    return t, g.as_double()

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
    self.max_iterations = 150
    self.weight = data_weight
    self.sites_cart_start = self.model.get_xray_structure().sites_cart()
    self.show(model=self.model, weight = data_weight)
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

  def get_shift(self, other):
    s1 = self.sites_cart_start
    s2 = other.sites_cart()
    return flex.mean(flex.sqrt((s1 - s2).dot()))

  def get_shift2(self, m1, m2):
    s1 = m1.get_sites_cart()
    s2 = m2.get_sites_cart()
    return flex.mean(flex.sqrt((s1 - s2).dot()))

  def geometry_is_good(self, stats):
    b, a = stats.bond().mean, stats.angle().mean
    #return round(b, 3) <= self.max_bond_rmsd and round(a, 2) <= 1.5
    return b <= self.max_bond_rmsd and a <= 1.5

  def macro_cycle(self, weight):
    def stalled(x):
      for it in x:
        if x.count(it) > 3:
          return True
      return False
    up   = 0
    down = 0
    slow_down = 0
    previous_good = None
    bond_rmsds = []
    cntr = 0
    while True: # Endless loop!
      print("cycle:", cntr)
      cntr+=1
      stats = self.show(model = self.model, weight = weight, prefix="  start:")
      model = self.run_one(weight = weight)
      stats = self.show(model = model, weight = weight, prefix="  final:")
      b = stats.bond().mean
      bond_rmsds.append(round(b,3))
      if(stalled(bond_rmsds)):
        print("<<<<< weight optimization stalled: quitting >>>>>")
        break
      if(self.geometry_is_good(stats)):
        up += 1
        previous_good = weight
        if(b<0.03 or stalled(bond_rmsds)): weight = weight*2
        else:
          weight = weight + 0.3*weight
        self.model.set_sites_cart(sites_cart = model.get_sites_cart())
      else:
        down += 1
        if(b>0.03):
          weight = weight/2
        else:
          slow_down += 1
          if(slow_down<3):
            weight = weight - 0.3*weight
          else:
            weight = weight/2
            slow_down = 0
      if(up>0 and down>0):
        print("<<<<< weight optimization oscillates: quitting >>>>>")
        break
    ####
    print()
    if(previous_good is None):
      print("Using last weight to continue:", weight)
      previous_good = weight
    ####
    for it in [1,2,3,4,5]:
      stats = self.show(model = self.model, weight = previous_good, prefix="  start:")
      model = self.run_one(weight = previous_good)
      stats = self.show(model = model, weight = previous_good, prefix="  final:")
      if(self.geometry_is_good(stats)):
        shift = self.get_shift2(model, self.model)
        self.model.set_sites_cart(sites_cart = model.get_sites_cart())
        if(shift<0.01):
          print("<<<<< shift fell below 0.01: quitting >>>>>")
          break
      else:
        print("<<<<< geometry got worse: quitting >>>>>")
        break

  def show(self, model, weight, prefix=""):
    stats = model.geometry_statistics(use_hydrogens=False)
    s = stats.show_short()
    s = s.split()
    s = " ".join(s)
    dist = self.get_shift(other=model.get_xray_structure())
    if(weight is not None): w = "%8.4f"%weight
    else:                   w = "%5s"%str(None)
    cc_mask = refine.show_cc(
      map_data=self.map_data, xray_structure=model.get_xray_structure())
    print(prefix, "weight=%s"%w, s, "shift=%6.4f"%dist, \
      "cc_mask=%6.4f"%cc_mask)
    return stats

  def get_weight(self):
    o = maptbx.target_and_gradients_simple(
      unit_cell     = self.model.crystal_symmetry().unit_cell(),
      map_target    = self.map_data,
      sites_cart    = self.model.get_sites_cart(),
      selection     = flex.bool(self.model.size(), True),
      interpolation = "tricubic"
      )
    hd_sel = self.model.get_xray_structure().hd_selection()
    g_map = o.gradients().select(~hd_sel)
    _, g_geo = self.restraints_manager.target_and_gradients(
      sites_cart = self.model.get_sites_cart())
    g_geo = g_geo.select(~hd_sel)
    #
    g_map_d = flex.sqrt(g_map.dot())
    sel = g_map_d>flex.mean(g_map_d)*3
    y = g_map.select(~sel).norm()
    #-
    g_geo_d = flex.sqrt(g_geo.dot())
    sel = g_geo_d>flex.mean(g_geo_d)*3
    x = g_geo.select(~sel).norm()
    #
    #y = g_map.norm()
    #
    #x = g_geo.norm()
    ################
    if(y != 0.0): return x/y
    else:         return 1.0 # ad hoc default fallback

  def run(self):
    weight = self.get_weight()
    print("Initial weight estimate from ratio of grad norms:", weight)
    self.macro_cycle(weight = weight)
    return self.model

  def run_one(self, weight):
    model = self.model.deep_copy()
    xrs = model.get_xray_structure()
    uc = xrs.crystal_symmetry().unit_cell()
    not_hd_sel = ~model.selection(string = "element H or element D")
    refined = cctbx.maptbx.real_space_refinement_simple.lbfgs(
      unit_cell                       = uc,
      gradients_method                = "tricubic",
      sites_cart                      = xrs.sites_cart(),
      density_map                     = self.map_data,
      geometry_restraints_manager     = self.restraints_manager,
      selection_variable_real_space   = not_hd_sel,
      real_space_target_weight        = weight,
      real_space_gradients_delta      = 0.25,
      gradient_only                   = self.gradient_only,
      line_search                     = self.line_search,
      lbfgs_core_params               = self.lbfgs_core_params,
      lbfgs_termination_params        = self.lbfgs_termination_params,
      lbfgs_exception_handling_params = self.lbfgs_exception_handling_params)
    model.set_sites_cart(sites_cart=refined.sites_cart)
    return model
