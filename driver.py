from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import pickle
import scitbx.lbfgs
from libtbx import easy_pickle
from libtbx import adopt_init_args
from cctbx import xray
from scitbx.array_family import flex
from . import calculator as calculator_module
from libtbx import group_args
from ase.optimize.lbfgs import LBFGS
import numpy
from libtbx.test_utils import approx_equal
from scitbx import minimizers
import sys, math
from libtbx.utils import Sorry

log = sys.stdout

class minimizer(object):
  def __init__(self,
        stpmax,
        log_switch,
        calculator,
        max_iterations,
        max_bond_rmsd,
        gradient_only,
        line_search,
        mode,
        geometry_rmsd_manager=None,
        preopt_params = None):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.x_previous = None
    self.number_of_function_and_gradients_evaluations=0
    self.lbfgs_core_params = scitbx.lbfgs.core_parameters(
      stpmin = 1.e-9,
      stpmax = stpmax)
    self.nstep = 0
    if preopt_params is not None:
       self.minimzer  = self.prcg_min (
       params = preopt_params
      )
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      gradient_only=gradient_only,
      line_search=line_search,
      log=log_switch,
      core_params=self.lbfgs_core_params,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations,
        min_iterations = max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_step_at_upper_bound=False,#True,
        ignore_line_search_failed_maxfev=True,
        ignore_line_search_failed_xtol=True,
        ignore_search_direction_not_descent=True
        )
        )

  def _get_bond_rmsd(self):
    b_mean = None
    if(self.geometry_rmsd_manager is not None):
      b_mean = self.geometry_rmsd_manager.bond_rmsd(
        sites_cart = flex.vec3_double(self.x))
    return b_mean

  # Enabling this makes it faster
  # TMP disabled to see the effect on refinement with ANI
  #
  #def callback_after_step(self, minimizer=None):
  #  if(self.geometry_rmsd_manager is not None and self.mode=="weight"):
  #    assert self.max_bond_rmsd is not None
  #    b_mean = self._get_bond_rmsd()
  #    print("    self.max_bond_rmsd", self.max_bond_rmsd, b_mean,
  #       self.number_of_function_and_gradients_evaluations)
  #    if(b_mean>self.max_bond_rmsd and
  #       self.number_of_function_and_gradients_evaluations-3>20):
  #      return True
  #  if(self.calculator.converged()): return True

  def compute_functional_and_gradients(self):
    self.number_of_function_and_gradients_evaluations += 1
    #print("  step: %3d bond rmsd: %8.6f"%(
    #  self.number_of_function_and_gradients_evaluations, self._get_bond_rmsd()),
    #  "%6.4f"%self.calculator.fmodel.r_work(),
    #  self.mode)
    return self.calculator.target_and_gradients(x = self.x)

  def prcg_min(self,params):
    import numpy as np
    import os
    """
    polak-ribiere conjugate gradient minimizer
    with simple step scaling and steepest decent steps
    """
    max_iter = params["maxiter"]
    max_step= params["stpmax"]
    switch_step= params["iswitch"] - 1
    gconv = params["gconv"]

    dim=len(self.x)
    x_new=self.x.as_numpy_array()
    x_old=x_new
    g_old=np.zeros(dim)
    gg=np.zeros(dim)
    xx=np.zeros(dim)
    step=np.zeros(dim)
    gg=np.zeros(dim)
    conv=False
    step_old=np.zeros(dim)

    self.x=flex.double(x_new.tolist())
    for iter in range(max_iter):
      self.eg=self.calculator.target_and_gradients(x = self.x )
      g=np.array(list(self.eg[1]))
      e=self.eg[0]
      gnorm=np.linalg.norm(g)
      # gnorm=self.eg[1].norm()/23.0605480121
      #self.number_of_function_and_gradients_evaluations += 1
      print('iter= %i E=%12.5E  G=%0.2f' % (iter+1,e,gnorm))
      #
      if gnorm <=gconv and iter >1:
        print('gnorm pre-convergenced!')
        self.x=flex.double(x_new.tolist())
        break
      #
      if iter<=switch_step:
        print('SD step')
        step=-g
      else:
        print('CG step')
        gg=g-g_old
        gdgg=np.vdot(g,gg)
        gdg=np.vdot(g_old,g_old)
        sdgg=np.vdot(step_old,gg)
        sds=np.vdot(step_old,step_old)
        alpha=sds/sdgg
        beta=gdgg/gdg
        step=-g + beta*step_old
        step*=alpha

      snorm=np.linalg.norm(step)
      if snorm>=max_step:
        step*=max_step/snorm
        # print 'step norm',snorm
      x_new =x_old + step
      #
      self.x=flex.double(x_new.tolist())
      e_old=e
      g_old=g
      x_old=x_new
      step_old=step

class clustering_update(object):
  def __init__(self, pre_sites_cart, log, rmsd_tolerance):
    self.pre_sites_cart = pre_sites_cart
    self.log = log
    self.rmsd_tolerance = rmsd_tolerance

  def re_clustering_check(self, sites_cart):
    rmsd_diff = pre_sites_cart.rms_difference(sites_cart)
    if(rmsd_diff < self.rmsd_tolerance):
      self.redo_clustering = False
    else:
      print(" rmsd_diff: ", rmsd_diff, "--> need to redo clustering", file=self.log)
      self.redo_clustering = True
      self.pre_sites_cart = sites_cart

  def re_clustering(self, calculator):
    try:    sites_cart = calculator.fmodel.xray_structure.sites_cart()
    except: sites_cart = calculator.model.get_sites_cart()
    rmsd_diff = self.pre_sites_cart.rms_difference(sites_cart)
    print(rmsd_diff, self.rmsd_tolerance)
    if(rmsd_diff > self.rmsd_tolerance):
      print(" rmsd_diff: ", rmsd_diff, "--> need to redo clustering", file=self.log)
      calculator.restraints_manager.fragment_manager.set_up_cluster_qm()
      self.pre_sites_cart = sites_cart

class minimizer_ase(object):
  def __init__(self, calculator, params, max_iterations, geometry_rmsd_manager):
    self.params = params
    self.max_iterations = max_iterations
    self.calculator = calculator
    self.geometry_rmsd_manager = geometry_rmsd_manager
    self.ase_atoms = calculator.ase_atoms
    self.ase_atoms.set_positions(flex.vec3_double(self.calculator.x))
    self.opt = LBFGS(atoms = self.ase_atoms)
    self.number_of_function_and_gradients_evaluations = 0
    self.b_rmsd = self._get_bond_rmsd(
      sites_cart = flex.vec3_double(self.calculator.x))
    print("  step: %3d bond rmsd: %8.6f"%(
      self.number_of_function_and_gradients_evaluations, self.b_rmsd))
    self.run(nstep = max_iterations)
    # Syncing and cross-checking begin
    e = 1.e-4
    assert approx_equal(
      self.ase_atoms.get_positions(), self.opt.atoms.get_positions(), e)
    self.calculator.update(x = self.opt.atoms.get_positions())
    assert approx_equal(self.calculator.x, self.opt.atoms.get_positions(), e)
    if(params.refine.mode=="refine"):
      assert approx_equal(flex.vec3_double(self.calculator.x),
        self.calculator.fmodel.xray_structure.sites_cart(), e)
    else:
      assert approx_equal(flex.vec3_double(self.calculator.x),
        self.calculator.xray_structure.sites_cart(), e)
    b_rmsd = self._get_bond_rmsd(sites_cart = flex.vec3_double(self.calculator.x))
    assert approx_equal(self.b_rmsd, b_rmsd, e)
    # Syncing and cross-checking end

  def _get_bond_rmsd(self, sites_cart):
    b_mean = None
    if(self.geometry_rmsd_manager is not None):
      energies_sites = \
        self.geometry_rmsd_manager.geometry.energies_sites(
          sites_cart        = sites_cart.select(s),
          compute_gradients = False)
      b_mean = energies_sites.bond_deviations()[2]
    return b_mean

  def step(self):
    sites_cart = flex.vec3_double(self.opt.atoms.get_positions())
    t,g = self.calculator.target_and_gradients(x = sites_cart)
    forces = numpy.array(g) * (-1)
    self.opt.step(forces)
    self.number_of_function_and_gradients_evaluations += 1
    self.calculator.update(x = self.opt.atoms.get_positions())
    #
    self.b_rmsd = self._get_bond_rmsd(
      sites_cart = flex.vec3_double(self.opt.atoms.get_positions()))
    print("  step: %3d bond rmsd: %8.6f"%(
      self.number_of_function_and_gradients_evaluations, self.b_rmsd))
    if(self.params.refine.mode=="refine" and
       self.b_rmsd>self.params.refine.max_bond_rmsd and
       self.number_of_function_and_gradients_evaluations>20):
      return False
    #
    return True

  def run(self, nstep):
    for i in range(nstep):
      v = self.step()
      if(not v): return

def run_minimize(calculator, params, geometry_rmsd_manager, mode):
  minimized = None
  if  (mode == "weight"):
    max_iterations = params.refine.max_iterations_weight
    max_bond_rmsd  = params.refine.max_bond_rmsd
  elif(mode == "refine"):
    max_iterations = params.refine.max_iterations_refine
    max_bond_rmsd  = None
  if(params.refine.use_ase_lbfgs):
    minimized = minimizer_ase(
      calculator            = calculator,
      params                = params,
      max_iterations        = max_iterations,
      geometry_rmsd_manager = geometry_rmsd_manager)
  else:
    log_switch = None
    preopt = None
    if (params.refine.pre_opt):
      preopt={
      'stpmax'  :params.refine.pre_opt_stpmax,
      'maxiter' :params.refine.pre_opt_iter,
      'iswitch' :params.refine.pre_opt_switch,
      'gconv'   :params.refine.pre_opt_gconv,
    }
    if (params.refine.opt_log or params.debug): log_switch=log
    if(max_iterations > 0):
      minimized = minimizer(
        log_switch            = log_switch,
        calculator            = calculator,
        stpmax                = params.refine.stpmax,
        gradient_only         = params.refine.gradient_only,
        line_search           = params.refine.line_search,
        max_iterations        = max_iterations,
        max_bond_rmsd         = max_bond_rmsd,
        mode                  = mode,
        geometry_rmsd_manager = geometry_rmsd_manager,
        preopt_params         = preopt)
  return minimized

def run_collect(n_fev, results, fmodel, geometry_rmsd_manager, calculator):
  cctbx_rm_bonds_rmsd = calculator_module.get_bonds_rmsd(
    restraints_manager = geometry_rmsd_manager.geometry,
    xrs                = fmodel.xray_structure)
  results.update(
    r_work                  = fmodel.r_work(),
    r_free                  = fmodel.r_free(),
    b                       = cctbx_rm_bonds_rmsd,
    xrs                     = fmodel.xray_structure,
    restraints_weight_scale = calculator.restraints_weight_scale,
    n_fev                   = n_fev)

def refine(fmodel,
           params,
           monitor,
           calculator,
           geometry_rmsd_manager,
           STOP_ONCE_REACHED=True):
  if(not params.refine.refine_sites): return
  # Ugly!
  try:    clustering = calculator.restraints_manager.clustering
  except: clustering = False
  if(clustering):
    cluster_qm_update = clustering_update(
      pre_sites_cart = calculator.fmodel.xray_structure.sites_cart(),
      log            = log,
      rmsd_tolerance = params.refine.rmsd_tolerance * 100,
      verbose        = params.debug)
  print("Start:", file=log)
  monitor.show(prefix="  ")

  if(not params.refine.skip_weight_search):
    print("Optimal weight search:", file=log)
    fmodel_copy = calculator.fmodel.deep_copy()
    data_weight = calculator_module.compute_weight(
      fmodel             = fmodel_copy,
      hdm                = calculator.hdm,
      restraints_manager = calculator.restraints_manager)
    calculator.setw(
      data_weight             = data_weight,
      restraints_weight_scale = 1.,
      restraints_weight       = 1.)
    # Weight control stuff
    r_free_best = fmodel_copy.r_free()
    up   = 0
    down = 0
    r_frees    = flex.double()
    b_rmsds    = flex.double()
    a_rmsds    = flex.double()
    sites_cart = []
    restraints_weight_scale = flex.double()
    print("Data weight (initial):", data_weight, file=log)
    #
    # Loop over weight search cycles
    for weight_cycle in range(params.refine.number_of_weight_search_cycles):
      fmodel = fmodel_copy.deep_copy() # Always use initial unchanged fmodel
      calculator.reset_fmodel(fmodel = fmodel)
      calculator.update_fmodel()
      if(clustering):
        cluster_qm_update.re_clustering(calculator)
      # Run minimization with given weight
      n_fev = 0
      for mc in range(1): # Just do once. Why one needs to do it more than once!
        minimized = run_minimize(
          calculator            = calculator,
          params                = params,
          geometry_rmsd_manager = geometry_rmsd_manager,
          mode                  = "weight")
        if(minimized is not None):
          calculator.reset_fmodel(fmodel = fmodel)
          calculator.update_fmodel()
          n_fev += minimized.number_of_function_and_gradients_evaluations
          break
      if(minimized is None): continue
      monitor.update(fmodel = fmodel)
      rws = calculator.restraints_weight_scale
      b_rmsds.append(round(monitor.b_rmsd,3))
      a_rmsds.append(round(monitor.a_rmsd,2))
      # Show
      msg = "rws(prev): %6.3f rws: %6.3f n_fev: %d"%(
        rws, calculator.restraints_weight_scale, n_fev)
      monitor.show(prefix="%d:"%weight_cycle, suffix = msg)
      # Sanity check:
      assert approx_equal(fmodel.r_work(), calculator.fmodel.r_work(), 1.e-4)
      #
      # DEFINE STOPPING RULE
      GOOD = monitor.b_rmsd <= params.refine.max_bond_rmsd and \
             monitor.a_rmsd <= params.refine.max_angle_rmsd
      #
      if(GOOD):
        down += 1
        restraints_weight_scale.append(rws)
        calculator.scale_restraints_weight_down(scale=1.5)
        r_frees   .append(monitor.r_work)
        sites_cart.append(monitor.sites_cart)
        if(params.refine.stop_one_found_first_good_weight):
          print("Optimal weight found. Stopping weight search.", file=log)
          break
      else:
        up += 1
        calculator.scale_restraints_weight_up(scale=1.5)
      # Show
      msg = "rws(prev): %6.3f rws: %6.3f n_fev: %d (NEXT)"%(
        rws, calculator.restraints_weight_scale, n_fev)
      monitor.show(prefix="%d:"%weight_cycle, suffix = msg)
      #
      # Ready to stop?
      if(up>0 and down>0):
        print("Flip happened. Stopping weight search.", file=log)
        break

      if(b_rmsds.size()>3):
        v1 = list(set(b_rmsds[-2:]))
        v2 = list(set(a_rmsds[-2:]))
        if(len(v1)==1 and len(v2)==1): # XXX See if this tighter!
          print("Bond rmsd stalled. Stopping weight search.", file=log)
          break
    # Ok, done with weights.
    if(r_frees.size()==0): # Fallback
      print("Weight search yields no result. Using last result.", file=log)
      print("Taking result with Rwork, Rfree: %6.4f %6.4f"%(
        fmodel.r_work(), fmodel.r_free()), file=log)
    else: # Choose best result
      s = flex.sort_permutation(r_frees)
      index = s[0]
      calculator.setw(restraints_weight_scale = restraints_weight_scale[index])
      print("Best Rfree from above: %6.4f"%r_frees[index], file=log)
      print("Best restraints scale:", round(restraints_weight_scale[index],3), file=log)
      xrs = fmodel.xray_structure
      xrs.set_sites_cart(sites_cart[index])
      fmodel.update_xray_structure(
        xray_structure = xrs,
        update_f_calc  = True,
        update_f_mask  = True)
    fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
    print("Best Rwork, Rfree (at refinement start): %6.4f %6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)
    monitor.update(fmodel = fmodel)
    print("At refinement start", file = log)
    monitor.show()
  #
  # Done with weights. Now let's refine!
  #
  print("Refinement:", file=log)
  calculator.reset_fmodel(fmodel = fmodel)
  calculator.update_fmodel()
  assert not calculator.converged()
  for refine_cycle in range(params.refine.number_of_refine_cycles):
    calculator.reset_fmodel(fmodel=fmodel)
    if(clustering):
      cluster_qm_update.re_clustering(calculator)
    #
    n_fev = 0
    for mc in range(params.refine.number_of_macro_cycles):
      minimized = run_minimize(
        calculator            = calculator,
        params                = params,
        geometry_rmsd_manager = geometry_rmsd_manager,
        mode                  = "refine")
      if(minimized is not None):
        calculator.reset_fmodel(fmodel = fmodel)
        calculator.update_fmodel()
        n_fev += minimized.number_of_function_and_gradients_evaluations
        break
    if(minimized is None): continue
    monitor.update(fmodel = fmodel)
    monitor.show(prefix="%d:"%refine_cycle, suffix="n_fev: %d"%n_fev)
    monitor.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(refine_cycle)+"_refine_cycle.pdb")
    if(calculator.converged()):
      print("Refinement converged. Stopping now.", file=log)
      break
  #print("At end of further refinement:", file=log)
  #print("calculator(refine), total_time (target_and_gradients)", calculator.total_time)
  #print("calculator(refine), number_of_target_and_gradients_calls (target_and_gradients)",
  #  calculator.number_of_target_and_gradients_calls)
  print("At refinement end:", file = log)
  monitor.show()
  assert approx_equal(fmodel.r_work(), calculator.fmodel.r_work(), 1.e-4)
  print("Best Rwork, Rfree (after refinement): %6.4f %6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)

def opt(model, params, monitor, calculator):
  assert model == calculator.model
  log_switch = None
  if (params.refine.opt_log or params.debug): log_switch=log
  monitor.show(prefix="start:")
  if(params.cluster.clustering):
    cluster_qm_update = clustering_update(
      model.get_sites_cart(), log,
      params.cluster.re_calculate_rmsd_tolerance)
    print("\ninteracting pairs number:  ",\
      len(calculator.restraints_manager.fragment_manager.interaction_list), file=log)
  for micro_cycle in range(params.refine.number_of_micro_cycles):
    assert model == calculator.model
    if(params.cluster.clustering):
      cluster_qm_update.re_clustering(calculator)
    if(params.refine.minimizer == "lbfgs"):
      minimized = minimizers.lbfgs(
        calculator     = calculator,
        max_iterations = params.refine.max_iterations_refine,
        gradient_only  = params.refine.gradient_only,
        stpmax         = params.refine.stpmax)
    else:
      assert not params.refine.gradient_only
      minimized = minimizers.lbfgsb(
        calculator     = calculator,
        max_iterations = params.refine.max_iterations_refine)
    minimized.show(log = log_switch, prefix="")
    calculator.apply_x()

    monitor.update(model = model)
    monitor.show(prefix="cycle %d: "%micro_cycle)
    monitor.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(micro_cycle)+"_opt_cycle.pdb")

    if(calculator.converged()):
      print("Convergence reached. Stopping now.", file=log)
      break

  print("calculator(opt), total_time (target_and_gradients)", calculator.total_time)
  print("calculator(opt), number_of_target_and_gradients_calls (target_and_gradients)",
    calculator.number_of_target_and_gradients_calls)
