from __future__ import division
import os
import pickle
import scitbx.lbfgs
from libtbx import easy_pickle
from libtbx import adopt_init_args
from cctbx import xray
from scitbx.array_family import flex
import calculator as calculator_module
from libtbx.utils import Sorry
from libtbx import group_args
from ase.optimize.lbfgs import LBFGS
import numpy
from libtbx.test_utils import approx_equal

class convergence(object):
  def __init__(self, params, fmodel=None, xray_structure=None):
    self.r_start=None
    if(fmodel is not None):
      self.r_start = fmodel.r_work()
      self.sites_cart_start = fmodel.xray_structure.sites_cart()
    else:
      self.sites_cart_start = xray_structure.sites_cart()
    self.r_tolerance=params.refine.r_tolerance
    self.max_bond_rmsd=params.refine.max_bond_rmsd
    self.rmsd_tolerance=params.refine.rmsd_tolerance
    self.use_convergence_test = params.refine.use_convergence_test
    self.number_of_convergence_occurances=0
    #
    self.rws = flex.double()
    self.rfs = flex.double()
    self.gaps = flex.double()
    self.restraints_weight_scales = flex.double()

  def is_converged(self, fmodel, bond_rmsd=None, restraints_weight_scale=None):
    #
    if not self.use_convergence_test: return False
    #
    rw = fmodel.r_work()
    rf = fmodel.r_free()
    gap = rf-rw
    self.rws                     .append(rw)
    self.rfs                     .append(rf)
    self.gaps                    .append(gap)
    if(restraints_weight_scale is not None):
      self.restraints_weight_scales.append(restraints_weight_scale)
    #
    if(restraints_weight_scale is not None):
      rwc = self.restraints_weight_scales
      i_last = self.rws.size()-1
      if(i_last>3):
        rwc_3 = rwc[i_last]
        rwc_2 = rwc[i_last-1]
        rwc_1 = rwc[i_last-2]
        rws123 = [rwc_1, rwc_2, rwc_3]
        for rwc_i in rws123:
          if(rws123.count(rwc_i)>1):
            return True
    #
    sites_cart = fmodel.xray_structure.sites_cart()
    r_diff=abs(self.r_start-rw)
    rmsd_diff=self.sites_cart_start.rms_difference(sites_cart)
    self.sites_cart_start = sites_cart
    self.r_start=rw
    if(r_diff<self.r_tolerance and rmsd_diff<self.rmsd_tolerance and
       (bond_rmsd is not None and bond_rmsd<self.max_bond_rmsd)):
      self.number_of_convergence_occurances+=2
    if(self.number_of_convergence_occurances==2 or rw<0.005):
      return True and self.use_convergence_test
    else: return False and self.use_convergence_test

  def is_geometry_converged(self, sites_cart):
    if not self.use_convergence_test: return False
    rmsd_diff=self.sites_cart_start.rms_difference(sites_cart)
    if(rmsd_diff<self.rmsd_tolerance):
      return True

  def geometry_exploded(self, fmodel, geometry_rmsd_manager):
    result = False
    cctbx_rm_bonds_rmsd = calculator_module.get_bonds_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = fmodel.xray_structure)
    if(cctbx_rm_bonds_rmsd>self.max_bond_rmsd*2.0):
      result = True
    return result

  def is_weight_scale_converged(self, restraints_weight_scale):
    return restraints_weight_scale in self.restraints_weight_scales

class minimizer(object):
  def __init__(self,
        stpmax,
        log_switch,
        calculator,
        max_iterations,
        max_bond_rmsd,
        gradient_only,
        line_search,
        results,
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
        max_iterations=max_iterations),
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
      s = self.calculator.not_hd_selection
      energies_sites = \
        self.geometry_rmsd_manager.geometry.select(s).energies_sites(
          sites_cart        = flex.vec3_double(self.x).select(s),
          compute_gradients = False)
      b_mean = energies_sites.bond_deviations()[2]
    return b_mean

  def callback_after_step(self, minimizer=None):
    if(self.geometry_rmsd_manager is not None or self.mode=="refine"):
      b_mean = self._get_bond_rmsd()
      if(self.mode=="refine" and
         b_mean>self.max_bond_rmsd and
         self.number_of_function_and_gradients_evaluations-3>20):
        return True

  def compute_functional_and_gradients(self):
    self.number_of_function_and_gradients_evaluations += 1
    # Ad hoc damping shifts; note arbitrary 1.0 below
    #x_current = self.x
    #if(self.x_previous is None):
    #  self.x_previous = x_current.deep_copy()
    #else:
    #  xray.ext.damp_shifts(previous=self.x_previous, current=x_current,
    #    max_value=1.0)
    #  self.x_previous = x_current.deep_copy()
    #
    print "  step: %3d bond rmsd: %8.6f"%(
      self.number_of_function_and_gradients_evaluations, self._get_bond_rmsd())
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
      print 'iter= %i E=%12.5E  G=%0.2f' % (iter+1,e,gnorm)
      #
      if gnorm <=gconv and iter >1:
        print 'gnorm pre-convergenced!'
        self.x=flex.double(x_new.tolist())
        break
      #
      if iter<=switch_step:
        print 'SD step'
        step=-g
      else:
        print 'CG step'
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
      print >> self.log, " rmsd_diff: ", rmsd_diff, "--> need to redo clustering"
      self.redo_clustering = True
      self.pre_sites_cart = sites_cart

  def re_clustering(self, calculator):
    sites_cart = calculator.fmodel.xray_structure.sites_cart()
    rmsd_diff = self.pre_sites_cart.rms_difference(sites_cart)
    if(rmsd_diff > self.rmsd_tolerance):
      print >> self.log, " rmsd_diff: ", rmsd_diff, "--> need to redo clustering"
      calculator.restraints_manager.fragments.set_up_cluster_qm()
      print >> self.log, " interacting pairs number:  ", \
        calculator.restraints_manager.fragments.interacting_pairs
      self.pre_sites_cart = sites_cart

class restart_data(object):
  def __init__(self, geometry_rmsd_manager, fmodel=None, xray_structure=None):
    assert [xray_structure, fmodel].count(None) == 1
    rst_data = {}
    if(fmodel is not None): rst_data["fmodel"] = fmodel
    else:                   rst_data["xrs"] = xray_structure
    rst_data["geometry_rmsd_manager"] = geometry_rmsd_manager
    self.rst_data = rst_data

  def write_rst_file(self, rst_file, weight_cycle = None, refine_cycle = None,
                     micro_cycle = None, fmodel = None, weights = None,
                     conv_test = None, results = None, xray_structure = None):
    self.rst_data["weight_cycle"] = weight_cycle
    self.rst_data["refine_cycle"] = refine_cycle
    self.rst_data["micro_cycle"] = micro_cycle
    self.rst_data["rst_fmodel"] = fmodel
    self.rst_data["rst_xray_structure"] = xray_structure
    self.rst_data["weights"] = weights
    self.rst_data["conv_test"] = conv_test
    self.rst_data["results"] = results
    easy_pickle.dump(file_name=rst_file, obj=self.rst_data)

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
    print "  step: %3d bond rmsd: %8.6f"%(
      self.number_of_function_and_gradients_evaluations, self.b_rmsd)
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
      s = self.calculator.not_hd_selection
      energies_sites = \
        self.geometry_rmsd_manager.geometry.select(s).energies_sites(
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
    print "  step: %3d bond rmsd: %8.6f"%(
      self.number_of_function_and_gradients_evaluations, self.b_rmsd)
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

def run_minimize(calculator, params, results, geometry_rmsd_manager, mode):
  assert mode in ["weight", "refine"]
  result = None
  try:
    result = run_minimize_(
      calculator            = calculator,
      params                = params,
      results               = results,
      geometry_rmsd_manager = geometry_rmsd_manager,
      mode                  = mode)
  except Exception as e:
    print "minimization failed:", e
    result = None
  return result

def run_minimize_(calculator, params, results, geometry_rmsd_manager, mode):
  minimized = None
  if  (mode == "weight"): max_iterations = params.refine.max_iterations_weight
  elif(mode == "refine"): max_iterations = params.refine.max_iterations_refine
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
    if (params.refine.opt_log or params.debug): log_switch=results.log
    if(max_iterations > 0):
      minimized = minimizer(
        log_switch            = log_switch,
        calculator            = calculator,
        stpmax                = params.refine.stpmax,
        gradient_only         = params.refine.gradient_only,
        line_search           = params.refine.line_search,
        max_iterations        = max_iterations,
        max_bond_rmsd         = params.refine.max_bond_rmsd,
        results               = results,
        mode                  = params.refine.mode,
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
    restraints_weight_scale = calculator.weights.restraints_weight_scale,
    n_fev                   = n_fev)

def refine(fmodel,
           params,
           results,
           calculator,
           geometry_rmsd_manager):
  if(not params.refine.refine_sites): return
  rst_file = params.rst_file
  rst_data = restart_data(fmodel=fmodel,
    geometry_rmsd_manager=geometry_rmsd_manager)
  if(os.path.isfile(rst_file)):
    with open(rst_file, 'rb') as handle:
      rst_file_data = pickle.load(handle)
      weight_cycle_start = rst_file_data["weight_cycle"]
      refine_cycle_start = rst_file_data["refine_cycle"]
      print >> results.log
      print >> results.log, "*"*50
      print >> results.log, "restarts from weight_cycle: %d, refine_cycle: %s"%(
        weight_cycle_start, refine_cycle_start)
      print >> results.log, "*"*50
      print >> results.log
  else:
    weight_cycle_start = 1
    refine_cycle_start = None
  if(weight_cycle_start==1):
    conv_test = convergence(fmodel = fmodel, params = params)
  else:
    conv_test = rst_file_data["conv_test"]
  try:
    clustering = calculator.restraints_manager.clustering
  except :
    clustering = False
  if(clustering):
    cluster_qm_update = clustering_update(
      pre_sites_cart = calculator.fmodel.xray_structure.sites_cart(),
      log            = results.log,
      rmsd_tolerance = params.refine.rmsd_tolerance * 100,
      verbose=params.debug,
      )
    print >> results.log, "\ninteracting pairs number:  ", \
      calculator.restraints_manager.fragments.interacting_pairs
  weight_cycle = weight_cycle_start
  print >> results.log, "Start:"
  results.show(prefix="  ")
  if(refine_cycle_start is not None):
    print >> results.log, \
     "Found optimal weight. Proceed to further refinement with this weight."
    fmodel = calculator.fmodel.deep_copy()
  elif(not params.refine.skip_initial_weight_optimization):
    print >> results.log, "Optimal weight search:"
    fmodel_copy = calculator.fmodel.deep_copy()
    for weight_cycle in xrange(weight_cycle_start,
                               params.refine.number_of_weight_search_cycles+1):
      if((weight_cycle!=1 and weight_cycle==weight_cycle_start)):
        fmodel = calculator.fmodel.deep_copy()
        if params.debug: print '>>> Using calculator fmodel'
      else:
        fmodel = fmodel_copy.deep_copy()
        if params.debug: print '>>> Using fmodel_copy fmodel'
      calculator.reset_fmodel(fmodel = fmodel)
      if(clustering):
        cluster_qm_update.re_clustering(calculator)
      # Calculate weight
      calculator.calculate_weight(verbose=params.debug)
      # Collect state
      rst_data.write_rst_file(
        rst_file     = rst_file,
        refine_cycle = None,
        weight_cycle = weight_cycle,
        fmodel       = fmodel,
        weights      = calculator.weights,
        conv_test    = conv_test,
        results      = results)
      # Run minimization
      n_fev = 0
      for mc in xrange(params.refine.number_of_macro_cycles):
        minimized = run_minimize(
          calculator            = calculator,
          params                = params,
          results               = results,
          geometry_rmsd_manager = geometry_rmsd_manager,
          mode                  = "weight")
        if(minimized is not None):
          calculator.reset_fmodel(fmodel = fmodel)
          calculator.update_fmodel()
          n_fev += minimized.number_of_function_and_gradients_evaluations
          break
      if(minimized is None): continue
      # collect
      run_collect(
        n_fev                 = n_fev,
        results               = results,
        fmodel                = fmodel,
        calculator            = calculator,
        geometry_rmsd_manager = geometry_rmsd_manager)
      # Dump model at this step
      results.write_pdb_file(
        output_folder_name = params.output_folder_name,
        output_file_name   = str(weight_cycle)+"_weight_cycle.pdb")
      # show this step
      results.show(prefix="  ")
      # Converged?
      is_converged = conv_test.is_converged(
        bond_rmsd = results.bs[len(results.bs)-1],
        fmodel                  = fmodel,
        restraints_weight_scale = calculator.weights.restraints_weight_scale)
      if(is_converged):
        print >> results.log, "Converged (model)."
        break
      calculator.weights.adjust_restraints_weight_scale(
        fmodel                = fmodel,
        geometry_rmsd_manager = geometry_rmsd_manager,
        max_bond_rmsd         = params.refine.max_bond_rmsd,
        scale                 = params.refine.adjust_restraints_weight_scale_value)
      # Converged?
      is_converged = conv_test.is_weight_scale_converged(
        restraints_weight_scale = calculator.weights.restraints_weight_scale)
      if(is_converged):
        print >> results.log, "Converged (weight scale)."
        break
      calculator.weights.\
          add_restraints_weight_scale_to_restraints_weight_scales()
    print >> results.log, "At end of weight search:"
    results.show(prefix="  ")
    #
    # Done with weight search. Now get best result and refine it further.
    #
    xrs_best, dummy, dummy, wsc_best = results.choose_best(
      use_r_work = params.refine.choose_best_use_r_work)
    fmodel.update_xray_structure(
      xray_structure = xrs_best,
      update_f_calc  = True,
      update_f_mask  = True)
    fmodel.update_all_scales(remove_outliers=False)
    calculator.update_restraints_weight_scale(restraints_weight_scale=wsc_best)
    ####
    results.reset_custom()
    run_collect(
      n_fev                 = 0,
      results               = results,
      fmodel                = fmodel,
      calculator            = calculator,
      geometry_rmsd_manager = geometry_rmsd_manager)
    ####
    # show best
    print >> results.log, "At start of further refinement:"
    results.show(prefix="  ")
    print >> results.log, "Start further refinement:"
    refine_cycle_start = 1
  if(refine_cycle_start is None): refine_cycle_start=1
  #
  if(params.refine.skip_initial_weight_optimization):
    calculator.calculate_weight(verbose=params.debug)
  #
  for refine_cycle in xrange(refine_cycle_start,
                      params.refine.number_of_refine_cycles+refine_cycle_start):
    calculator.reset_fmodel(fmodel=fmodel)
    if(clustering):
      cluster_qm_update.re_clustering(calculator)
    rst_data.write_rst_file(
      rst_file     = rst_file,
      refine_cycle = refine_cycle,
      weight_cycle = weight_cycle,
      fmodel       = fmodel,
      weights      = calculator.weights,
      conv_test    = conv_test,
      results      = results)
    #
    n_fev = 0
    for mc in xrange(params.refine.number_of_macro_cycles):
      minimized = run_minimize(
        calculator            = calculator,
        params                = params,
        results               = results,
        geometry_rmsd_manager = geometry_rmsd_manager,
        mode                  = "refine")
      if(minimized is not None):
        calculator.reset_fmodel(fmodel = fmodel)
        calculator.update_fmodel()
        n_fev += minimized.number_of_function_and_gradients_evaluations
        break
    if(minimized is None): continue
    run_collect(
      n_fev                 = n_fev,
      results               = results,
      fmodel                = fmodel,
      calculator            = calculator,
      geometry_rmsd_manager = geometry_rmsd_manager)
    results.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(refine_cycle)+"_refine_cycle.pdb")
    results.show(prefix="  ")
    if(conv_test.is_converged(fmodel=fmodel)):
      print >> results.log, "Converged (model)."
      break
    calculator.weights.adjust_restraints_weight_scale(
      fmodel                = fmodel,
      geometry_rmsd_manager = geometry_rmsd_manager,
      max_bond_rmsd         = params.refine.max_bond_rmsd,
      scale                 = 1.2)
  rst_data.write_rst_file(
      rst_file     = rst_file,
      refine_cycle = refine_cycle+1,
      weight_cycle = weight_cycle,
      fmodel       = fmodel,
      weights      = calculator.weights,
      conv_test    = conv_test,
      results      = results)
  print >> results.log, "At end of further refinement:"
  results.show(prefix="  ")

def opt(xray_structure,
        params,
        results,
        calculator,
        geometry_rmsd_manager):
  log_switch = None
  if (params.refine.opt_log or params.debug): log_switch=results.log
  rst_file = params.rst_file
  rst_data = restart_data(xray_structure, geometry_rmsd_manager)
  if(os.path.isfile(rst_file)):
    with open(rst_file, 'rb') as handle:
      rst_file_data = pickle.load(handle)
      micro_cycle_start = rst_file_data["micro_cycle"]
      print >> results.log, "\n***********************************************************"
      print >> results.log, "restarts from micro_cycle: %d"%(micro_cycle_start)
      print >> results.log, "***********************************************************\n"
      ## check the restart fmodel
      xray_structure = calculator.xray_structure
  else:
    micro_cycle_start = 1
  try:
    clustering = calculator.restraints_manager.clustering
  except :
    clustering = False
  if(clustering):
    cluster_qm_update = clustering_update(
      calculator.xray_structure.sites_cart(), results.log, \
      params.rmsd_tolerance * 100)
    print >> results.log, "\ninteracting pairs number:  ",\
      calculator.restraints_manager.fragments.interacting_pairs
  results.show(prefix="start")
  for micro_cycle in xrange(micro_cycle_start,
                        params.refine.number_of_micro_cycles+micro_cycle_start):
    if(clustering):
      cluster_qm_update.re_clustering(calculator)
    conv_test = convergence(
      xray_structure=calculator.xray_structure, params=params)
    rst_data.write_rst_file(rst_file, micro_cycle=micro_cycle,
      xray_structure=xray_structure,
      results=results)
    minimized = run_minimize(
      calculator            = calculator,
      params                = params,
      results               = results,
      geometry_rmsd_manager = geometry_rmsd_manager,
      mode                  = "refine")
    #if minimized is None: break
    calculator.update_xray_structure()
    cctbx_rm_bonds_rmsd = calculator_module.get_bonds_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = xray_structure)
    results.update(
      b      = cctbx_rm_bonds_rmsd,
      xrs    = xray_structure,
      n_fev  = minimized.number_of_function_and_gradients_evaluations)
    results.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(micro_cycle)+"_opt_cycle.pdb")
    results.show(prefix="micro_cycle")
    if(conv_test.is_geometry_converged(
       sites_cart = xray_structure.sites_cart())):
      print >> results.log, " Convergence at micro_cycle:", micro_cycle
      break
    params.refine.pre_opt = None # no pre-optimizations after first macrocycle
  rst_data.write_rst_file(rst_file, micro_cycle=micro_cycle+1,
    xray_structure=xray_structure, results=results)


def run_gradient(calculator):
    eg=calculator.target_and_gradients(x = calculator.x)
    return list(eg[1])

