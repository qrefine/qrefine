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

class convergence(object):
  def __init__(self, fmodel, params):
    self.r_start = fmodel.r_work()
    self.sites_cart_start = fmodel.xray_structure.sites_cart()
    self.r_tolerance=params.refine.r_tolerance
    self.max_bond_rmsd=params.refine.max_bond_rmsd
    self.rmsd_tolerance=params.refine.rmsd_tolerance
    self.use_convergence_test = params.refine.use_convergence_test
    self.number_of_convergence_occurances=0
    #
    self.rws = flex.double()
    self.rfs = flex.double()
    self.gaps = flex.double()
    self.restraints_weight_scale = flex.double()

  def is_converged(self, fmodel, restraints_weight_scale=None):
    #
    rw = fmodel.r_work()
    rf = fmodel.r_free()
    gap = rf-rw
    self.rws                     .append(rw)
    self.rfs                     .append(rf)
    self.gaps                    .append(gap)
    if(restraints_weight_scale is not None):
      self.restraints_weight_scale.append(restraints_weight_scale)
    #
    if(restraints_weight_scale is not None):
      rwc = self.restraints_weight_scale
      i_last = self.rws.size()-1
      if(i_last>3):
        rwc_3 = rwc[i_last]
        rwc_2 = rwc[i_last-1]
        rwc_1 = rwc[i_last-2]
        rws123 = [rwc_1, rwc_2, rwc_3]
        for rwc_i in rws123:
          if(rws123.count(rwc_i)>1): return True
    #
    r = rw
    sites_cart = fmodel.xray_structure.sites_cart()
    r_diff=abs(self.r_start-r)
    rmsd_diff=self.sites_cart_start.rms_difference(sites_cart)
    self.sites_cart_start = sites_cart
    self.r_start=r
    if(r_diff<self.r_tolerance and rmsd_diff<self.rmsd_tolerance):
      self.number_of_convergence_occurances+=2
    if(self.number_of_convergence_occurances==2 or r<0.005):
      return True and self.use_convergence_test
    else: return False and self.use_convergence_test

  def is_geometry_converged(self, sites_cart):
    rmsd_diff=self.sites_cart_start.rms_difference(sites_cart)
    if(rmsd_diff<self.rmsd_tolerance):
      return True

  def geometry_exploded(self, fmodel, geometry_rmsd_manager):
    result = False
    cctbx_rm_bonds_rmsd = calculator_module.get_bonds_angles_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = fmodel.xray_structure)
    if(cctbx_rm_bonds_rmsd>self.max_bond_rmsd*2.0):
      result = True
    return result

  def is_weight_scale_converged(self, restraints_weight_scale,
        restraints_weight_scales):
    return restraints_weight_scale in restraints_weight_scales

class minimizer(object):
  def __init__(self,
        stpmax,
        calculator,
        max_iterations,
        gradient_only,
        line_search):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.x_previous = None
    self.number_of_function_and_gradients_evaluations=0
    self.lbfgs_core_params = scitbx.lbfgs.core_parameters(
      stpmin = 1.e-9,
      stpmax = stpmax)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      gradient_only=gradient_only,
      line_search=line_search,
      core_params=self.lbfgs_core_params,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))

  def compute_functional_and_gradients(self):
    self.number_of_function_and_gradients_evaluations += 1
    # Ad hoc damping shifts; note arbitrary 1.0 below
    x_current = self.x
    if(self.x_previous is None):
      self.x_previous = x_current.deep_copy()
    else:
      xray.ext.damp_shifts(previous=self.x_previous, current=x_current,
        max_value=1.)
      self.x_previous = x_current.deep_copy()
    #
    return self.calculator.target_and_gradients(x = self.x)

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
  def __init__(self, fmodel, geometry_rmsd_manager):
    rst_data = {}
    rst_data["fmodel"] = fmodel
    rst_data["geometry_rmsd_manager"] = geometry_rmsd_manager
    self.rst_data = rst_data

  def write_rst_file(self, rst_file, weight_cycle = None, refine_cycle = None,
                     micro_cycle = None, fmodel = None, weights = None,
                     conv_test = None, results = None):
    self.rst_data["weight_cycle"] = weight_cycle
    self.rst_data["refine_cycle"] = refine_cycle
    self.rst_data["micro_cycle"] = micro_cycle
    self.rst_data["rst_fmodel"] = fmodel
    self.rst_data["weights"] = weights
    self.rst_data["conv_test"] = conv_test
    self.rst_data["results"] = results
    easy_pickle.dump(file_name=rst_file, obj=self.rst_data)

def run_minimize(calculator, params, results):
  minimized = None
  try:
    if(params.refine.max_iterations > 0):
      minimized = minimizer(
        calculator     = calculator,
        stpmax         = params.refine.stpmax,
        gradient_only  = params.refine.gradient_only,
        line_search    = params.refine.line_search,
        max_iterations = params.refine.max_iterations)
  except Exception as e:
    print >> results.log, e
    print >> results.log, "  minimizer failed"
    results.log.flush()
  return minimized

def run_collect(n_fev, results, fmodel, geometry_rmsd_manager, calculator):
  cctbx_rm_bonds_rmsd = calculator_module.get_bonds_angles_rmsd(
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
  rst_data = restart_data(fmodel, geometry_rmsd_manager)
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
      rmsd_tolerance = params.refine.rmsd_tolerance * 100)
    print >> results.log, "\ninteracting pairs number:  ", \
      calculator.restraints_manager.fragments.interacting_pairs
  #
  # Optimal weight search loop
  #
  print >> results.log, "Start:"
  results.show(prefix="  ")
  print >> results.log, "Optimal weight search:"
  fmodel_copy = fmodel.deep_copy()
  for weight_cycle in xrange(weight_cycle_start,
                             params.refine.number_of_weight_search_cycles+1):
    if((weight_cycle!=1 and weight_cycle==weight_cycle_start) or
       (refine_cycle_start is not None)):
      fmodel = calculator.fmodel.deep_copy()
    else:
      fmodel = fmodel_copy.deep_copy()
    if(refine_cycle_start is not None):
      print >> results.log, \
        "Found the optimal weight. Proceed to further refinement with this weight."
      break
    calculator.reset_fmodel(fmodel = fmodel)
    if(clustering):
      cluster_qm_update.re_clustering(calculator)
    # Calculate weight
    try:
      calculator.calculate_weight()
    except Exception as e:
      print >> results.log, e
      raise Sorry("Failed to get weight")
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
      minimized = run_minimize(calculator=calculator, params=params,
        results=results)
      calculator.reset_fmodel(fmodel = fmodel)
      calculator.update_fmodel()
      n_fev += minimized.number_of_function_and_gradients_evaluations
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
    #
    if(conv_test.is_converged(fmodel=fmodel, restraints_weight_scale =
       calculator.weights.restraints_weight_scale)):
      break
    calculator.weights.adjust_restraints_weight_scale(
      fmodel                = fmodel,
      geometry_rmsd_manager = geometry_rmsd_manager,
      max_bond_rmsd         = params.refine.max_bond_rmsd)
    results.show(prefix="  ")
    if(params.refine.mode == "refine"):
      calculator.weights.\
        add_restraints_weight_scale_to_restraints_weight_scales()
  print >> results.log, "At end of weight search:"
  results.show(prefix="  ")
  #
  # Done with weight search. Now get best result and refine it further.
  #
  xrs_best, dummy, dummy, wsc_best = results.choose_best()
  fmodel.update_xray_structure(
    xray_structure = xrs_best,
    update_f_calc  = True,
    update_f_mask  = True)
  fmodel.update_all_scales(remove_outliers=False)
  calculator.update_restraints_weight_scale(restraints_weight_scale=wsc_best)
  ####
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
  if(refine_cycle_start is None): refine_cycle_start = 1
  for refine_cycle in xrange(refine_cycle_start, 5+1):
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
      minimized = run_minimize(calculator=calculator, params=params,
        results=results)
      calculator.reset_fmodel(fmodel = fmodel)
      calculator.update_fmodel()
      n_fev += minimized.number_of_function_and_gradients_evaluations
    #for micro_cycle in xrange(params.number_of_macro_cycles):
    #  minimized = run_minimize(calculator=calculator, params=params,
    #    results=results)
    #  calculator.update_fmodel()
    #
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
      break
    calculator.weights.adjust_restraints_weight_scale(
        fmodel                = fmodel,
        geometry_rmsd_manager = geometry_rmsd_manager,
        max_bond_rmsd         = params.refine.max_bond_rmsd)
  print >> results.log, "At end of further refinement:"
  results.show(prefix="  ")

def opt(fmodel,
        params,
        results,
        calculator,
        geometry_rmsd_manager):
  rst_file = params.rst_file
  rst_data = restart_data(fmodel, geometry_rmsd_manager)
  if(os.path.isfile(rst_file)):
    with open(rst_file, 'rb') as handle:
      rst_file_data = pickle.load(handle)
      micro_cycle_start = rst_file_data["micro_cycle"]
      print >> results.log, "\n***********************************************************"
      print >> results.log, "restarts from micro_cycle: %d"%( micro_cycle_start)
      print >> results.log, "***********************************************************\n"
  else:
    micro_cycle_start = 1
  try:
    clustering = calculator.restraints_manager.clustering
  except :
    clustering = False
  if(clustering):
    cluster_qm_update = clustering_update(
      calculator.fmodel.xray_structure.sites_cart(), results.log, \
      params.rmsd_tolerance * 100)
    print >> results.log, "\ninteracting pairs number:  ",\
      calculator.restraints_manager.fragments.interacting_pairs
  results.show(prefix="start")
  for micro_cycle in xrange(micro_cycle_start, params.refine.number_of_micro_cycles+1):
    if(clustering):
      cluster_qm_update.re_clustering(calculator)
    conv_test = convergence(fmodel=calculator.fmodel, params=params)
    rst_data.write_rst_file(rst_file, micro_cycle=micro_cycle, fmodel=fmodel,
      results=results)
    minimized = minimizer(calculator = calculator,
                          stpmax = params.refine.stpmax,
                          gradient_only = params.refine.gradient_only,
                          line_search = params.refine.line_search,
                          max_iterations = params.refine.max_iterations)
    calculator.update_fmodel_opt()
    cctbx_rm_bonds_rmsd = calculator_module.get_bonds_angles_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = fmodel.xray_structure)
    results.update(
      r_work = fmodel.r_work(),
      r_free = fmodel.r_free(),
      b      = cctbx_rm_bonds_rmsd,
      xrs    = fmodel.xray_structure,
      n_fev  = minimized.number_of_function_and_gradients_evaluations)
    results.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(micro_cycle)+"_opt_cycle.pdb")
    results.show(prefix="micro_cycle")
    if(conv_test.is_geometry_converged(
       sites_cart = fmodel.xray_structure.sites_cart())):
      print >> results.log, " Convergence at micro_cycle:", micro_cycle
      break
