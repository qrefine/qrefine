from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import scitbx.lbfgs
from libtbx import easy_pickle
from . import calculator as calculator_module
from . import weights
from libtbx import group_args
from libtbx.test_utils import approx_equal
from scitbx import minimizers
import os

class clustering_update(object):
  """
  Re-clusters (fragment_manager.set_up_cluster_qm) once the model moved by
  more than rmsd_tolerance from the coordinates of the last clustering.
  """
  def __init__(self, pre_sites_cart, log, rmsd_tolerance):
    self.pre_sites_cart = pre_sites_cart
    self.log = log
    self.rmsd_tolerance = rmsd_tolerance

  def re_clustering(self, calculator, sites_cart):
    rmsd_diff = self.pre_sites_cart.rms_difference(sites_cart)
    if(rmsd_diff > self.rmsd_tolerance):
      print(" rmsd_diff: ", rmsd_diff, "--> need to redo clustering", file=self.log)
      calculator.restraints_manager.fragment_manager.set_up_cluster_qm()
      self.pre_sites_cart = sites_cart.deep_copy()

def stop_reason(minimizer, n_iterations, max_iterations):
  """
  Why a minimization ended: "error" when the minimizer swallowed a
  line-search error (its error attribute is set), "cap" when it used its
  iteration budget, "conv" when its own convergence test stopped it.
  """
  if(getattr(minimizer, "error", None) is not None): return "error"
  if(n_iterations >= max_iterations): return "cap"
  return "conv"

def run_minimize(calculator, params, mode=None):
  max_iterations = params.refine.max_iterations_refine
  if(mode == "weight"):
    max_iterations = params.refine.max_iterations_weight
  assert params.refine.minimizer in ["lbfgs", "lbfgsb"]
  if(params.refine.minimizer == "lbfgs"):
      core_params = scitbx.lbfgs.core_parameters(
        max_shift            = params.refine.stpmax,
        first_step_max_shift = calculator_module.FIRST_STEP_MAX_SHIFT)
      minimized = minimizers.lbfgs(
        calculator     = calculator,
        mode           = "lbfgs",
        core_params    = core_params,
        max_iterations = max_iterations,
        gradient_only  = params.refine.gradient_only)
  else:
    assert not params.refine.gradient_only
    minimized = minimizers.lbfgs(
      mode           = "lbfgsb",
      calculator     = calculator,
      max_iterations = max_iterations)
  calculator.apply_x()
  # There is no early stop in the calculator: every trial of the weight
  # search gets the same budget.
  if(params.refine.minimizer == "lbfgs"):
    n_iterations = minimized.minimizer.iter()
  else:
    n_iterations = minimized.minimizer.n_calls
  return group_args(
    minimized   = minimized,
    stop_reason = stop_reason(minimized.minimizer, n_iterations,
      max_iterations))

# Cycles added to refine.number_of_refine_cycles while the last cycle ends
# with the bond rmsd above refine.max_bond_rmsd.
MAX_EXTRA_CYCLES = 3

def bond_rmsd_weight_factor(bond_rmsd, max_bond_rmsd):
  """
  The factor (> 1) for the restraints weight scale after a refinement cycle
  that ended with bond_rmsd above max_bond_rmsd. The bond rmsd of a
  minimized model falls roughly as scale**-0.8 (fixed-weight scans, 23
  models at 3-4 A): the exponent 1.5 aims a little below the limit; the
  factor is kept between 1.1 and 2 per cycle.
  """
  assert bond_rmsd > max_bond_rmsd > 0
  return min(2.0, max(1.1, (bond_rmsd / max_bond_rmsd)**1.5))

def refine(fmodel,
           params,
           monitor,
           calculator,
           geometry_rmsd_manager,
           log,
           intermediate_folder=None,
           bulk_solvent_params=None):
  """
  Weight search (unless refine.skip_weight_search or refine.data_weight) and
  refinement cycles. Every print goes to log; the model of
  each refinement cycle is written into intermediate_folder unless it is
  None; bulk_solvent_params goes to every update_all_scales call. Returns
  group_args(n_cycles, n_fev): refinement cycles run and the target and
  gradient evaluations in them (the weight search not counted).
  """
  if(not params.refine.refine_sites):
    return group_args(n_cycles = 0, n_fev = 0)
  clustering = params.cluster.clustering
  if(clustering):
    cluster_qm_update = clustering_update(
      pre_sites_cart = calculator.sites_cart_nuclear.deep_copy(),
      log            = log,
      rmsd_tolerance = params.cluster.re_calculate_rmsd_tolerance)
  print("Start:", file=log)
  monitor.show(prefix="  ", full=True)

  #
  # Weights. refine.data_weight (when given) is the data weight and by-passes
  # the search, else the data weight is the ratio of the gradient norms;
  # refine.restraints_weight_scale is the restraints weight scale without the
  # search and its start value with it.
  #
  skip_weight_search = (params.refine.skip_weight_search or
                        params.refine.data_weight is not None)
  data_weight = params.refine.data_weight
  if(data_weight is None):
    data_weight = calculator_module.compute_weight(
      fmodel             = calculator.fmodel.deep_copy(),
      hdm                = calculator.hdm,
      exclude_selection  = calculator.exclude_selection,
      restraints_manager = calculator.restraints_manager,
      sites_cart         = calculator.sites_cart_nuclear,
      mapper             = calculator.mapper)
  calculator.setw(
    data_weight             = data_weight,
    restraints_weight_scale = params.refine.restraints_weight_scale,
    restraints_weight       = 1.)
  if(skip_weight_search):
    print("No weight search: data weight %s (%s), restraints weight scale %s"%(
      data_weight,
      "refine.data_weight" if params.refine.data_weight is not None
        else "gradient norms",
      params.refine.restraints_weight_scale), file=log)
  else:
    print("Optimal weight search:", file=log)
    fmodel_copy = calculator.fmodel.deep_copy()
    c = weights.controller(
      start_scale = params.refine.restraints_weight_scale,
      factor      = params.refine.adjust_restraints_weight_scale_value,
      max_trials  = params.refine.number_of_weight_search_cycles,
      r_tolerance = params.refine.r_tolerance)
    print("Data weight (initial):", data_weight, file=log)
    # Nuclear start sites: fmodel_copy holds their scattering view
    sites_cart_start = calculator.sites_cart_nuclear.deep_copy()
    #
    # Loop over weight search trials (qrefine.weights: next scale, stop,
    # choice)
    while(not c.done()):
      rws = c.scale()
      fmodel = fmodel_copy.deep_copy() # Always use initial unchanged fmodel
      calculator.reset_fmodel(fmodel = fmodel, sites_cart = sites_cart_start)
      calculator.setw(restraints_weight_scale = rws)
      if(clustering):
        cluster_qm_update.re_clustering(calculator,
          sites_cart = calculator.sites_cart_nuclear)
      # Run minimization with given weight
      m = run_minimize(
        calculator = calculator,
        params     = params,
        mode       = "weight")
      # The trial is scored like a refinement cycle: with the bulk-solvent
      # mask and the scales of its own model (those of the start model
      # shift R-free by more than the trials differ).
      fmodel.update_xray_structure(
        update_f_calc       = True,
        update_f_mask       = True,
        force_update_f_mask = True)
      fmodel.update_all_scales(remove_outliers=False,
        refine_hd_scattering=False, params=bulk_solvent_params)
      monitor.update(fmodel = fmodel, sites_cart = calculator.sites_cart_nuclear)
      # Show
      msg = "rws: %.3g n_fev: %d stop: %s"%(
        rws, calculator.n_tg_calls, m.stop_reason)
      monitor.show(prefix="%d:"%len(c.trials), suffix = msg)
      # Sanity check:
      assert approx_equal(fmodel.r_work(), calculator.fmodel.r_work(), 1.e-4)
      c.record(
        scale      = rws,
        good       = monitor.b_rmsd <= params.refine.max_bond_rmsd,
        r_free     = monitor.r_free,
        sites_cart = monitor.sites_cart)
    # Ok, done with weights.
    ch = c.choice()
    print("Weight search: %d trials, chosen rws %.3f (%s)"%(
      len(c.trials), ch.scale, ch.reason), file=log)
    if(c.stopped_by_cap()):
      print("Weight search stopped by number_of_weight_search_cycles=%d "
        "before a stopping rule applied (R-free may still improve): "
        "consider raising it." % params.refine.number_of_weight_search_cycles,
        file=log)
    if(ch.reason == "fallback"):
      print("Weight search yields no result. Using last result.", file=log)
    else:
      print("Best Rfree from above: %6.4f"%ch.r_free, file=log)
      print("Best restraints scale:", round(ch.scale, 3), file=log)
    calculator.setw(restraints_weight_scale = ch.scale)
    calculator.set_sites(ch.sites_cart)
    fmodel.update_xray_structure(
      update_f_calc       = True,
      update_f_mask       = True,
      force_update_f_mask = True)
    fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False,
      params=bulk_solvent_params)
    print("Best Rwork, Rfree (at refinement start): %6.4f %6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)
    monitor.update(fmodel = fmodel, sites_cart = calculator.sites_cart_nuclear)
    print("At refinement start", file = log)
    monitor.show(full=True)
  #
  # Done with weights. Now let's refine!
  #
  print("Refinement:", file=log)
  calculator.reset_fmodel(fmodel = fmodel)
  assert not calculator.converged()
  n_cycles, n_fev_total = 0, 0
  #
  # Refinement cycles. The weight is the one of the search, found with
  # shorter minimizations: a cycle that ends with the bond rmsd above
  # refine.max_bond_rmsd scales the restraints weight up for the next one
  # (bond_rmsd_weight_factor), and up to MAX_EXTRA_CYCLES cycles are added
  # while the last cycle is above the limit.
  #
  max_bond_rmsd = params.refine.max_bond_rmsd
  n_extra = 0
  refine_cycle = 0
  while(refine_cycle < params.refine.number_of_refine_cycles + n_extra):
    if(clustering):
      cluster_qm_update.re_clustering(calculator,
        sites_cart = calculator.sites_cart_nuclear)
    m = run_minimize(
      calculator = calculator,
      params     = params,
      mode       = "refine")
    n_fev = calculator.n_tg_calls
    n_cycles += 1
    n_fev_total += n_fev
    # The bulk-solvent mask of the moved model, so that the R factors of
    # this cycle are those of the model as written (the minimizer updates
    # f_calc only).
    fmodel.update_xray_structure(
      update_f_calc       = True,
      update_f_mask       = True,
      force_update_f_mask = True)
    fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False,
      params=bulk_solvent_params)
    calculator.reset_fmodel(fmodel = fmodel)
    monitor.update(fmodel = fmodel, sites_cart = calculator.sites_cart_nuclear)
    msg = "rws: %.3g n_fev: %d stop: %s"%(
      calculator.restraints_weight_scale, n_fev, m.stop_reason)
    monitor.show(prefix="%d:"%refine_cycle, suffix=msg)
    if(intermediate_folder is not None):
      monitor.write_pdb_file(
        output_folder_name = intermediate_folder,
        output_file_name   = str(refine_cycle)+"_refine_cycle")
    if(calculator.converged()):
      print("Refinement converged. Stopping now.", file=log)
      break
    refine_cycle += 1
    if(max_bond_rmsd is not None and monitor.b_rmsd > max_bond_rmsd):
      rws = calculator.restraints_weight_scale
      rws_new = rws * bond_rmsd_weight_factor(monitor.b_rmsd, max_bond_rmsd)
      is_last = refine_cycle >= params.refine.number_of_refine_cycles + n_extra
      if(is_last and n_extra >= MAX_EXTRA_CYCLES):
        print("  bond rmsd %.4f is above max_bond_rmsd=%g after %d extra "
          "cycles: stopping."%(monitor.b_rmsd, max_bond_rmsd, n_extra),
          file=log)
        break
      calculator.setw(restraints_weight_scale = rws_new)
      print("  bond rmsd %.4f > max_bond_rmsd=%g: restraints weight scale "
        "%.3g -> %.3g%s"%(monitor.b_rmsd, max_bond_rmsd, rws, rws_new,
        " (extra cycle)" if is_last else ""), file=log)
      if(is_last): n_extra += 1
  print("At refinement end:", file = log)
  monitor.show(full=True)
  assert approx_equal(fmodel.r_work(), calculator.fmodel.r_work(), 1.e-4)
  print("Best Rwork, Rfree (after refinement): %6.4f %6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)
  return group_args(n_cycles = n_cycles, n_fev = n_fev_total)

def opt(params, monitor, calculator, log, working_folder,
        intermediate_folder=None):
  """
  Geometry optimization cycles. Every print goes to log; the model of each
  cycle is written into intermediate_folder unless it is None; with
  calculator.debug the targets and gradients are pickled into
  working_folder (not with working_folder None). Returns group_args(n_cycles,
  n_fev): optimization cycles run and target and gradient evaluations.
  """
  assert monitor.model == calculator.model
  monitor.show(prefix="start:", full=True)
  if(params.cluster.clustering):
    cluster_qm_update = clustering_update(
      calculator.model.get_sites_cart(), log,
      params.cluster.re_calculate_rmsd_tolerance)
    print("\ninteracting pairs number:  ",\
      len(calculator.restraints_manager.fragment_manager.interaction_list), file=log)
  n_cycles = 0
  for micro_cycle in range(params.refine.number_of_micro_cycles):
    n_cycles += 1
    assert monitor.model == calculator.model
    if(params.cluster.clustering):
      cluster_qm_update.re_clustering(calculator,
        sites_cart = calculator.model.get_sites_cart())
    run_minimize(calculator = calculator, params = params)
    monitor.update(model = calculator.model)
    monitor.show(prefix="cycle %3d: f=%.9g nfev: %3d | "%(
      micro_cycle, calculator.f, calculator.number_of_target_and_gradients_calls))
    if(intermediate_folder is not None):
      monitor.write_pdb_file(
        output_folder_name = intermediate_folder,
        output_file_name   = str(micro_cycle)+"_opt_cycle")
    if(calculator.converged()):
      print("Convergence reached. Stopping now.", file=log)
      break
  print("calculator(opt), total_time (target_and_gradients)",
    calculator.total_time, file=log)
  print("calculator(opt), number_of_target_and_gradients_calls (target_and_gradients)",
    calculator.number_of_target_and_gradients_calls, file=log)
  if(calculator.debug and working_folder is not None):
    o = group_args(
      targets   = calculator.all_targets,
      gradients = calculator.all_gradients)
    os.makedirs(working_folder, exist_ok=True)
    easy_pickle.dump(os.path.join(working_folder,
      "debug_optimization_clustering%s.pkl"%str(params.cluster.clustering)), o)
  return group_args(n_cycles = n_cycles,
    n_fev = calculator.number_of_target_and_gradients_calls)
