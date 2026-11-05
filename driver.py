from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import scitbx.lbfgs
from libtbx import easy_pickle
from scitbx.array_family import flex
from . import calculator as calculator_module
from libtbx import group_args
from libtbx.test_utils import approx_equal
from scitbx import minimizers
import sys

log = sys.stdout

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

def run_minimize(calculator, params, mode=None):
  max_iterations = params.refine.max_iterations_refine
  if(mode == "weight"):
    max_iterations = params.refine.max_iterations_weight
  assert params.refine.minimizer in ["lbfgs", "lbfgsb"]
  if(params.refine.minimizer == "lbfgs"):
      core_params = scitbx.lbfgs.core_parameters(
        stpmin = 1.e-9,
        stpmax = params.refine.stpmax)
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
  return minimized

def refine(fmodel,
           params,
           monitor,
           calculator,
           geometry_rmsd_manager):
  if(not params.refine.refine_sites): return
  clustering = params.cluster.clustering
  if(clustering):
    cluster_qm_update = clustering_update(
      pre_sites_cart = calculator.sites_cart_nuclear.deep_copy(),
      log            = log,
      rmsd_tolerance = params.cluster.re_calculate_rmsd_tolerance)
  print("Start:", file=log)
  monitor.show(prefix="  ")

  if(not params.refine.skip_weight_search):
    print("Optimal weight search:", file=log)
    fmodel_copy = calculator.fmodel.deep_copy()
    data_weight = calculator_module.compute_weight(
      fmodel             = fmodel_copy,
      hdm                = calculator.hdm,
      exclude_selection  = calculator.exclude_selection,
      restraints_manager = calculator.restraints_manager,
      sites_cart         = calculator.sites_cart_nuclear,
      mapper             = calculator.mapper)
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
    # Nuclear start sites: fmodel_copy holds their scattering view
    sites_cart_start = calculator.sites_cart_nuclear.deep_copy()
    #
    # Loop over weight search cycles
    for weight_cycle in range(params.refine.number_of_weight_search_cycles):
      fmodel = fmodel_copy.deep_copy() # Always use initial unchanged fmodel
      calculator.reset_fmodel(fmodel = fmodel, sites_cart = sites_cart_start)
      if(clustering):
        cluster_qm_update.re_clustering(calculator,
          sites_cart = calculator.sites_cart_nuclear)
      # Run minimization with given weight
      run_minimize(
        calculator = calculator,
        params     = params,
        mode       = "weight")
      monitor.update(fmodel = fmodel, sites_cart = calculator.sites_cart_nuclear)
      b_rmsds.append(round(monitor.b_rmsd,3))
      a_rmsds.append(round(monitor.a_rmsd,2))
      rws = calculator.restraints_weight_scale
      # Show
      msg = "rws: %6.3f n_fev: %d"%(rws, calculator.n_tg_calls)
      monitor.show(prefix="%d:"%weight_cycle, suffix = msg)
      # Sanity check:
      assert approx_equal(fmodel.r_work(), calculator.fmodel.r_work(), 1.e-4)
      #
      # DEFINE STOPPING RULE
      GOOD = monitor.b_rmsd <= params.refine.max_bond_rmsd #and \
             #monitor.a_rmsd <= params.refine.max_angle_rmsd
      #
      if(GOOD):
        down += 1
        restraints_weight_scale.append(rws)
        calculator.scale_restraints_weight_down(scale=1.5)
        r_frees   .append(monitor.r_free)
        sites_cart.append(monitor.sites_cart)
        #if(params.refine.stop_one_found_first_good_weight):
        #  print("Optimal weight found. Stopping weight search.", file=log)
        #  break
      else:
        up += 1
        calculator.scale_restraints_weight_up(scale=1.5)
      # Show
      msg = "rws: %6.3f n_fev: %d (NEXT)"%(
        calculator.restraints_weight_scale, 0)
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
      calculator.set_sites(sites_cart[index])
      fmodel.update_xray_structure(
        update_f_calc  = True,
        update_f_mask  = True
        )
    fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
    print("Best Rwork, Rfree (at refinement start): %6.4f %6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)
    monitor.update(fmodel = fmodel, sites_cart = calculator.sites_cart_nuclear)
    print("At refinement start", file = log)
    monitor.show()
  #
  # Done with weights. Now let's refine!
  #
  print("Refinement:", file=log)
  calculator.reset_fmodel(fmodel = fmodel)
  assert not calculator.converged()
  for refine_cycle in range(params.refine.number_of_refine_cycles):
    if(clustering):
      cluster_qm_update.re_clustering(calculator,
        sites_cart = calculator.sites_cart_nuclear)
    run_minimize(
      calculator = calculator,
      params     = params,
      mode       = "refine")
    n_fev = calculator.n_tg_calls
    fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
    calculator.reset_fmodel(fmodel = fmodel)
    monitor.update(fmodel = fmodel, sites_cart = calculator.sites_cart_nuclear)
    msg = "rws: %6.3f n_fev: %d"%(
      calculator.restraints_weight_scale, n_fev)
    monitor.show(prefix="%d:"%refine_cycle, suffix=msg)
    monitor.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(refine_cycle)+"_refine_cycle.pdb")
    if(calculator.converged()):
      print("Refinement converged. Stopping now.", file=log)
      break
  print("At refinement end:", file = log)
  monitor.show()
  assert approx_equal(fmodel.r_work(), calculator.fmodel.r_work(), 1.e-4)
  print("Best Rwork, Rfree (after refinement): %6.4f %6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)

def opt(params, monitor, calculator):
  assert monitor.model == calculator.model
  monitor.show(prefix="start:")
  if(params.cluster.clustering):
    cluster_qm_update = clustering_update(
      calculator.model.get_sites_cart(), log,
      params.cluster.re_calculate_rmsd_tolerance)
    print("\ninteracting pairs number:  ",\
      len(calculator.restraints_manager.fragment_manager.interaction_list), file=log)
  for micro_cycle in range(params.refine.number_of_micro_cycles):
    assert monitor.model == calculator.model
    if(params.cluster.clustering):
      cluster_qm_update.re_clustering(calculator,
        sites_cart = calculator.model.get_sites_cart())
    run_minimize(calculator = calculator, params = params)
    monitor.update(model = calculator.model)
    monitor.show(prefix="cycle %3d: f=%.9g nfev: %3d | "%(
      micro_cycle, calculator.f, calculator.number_of_target_and_gradients_calls))
    monitor.write_pdb_file(
      output_folder_name = params.output_folder_name,
      output_file_name   = str(micro_cycle)+"_opt_cycle.pdb")
    if(calculator.converged()):
      print("Convergence reached. Stopping now.", file=log)
      break
  print("calculator(opt), total_time (target_and_gradients)",
    calculator.total_time, file=log)
  print("calculator(opt), number_of_target_and_gradients_calls (target_and_gradients)",
    calculator.number_of_target_and_gradients_calls, file=log)
  if calculator.debug:
    o = group_args(
      targets   = calculator.all_targets,
      gradients = calculator.all_gradients)
    easy_pickle.dump(
      "debug_optimization_clustering%s.pkl"%str(params.cluster.clustering), o)
