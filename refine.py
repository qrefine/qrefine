"""This is the entry point into Q|R. It takes user inputs,
   and then constructs all of the objects
   needed to carry out the quantum refinement.

   - fmodel (crystallographic information)
  - calculator (composite object)
  - restraints_manager (computes energy and gradients using either qm codes or
    cctbx (standard))
  - geometry_restraints_manager (analyses geometry e.g. bond RMSDs)
  - weights (scale factors needed to scale up or down data versus restraints
    contributions)
  Then we process them by the refinement/optimization engine, driver.py:

   driver.refine(params   = params,
               fmodel     = fmodel,
               calculator = calculator_manager,
               results    = results)
   results_manager (store all reportable information, and write it out as a log,
   and also write our final pdb structure.)
   """
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import time
import pickle

from numpy.lib.function_base import select
import mmtbx.command_line
import mmtbx.f_model
import mmtbx.utils
from libtbx.utils import Sorry
from libtbx import easy_pickle
from libtbx import group_args
from libtbx.utils import null_out
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from qrefine.fragment import fragments

from qrefine import calculator
from qrefine import driver
from qrefine import restraints
from qrefine import cluster_restraints
from qrefine import results

from qrefine.super_cell import expand
import mmtbx.model.statistics
from libtbx import Auto
from ase.io import read as ase_io_read
from cctbx import maptbx
from scitbx.array_family import flex

def show_cc(map_data, xray_structure, log=None):
  import mmtbx.maps.mtriage
  from mmtbx.maps.correlation import five_cc
  xrs = xray_structure
  xrs.scattering_type_registry(table="electron")
  d99 = maptbx.d99(
    map = map_data, crystal_symmetry = xrs.crystal_symmetry()).result.d99
  if(log is not None): print("Resolution of map is: %6.4f" %d99, file=log)
  result = five_cc(
    map               = map_data,
    xray_structure    = xrs,
    d_min             = d99,
    compute_cc_box    = False,
    compute_cc_mask   = True,
    compute_cc_peaks  = False,
    compute_cc_volume = False).result
  if(log is not None):
    print("Map-model correlation coefficient (CC)", file=log)
    print("  CC_mask  : %6.4f" %result.cc_mask, file=log)
  return result.cc_mask

def create_fragment_manager(
      model,
      params,
      file_name = os.path.join("ase","tmp_ase.pdb")):
  if(not params.cluster.clustering): return None
  return fragments(
    cif_objects                = model._restraint_objects,
    working_folder             = os.path.split(file_name)[0]+ "/",
    clustering_method          = params.cluster.clustering_method,
    altloc_method              = params.cluster.altloc_method,
    maxnum_residues_in_cluster = params.cluster.maxnum_residues_in_cluster,
    charge_embedding           = params.cluster.charge_embedding,
    two_buffers                = params.cluster.two_buffers,
    pdb_hierarchy              = model.get_hierarchy(),
    qm_engine_name             = params.quantum.engine_name,
    crystal_symmetry           = model.crystal_symmetry(),
    debug                      = params.debug,
    charge_cutoff              = params.cluster.charge_cutoff,
    save_clusters              = params.cluster.save_clusters,
    select_within_radius       = params.cluster.select_within_radius,
    bond_with_altloc_flag      = params.cluster.bond_with_altloc)

class hd_mapper(object):
  """
  Utility function that allows to use a model with exhangeable H/D sites as
  single protonated model by re-mapping 'short' gradients into the full-size
  gradients.
  """

  def __init__(self, model):
    self.model = model
    self.ehd = model.get_hierarchy().exchangeable_hd_selections()
    self.elements = [e.strip().upper() for e in
      model.get_hierarchy().atoms().extract_element()]
    self.size = model.size()
    self.keep = flex.bool(self.size, True)
    for pair in self.ehd:
      i,j = pair[0][0], pair[1][0]
      for it in [i,j]:
        if self.elements[it]=="D": self.keep[it] = False
    self.n_keep = self.keep.count(True)
    self._model = model.select(self.keep)
    self._model.de_deuterate()

  def get_single_model(self):
    return self._model

  def shrink(self, array):
    return array.select(self.keep)

  def check(self, array):
    assert array.size() == self.size
    for pair in self.ehd:
      i,j = pair[0][0], pair[1][0]
      assert array[i] == array[j]

  def average(self, array):
    for pair in self.ehd:
      i,j = pair[0][0], pair[1][0]
      ai = array[i]
      ai = flex.vec3_double([array[i],])
      aj = flex.vec3_double([array[j],])
      aa = (ai+aj)*0.5
      array[i] = aa[0]
      array[j] = aa[0]
    return array

  def expand(self, g_short):
    assert self.n_keep == g_short.size()
    g_all = flex.vec3_double(self.size, [0,0,0])
    g_all = g_all.set_selected(self.keep, g_short)
    for pair in self.ehd:
      i,j = pair[0][0], pair[1][0]
      if   g_all[i] == (0,0,0): g_all[i] = g_all[j]
      elif g_all[j] == (0,0,0): g_all[j] = g_all[i]
      else: assert 0, [g_all[j] , g_all[i]]
      #if self.elements[i]=="H": g_all[j] = g_all[i]
      #else:                     g_all[i] = g_all[j]
    for g in g_all:
      assert g != (0,0,0)
    return g_all

def create_restraints_manager(params, model, altlocs_present):
  restraints_source = restraints.restraints(params = params, model = model)
  #
  # General case of altlocs. Developmenal code. NOT IN PRODUCTION.
  #
  if(altlocs_present):
    print("altlocs_present", altlocs_present)
    assert not model.altlocs_present_only_hd()
    return restraints.from_altlocs2(
      model  = model,
      params = params,
      method = params.cluster.altloc_method)

  # No altlocs of any kind with expansion.
  #
  if(params.expansion):
    assert not params.cluster.clustering
    #
    # TMP development insert to investigate Rss BEGIN
    #
    if(params.debug_rss):
      from collections import OrderedDict
      import time
      from libtbx import easy_pickle
      result = OrderedDict()
      RSS = [10,9,8,7,6,5,4,3,2,1]
      for i, rss in enumerate(RSS):
        t0 = time.time()
        params.cluster.select_within_radius=rss
        FE = restraints.from_expansion(
          params            = params,
          restraints_source = restraints_source)
        _, g = FE.target_and_gradients(sites_cart = model.get_sites_cart())
        result[rss] = g
        print("Trying rss", rss, time.time()-t0)
      easy_pickle.dump(
        file_name="%s_gradients.pkl"%params.output_file_name_prefix, obj = result)
      STOP()
    #
    # TMP development insert to investigate Rss END
    #
    print("Restraints are from restraints.from_expansion")
    return restraints.from_expansion(
      params            = params,
      restraints_source = restraints_source)
  #
  # No altlocs of any kind, no expansion
  #
  else:
    if(params.cluster.clustering):
      fragment_manager = create_fragment_manager(params = params, model = model)
      return cluster_restraints.from_cluster(
        restraints_manager = restraints_source.restraints_manager,
        fragment_manager   = fragment_manager,
        parallel_params    = params.parallel)
    else:
      # restraints=cctbx clustering=false expansion=false
      return restraints_source.restraints_manager

def create_calculator(params, restraints_manager, model, fmodel=None, hdm=None,
                      exclude_selection=None):
  if(params.refine.refine_sites):
    if(params.refine.mode == "refine"):
      assert model is not None
      return calculator.sites(
        fmodel             = fmodel,
        hdm                = hdm,
        exclude_selection  = exclude_selection,
        restraints_manager = restraints_manager,
        dump_gradients     = params.dump_gradients)
    else:
      # XXX
      # XXX exclude_selection and hdm are not used XXX
      # XXX
      return calculator.sites_opt(
        restraints_manager = restraints_manager,
        model              = model,
        dump_gradients     = params.dump_gradients,
        max_shift          = params.refine.stpmax,
        shift_eval         = params.refine.shift_evaluation)

def set_qm_defaults(params, log):
  outl = ''
  if params.quantum.engine_name=='mopac':
    if params.quantum.method==Auto:
      params.quantum.method='PM7'
      outl += '  Setting QM method to PM7 (mopac)\n'
    if params.quantum.basis==Auto:
      params.quantum.basis=''
  if params.quantum.engine_name=='xtb':
    if params.quantum.method==Auto:
      params.quantum.method=' --gfn 2 --etemp 500 --acc 0.1 --gbsa h2o'
      print(' Default method for xtb is %s' % (
          params.quantum.method,
          ), file=log)
      params.quantum.basis = ''
  if params.quantum.engine_name=='aimnet2':
    print(file=log)
    print("Setting aimnet2 specific defaults:", file=log)
    if params.quantum.method==Auto:
      params.quantum.method="aimnet2-qr"
      print(f"  Default Aimnet2 model: {params.quantum.method}", file=log)
      try:
        from aimnet2calc import AIMNet2ASE
      except ModuleNotFoundError as e:
        print(str(e), file=log)
        raise Sorry("AIMNet2 does not seem to be installed or available.")
  else:
    if params.quantum.method==Auto:
      params.quantum.method='HF'
      outl += '  Setting QM method to HF\n'
    if params.quantum.basis==Auto:
      params.quantum.basis='6-31g'
      outl += '  Setting QM basis to 6-31g\n'
  if outl:
    print(' Setting QM defaults', file=log)
    print(outl, file=log)

  if params.quantum.engine_name=='mopac':
    if params.quantum.basis:
      print('  Because engine is %s basis set %s ignored' % (
        params.quantum.engine_name,
        params.quantum.basis,
        ), file=log)
      params.quantum.basis = ''
    if params.quantum.method=='hf': # default
      print(' Default method set as PM7 (mopac)', file=log)
      params.quantum.method='PM7'


def run(model, fmodel, map_data, params, rst_file, prefix, log):
  set_qm_defaults(params,log)
  if(params.cluster.clustering):
    params.refine.gradient_only = True
    print(" params.gradient_only", params.refine.gradient_only, file=log)
  if (model.size() > params.max_atoms):
    print("Too many atoms. Can take forever or crash."*50, file=log)

  geometry_rmsd_manager = calculator.cctbx_geometry(model=model)

  #
  if(params.refine.dry_run): return
  #
  monitor = results.manager(
    model                 = model,
    geometry_rmsd_manager = geometry_rmsd_manager,
    log                   = log)
  monitor.update(fmodel = fmodel, model = model)
  if(params.rst_file is None):
    if(params.output_file_name_prefix is None):
      params.output_file_name_prefix = prefix
    if(os.path.exists(params.output_folder_name) is False):
      os.mkdir(params.output_folder_name)
    params.output_file_name_prefix = os.path.basename(params.output_file_name_prefix)
    params.rst_file = os.path.abspath(params.output_folder_name + "/" + \
      params.output_file_name_prefix + ".rst.pickle")
  if os.path.isfile(params.rst_file):
    os.remove(params.rst_file)
  start_fmodel = fmodel
  start_ph = None # is it used anywhere? I don't see where it is used!

  # XXX This needs to be consolidated !
  altlocs_present = model.altlocs_present()

  hdm = None
  if model.altlocs_present_only_hd():
    print("Using hd_mapper")
    hdm = hd_mapper(model = model)
    altlocs_present = False

  exclude_selection = None
  if params.refine.exclude is not None:
    exclude_selection = model.selection(string = params.refine.exclude)

  if [exclude_selection, hdm].count(None)==0:
    raise Sorry("exclude_selection is incompatible with presence of H/D.")
  if params.cluster.clustering:
    if exclude_selection is not None:
      raise Sorry("exclude_selection is incompatible with clustering=True")

  if hdm is not None:
    print("Using hd_mapper again!")
    model_for_restraints = hdm.get_single_model()
  elif exclude_selection is not None:
    model_for_restraints = model.select(~exclude_selection)
  else:
    model_for_restraints = model

  restraints_manager = create_restraints_manager(
    params=params, model=model_for_restraints, altlocs_present=altlocs_present)

  # XXX This needs to be consolidated !

  if(map_data is not None and params.refine.mode == "refine"):
    model.geometry_statistics(use_hydrogens=False).show()
    show_cc(
      map_data        = map_data,
      xray_structure  = model.get_xray_structure(),
      log             = log)
    O = calculator.sites_real_space(
      model                 = model,
      geometry_rmsd_manager = geometry_rmsd_manager,
      max_bond_rmsd         = params.refine.max_bond_rmsd,
      map_data              = map_data,
      data_weight           = params.refine.data_weight,
      refine_cycles         = params.refine.number_of_refine_cycles,
      skip_weight_search    = params.refine.skip_weight_search,
      stpmax                = params.refine.stpmax,
      gradient_only         = params.refine.gradient_only,
      line_search           = params.refine.line_search,
      restraints_manager    = restraints_manager,
      max_iterations        = params.refine.max_iterations_refine,
      log                   = log)
    model = O.run()
    of = open("real_space_refined.pdb", "w")
    print(model.model_as_pdb(output_cs=True), file=of)
    of.close()
    model.geometry_statistics(use_hydrogens=False).show()
    show_cc(
      map_data=map_data,
      xray_structure=model.get_xray_structure(),
      log=log)
    return
  else:
    calculator_manager = create_calculator(
      fmodel             = start_fmodel,
      model              = model,
      params             = params,
      restraints_manager = restraints_manager,
      hdm                = hdm,
      exclude_selection  = exclude_selection)
    if(params.refine.mode == "refine"):
      #
      # Optimize H
      #
      if fmodel.f_obs().d_min()<1.2:
        import mmtbx.hydrogens
        print("Optimize H with d.o.f.", file=log)
        print ("  start: r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()))
        mmtbx.hydrogens.fit_rotatable2(model=model, fmodel=fmodel)
        print ("  final: r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()))
        monitor.update(fmodel = fmodel, model = model)
      #
      # PRE-OPTIMIZATION (CODITIONAL)
      #
      stats = monitor.get_stats()
      if(stats.rama_z_score().whole.value<-2.5 or
         stats.clash().score>10 or
         fmodel.f_obs().d_min()>3):
        print("Initial geometry regularization:")
        print("  r_work: %6.4f r_free: %6.4f"%(
          fmodel.r_work(), fmodel.r_free()))
        calculator_sites_opt_pre = calculator.sites_opt(
          restraints_manager      = restraints_manager,
          model                   = model,
          dump_gradients          = False,
          max_shift               = params.refine.stpmax,
          shift_eval              = params.refine.shift_evaluation,
          use_callback_after_step = True,
          convergence_threshold   = 0.3)
        driver.opt(
          params     = params,
          calculator = calculator_sites_opt_pre,
          monitor    = monitor)
        fmodel.xray_structure.set_sites_cart(
          sites_cart = monitor.model.get_sites_cart())
        fmodel.update_xray_structure(update_f_calc  = True)
        print("  r_work: %6.4f r_free: %6.4f"%(
          fmodel.r_work(), fmodel.r_free()))
      # PRE-OPTIMIZATION (CODITIONAL) END
      driver.refine(
        params                = params,
        fmodel                = fmodel,
        geometry_rmsd_manager = geometry_rmsd_manager,
        calculator            = calculator_manager,
        monitor               = monitor)
    else:
      driver.opt(
        params     = params,
        model      = model,
        calculator = calculator_manager,
        monitor    = monitor)
    monitor.finalize(
      input_file_name_prefix  = prefix,
      output_file_name_prefix = params.output_file_name_prefix,
      output_folder_name      = params.output_folder_name)
