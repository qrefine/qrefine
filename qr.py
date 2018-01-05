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
   results_manager (store all reportable infomation, and write it out as a log,
   and also write our final pdb structure.)
   """
from __future__ import division

import os
import sys
import time
import pickle
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
from fragment import fragments
import calculator
import driver
import restraints
import cluster_restraints
import results
from qrefine.super_cell import expand
import mmtbx.model.statistics

master_params_str ="""

max_atoms = 15000
  .type = int
debug = False
    .type = bool
cluster{
  charge_cutoff = 8.0
    .type = float
  clustering = False
    .type = bool
  charge_embedding = False
    .type = bool
  two_buffers = False
    .type = bool
  maxnum_residues_in_cluster = 15
    .type = int
  clustering_method = gnc  *bcc
    .type = choice(multi=False)
  altloc_method = *subtract average
    .type = choice(multi=False)
}

restraints = cctbx *qm
.type = choice(multi=False)
qm_engine_name = *mopac terachem turbomole pyscf orca gaussian
.type = choice(multi=False)
charge= None
.type = int
basis = "sto-3g"
.type = str

refine {
  dry_run=False
    .type = bool
  sf_algorithm = *direct fft
    .type = choice(multi=False)
  refinement_target_name = *ml ls_wunit_k1
    .type = choice
  mode = opt *refine
    .type = choice(multi=False)
  number_of_macro_cycles=1
    .type = int
  number_of_weight_search_cycles=50
    .type = int
  number_of_micro_cycles=50
    .type = int
  data_weight=None
    .type = float
  max_iterations = 50
    .type = int
  line_search = True
    .type = bool
  stpmax = 1.e9
    .type = float
  gradient_only = False
    .type = bool
  update_all_scales = True
    .type = bool
  refine_sites = True
    .type = bool
  refine_adp = False
    .type = bool
  restraints_weight_scale = 1.0
    .type = float
  shake_sites = False
    .type = bool
  use_convergence_test = True
    .type = bool
  max_bond_rmsd = 0.03
    .type = float
  max_r_work_r_free_gap = 5.0
    .type = float
  r_tolerance = 0.001
    .type = float
  rmsd_tolerance = 0.01
    .type = float
}

parallel_params {
method = *multiprocessing slurm pbs sge lsf threading
.type = choice(multi=False)
nproc = None
.type = int
qsub_command = None
.type = str
}

output_file_name_prefix = None
.type = str
output_folder_name = "pdb"
.type = str
shared_disk = True
.type = bool
rst_file = None
.type = str

dump_gradients=None
  .type = str
"""

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=master_params_str)

def create_fmodel(cmdline, log):
  fmodel = mmtbx.f_model.manager(
    f_obs          = cmdline.f_obs,
    r_free_flags   = cmdline.r_free_flags,
    xray_structure = cmdline.xray_structure,
    target_name    = cmdline.params.refine.refinement_target_name)
  if(cmdline.params.refine.update_all_scales):
    fmodel.update_all_scales(remove_outliers=False)
    fmodel.show(show_header=False, show_approx=False)
  print >> log, "Initial r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), fmodel.r_free())
  log.flush()
  return fmodel

#def process_model_file(pdb_file_name, cif_objects, crystal_symmetry):
#  params = monomer_library.pdb_interpretation.master_params.extract()
#  params.use_neutron_distances = True
#  params.restraints_library.cdl = False
#  params.sort_atoms = False
#  processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
#    pdb_interpretation_params = params,
#    stop_for_unknowns         = True,
#    log                       = null_out(),
#    crystal_symmetry          = crystal_symmetry,
#    cif_objects               = cif_objects,
#    use_neutron_distances     = True)
#  processed_pdb_file, junk = processed_pdb_files_srv.\
#    process_pdb_files(pdb_file_names = [pdb_file_name])
#  xray_structure = processed_pdb_file.xray_structure()
#  sctr_keys = \
#    xray_structure.scattering_type_registry().type_count_dict().keys()
#  has_hd = "H" in sctr_keys or "D" in sctr_keys
#  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
#  #super_cell = expand( # XXX needs to be optional?
#  #  pdb_hierarchy        = pdb_hierarchy,
#  #  crystal_symmetry     = crystal_symmetry,
#  #   select_within_radius = 15) # XXX needs to be a parameter?
#  return group_args(
#    processed_pdb_file = processed_pdb_file,
#    pdb_hierarchy      = pdb_hierarchy,
#    xray_structure     = xray_structure,
#    #super_cell         = super_cell,
#    has_hd             = has_hd)

def process_model_file(pdb_file_name, cif_objects, crystal_symmetry):
  import iotbx.pdb
  import mmtbx
  import mmtbx.monomer_library.server
  import mmtbx.monomer_library.pdb_interpretation
  import mmtbx.restraints
  from mmtbx import monomer_library

  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib    = mmtbx.monomer_library.server.ener_lib()
  params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  params.sort_atoms = False

  pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
  ph = iotbx.pdb.input(file_name=pdb_file_name).construct_hierarchy()
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
       mon_lib_srv              = mon_lib_srv,
       ener_lib                 = ener_lib,
       params                   = params,
       pdb_hierarchy            = ph,
       strict_conflict_handling = False,
       crystal_symmetry         = pdb_inp.crystal_symmetry(), #XXX I hope this is correct
       force_symmetry           = True,
       log                      = null_out())
  xrs = processed_pdb_file.xray_structure()
  sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  return group_args(
      processed_pdb_file = processed_pdb_file,
      pdb_hierarchy      = ph,
       xray_structure    = xrs,
      has_hd             = has_hd)

def create_fragment_manager(
      cif_objects,
      pdb_hierarchy,
      crystal_symmetry,
      params,
      file_name      = os.path.join("ase","tmp_ase.pdb")):
  if(not params.cluster.clustering): return None
  return fragments(
    cif_objects                = cif_objects,
    working_folder             = os.path.split(file_name)[0]+ "/",
    clustering_method          = params.cluster.clustering_method,
    altloc_method              = params.cluster.altloc_method,
    maxnum_residues_in_cluster = params.cluster.maxnum_residues_in_cluster,
    charge_embedding           = params.cluster.charge_embedding,
    two_buffers                = params.cluster.two_buffers,
    pdb_hierarchy              = pdb_hierarchy,
    qm_engine_name             = params.qm_engine_name,
    crystal_symmetry           = crystal_symmetry,
    debug                      = params.debug,
    charge_cutoff              = params.cluster.charge_cutoff)

def create_restraints_manager(
      cmdline,
      model,
      fragment_manager=None):
  if(cmdline.params.restraints == "cctbx"):
    assert model.processed_pdb_file is not None
    restraints_manager = restraints.from_cctbx(
      processed_pdb_file = model.processed_pdb_file,
      has_hd             = model.has_hd)
  else:
    assert cmdline.params.restraints == "qm"
    restraints_manager = restraints.from_qm(
      cif_objects                = cmdline.cif_objects,
      basis                      = cmdline.params.basis,
      pdb_hierarchy              = model.pdb_hierarchy,
      charge                     = cmdline.params.charge,
      qm_engine_name             = cmdline.params.qm_engine_name,
      crystal_symmetry           = cmdline.crystal_symmetry,
      clustering                 = cmdline.params.cluster.clustering)
  return restraints_manager

def create_calculator(weights, fmodel, params, restraints_manager):
  if(weights is None):
    weights = calculator.weights(
      shake_sites             = params.refine.shake_sites ,
      restraints_weight       = 1.0,
      data_weight             = params.refine.data_weight,
      restraints_weight_scale = params.refine.restraints_weight_scale)
  if(params.refine.refine_sites):
    if(params.refine.mode == "refine"):
      return calculator.sites(
        fmodel             = fmodel,
        restraints_manager = restraints_manager,
        weights            = weights,
        dump_gradients     = params.dump_gradients)
    else:
      return calculator.sites_opt(
        restraints_manager = restraints_manager,
        fmodel             = fmodel,
        dump_gradients     = params.dump_gradients)
  if(params.refine.refine_adp):
    return calculator.adp(
      fmodel             = fmodel,
      restraints_manager = restraints_manager,
      weights            = weights)

def run(cmdline, log):
  params = cmdline.params
  model = process_model_file(
    pdb_file_name    = cmdline.pdb_file_names[0],
    cif_objects      = cmdline.cif_objects,
    crystal_symmetry = cmdline.crystal_symmetry)
  if(params.cluster.clustering):
    params.refine.gradient_only = True
    print >> log, " params.gradient_only", params.refine.gradient_only
  # RESTART
  if(os.path.isfile(str(cmdline.params.rst_file))):
    print >> log, "restart info is loaded from %s" % cmdline.params.rst_file
    rst_data = easy_pickle.load(cmdline.params.rst_file)
    fmodel = rst_data["fmodel"]
    results_manager = rst_data["results"]
    results_manager.log = log
    weights = rst_data["weights"]
    geometry_rmsd_manager = rst_data["geometry_rmsd_manager"]
    start_fmodel = rst_data["rst_fmodel"]
    start_ph = cmdline.pdb_hierarchy.deep_copy().adopt_xray_structure(
      start_fmodel.xray_structure)
  else:
    weights = None
    if (cmdline.pdb_hierarchy.atoms().size() > params.max_atoms):
      raise Sorry("Too many atoms.")
    cmdline.working_phil.show(out=log, prefix="   ")
    fmodel = create_fmodel(cmdline=cmdline, log=log)
    geometry_rmsd_manager = restraints.from_cctbx(
      processed_pdb_file = model.processed_pdb_file,
      has_hd             = model.has_hd).geometry_restraints_manager
    cctbx_rm_bonds_rmsd = calculator.get_bonds_angles_rmsd(
      restraints_manager = geometry_rmsd_manager.geometry,
      xrs                = fmodel.xray_structure)
    # Show full model statistics
    if(0):
      # block this:
      mso = mmtbx.model.statistics.geometry(
        pdb_hierarchy               = model.pdb_hierarchy.select(sel),
        geometry_restraints_manager = geometry_rmsd_manager,
        use_hydrogens               = False)
      mso.show(lowercase=True, out=log)
    #
    if(params.refine.dry_run): return
    #
    results_manager = results.manager(
      r_work                  = fmodel.r_work(),
      r_free                  = fmodel.r_free(),
      b                       = cctbx_rm_bonds_rmsd,
      xrs                     = fmodel.xray_structure,
      max_bond_rmsd           = params.refine.max_bond_rmsd,
      max_r_work_r_free_gap   = params.refine.max_r_work_r_free_gap,
      pdb_hierarchy           = cmdline.pdb_hierarchy,
      mode                    = params.refine.mode,
      log                     = log,
      restraints_weight_scale = params.refine.restraints_weight_scale)
    if(params.rst_file is None):
      if(params.output_file_name_prefix is None):
        params.output_file_name_prefix = \
          os.path.basename(cmdline.pdb_file_names[0])[:-4]
      if(os.path.exists(params.output_folder_name) is False):
        os.mkdir(params.output_folder_name)
      params.rst_file = params.output_folder_name + "/" + \
        params.output_file_name_prefix + ".rst.pickle"
    if os.path.isfile(params.rst_file):
      os.remove(params.rst_file)
    print >> log, "\n***********************************************************"
    print >> log, "restart info will be stored in %s" % params.rst_file
    print >> log, "***********************************************************\n"
    start_fmodel = fmodel
    start_ph = None # is it used anywhere? I don't see where it is used!
  fragment_manager = create_fragment_manager(
    params           = params,
    pdb_hierarchy    = model.pdb_hierarchy,
    cif_objects      = cmdline.cif_objects,
    crystal_symmetry = cmdline.crystal_symmetry)
  restraints_manager = create_restraints_manager(
    cmdline          = cmdline,
    model            = model)
  if(fragment_manager is not None):
    cluster_restraints_manager = cluster_restraints.from_cluster(
      restraints_manager = restraints_manager,
      fragment_manager   = fragment_manager,
      parallel_params    = params.parallel_params)
  rm = restraints_manager
  if(fragment_manager is not None):
    rm = cluster_restraints_manager
  calculator_manager = create_calculator(
    weights            = weights,
    fmodel             = start_fmodel,
    params             = params,
    restraints_manager = rm)
  if(params.refine.mode == "refine"):
    driver.refine(
      params                = params,
      fmodel                = fmodel,
      geometry_rmsd_manager = geometry_rmsd_manager,
      calculator            = calculator_manager,
      results               = results_manager)
  else:
    driver.opt(
      params                = params,
      fmodel                = fmodel,
      geometry_rmsd_manager = geometry_rmsd_manager,
      calculator            = calculator_manager,
      results               = results_manager)
  xrs_best = results_manager.finalize(
    input_file_name_prefix  = os.path.basename(cmdline.pdb_file_names[0])[:-4],
    output_file_name_prefix = params.output_file_name_prefix,
    output_folder_name      = params.output_folder_name)

if (__name__ == "__main__"):
  t0 = time.time()
  log = sys.stdout
  args = sys.argv[1:]
  cmdline = mmtbx.command_line.load_model_and_data(
      args          = args,
      master_phil   = get_master_phil(),
      create_fmodel = False,
      out           = log)
  run(cmdline=cmdline, log = log)
  print >> log, "Time: %6.4f"%(time.time()-t0)
