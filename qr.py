from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
import iotbx.phil
import os
from libtbx.utils import null_out
from cctbx.array_family import flex
import mmtbx.model
from libtbx import group_args

from qrefine import conformers, api


master_phil_str = """
auto_cust = True
  .type = bool
  .help = Allow to overwrite parameters inside the program: some user-defined \
          parameters will be ignored.
debug = False
  .type = bool
  .help = flag to control verbosity of output for debugging problematic code.
scattering_table = wk1995  it1992  *n_gaussian electron neutron
  .type = choice
  .help = Choices of scattering table for structure factors calculations
cluster{
  charge_cutoff = 8.0
    .type = float
    .help = distance for point charge cutoff
  clustering = True
    .type = bool
    .help = enable/disable clustering
  charge_embedding = False
    .type = bool
    .help = point charge embedding
  two_buffers = False
    .type = bool
    .help = two-buffers are used when gradients for the whole system do not match the joined gradients from fragments \
            when only as single buffer was used to surround each clusters.
  distance_buffer = False
    .type = bool
    .help = In addition to the interaction-based buffer, add to each fragment \
            every residue of the expanded system that has an atom within \
            distance_buffer_radius of the cluster. Captures neighbors that \
            are too far for the interaction criterion but still affect the \
            gradients (e.g. charged or polar groups 4-6 A away).
  distance_buffer_radius = 5.0
    .type = float
    .help = Radius (Angstrom) for distance_buffer. Larger values improve \
            agreement with the non-clustered gradients at the cost of larger \
            fragments.
  maxnum_residues_in_cluster = 15
    .type = int
    .help = maximum number of residues in a cluster
  bcc_threshold = 9
    .type = int
    .help = threshold value for bcc clustering
  select_within_radius = 10
    .type = int
    .help = supersphere expansion radius
  clustering_method = *bcc
    .type = choice(multi=False)
    .help = clustering algorithm (betweenness centrality clustering)
  save_clusters = True
    .type = bool
    .help = save currently used fragments and clusters to disk as PDBs.
  re_calculate_rmsd_tolerance = 0.5
    .type = float
    .help = Re-calculate clusters once the model shifted by more than \
            re_calculate_rmsd_tolerance from initial
}

restraints = cctbx *qm
  .type = choice(multi=False)
  .help = Choice of restraints: cctbx (fast) or any of available QM engines.
altlocs {
  handling = *conformers native
    .type = choice(multi=False)
    .help = conformers: split the model into single-conformer models, compute \
            gradients for each and combine them (works with any restraints \
            source). native: hand the model with altlocs to the restraints \
            source as is. Only restraints=cctbx can do this.
  reconcile = True
    .type = bool
    .help = Relabel altlocs of interacting residues so that the same label \
            means the same conformer (iotbx.pdb.hierarchy.reconcile_altlocs) \
            before restraints are built.
}
expansion = False
  .type = bool
  .help = Expand input model into super-sphere
quantum {
  engine_name = *mopac fairchem aimnet2 aimnet2-old torchani terachem turbomole pyscf orca gaussian xtb server
    .type = choice(multi=False)
    .help = choose the QM program
  basis = Auto
    .type = str
    .help = pre-defined defaults
  charge= None
    .type = int
    .help = The formal charge of the entire molecule
  method = Auto
    .type = str
    .help = Defaults to HF for all but MOPAC (PM7), xTB (GFN2) and TorchANI (ani-1x_8x)
  memory = None
    .type = str
    .help = memory for the QM program
  nproc = 1
    .type = int
    .help = number of parallel processes for the QM program
  server_url =  http://127.0.0.1:8000
    .type = str
    .help = address (http://address:port) of the server API if engine_name=server
  qm_addon = gcp dftd3 gcp-d3
    .type = choice(multi=False)
    .help = allows additional calculations of the gCP and/or DFT-D3 corrections using their stand-alone programs
  qm_addon_method = None
    .type = str
    .help = specifies flags for the qm_addon. See manual for details.
}

refine {
  dry_run=False
    .type = bool
    .help = do not perform calculations, only setup steps
  sf_algorithm = *direct fft
    .type = choice(multi=False)
    .help = algorithm used to compute structure factors, either the direct method or fast fourier transform.
  refinement_target_name = *Auto ml mlhl ls_wunit_k1
    .type = choice
    .help = Data target type: least-squares or maximum-likelihood. Auto: ml, \
            or mlhl when the data carry experimental phases.
  mode = opt *refine gtest
    .type = choice(multi=False)
    .help = choose between refinement, geometry optimization or gradient test
  number_of_macro_cycles=1
    .type = int
    .help = number of macro cycles used in the refinement procedure
  number_of_weight_search_cycles=12
    .type = int
    .help = Maximum number of weight search trials. In reciprocal space the \
            search usually stops earlier, at the R-free plateau or the \
            flip. The real-space weight search uses the same cap.
  number_of_refine_cycles=5
    .type = int
    .help = maximum number of refinement cycles
  number_of_micro_cycles=50
    .type = int
    .help = maximum number of micro cycles used in refinement
  data_weight=None
    .type = float
    .help = Allows to specify data/restraints weight (by-pass weight search)
  skip_weight_search = False
    .type = bool
    .help = Calculate weight based on gradient norms (skip weight optimization)
  adjust_restraints_weight_scale_value = 2
    .type = float
    .help = Factor by which the restraints weight scale changes between \
            weight search trials
  max_iterations_weight = 50
    .type = int
    .help = Max number of trial refinement iterations for weight search
  max_iterations_refine = 50
    .type = int
    .help = Max number of actual refinement iterations
  use_ase_lbfgs = False
    .type = bool
    .help = used for debugging the lbfgs minimizer from cctbx.
  line_search = True
    .type = bool
    .help = flag to use a line search in minimizer.
  stpmax = 0.2
    .type = float
    .help = maximum step length, empirically we find 3 for cctbx, but 0.2 is better for QM methods.
  gradient_only = False
    .type = bool
    .help = use the gradient only line search according to JA Snyman 2005.
  refine_sites = True
    .type = bool
    .help = only refine the cartesian coordinates of the molecular system.
  refine_adp = False
    .type = bool
    .help = adp refinement are not currently supported.
  restraints_weight_scale = 1.0
    .type = float
    .help = Scale factor for restraints term
  shake_sites = False
    .type = bool
    .help = Randomize coordinates prior weight calculation
  use_convergence_test = True
    .type = bool
    .help = Check if refinement converged (for earlier termination)
  max_bond_rmsd = 0.03
    .type = float
    .help = Max bond RMSD to accept optimized weight
  max_angle_rmsd = 1.7
    .type = float
    .help = Max angle RMSD to accept optimized weight
  max_r_work_r_free_gap = 5.0
    .type = float
    .help = Max difference between Rfree and Rwork to accept optimized weight
  r_tolerance = 0.001
    .type = float
    .help = R factor delta between refinement cycles to determine convergence, \
            in the weight search, R-free changes within it count as a plateau
  rmsd_tolerance = 0.01
    .type = float
    .help = maximum acceptable tolerance for the rmsd.
  opt_log = False
    .type = bool
    .help = additional output of the L-BFGS optimizer
  pre_opt = False
    .type = bool
    .help = pre-optimization using steepest decent (SD) and conjugate gradient (CG) techniques w/o line search
  pre_opt_stpmax = 0.1
    .type = float
    .help = step size
  pre_opt_iter= 10
    .type = int
    .help = max. iterations for pre-optimizer
  pre_opt_switch = 2
    .type = int
    .help = max. iterations before switching from SD to CG
  pre_opt_gconv = 3000
    .type = float
    .help = gradient norm convergence threshold for pre-optimizer
  minimizer = lbfgs *lbfgsb
    .type = choice(multi=False)
    .help = Choice of minimizer to use. LBFGS-B requires target and gradients.\
            LBFGS can use either gradients only or both, target and gradients.
  shift_evaluation = max *mean
    .type = choice(multi=False)
    .help = Use mean or max coordinate shift to decide on convergence
  exclude = None
    .type = str
    .help = Atom selection string to select atoms to be exluded from refinement
  freeze {
    selection = None
      .type = str
    r_min = 4
      .type = float
    level = 0 1 2 *3
      .type=choice(multi=False)
    passes = 1
      .type = int
  }
}

parallel {
  method = *multiprocessing slurm pbs sge lsf threading
    .type = choice(multi=False)
    .help = type of parallel mode and efficient method of processes on the \
            current computer. The others are queueing protocols with the \
            expection of threading which is not a safe choice.
  nproc = 1
    .type = int
    .help = Number of processes to use
  qsub_command = None
    .type = str
    .help = Specific command to use on the queue system
}

output_file_name_prefix = None
  .type = str
  .help = Output file name prefix
output_folder_name = "pdb"
  .type = str
  .help = Output folder name
write_nuclear_h_positions = False
  .type = bool
  .help = Write the refined nuclear H positions. Default: H positions follow \
          the scattering-table convention, X-ray X-H lengths for X-ray \
          tables, nuclear for neutron and electron.
dump_gradients=None
  .type = str
  .help = used for debugging gradients when clustering.
"""

def get_default_params():
  return iotbx.phil.parse(
    input_string = master_phil_str,
    process_includes = True
  ).extract()

def _header(log, text):
  """The section header of ProgramTemplate.header, printed into log."""
  print("-"*79, file=log)
  print(text, file=log)
  print("*"*len(text), file=log)

def prepare_inputs(data_manager, params, log):
  """
  The objects qr.refine refines, from the inputs of data_manager (one
  model; at most one reflection file or one map): the fmodel (scaled, its
  structure a copy of the model's, in the model's atom order), the
  normalized map, the model (altlocs reconciled, electron table
  neutralised, processed with neutron distances, CDL off, input atom order,
  scattering dictionaries at the data or map resolution). Prints into log.
  Returns group_args(model, fmodel, map_data, d_min, prefix): fmodel and
  map_data None when absent; d_min the resolution of the reflection data
  (None without them), as api.finalize_params takes it; prefix the model
  file name without extension.
  """
  has_ma = data_manager.has_miller_arrays(
    expected_n=1,
    exact_count=True,
    raise_sorry=False)
  has_map = data_manager.has_real_maps(
    expected_n=1,
    exact_count=True,
    raise_sorry=False)
  _header(log, "Refinement start")
  data_manager.update_all_defaults(
    data_type = params.scattering_table)
  # fmodel stuff
  fmodel=None
  if(has_ma):
    _header(log, "Extracting fmodel")
    fmodel = data_manager.get_fmodel(
      scattering_table = params.scattering_table)
  # real map
  map_data = None
  if(has_map):
    map_data = data_manager.get_real_map().map_data()
    map_data = map_data - flex.mean(map_data)
    sd = map_data.sample_standard_deviation()
    if sd is not None and sd != 0:
      map_data = map_data/sd
  #
  # model stuff
  _header(log, "Extracting model")
  model_file_name = data_manager.get_model_names()[0]
  model = data_manager.get_model(filename = model_file_name)
  if(params.altlocs.reconcile and
     conformers.n_conformers(model.get_hierarchy()) > 1):
    _header(log, "Reconciling altlocs")
    # Relabel a fresh hierarchy in input-file order: the hierarchy of
    # model is built with sort_atoms=True, and process() below rebuilds
    # the hierarchy from the input (sort_atoms=False) anyway.
    model_input = model.get_model_input()
    if(model_input is not None):
      hierarchy = model_input.construct_hierarchy(
        sort_atoms=False).deep_copy()
    else:
      hierarchy = model.get_hierarchy().deep_copy()
    r = hierarchy.reconcile_altlocs(
      crystal_symmetry = model.crystal_symmetry(), log = log)
    print("  %d residues relabeled" % len(r.swaps), file=log)
    if(len(r.swaps) > 0):
      # Model without input: process() would rebuild the hierarchy from the
      # input, undoing the relabel. The HELIX/SHEET annotation is carried
      # over; links are found automatically by pdb_interpretation (it never
      # reads the file's LINK/SSBOND records); TER positions are lost.
      hierarchy.atoms().reset_i_seq()
      hierarchy.atoms_reset_serial()
      ss_annotation = model.get_ss_annotation()
      model = mmtbx.model.manager(
        model_input       = None,
        pdb_hierarchy     = hierarchy,
        crystal_symmetry  = model.crystal_symmetry(),
        restraint_objects = model.get_restraint_objects(),
        log               = null_out())
      if(ss_annotation is not None):
        model.set_ss_annotation(ann = ss_annotation)
      print("  model rebuilt from its hierarchy (secondary structure "
        "annotation carried over, links found automatically)", file=log)
    if(r.unresolved):
      print("  WARNING: altlocs of %d groups of residues cannot be made "
        "consistent; refinement continues with the labels as they are." %
        len(r.unresolved), file=log)
  if(params.scattering_table == "electron"):
    model.neutralize_scatterers()
  pi_params = api.pdb_interpretation_params()
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = pi_params)
  # Resolution of the dictionaries as in phenix.refine: of the data, of the
  # map, else 1 A.
  d_min = api.dictionaries_d_min(model=model, fmodel=fmodel,
    map_data=map_data)
  model.setup_scattering_dictionaries(
    scattering_table  = params.scattering_table,
    d_min             = d_min,
    log               = log)
  # fmodel holds the structure of the data manager's (sorted) model; the
  # refined model may differ in atom order (altloc reconciliation rebuilds
  # it from the input-order hierarchy). Give fmodel a copy of the model's
  # structure (never the same object: the X-H mapper keeps nuclear and
  # scattering sites apart), then scale.
  if(fmodel is not None):
    fmodel.update_xray_structure(
      xray_structure = model.get_xray_structure().deep_copy_scatterers(),
      update_f_calc  = True,
      update_f_mask  = True)
    fmodel.update_all_scales(log=log, refine_hd_scattering=False)
    fmodel.show(log=log, show_header=False)
    _header(log, "Starting r-factors:")
    print("r_work=%6.4f r_free=%6.4f"%(
      fmodel.r_work(), fmodel.r_free()), file=log)
  return group_args(
    model    = model,
    fmodel   = fmodel,
    map_data = map_data,
    d_min    = d_min if fmodel is not None else None,
    prefix   = os.path.splitext(model_file_name)[0])

def write_output(result, params, prefix, log):
  """
  Write the refined model of result (api.refine) and return the file name:
  pdb/<prefix>_refined.<ext> (folder params.output_folder_name, prefix
  params.output_file_name_prefix if set; directory part dropped) with the
  sites result.sites_cart, or result.sites_cart_nuclear with
  write_nuclear_h_positions=True; real-space mode writes
  real_space_refined.<ext> into the current directory. <ext> is .pdb, or
  .cif when the identifiers do not fit the PDB format. result.model is
  left with its (nuclear) sites.
  """
  model = result.model
  ext = ".pdb" if model.can_be_output_as_pdb() else ".cif"
  if(result.mode == "real_space"):
    file_name = "real_space_refined" + ext
    with open(file_name, "w") as of:
      print(model.as_pdb_or_mmcif_string(target_format="pdb"), file=of)
    return file_name
  output_file_name_prefix = os.path.basename(
    prefix if params.output_file_name_prefix is None
    else params.output_file_name_prefix)
  base_name = "%s_refined%s" % (output_file_name_prefix, ext)
  folder = params.output_folder_name
  os.makedirs(folder, exist_ok=True)
  file_name = os.path.join(folder, base_name)
  sites_cart = result.sites_cart
  if(params.write_nuclear_h_positions):
    sites_cart = result.sites_cart_nuclear
  nuclear = model.get_sites_cart()
  model.set_sites_cart(sites_cart)
  try:
    with open(file_name, "w") as fo:
      fo.write(model.as_pdb_or_mmcif_string(target_format="pdb"))
  finally:
    model.set_sites_cart(nuclear)
  print("See %s in %s folder." % (base_name, folder), file=log)
  return file_name

class Program(ProgramTemplate):

  description = """
qr.refine is an open-source module that carries out refinement of
bio-macromolecules utilizing chemical restraints QM calculations.

Example:
qr.refine model.pdb model.mtz [<param_name>=<param_value>] ...
"""

  datatypes = ['model', 'phil', 'miller_array', 'restraint', 'real_map']

  master_phil_str = master_phil_str

  def validate(self):
    print('Validate inputs:', file=self.logger)
    self.data_manager.has_models(
      expected_n=1,
      exact_count=True,
      raise_sorry=True)
    self.has_ma = self.data_manager.has_miller_arrays(
      expected_n=1,
      exact_count=True,
      raise_sorry=False)
    self.has_map = self.data_manager.has_real_maps(
      expected_n=1,
      exact_count=True,
      raise_sorry=False)
    self.has_data = self.has_ma or self.has_map

  def run(self):
    inputs = prepare_inputs(data_manager = self.data_manager,
      params = self.params, log = self.logger)
    has_map = inputs.map_data is not None
    # All parameter rules, once.
    api.finalize_params(
      params           = self.params,
      d_min            = inputs.d_min,
      scattering_table = self.params.scattering_table,
      has_map          = has_map,
      has_data         = inputs.fmodel is not None or has_map,
      log              = self.logger,
      n_conformers     = conformers.n_conformers(inputs.model.get_hierarchy()))
    # Fail fast: engine availability; model readiness (altloc completion,
    # which is the job of mmtbx.reduce2 qmready=true; caps, H, elements,
    # crystal neighbours with restraints=qm).
    api.check_engine(self.params)
    api.check_model(inputs.model, self.params, has_map = has_map)
    self._print(self.get_program_phil_str())
    if(inputs.fmodel is not None):
      self.header("Setting refinement target")
      target_name = api.resolve_target_name(self.params, inputs.fmodel)
      print(target_name, file=self.logger)
      inputs.fmodel.set_target_name(target_name=target_name)
    #
    # run qrefine
    result = api.refine(
      model               = inputs.model,
      params              = self.params,
      fmodel              = inputs.fmodel,
      map_data            = inputs.map_data,
      working_folder      = "ase",
      intermediate_folder = self.params.output_folder_name,
      log                 = self.logger)
    if(result is None): return # refine.dry_run
    write_output(result = result, params = self.params,
      prefix = inputs.prefix, log = self.logger)
