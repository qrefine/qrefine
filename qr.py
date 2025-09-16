from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
from libtbx import group_args
import iotbx.pdb
import iotbx.pdb.fetch
import os
import sys
from libtbx import easy_pickle
from libtbx.utils import Sorry
from cctbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out

from qrefine import refine


master_phil_str = """
auto_cust = True
  .type = bool
  .help = Allow to overwrite parameters inside the program: some user-defined \
          parameters will be ignored.
max_atoms = 50000
  .type = int
  .help = maximum number of atoms
debug_rss = False
  .type=bool
debug = False
  .type = bool
  .help = flag to control verbosity of output for debugging problematic code.
scattering_table = wk1995  it1992  *n_gaussian electron neutron
  .type = choice
  .help = Choices of scattering table for structure factors calculations
use_reduce = True
  .type = bool
  .help = Use Reduce/ReadySet in completion/capping
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
  maxnum_residues_in_cluster = 15
    .type = int
    .help = maximum number of residues in a cluster
  bcc_threshold = 9
    .type = int
    .help = threshold value for bcc clustering
  select_within_radius = 10
    .type = int
    .help = supersphere expansion radius
  clustering_method = gnc  *bcc
    .type = choice(multi=False)
    .help = type of clustering algorithm
  altloc_method = average *subtract
    .type = choice(multi=False)
    .help = two strategies on how to join energies from multiple energy and \
            gradient calculations are performed for alternate locations.\
            altloc_method=average does not work! (known issue).
  save_clusters = True
    .type = bool
    .help = save currently used fragments and clusters to disk as PDBs.
  bond_with_altloc = True
    .type = bool
    .help = Disable bond_with_altloc check
  re_calculate_rmsd_tolerance = 0.5
    .type = float
    .help = Re-calculate clusters once the model shifted by more than \
            re_calculate_rmsd_tolerance from initial
}

restraints = cctbx *qm
  .type = choice(multi=False)
  .help = Choice of restraints: cctbx (fast) or any of available QM engines.
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
  refinement_target_name = *ml ls_wunit_k1
    .type = choice
    .help = Data target type: least-squares or maximum-likelihood.
  mode = opt *refine gtest
    .type = choice(multi=False)
    .help = choose between refinement, geometry optimization or gradient test
  number_of_macro_cycles=1
    .type = int
    .help = number of macro cycles used in the refinement procedure
  number_of_weight_search_cycles=50
    .type = int
    .help = Number of attempts to find optimal weight
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
  stop_one_found_first_good_weight = True
    .type = bool
    .help = Set True for low-res (or if use regularization), and False for high-res
  max_r_work_r_free_gap = 5.0
    .type = float
    .help = Max difference between Rfree and Rwork to accept optimized weight
  r_tolerance = 0.001
    .type = float
    .help = R factor delta beteen refinement cycles to determine convergence
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
  minimizer = *lbfgs lbfgsb
    .type = choice(multi=False)
    .help = Choice of minimizer to use. LBFGS-B requires target and gradients.\
            LBFGS can use either gradients only or both, target and gradients.
  shift_evaluation = max *mean
    .type = choice(multi=False)
    .help = Use mean or max coordinate shift to decide on convergence
  exclude = None
    .type = str
    .help = Atom selection string to select atoms to be exluded from refinement
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
rst_file = None
  .type = str
  .help = Restart file to use for determining location in run. Loads previous \
          results of weight calculations.
dump_gradients=None
  .type = str
  .help = used for debugging gradients when clustering.
"""

def get_default_params():
  return iotbx.phil.parse(
    input_string = master_phil_str,
    process_includes = True
  ).extract()

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
    #
    assert [self.params.expansion,
            self.params.cluster.clustering].count(True) != 2
    #
    if self.params.refine.minimizer == "lbfgsb":
      if self.params.refine.gradient_only:
        raise Sorry("gradient_only must be False for lbfgsb.")
    #
    if self.params.expansion or self.params.cluster.clustering:
      self.params.refine.minimizer="lbfgs"
      self.params.refine.gradient_only=True
      print("\n  expansion/clustering require minimizer=lbfgs and gradient_only=True \n",
        file=self.logger)

  def auto_cust(self):
    #
    # XXX Temporary work-around. Future: use the same machinery as phenix.refine
    #
    if(not self.params.auto_cust): return
    #
    # AIMNet2 specific settings (as used in tests for the paper).
    #
    if(self.params.refine.mode=="refine" and
       self.params.quantum.engine_name=="aimnet2"):
      msg="""
The following settings have been auto-set to match refine.mode=refine
and quantum.engine_name=aimnet2:
"""
      self.params.restraints="qm"
      self.params.cluster.clustering=False
      self.params.refine.number_of_weight_search_cycles=10
      self.params.refine.number_of_refine_cycles=5
      self.params.use_reduce=False
      self.params.cluster.select_within_radius=7
      self.params.refine.max_iterations_weight=100
      if self.fmodel is not None:
        self.params.expansion=True
        self.params.refine.minimizer="lbfgs"
        self.params.refine.gradient_only=True

        d_min = self.fmodel.f_obs().d_min()

        if(d_min<1.2):
          self.params.refine.max_bond_rmsd=0.025
          self.params.refine.max_angle_rmsd=2.5
        elif(d_min>=1.2 and d_min<3):
          self.params.refine.max_bond_rmsd=0.02
          self.params.refine.max_angle_rmsd=2.0
        else:
          self.params.refine.max_bond_rmsd=0.01
          self.params.refine.max_angle_rmsd=1.7
      elif self.map_data is not None:
        self.params.refine.max_bond_rmsd=0.009
        self.params.expansion=False
        self.params.refine.minimizer="lbfgsb"
        self.params.refine.gradient_only=False
      else: assert 0
      self.params.refine.stop_one_found_first_good_weight=False
      print(msg, file=self.logger)
      print("  restraints                              ", self.params.restraints, file=self.logger)
      print("  cluster.clustering                      ", self.params.cluster.clustering, file=self.logger)
      print("  refine.minimizer                        ", self.params.refine.minimizer, file=self.logger)
      print("  refine.number_of_weight_search_cycles   ", self.params.refine.number_of_weight_search_cycles, file=self.logger)
      print("  refine.number_of_refine_cycles          ", self.params.refine.number_of_refine_cycles, file=self.logger)
      print("  expansion                               ", self.params.expansion, file=self.logger)
      print("  use_reduce                              ", self.params.use_reduce, file=self.logger)
      print("  cluster.select_within_radius            ", self.params.cluster.select_within_radius, file=self.logger)
      print("  refine.max_iterations_weight            ", self.params.refine.max_iterations_weight, file=self.logger)
      print("  refine.max_bond_rmsd                    ", self.params.refine.max_bond_rmsd, file=self.logger)
      print("  refine.max_angle_rmsd                   ", self.params.refine.max_angle_rmsd, file=self.logger)
      print("  refine.stop_one_found_first_good_weight ", self.params.refine.stop_one_found_first_good_weight, file=self.logger)
      print("  refine.gradient_only                    ", self.params.refine.gradient_only, file=self.logger)
      #
      print("\n***", file=self.logger)
      print("To disable auto-custing of parameters use auto_cust=False", file=self.logger)
      print("***\n", file=self.logger)
      if(self.model.altlocs_present() and not
         self.model.altlocs_present_only_hd()):
        raise Sorry("Alternative conformations are not supported with AIMNet2.")

  def run(self):
    self.header("Refinement start")
    self.data_manager.update_all_defaults(
      data_type = self.params.scattering_table)
    # fmodel stuff
    self.fmodel=None
    if(self.has_ma):
      self.header("Extracting fmodel")
      self.fmodel = self.data_manager.get_fmodel(
        scattering_table = self.params.scattering_table)
      self.fmodel.update_all_scales(log=self.logger, refine_hd_scattering=False)
      self.fmodel.show(log=self.logger, show_header=False)
      self.header("Starting r-factors:")
      print("r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.logger)
      self.header("Setting refinement target")
      print(self.params.refine.refinement_target_name, file=self.logger)
      self.fmodel.set_target_name(
        target_name=self.params.refine.refinement_target_name)
    # real map
    self.map_data = None
    if(self.has_map):
      self.map_data = self.data_manager.get_real_map().map_data()
      self.map_data = self.map_data - flex.mean(self.map_data)
      sd = self.map_data.sample_standard_deviation()
      if sd is not None and sd != 0:
        self.map_data = self.map_data/sd
      if(self.map_data is not None):
        self.params.refine.mode="refine"
    #
    if(not self.has_data):
      self.params.refine.mode="opt"
    # model stuff
    self._print(self.get_program_phil_str())
    self.header("Extracting model")
    model_file_name = self.data_manager.get_model_names()[0]
    model_names = self.data_manager.get_model_names()
    self.model = self.data_manager.get_model(filename = model_file_name)
    if(self.params.scattering_table == "electron"):
      self.model.neutralize_scatterers()
    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
    params.pdb_interpretation.use_neutron_distances = True
    params.pdb_interpretation.restraints_library.cdl = False
    params.pdb_interpretation.sort_atoms = False
    self.model.process(make_restraints=True, grm_normalization=False,
      pdb_interpretation_params = params)
    self.model.setup_scattering_dictionaries(
      scattering_table  = self.params.scattering_table,
      d_min             = 1.0,
      log               = self.logger)
    # Auto-cust (set) parameters
    #
    self.auto_cust()
    #
    # run qrefine
    refine.run(
      model    = self.model,
      fmodel   = self.fmodel,
      map_data = self.map_data,
      params   = self.params,
      rst_file = self.params.rst_file,
      prefix   = model_file_name[:-4],
      log      = self.logger)
