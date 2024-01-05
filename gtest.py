import sys, time, os
from qrefine import __version__
import libtbx
import mmtbx
from qrefine.utils import hierarchy_utils
from iotbx.cli_parser import run_program
from qrefine import gtest
from qrefine import qr
from libtbx.program_template import ProgramTemplate
from .refine import set_qm_defaults

def get_default_params():
    return iotbx.phil.parse(input_string=master_phil_str, process_includes=True).extract()


local_phil = """
    gtest{
    g_scan  = 10 15 20
        .type = str
        .help = sequence of numbers specifying maxnum_residues_in_cluster for gradient convergence test (mode=gtest), treat as string!
      g_mode = None
        .type = int
        .help = manual control over gradient test loops (0=expansion ,1=standard, 2=standard + point-charges, 3=two buffer 4=two_buffer+ point-charges)
        }
    """

def run_g_test(params, model, weights, start_fmodel, log):
  import numpy as np
  import time
  from .fragment import write_cluster_and_fragments_pdbs
  from .refine import create_fragment_manager, create_restraints_manager,create_calculator

  # determine what kind of buffer to calculate
  g_mode=[]
  if params.gtest.g_mode is None:
    g_mode.append(1)
    if params.cluster.charge_embedding:
      g_mode.append(2)
    if params.cluster.two_buffers:
      g_mode.append(3)
    if params.cluster.two_buffers and params.cluster.charge_embedding:
      g_mode.append(4)
  else:
    g_mode.append(params.gtest.g_mode)


  # reset flags
  params.cluster.clustering=True
  params.cluster.save_clusters=True
  params.cluster.charge_embedding=False
  params.cluster.two_buffers=False
  grad=[]
  idx=0
  idl=[]

  # input for cluster size
  cluster_scan=sorted([int(x) for x in params.gtest.g_scan.split()])

  if g_mode[0]==0:
    print('warning: supersphere calculation!', file=log)
    params.cluster.clustering=False
    params.expansion=True
    cluster_scan=[0]
    clusters=[]

  n_grad=len(cluster_scan)*len(g_mode)
  print('Calculating %3i gradients \n' % (n_grad), file=log)
  print('Starting loop over different fragment sizes', file=log)
  for ig in g_mode:

    print('loop for g_mode = %i ' % (ig), file=log)
    if ig == 2:
      print('pc on', file=log)
      params.cluster.charge_embedding=True
    if ig == 3:
      print('two_buffers on, pc off', file=log)
      params.cluster.charge_embedding=False
      params.cluster.two_buffers=True
    if ig == 4:
      print('two_buffers on, pc on', file=log)
      params.cluster.charge_embedding=True
      params.cluster.two_buffers=True

    for max_cluster in cluster_scan:
      idl.append([ig,max_cluster])
      print('g_mode: %s' % (" - ".join(map(str,idl[idx]))), file=log)
      t0 = time.time()
      print("~max cluster size ",max_cluster, file=log)
      params.cluster.maxnum_residues_in_cluster=max_cluster
      fragment_manager = create_fragment_manager(params = params, model = model)
      restraints_manager = create_restraints_manager(params, model)
      rm = restraints_manager
      if(fragment_manager is not None):
        print(f"time taken for fragments {time.time() - t0:.2f}")
        frags=fragment_manager
        print('~  # clusters  : ',len(frags.clusters), file=log)
        print('~  list of atoms per cluster:', file=log)
        print('~   ',[len(x) for x in frags.cluster_atoms], file=log)
        print('~  list of atoms per fragment:', file=log)
        print('~   ',[len(x) for x in frags.fragment_super_atoms], file=log)

        # save fragment data. below works
        # better way is to make a single PDB file with chain IDs
        label="-".join(map(str,idl[idx]))
        write_cluster_and_fragments_pdbs(
          fragments=frags.get_fragment_extracts(),directory=label)

      calculator_manager = create_calculator(
        fmodel             = start_fmodel,
        model              = model,
        params             = params,
        restraints_manager = rm)
      grad=list(calculator_manager.target_and_gradients())[1]
      print(f"~   gnorm {np.linalg.norm(grad):.4f}", file=log)
      print(f"~   max_g {max(abs(i) for i in grad):.4f}  min_g {min(abs(i) for i in grad):.4f}", file=log)
      name="-".join(map(str,idl[idx]))
      np.save(name,grad)
      idx+=1
      print(f"total time for gradient {time.time() - t0:.2f} \n\n", file=log)

  print('ready to run qr.granalyse!', file=log)


class Program(ProgramTemplate):
    description = """
    qr.gtest is an open-source module that carries out refinement of
    bio-macromolecules utilizing chemical restraints QM calculations.

    Example:
    qr.refine model.pdb model.mtz [<param_name>=<param_value>] ...
    """

    datatypes = ["model", "phil"]

    master_phil_str = qr.master_phil_str + local_phil


    def validate(self):
        self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
        self.has_data = False

    def run(self):
        set_qm_defaults(self.params, self.logger)
        self._print(self.get_program_phil_str())
        model_file_name = self.data_manager.get_model_names()[0]
        model_names = self.data_manager.get_model_names()
        self.model = self.data_manager.get_model(filename=model_file_name)
        self.params.refine.mode = "gtest"
        run_g_test(params=self.params, model=self.model, weights=None, start_fmodel=None, log=self.logger)
