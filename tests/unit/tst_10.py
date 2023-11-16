from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import random
from ase.units import Hartree, Bohr, mol, kcal
import iotbx.pdb
import mmtbx.command_line
import libtbx.load_env
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from qrefine.cluster_restraints import from_cluster
from qrefine.restraints import from_qm, from_cctbx
from qrefine.fragment import fragments
from qrefine.clustering import betweenness_centrality_clustering
from qrefine.tests.unit import run_tests
import mmtbx.model
from qrefine import qr
from qrefine.utils import hierarchy_utils

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

# fix random seed in this script
random.seed(0)
flex.set_random_seed(0)


master_params_str ="""
method = *multiprocessing pbs sge lsf threading
.type = choice(multi=False)
nproc = 1
.type = int
qsub_command = None
.type = str
"""

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=master_params_str)

def get_model():
  file_name = os.path.join(qr_unit_tests,"data_files","helix.pdb")
  pdb_inp = iotbx.pdb.input(file_name)
  model = hierarchy_utils.process_model_file(
    pdb_file_name = file_name,
    cif_objects = None,
    crystal_symmetry=pdb_inp.crystal_symmetry()).model
  return model

def run(prefix):
  """
  Exercise combined energy and gradients from cluster qm.
  
  XXX TEST FAILS, BOTH: ["cctbx","qm"]. Is this related to tst_02 failing???? XXX
  
  """
  for restraints in ["cctbx","qm"]:
    if 0:
      print("Using restraints:", restraints)
    result = []
    for clustering in [True, False]:
      if 0:
        print("  clustering", clustering, "-"*30)
      model = get_model()
      if(restraints=="qm"):
        fq = from_qm(
          pdb_hierarchy    = model.get_hierarchy(),
          qm_engine_name   = "mopac",
          method           = "PM3",
          crystal_symmetry = model.crystal_symmetry(),
          clustering       = clustering)
      else:
        fq = from_cctbx(restraints_manager = model.get_restraints_manager())
      if(clustering):
        fm = fragments(
          working_folder             = os.path.split("./ase/tmp_ase.pdb")[0]+ "/",
          clustering_method          = betweenness_centrality_clustering,
          maxnum_residues_in_cluster = 8,
          charge_embedding           = False,
          two_buffers                = False,
          fast_interaction           = True,
          pdb_hierarchy              = model.get_hierarchy().deep_copy(), # deep copy just in case
          qm_engine_name             = "mopac",
          crystal_symmetry           = model.crystal_symmetry())
        fc = from_cluster(
          restraints_manager = fq,
          fragment_manager   = fm,
          parallel_params    =get_master_phil().extract())
      else:
        fc = fq
      energy, gradients = fc.target_and_gradients(sites_cart=model.get_sites_cart())
      if(restraints=="qm"):
        energy = energy*(kcal/mol)*(kcal/mol)/Hartree
        gradients = gradients*(kcal/mol)*(kcal/mol)*(Bohr/Hartree)
      gradients = gradients.as_double()
      result.append(gradients.deep_copy())
    #
    diff = flex.abs(result[0] - result[1])
    for d in diff:
      print(d)
    max_diff = flex.max(diff)
    print("  max(diff_grad):", max_diff)
    assert max_diff < 1.e-9

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=True)
