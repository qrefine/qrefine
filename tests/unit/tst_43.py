from __future__ import division

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
import run_tests
import mmtbx.model
from libtbx.utils import null_out
from qrefine import qr

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
  file_name = os.path.join(qr_unit_tests,"data_files","h_altconf_complete.pdb")
  pdb_inp = iotbx.pdb.input(file_name)
  model = qr.process_model_file(
    pdb_file_name = file_name,
    cif_objects = None,
    crystal_symmetry=pdb_inp.crystal_symmetry()).model
  return model

def run(maxnum_residues_in_cluster):
  result = []
  for clustering in [True, False]:
    print "  clustering", clustering, "-"*30
    model = get_model()
    fq = from_cctbx(restraints_manager = model.get_restraints_manager())
    if(clustering):
      fm = fragments(
       working_folder             = os.path.split("./ase/tmp_ase.pdb")[0]+ "/",
       clustering_method          = betweenness_centrality_clustering,
       maxnum_residues_in_cluster = maxnum_residues_in_cluster,
       altloc_method              = "subtract",
       charge_embedding           = False,
       two_buffers                = False,
       clustering                 = clustering,
       pdb_hierarchy              = model.get_hierarchy().deep_copy(),
       qm_engine_name             = "mopac",
       fast_interaction           = True,
       crystal_symmetry           = model.crystal_symmetry())
    else:
      fc = fq
    fc = from_cluster(
      restraints_manager = fq,
      fragment_manager   = fm,
      parallel_params    = get_master_phil().extract())
    energy, gradients = fc.target_and_gradients(sites_cart=model.get_sites_cart())
    gradients = gradients.as_double()
    result.append(gradients.deep_copy())
  diff = flex.abs(result[0] - result[1])
  max_diff = flex.max(diff)
  #print "  max(diff_grad):", max_diff
  assert max_diff < 1.e-9

if(__name__ == "__main__"):
  for maxnum_residues_in_cluster in [2, 15]:
    print "Using maxnum_residues_in_cluster:", maxnum_residues_in_cluster
    run(maxnum_residues_in_cluster=maxnum_residues_in_cluster)
