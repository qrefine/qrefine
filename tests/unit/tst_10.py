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
from libtbx.utils import null_out

import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.restraints
from mmtbx import monomer_library

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
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  h = pdb_inp.construct_hierarchy()
  #
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib    = mmtbx.monomer_library.server.ener_lib(
    use_neutron_distances = True)
  #
  params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  params.sort_atoms = False
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = ener_lib,
    params                   = params,
    pdb_hierarchy            = h,
    strict_conflict_handling = False,
    crystal_symmetry         = pdb_inp.crystal_symmetry(),
    force_symmetry           = True,
    log                      = null_out())
  xrs = processed_pdb_file.xray_structure()
  sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    assume_hydrogens_all_missing = not has_hd,
    plain_pairs_radius           = 5.0)
  return mmtbx.restraints.manager(
     geometry = geometry, normalization = False), pdb_inp.crystal_symmetry(),h

def run(prefix, verbose=False):
  """
  Exercise combined energy and gradients from cluster qm.
  """
  for restraints in ["cctbx","qm"]:
    if verbose: print("Using restraints:", restraints)
    result = []
    for clustering in [True, False]:
      if verbose: print("  clustering", clustering, "-"*30)
      rm, cs, h = get_model()
      if(restraints=="qm"):
        fq = from_qm(
          pdb_hierarchy    = h,
          qm_engine_name   = "mopac",
          method           = "PM3",
          crystal_symmetry = cs,
          clustering       = clustering)
      else:
        fq = from_cctbx(restraints_manager = rm)
      if(clustering):
        fm = fragments(
          working_folder             = os.path.split("./ase/tmp_ase.pdb")[0]+ "/",
          clustering_method          = betweenness_centrality_clustering,
          maxnum_residues_in_cluster = 2,
          charge_embedding           = False,
          two_buffers                = False,
          fast_interaction           = True,
          pdb_hierarchy              = h.deep_copy(), # deep copy just in case
          qm_engine_name             = "mopac",
          crystal_symmetry           = cs)
        fragment_sizes = [f.iselection().size() for f in fm.fragment_selections]
        fragment_sizes.sort()
        assert fragment_sizes == [44, 57, 65, 85]
        fc = from_cluster(
          restraints_manager = fq,
          fragment_manager   = fm,
          parallel_params    = get_master_phil().extract())
      else:
        fc = fq
      energy, gradients = fc.target_and_gradients(sites_cart=h.atoms().extract_xyz())
      if(restraints=="qm"):
        energy = energy*(kcal/mol)*(kcal/mol)/Hartree
        gradients = gradients*(kcal/mol)*(kcal/mol)*(Bohr/Hartree)
      gradients = gradients.as_double()
      result.append(gradients.deep_copy())
    #
    diff = flex.abs(result[0] - result[1])
    max_diff = flex.max(diff)
    if verbose:  print("  max(diff_grad):", max_diff)
    if(restraints=="cctbx"):
      assert max_diff < 1.e-9
    else:
      assert max_diff < 5.e-5

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
