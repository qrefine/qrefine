from __future__ import division
from __future__ import absolute_import

import os
import time
from qrefine.tests.unit import run_tests
import libtbx.load_env
from mmtbx import monomer_library
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation
from scitbx.array_family import flex

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def get_bond_rmsd(file_name):
  mon_lib_srv = monomer_library.server.server()
  ener_lib    = monomer_library.server.ener_lib(use_neutron_distances=True)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = file_name,
    params         = params,
    force_symmetry = True)
  sites_cart = processed_pdb_file.xray_structure().sites_cart()
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    assume_hydrogens_all_missing = False,
    plain_pairs_radius = 5.0)
  es = geometry.energies_sites(sites_cart = sites_cart)
  return es.bond_deviations()[2], sites_cart

def run(prefix):
  """
  Exercise standard (cctbx-based restraints) optimization.
  Assert mode=opt is approximately equivalent to data_weight=0.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  # Run optimization
  print("setup_helix_done")
  run_tests.run_cmd(prefix,args = ["restraints=cctbx","mode=opt",
    "clustering=False","minimizer=lbfgsb", "number_of_micro_cycles=3",
    "max_iterations_refine=100"])
  print(prefix)
  print("run_test_done")
  assert get_bond_rmsd(file_name=os.path.join(qr_unit_tests,"data_files","m00_poor.pdb"))[0] > 0.1
  result1, sc1 = get_bond_rmsd(file_name=os.path.join(prefix,"m00_poor_refined.pdb"))
  assert result1 < 0.001, result1
  #Run refinement without data term
  run_tests.run_cmd(prefix,args = ["restraints=cctbx","data_weight=0",
    "clustering=False","minimizer=lbfgsb", "number_of_micro_cycles=3",
    "max_iterations_refine=100", "stpmax=999"])
  result2, sc2 = get_bond_rmsd(file_name=os.path.join(prefix,"m00_poor_refined.pdb"))
  assert result2 < 0.01, result2
  #
  dist = flex.mean(flex.sqrt((sc1 - sc2).dot()))
  assert dist < 0.15

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)

