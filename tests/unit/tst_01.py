from __future__ import division

import os
import time
import run_tests
import libtbx.load_env
from mmtbx import monomer_library
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation

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
  return es.bond_deviations()[2]

def run(prefix = "tst_01"):
  """
  Exercise standard (cctbx-based restraints) optimization.
  Assert mode=opt is approximately equivalent to data_weight=0.
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  # Run optimization
  run_tests.run_cmd(prefix,args = ["restraints=cctbx","mode=opt"])
  assert get_bond_rmsd(file_name=os.path.join(qr_unit_tests,"data_files","m00_poor.pdb")) > 0.1
  assert get_bond_rmsd(file_name=os.path.join(prefix,"m00_poor_refined.pdb")) < 0.0009
  #Run refinement without data term
  run_tests.run_cmd(prefix,args = ["restraints=cctbx","data_weight=0"])
  assert get_bond_rmsd(file_name=os.path.join(prefix,"m00_poor_refined.pdb")) < 0.0009
  # Cleanup
  run_tests.clean_up(prefix)

  return 0

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_01"
  rc = run(prefix)
  print prefix +":  OK  " + "Time: %6.2f (s)"%(time.time()-t0)
  run_tests.clean_up(prefix)

