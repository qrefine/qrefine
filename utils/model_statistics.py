from __future__ import print_function
import sys
import iotbx.pdb
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import model_statistics
from libtbx.utils import null_out

def get_model_stat(pdb_file_name=None, pdb_hierarchy=None,
                   crystal_symmetry=None, show=True):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib    = mmtbx.monomer_library.server.ener_lib()
  params = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  params.sort_atoms = False
  params.clash_guard.nonbonded_distance_threshold=None
  #
  assert [pdb_file_name, pdb_hierarchy].count(None)==1
  assert [crystal_symmetry, pdb_hierarchy].count(None) in [0,2]
  if(pdb_file_name is not None):
    pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
    crystal_symmetry = pdb_inp.crystal_symmetry()
    pdb_hierarchy = \
      iotbx.pdb.input(file_name=pdb_file_name).construct_hierarchy()
  # XXX centralize. 3rd copy-paste!
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.clash_guard.nonbonded_distance_threshold=None
  params.disable_uc_volume_vs_n_atoms_check=False
  params.use_neutron_distances = True
  params.restraints_library.cdl = False
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
       mon_lib_srv              = mon_lib_srv,
       ener_lib                 = ener_lib,
       params                   = params,
       pdb_inp                  = pdb_hierarchy.as_pdb_input(),
       strict_conflict_handling = False,
       crystal_symmetry         = crystal_symmetry,
       force_symmetry           = True,
       log                      = null_out())
  xrs = processed_pdb_file.xray_structure()
  sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  #
  restraints_manager = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    assume_hydrogens_all_missing = not has_hd,
    plain_pairs_radius           = 5.0)
  #
  stats = model_statistics.geometry(
    pdb_hierarchy      = pdb_hierarchy,
    restraints_manager = restraints_manager,
    molprobity_scores  = True)
  if(show):
    print("-"*79)
    stats.show(out=sys.stdout)
    print("-"*79)
  #
  return stats
