from __future__ import print_function
import os, sys
from io import StringIO
import iotbx
import mmtbx.model
from libtbx.utils import Sorry
from libtbx.utils import null_out
from mmtbx import monomer_library
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation

n_terminal_amino_acid_codes = ['FVA']
c_terminal_amino_acid_codes = []


mon_lib_server = server.server()
get_class = iotbx.pdb.common_residue_names_get_class

def display_residue_group(rg):
  return '  residue_group: resseq="%s" icode="%s"' % (rg.resseq, rg.icode)

def display_residue_group_from_hierarchy(hierarchy,
                                         resseq=None,
                                         ):
  for rg in hierarchy.residue_groups():
    if rg.resseq.strip()==str(resseq):
      for atom in rg.atoms():
        print(atom.quote())

def display_atom_group(ag, verbose=False):
  if verbose:
    outl = display_atom_group(ag)
    for atom in ag.atoms():
      outl += '\n  %s ' % atom.format_atom_record()
    return outl
  return '    atom_group: resname="%s" altloc="%s"' % (ag.resname, ag.altloc)

def generate_residue_groups(hierarchy,
                            assert_no_alt_loc=False,
                            exclude_water=False,
                            verbose=False,
                            ):
  for rg in hierarchy.residue_groups():
    if verbose: display_residue_group(rg)
    if exclude_water:
      ag=rg.atom_groups()[0]
      if get_class(ag.resname)=='common_water': continue
    if assert_no_alt_loc:
      if len(rg.atom_groups())>1:
        raise Sorry("Contains alt. locs.")
    yield rg

def get_pdb_interpretation_params():
  master_params = iotbx.phil.parse(
    input_string=pdb_interpretation.master_params_str,
    process_includes=True)
  params = master_params.extract()
  params.automatic_linking.link_metals=True
  return params

def get_processed_pdb(pdb_filename=None,
                      pdb_hierarchy=None,
                      crystal_symmetry=None,
                      raw_records=None,
                      pdb_inp=None,
                      params=None,
                      cif_objects=None,
                      ):
  if params is None:
    params = get_pdb_interpretation_params()
  mon_lib_srv = monomer_library.server.server()
  if cif_objects:
    for cif_object in cif_objects:
      mon_lib_srv.process_cif_object(cif_object[1])
  ener_lib = monomer_library.server.ener_lib()
  ppf = monomer_library.pdb_interpretation.process(
    mon_lib_srv           = mon_lib_srv,
    ener_lib              = ener_lib,
    params                = params,
    pdb_hierarchy=pdb_hierarchy,
    crystal_symmetry=crystal_symmetry,
    file_name             = pdb_filename,
    pdb_inp               = pdb_inp,
    raw_records           = raw_records,
    keep_monomer_mappings = True,
    log                   = null_out())
  return ppf

def get_raw_records(pdb_inp=None,
                    pdb_hierarchy=None,
                    crystal_symmetry=None,
                   ):
  if crystal_symmetry:
    lines = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
  else:
    lines = pdb_hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.crystal_symmetry())
  return lines

def is_n_terminal_atom_group(atom_group):
  if atom_group.resname in n_terminal_amino_acid_codes: return True
  return False

def process_model_file(pdb_file_name, cif_objects, crystal_symmetry):
  import mmtbx
  import iotbx.pdb
  from libtbx.utils import null_out
  from libtbx import group_args
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
  model = mmtbx.model.manager(
    model_input               = pdb_inp,
    crystal_symmetry          = crystal_symmetry,
    restraint_objects         = cif_objects,
    log                       = null_out())
  model.process(make_restraints=True,
                grm_normalization=False,
                pdb_interpretation_params = params)
  return group_args(
    model              = model,
    pdb_hierarchy      = model.get_hierarchy(),
    xray_structure     = model.get_xray_structure(),
    cif_objects        = cif_objects,
    crystal_symmetry   = model.get_xray_structure().crystal_symmetry()
    )

# def process_model_file(pdb_file_name, scattering_table, log=null_out(), cif_objects=None,
#                        crystal_symmetry=None):
#   import iotbx.pdb
#   params = mmtbx.model.manager.get_default_pdb_interpretation_params()
#   params.pdb_interpretation.use_neutron_distances = True
#   params.pdb_interpretation.restraints_library.cdl = False
#   params.pdb_interpretation.sort_atoms = False
#   pdb_inp = iotbx.pdb.input(file_name=pdb_file_name)
#   validate_model_via_hierarchy(pdb_inp.construct_hierarchy())
#   model = mmtbx.model.manager(
#     model_input       = pdb_inp,
#     crystal_symmetry  = crystal_symmetry,
#     restraint_objects = cif_objects,
#     log               = log)
#   model.process(make_restraints=True, grm_normalization=True,
#     pdb_interpretation_params = params)
#   model.setup_scattering_dictionaries(
#     scattering_table  = scattering_table,
#     d_min             = 1.0,
#     log               = log)
#   return model

def model_from_hierarchy(pdb_hierarchy, crystal_symmetry, cif_objects=None, log=None):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  params.pdb_interpretation.flip_symmetric_amino_acids = False
  params.pdb_interpretation.correct_hydrogens=False
  # params.pdb_interpretation.const_shrink_donor_acceptor=0.6 # XXX change back to old default?
  if log is None: log = null_out()
  model = mmtbx.model.manager(
    model_input       = None,
    pdb_hierarchy     = pdb_hierarchy,
    crystal_symmetry  = crystal_symmetry,
    restraint_objects = cif_objects,
    log               = log)
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params,)
  return model
