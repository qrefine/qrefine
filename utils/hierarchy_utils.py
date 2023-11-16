from __future__ import print_function
import os, sys
from io import StringIO
import iotbx
from libtbx.utils import Sorry
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
    file_name             = pdb_filename,
    pdb_inp               = pdb_inp,
    raw_records           = raw_records,
    keep_monomer_mappings = True)
  return ppf

def write_hierarchy(pdb_filename, pdb_inp, hierarchy, underscore, verbose=False):
  output = "%s_%s.pdb" % (pdb_filename[:-4],
                          underscore,
                          )
  output = os.path.basename(output)
  f=open(output, "w")
  f.write(hierarchy.as_pdb_string(
    crystal_symmetry=pdb_inp.crystal_symmetry()),
          )
  f.close()
  if verbose: print("\n  Output written to: %s" % output)
  return output

def get_raw_records(pdb_inp=None,
                    pdb_hierarchy=None,
                    crystal_symmetry=None,
                   ):
  if crystal_symmetry:
    lines = pdb_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
  else:
    lines = pdb_hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.crystal_symmetry())
  return lines

def add_hydrogens_using_reduce(pdb_filename):
  assert 0

def add_hydrogens_using_ReadySet(pdb_filename, ligand_cache_directory=None):
  from elbow.command_line.ready_set import run_though_all_the_options
  pdb_lines = open(pdb_filename, 'r').read()
  output_file_name=pdb_filename.replace('.pdb',
                                        '.updated.pdb',
                                        )
  sys_std = None
  if 1:
    sys_std = sys.stdout
    sys.stdout = StringIO()
    print('NOT'*20)
  rc = run_though_all_the_options(
    pdb_lines,
    [], # args
    hydrogens=True,
    reduce_args=['-BUILD'],
    ligands=True,
    ligand_cache_directory=ligand_cache_directory,
    add_h_to_water=True,
    metals=False,
    output_file_name=output_file_name, # needed
    )
  if sys_std:
    print('NOT'*20)
    sys.stdout = sys_std
  if 0:
    #print 'overwriting',pdb_filename
    f=open(pdb_filename, 'w')
    for line in rc['cryst1']:
      f.write('%s\n' % line)
    f.write(rc['model_hierarchy'].as_pdb_string())
    f.close()
  else:
    return rc['model_hierarchy']

def remove_atom_from_chain(chain, atom):
  for rg in chain.residue_groups():
    for ag in rg.atom_groups():
      for at in ag.atoms():
        if atom==at:
          ag.remove_atom(atom)
          break

class smart_add_atoms(list):
  def __init__(self): pass

  def append(self, item):
    for chain1 in item:
      remove = []
      for atom1 in chain1.atoms():
        for chain_list in self:
          for chain2 in chain_list:
            for atom2 in chain2.atoms():
              if atom1.quote()==atom2.quote():
                remove.append(atom1)
      if remove:
        for atom in remove:
          remove_atom_from_chain(chain1, atom)
    list.append(self, item)

def is_n_terminal_residue(residue_group):
  residues = []
  for atom_group in residue_group.atom_groups():
    if atom_group.resname not in residues: residues.append(atom_group.resname)
  assert len(residues)==1
  if residues[0] in n_terminal_amino_acid_codes: return True
  return False

def is_n_terminal_atom_group(atom_group):
  if atom_group.resname in n_terminal_amino_acid_codes: return True
  return False

def process_model_file(pdb_file_name, cif_objects, crystal_symmetry):
  import mmtbx
  import iotbx.pdb
  from libtbx.utils import null_out
  from libtbx import group_args
  from ase.io import read as ase_io_read
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
    ase_atoms          = ase_io_read(pdb_file_name), # To be able to use ASE LBFGS
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
