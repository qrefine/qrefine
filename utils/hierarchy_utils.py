import os, sys
import StringIO
import iotbx
from libtbx.utils import Sorry
from mmtbx import monomer_library
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation

from iotbx.pdb import amino_acid_codes as aac
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
        print atom.quote()

def display_atom_group(ag, verbose=False):
  if verbose:
    outl = display_atom_group(ag)
    for atom in ag.atoms():
      outl += '\n  %s ' % atom.format_atom_record()
    return outl
  return '    atom_group: resname="%s" altloc="%s"' % (ag.resname, ag.altloc)


def generate_residues_via_conformer(hierarchy,
                                    backbone_only=False,
                                    verbose=False,
                                    ):
  backbone_asc = hierarchy.atom_selection_cache()
  backbone_sel = backbone_asc.selection("name ca or name c or name n or name o or name cb")
  backbone_hierarchy = hierarchy.select(backbone_sel)
  get_class = iotbx.pdb.common_residue_names_get_class
  loop_hierarchy=hierarchy
  if backbone_only: loop_hierarchy=backbone_hierarchy
  for model in loop_hierarchy.models():
    if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      if verbose: print 'chain: "%s"' % chain.id
      for conformer in chain.conformers():
        if verbose: print '  conformer: altloc="%s"' % (
          conformer.altloc)
#        while threes: del threes[0]
#        threes.start=None
#        threes.end=None
#        list_of_threes = []
        for residue in conformer.residues():
          if verbose:
            if residue.resname not in ["HOH"]:
              print '    residue: resname="%s" resid="%s"' % (
                residue.resname, residue.resid())
          if verbose: print '      residue class : %s' % get_class(residue.resname)
          if get_class(residue.resname) not in ["common_amino_acid",
                                                'modified_amino_acid',
                                              ]:
            # this needs to be moved to cctbx get_class
            #'ETA', # COOH terminal - not in modified
            if residue.resname not in aac.three_letter_l_given_three_letter_d:
              continue
          yield residue

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

def generate_protein_fragments(hierarchy,
                               geometry,
                               backbone_only=False,
                               use_capping_hydrogens=False,
                               verbose=False,
                               ):
  from mmtbx.conformation_dependent_library.multi_residue_class import \
    ThreeProteinResidues, RestraintsRegistry
  registry = RestraintsRegistry()
  threes = ThreeProteinResidues(geometry, registry=registry)
  for residue in generate_residues_via_conformer(hierarchy,
                                                 backbone_only=backbone_only,
                                                 verbose=verbose,
                                                 ):
    list.append(threes, residue)
    if verbose: print 'THREE',threes
    sub_unit = threes.provide_second_sub_unit_if_unlinked()
    if verbose: print 'THREE, SUBUNIT',threes, sub_unit
    if sub_unit:
      threes.start = True
      threes.end = True
      yield threes
      threes = sub_unit
  threes.start = True
  threes.end = True
  yield threes

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

def write_hierarchy(pdb_filename, pdb_inp, hierarchy, underscore):
  output = "%s_%s.pdb" % (pdb_filename[:-4],
                          underscore,
                          )
  output = os.path.basename(output)
  f=file(output, "wb")
  f.write(hierarchy.as_pdb_string(
    crystal_symmetry=pdb_inp.crystal_symmetry()),
          )
  f.close()
  #print "\n  Output written to: %s" % output
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

def add_hydrogens_using_ReadySet(pdb_filename,
                                 ligand_cache_directory=None,
                                 ):
  from elbow.command_line.ready_set import run_though_all_the_options
  pdb_lines = open(pdb_filename, 'rb').read()
  output_file_name=pdb_filename.replace('.pdb',
                                        '.updated.pdb',
                                        )
  sys_std = None
  if 1:
    sys_std = sys.stdout
    sys.stdout = StringIO.StringIO()
    print 'NOT'*20
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
    print 'NOT'*20
    sys.stdout = sys_std
  if 0:
    #print 'overwriting',pdb_filename
    f=file(pdb_filename, 'wb')
    for line in rc['cryst1']:
      f.write('%s\n' % line)
    f.write(rc['model_hierarchy'].as_pdb_string())
    f.close()
  else:
    return rc['model_hierarchy']

def attempt_to_squash_alt_loc(hierarchy):
  indices = hierarchy.altloc_indices()
  altlocs = filter(None,indices)
  if len(altlocs)==0: return hierarchy
  elif len(altlocs)>1: return None
  squash_hierarchy = hierarchy.deep_copy()
  for rg in squash_hierarchy.residue_groups():
    if len(rg.atom_groups())==1: continue
    ags = rg.atom_groups()
    detached_ag = ags[1].detached_copy()
    for atom in detached_ag.atoms():
      ags[0].append_atom(atom.detached_copy())
    rg.remove_atom_group(ags[1])
  return squash_hierarchy

def merge_atoms_at_end_to_residues(hierarchy):
  residues = {}
  for ag in hierarchy.atom_groups():
    # complication with alt.loc.
    key = ag.id_str()
    previous_instance = residues.setdefault(key, None)
    if previous_instance:
      # move atoms from here to there
      for atom in ag.atoms():
        previous_instance.append_atom(atom.detached_copy())
        ag.remove_atom(atom)
      rg = ag.parent()
      rg.remove_atom_group(ag)
      chain = rg.parent()
      chain.remove_residue_group(rg)
    residues[key] = ag
  return hierarchy

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
