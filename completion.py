from __future__ import print_function
from __future__ import absolute_import
import sys

import iotbx
from mmtbx.monomer_library import server
from libtbx.utils import Sorry

from iotbx.pdb import amino_acid_codes as aac

mon_lib_server = server.server()
get_class = iotbx.pdb.common_residue_names_get_class

from qrefine.utils import hierarchy_utils
from mmtbx.hydrogens.specialised_hydrogen_atoms import conditional_add_cys_hg_to_atom_group
from mmtbx.hydrogens.specialised_hydrogen_atoms import conditional_remove_cys_hg_to_atom_group
from mmtbx.ligands.hierarchy_utils import get_bonds_as_dict
from mmtbx.ligands.ready_set_utils import add_n_terminal_hydrogens_to_residue_group
from mmtbx.ligands.ready_set_utils import add_c_terminal_oxygens_to_residue_group
from mmtbx.ligands.ready_set_utils import generate_protein_fragments

log = sys.stdout

from cctbx.eltbx import distance_based_connectivity

def from_simple_connectivity(
      hierarchy,
      geometry_restraints_manager,
      use_capping_hydrogens=False,
      append_to_end_of_model=False,
      selection=None):
  atoms = hierarchy.atoms()
  bonds=get_bonds_as_dict(geometry_restraints_manager)
  conn = distance_based_connectivity.build_edge_list(
    sites_cart=atoms.extract_xyz(), elements=atoms.extract_element())
  dconn = {}
  for p in conn:
    dconn.setdefault(p[0],[]).append(p[1])
    dconn.setdefault(p[1],[]).append(p[0])
  additional_hydrogens = hierarchy_utils.smart_add_atoms()
  for k, v in zip(dconn.keys(), dconn.values()):
    a  = atoms[k]
    if selection is not None and not selection[a.i_seq]: continue # XXX
    rg = a.parent().parent()
    skip=False
    for resname in rg.unique_resnames():
      if get_class(resname)!="common_amino_acid": skip=True
    if skip: continue
    if(a.name == ' C  ' and len(v)!=3):
      #print(a.name, [atoms[i].name for i in v], rg.resseq)
      rc = add_c_terminal_oxygens_to_residue_group(
        rg,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
    if(a.name == ' N  ' and len(v)!=3):
      #print(a.name, [atoms[i].name for i in v], rg.resseq)
      rc = add_n_terminal_hydrogens_to_residue_group(
        rg,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)

  for rg in hierarchy.residue_groups():
    if use_capping_hydrogens:

      skip=True
      for a in rg.atoms():
        if selection is not None and not selection[a.i_seq]:
          skip=True
          break
      if skip: continue

      rc = conditional_add_cys_hg_to_atom_group(
        geometry_restraints_manager,
        rg,
        append_to_end_of_model=append_to_end_of_model)
      if rc: additional_hydrogens.append(rc)
    else:
      conditional_remove_cys_hg_to_atom_group(geometry_restraints_manager, rg)
  return additional_hydrogens

def iterate_over_threes(hierarchy,
                        geometry_restraints_manager,
                        use_capping_hydrogens=False,
                        append_to_end_of_model=False,
                        verbose=False,
                        ):
  atoms = hierarchy.atoms()
  ###
  def get_residue_group(residue):
    for atom in residue.atoms():
      atom = atoms[atom.i_seq]
      break
    return atom.parent().parent()
  ###
  bonds=get_bonds_as_dict(geometry_restraints_manager)
  ###
  additional_hydrogens=hierarchy_utils.smart_add_atoms()
  for three in generate_protein_fragments(
    hierarchy,
    geometry_restraints_manager,
    backbone_only=False,
    use_capping_hydrogens=use_capping_hydrogens,
  ):
    if verbose: print(three)
    if not len(three): continue
    ptr=0
    assert three.are_linked()
    if use_capping_hydrogens:
      for i in range(len(three)):
        rg = get_residue_group(three[i])
        rc = conditional_add_cys_hg_to_atom_group(
          geometry_restraints_manager,
          rg,
          append_to_end_of_model=append_to_end_of_model,
          )
        if rc:
          additional_hydrogens.append(rc)
    else:
      for i in range(len(three)):
        rg = get_residue_group(three[i])
        conditional_remove_cys_hg_to_atom_group(geometry_restraints_manager,
                                                rg,
                                                )
    # check if N-term residue - FVA
    n_term_done = False
    if three[0].resname in ['FVA',
                            ]:
      n_term_done = True
      ptr+=1
      assert ptr==1, 'ptr (%d) is not 1' % ptr
    if three.start and not n_term_done:
      ptr+=1
      assert ptr==1, 'ptr (%d) is not 1' % ptr
      rg = get_residue_group(three[0])
      #print("add_n_terminal_hydrogens_to_residue_group",rg.resseq)
      rc = add_n_terminal_hydrogens_to_residue_group(
        rg,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
      #hierarchy.reset_i_seq_if_necessary()
    c_term_done = False
    if three[-1].resname in ['ETA',
                             ]:
      c_term_done = True
      ptr-=1
      assert ptr==0, 'ptr (%d) is not 0' % ptr
    if three.end and not c_term_done:
      ptr-=1
      assert ptr==0, 'ptr (%d) is not 0' % ptr
      rg = get_residue_group(three[-1])
      #print("add_c_terminal_oxygens_to_residue_group",rg.resseq)
      rc = add_c_terminal_oxygens_to_residue_group(
        rg,
        bonds=bonds,
        use_capping_hydrogens=use_capping_hydrogens,
        append_to_end_of_model=append_to_end_of_model,
      )
      if rc: additional_hydrogens.append(rc)
    else:
      pass
  return additional_hydrogens

def add_terminal_hydrogens_qr(
    hierarchy,
    geometry_restraints_manager,
    add_to_chain_breaks=False,
    use_capping_hydrogens=False,  # instead of terminal H
    append_to_end_of_model=False, # useful for Q|R
    original_hierarchy=None,
    verbose=False,
    selection=None,
    ):
  # add N terminal hydrogens because Reduce only does it to resseq=1
  # needs to be alt.loc. aware for non-quantum-refine
  if True:
    additional_hydrogens = from_simple_connectivity(
      hierarchy,
      geometry_restraints_manager,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
      selection=selection,
      )
  else:
    additional_hydrogens=iterate_over_threes(
      hierarchy,
      geometry_restraints_manager,
      use_capping_hydrogens=use_capping_hydrogens,
      append_to_end_of_model=append_to_end_of_model,
    )

  if append_to_end_of_model and additional_hydrogens:
    from mmtbx.ligands.ready_set_utils import _add_atoms_from_chains_to_end_of_hierarchy
    tmp = []
    for group in additional_hydrogens:
      for chain in group:
        tmp.append(chain)
    _add_atoms_from_chains_to_end_of_hierarchy(hierarchy, tmp)

def remove_acid_side_chain_hydrogens(hierarchy, selection=None):
  from mmtbx.ligands.ready_set_basics import get_proton_info
  proton_element, proton_name = get_proton_info(hierarchy)
  removes = {"GLU" : "%sE2" % proton_element,
             "ASP" : "%sD2" % proton_element,
             }
  for ag in hierarchy.atom_groups():
    r = removes.get(ag.resname, None)
    if r is None: continue
    atom = ag.get_atom(r)
    skip=False
    if selection is not None and not selection[atom.i_seq]: skip=True
    if skip:
      continue
    if atom:
      ag.remove_atom(atom)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def __HELPER1(crystal_symmetry, hierarchy, params):
  import mmtbx.model
  import iotbx.pdb
  from libtbx import group_args
  # 1. Provide an explicit model_input to avoid the model.manager bug
  # that destroys the hierarchy, and to prevent in-place mutation.
  raw_records = hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records)
  # 2. Use model.manager with stop_for_unknowns=False to mimic the more
  # permissive behavior of the older hierarchy_utils
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    crystal_symmetry  = crystal_symmetry,
    stop_for_unknowns = False)
  p = model.get_default_pdb_interpretation_params()
  p.pdb_interpretation = params
  model.process(pdb_interpretation_params = p, make_restraints=True)
  # 3. Restore the exact original Cartesian coordinates. Converting to a PDB
  # string truncated them to 3 decimal places. Setting the sites_cart
  # directly from the original hierarchy restores the exact floating-point values.
  sites_cart = hierarchy.atoms().extract_xyz()
  model.get_hierarchy().atoms().set_xyz(sites_cart)
  return model

def complete_pdb_hierarchy(hierarchy,
                           geometry_restraints_manager,
                           use_capping_hydrogens=False,
                           append_to_end_of_model=False,
                           crystal_symmetry=None,
                           original_pdb_filename=None,
                           verbose=False,
                           debug=False,
                           use_reduce=True,
                           selection=None,# ONLY APPLIES FOR CAPPING
                          ):
  #
  # some validations
  #
  if not hierarchy.is_hierarchy_altloc_consistent():
    hierarchy.is_hierarchy_altloc_consistent(verbose=True)
    raise Sorry('Altloc structure of model not consistent. Make each altloc the same depth or remove completely.')
  original_hierarchy = None
  params = hierarchy_utils.get_pdb_interpretation_params()
  params.restraints_library.cdl=False
  params.automatic_linking.link_residues=True
  if use_capping_hydrogens:
    params.link_distance_cutoff=1.8 # avoid linking across a single missing AA
    if original_pdb_filename:
      original_pdb_inp = iotbx.pdb.input(original_pdb_filename)
      original_hierarchy = original_pdb_inp.construct_hierarchy()
  #
  # assume model is heavy-atom complete
  #
  if not use_capping_hydrogens:
    model = __HELPER1(crystal_symmetry=crystal_symmetry,
      hierarchy=hierarchy, params=params)
  #
  # need to use Reduce/ReadySet! to add hydrogens
  #
  if not use_capping_hydrogens:
    from mmtbx.building import extend_sidechains
    n_changed = extend_sidechains.extend_protein_model(
      model.get_hierarchy(),
      mon_lib_server,
      add_hydrogens=False,
    )

    if(use_reduce):
      print("Using reduce to add hydrogens",file=log)
      hierarchy = hierarchy_utils.add_hydrogens_using_reduce(
      pdb_hierarchy    = model.get_hierarchy(),
      crystal_symmetry = crystal_symmetry
      )
    else:
      print("Using ReadySet to add hydrogens",file=log)
      hierarchy = hierarchy_utils.add_hydrogens_using_ReadySet(output)
  #
  # remove side chain acid hydrogens - maybe not required since recent changes
  #
  model = __HELPER1(crystal_symmetry=crystal_symmetry, hierarchy=hierarchy, params=params)

  if not use_capping_hydrogens:
    remove_acid_side_chain_hydrogens(
      hierarchy=model.get_hierarchy(),
      selection=selection)
  #
  # add hydrogens in special cases
  #  eg ETA
  #  eg N - H, H2
  #
  hierarchy = model.get_hierarchy()
  model = __HELPER1(crystal_symmetry=crystal_symmetry, hierarchy=hierarchy, params=params)
  #
  # maybe more to cctbx
  #
  add_terminal_hydrogens_qr( model.get_hierarchy(),
                             model.get_restraints_manager().geometry,
                             use_capping_hydrogens=use_capping_hydrogens,
                             append_to_end_of_model=append_to_end_of_model,
                             original_hierarchy=original_hierarchy,
                             verbose=verbose,
                             selection=selection,
                            ) # in place
  model.get_hierarchy().atoms().set_chemical_element_simple_if_necessary()
  model.get_hierarchy().sort_atoms_in_place()
  return model

def run(pdb_filename=None,
        pdb_hierarchy=None,
        crystal_symmetry=None,
        model_completion=True,
        original_pdb_filename=None,
        append_to_end_of_model=True,
        use_reduce=True,
        selection=None, # ONLY APPLIES FOR CAPPING
        ):
  #
  # function as be used in two main modes
  #   1. completing a model with hydrogens in a protein-like manner
  #   2. completing a cluster with hydrogens in a QM-sensible manner
  #
  # Validation
  #
  if pdb_hierarchy:
    assert crystal_symmetry
    assert pdb_filename is None
  #
  # output
  #
  if model_completion:
    use_capping_hydrogens=False
    fname = 'complete' # only this uses reduce/ReadySet
  else:
    use_capping_hydrogens=True
    fname = 'capping'
  #
  # adjust parameters
  #
  params=None
  if use_capping_hydrogens:
    params = hierarchy_utils.get_pdb_interpretation_params()
    params.link_distance_cutoff=1.8

  if(pdb_filename is not None):
    pi = iotbx.pdb.input(file_name = pdb_filename)
    pdb_hierarchy = pi.construct_hierarchy()
    crystal_symmetry = pi.crystal_symmetry()

  MH = __HELPER1(crystal_symmetry=crystal_symmetry, hierarchy=pdb_hierarchy, params=params)

  #
  # guts
  #
  model = complete_pdb_hierarchy(
    MH.get_hierarchy(),
    MH.get_restraints_manager(),
    use_capping_hydrogens=use_capping_hydrogens,
    append_to_end_of_model=append_to_end_of_model, # needed for clustering
                                                   # code and Molprobity
    crystal_symmetry=MH.crystal_symmetry(), # used in get_raw_records. why
    original_pdb_filename=original_pdb_filename,
    verbose=False,
    use_reduce=use_reduce,
    selection=selection, # ONLY APPLIES FOR CAPPING
  )
  if pdb_filename:
    output = hierarchy_utils.write_hierarchy(
      pdb_filename,
      model.crystal_symmetry(),
      model.get_hierarchy(),
      fname,
    )
  return model.get_hierarchy()

def display_hierarchy_atoms(hierarchy, n=5):
  print('-'*80)
  for i, atom in enumerate(hierarchy.atoms()):
    print(atom.quote())
    if i>n: break

if __name__=="__main__":
  def _fake_phil_parse(arg):
    def _boolean(s):
      if s.lower() in ['1', 'true']: return True
      elif s.lower() in ['0', 'false']: return False
      else: return s
    rc = {arg.split('=')[0] : _boolean(arg.split('=')[1])}
    return rc
  args = sys.argv[1:]
  del sys.argv[1:]
  kwds={}
  remove=[]
  for i, arg in enumerate(args):
    if arg.find('=')>-1:
      kwds.update(_fake_phil_parse(arg))
      remove.append(i)
  remove.reverse()
  for r in remove: del args[r]
  if 'test_from_clustering' in args:
    args.remove('test_from_clustering')
    ppf = hierarchy_utils.get_processed_pdb(args[0])
    sites_cart = ppf.all_chain_proxies.pdb_hierarchy.atoms().extract_xyz()
    sites_cart[0]=(4.123456789, 7.7, 1.5)
    ppf.all_chain_proxies.pdb_hierarchy.atoms().set_xyz(sites_cart)
    kwds['pdb_hierarchy'] = ppf.all_chain_proxies.pdb_hierarchy
    kwds['crystal_symmetry'] = ppf.all_chain_proxies.pdb_inp.crystal_symmetry()
    display_hierarchy_atoms(kwds['pdb_hierarchy'])
    rc = run(None, **kwds)
    display_hierarchy_atoms(rc)
    assert 0, 'FINISHED TESTING'
  run(*tuple(args), **kwds)
