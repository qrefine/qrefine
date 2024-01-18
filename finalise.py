from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.development.ready_set
import os
import sys

import iotbx
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_run

from qrefine import charges
from qrefine import completion
from qrefine.utils import hierarchy_utils
from qrefine.tests.unit.skip import skip
import mmtbx.model
from libtbx.utils import null_out

qrefine = libtbx.env.find_in_repositories("qrefine")

def remove_alternative_locations(hierarchy):
  # should use cctbx
  for rg in hierarchy.residue_groups():
    if len(rg.atom_groups())==1: continue
    resnames = []
    for ag in rg.atom_groups():
      if ag.resname not in resnames: resnames.append(ag.resname)
    assert len(resnames)==1
    bag = rg.atom_groups()[0]
    bag.altloc = ""
    for i, ag in enumerate(rg.atom_groups()):
      if i==0: continue
      for atom in ag.atoms():
        if bag.get_atom(atom.name.strip()): continue
        da = atom.detached_copy()
        bag.append_atom(da)
      rg.remove_atom_group(ag)
  hierarchy.atoms_reset_serial()
  hierarchy.atoms().reset_i_seq()
  return hierarchy

def run_fetch_pdb(code):
  cmd = 'iotbx.fetch_pdb %s' % code
  print('Fetching files')
  print(cmd)
  easy_run.call(cmd)

def loop_over_dir(d):
  i=0
  for filename in os.listdir(d):
    if not filename.endswith('.pdb'): continue
    i+=1
    print('%s\n %3d %s\n%s' % ('*'*80, i, os.path.join(d, filename), '*'*80))
    if filename in skip:
      print('skipping')
      continue
    if os.path.exists(filename.replace('.pdb', '.updated.pdb')):
      run(filename.replace('.pdb', '.updated.pdb'))
    else:
      run(os.path.join(d, filename))

def run(pdb_filename,
        model_completion=True,
        keep_alt_loc=False,
        skip_validation=False,
        calculate_charge=False,
        append_to_end_of_model=False,
        neutron_option=None,
        hydrogen_atom_occupancies=0.,
        use_reduce=True
        ):
  #
  # do all *.pdb in a directory
  #
  if os.path.isdir(pdb_filename):
    loop_over_dir(pdb_filename)
    return
  #
  # fetch a PDB code
  #
  if pdb_filename.find(".pdb")==-1:
    if not os.path.exists('%s.pdb' % pdb_filename):
      run_fetch_pdb(pdb_filename)
    pdb_filename = '%s.pdb' % pdb_filename

  if 1: # read file option
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=pdb_filename)
  elif 0: # raw records example
    pdb_inp = iotbx.pdb.input(pdb_filename)
    hierarchy = pdb_inp.construct_hierarchy()
    raw_records = []
    for atom in hierarchy.atoms():
      raw_records.append(atom.format_atom_record())
    ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records)
  else: # from hierarchy
    pdb_inp = iotbx.pdb.input(pdb_filename)
    hierarchy = pdb_inp.construct_hierarchy()
    ppf = hierarchy_utils.get_processed_pdb(pdb_inp=hierarchy.as_pdb_input())

  if not skip_validation:
    model = mmtbx.model.manager(
      model_input      = None,
      crystal_symmetry = ppf.all_chain_proxies.pdb_inp.crystal_symmetry(),
      pdb_hierarchy    = ppf.all_chain_proxies.pdb_hierarchy)
    initial_model_statistics = model.geometry_statistics()

  # should use cctbx
  if keep_alt_loc: pass
  else:
    hierarchy = remove_alternative_locations(
      ppf.all_chain_proxies.pdb_hierarchy
    )

  ppf.all_chain_proxies.pdb_hierarchy.shift_to_origin(
    ppf.all_chain_proxies.pdb_inp.crystal_symmetry())
  ppf.all_chain_proxies.pdb_hierarchy.remove_residue_groups_with_atoms_on_special_positions_selective(ppf.all_chain_proxies.pdb_inp.crystal_symmetry())
  if 0:
    output = hierarchy_utils.write_hierarchy(
      pdb_filename, # uses to get output filename
      ppf.all_chain_proxies.pdb_inp,
      ppf.all_chain_proxies.pdb_hierarchy,
      "temp0")
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=output)
  else:
    raw_records = hierarchy_utils.get_raw_records(
      ppf.all_chain_proxies.pdb_inp,
      ppf.all_chain_proxies.pdb_hierarchy,
    )
    ppf = hierarchy_utils.get_processed_pdb(raw_records=raw_records)
  hetero_charges = charges.get_hetero_charges(
    ppf.all_chain_proxies.pdb_inp,
    ppf.all_chain_proxies.pdb_hierarchy,
    )

  if not hetero_charges:
    # some defaults
    hetero_charges = charges.default_ion_charges
  inter_residue_bonds = charges.get_inter_residue_bonds(ppf)
  #
  # extends side chains and add hydrogens
  #
  if model_completion:
    use_capping_hydrogens=False
    fname = 'complete'
  else:
    use_capping_hydrogens=True
    fname = 'capping'
    #assert 0 # model has H
  ppf = completion.complete_pdb_hierarchy(
    ppf.all_chain_proxies.pdb_hierarchy,
    ppf.geometry_restraints_manager(),
    pdb_filename=pdb_filename,
    pdb_inp=ppf.all_chain_proxies.pdb_inp,
    use_capping_hydrogens=use_capping_hydrogens,
    append_to_end_of_model=append_to_end_of_model,
    use_reduce=use_reduce
  )

  if calculate_charge:
    new_pdb_filename = hierarchy_utils.write_hierarchy(
      pdb_filename, # uses to get output filename
      ppf.all_chain_proxies.pdb_inp,
      ppf.all_chain_proxies.pdb_hierarchy,
      fname)
    ## need now inter_residue_bonds because of added hydrogens
    ##  maybe update the bond table!!!
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=new_pdb_filename)
    inter_residue_bonds = charges.get_inter_residue_bonds(ppf, verbose=True)
    total_charge = charges.calculate_pdb_hierarchy_charge(
      ppf.all_chain_proxies.pdb_hierarchy,
      hetero_charges=hetero_charges,
      inter_residue_bonds=inter_residue_bonds,
    )
    print("total_charge",total_charge)

  # Idealize H as riding
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  h = ppf.all_chain_proxies.pdb_hierarchy
  asc = h.atom_selection_cache()
  sel = asc.selection("element H or element D")
  model = mmtbx.model.manager(
    model_input               = None,
    pdb_hierarchy             = ppf.all_chain_proxies.pdb_hierarchy,
    crystal_symmetry          = ppf.all_chain_proxies.pdb_inp.crystal_symmetry(),
    log                       = null_out())
  if 0:
    ph=model.get_hierarchy()
    ph.write_pdb_file('test.pdb')
  model.process(make_restraints=True, grm_normalization=True,
    pdb_interpretation_params = params)
  model.idealize_h_riding()
  if neutron_option=='all_d':
    from mmtbx.ligands.ready_set_utils import perdeuterate_model_ligands
    perdeuterate_model_ligands(ppf.all_chain_proxies.pdb_hierarchy)
  model.set_occupancies(hydrogen_atom_occupancies, selection=sel)

  ## after no error getting total charge, write the completed pdb file
  hierarchy_utils.write_hierarchy(pdb_filename, # uses to get output filename
                                  ppf.all_chain_proxies.pdb_inp,
                                  model.get_hierarchy(),
                                  fname)
  if not skip_validation:
    final_model_statistics = model.geometry_statistics()

if __name__=="__main__":
  def _fake_phil_parse(arg):
    def _boolean(s):
      if s.lower() in ['1', 'true']: return True
      elif s.lower() in ['0', 'false']: return False
      else: assert 0
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
  run(*tuple(args), **kwds)

#  args = sys.argv[1:]
#  del sys.argv[1:]
#  run(*tuple(args))
