# LIBTBX_SET_DISPATCHER_NAME phenix.development.ready_set
import os
import sys

import iotbx
import iotbx.pdb
from libtbx import easy_run

import charges
import completion
from utils import hierarchy_utils

from skip import skip

def remove_alt_loc(hierarchy):
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

def run_ready_set(pdb_filename):
  from StringIO import StringIO
  assert pdb_filename.find('.pdb')>-1, 'ReadySet! only works on PDB, not CIF'
  cmd = 'phenix.ready_set %s  remove_water=true ' % pdb_filename
  print 'Running ReadySet!'
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO()
  ero.show_stdout(std)
  for line in std.getvalue().splitlines():
    print line
  return os.path.basename(pdb_filename).replace('.pdb', '.updated.pdb')
  # maybe read from stdout

def run_fetch_pdb(code):
  cmd = 'phenix.fetch_pdb %s' % code
  print 'Fetching files'
  print cmd
  easy_run.call(cmd)

def loop_over_dir(d):
  i=0
  for filename in os.listdir(d):
    if not filename.endswith('.pdb'): continue
    i+=1
    print '%s\n %3d %s\n%s' % ('*'*80, i, os.path.join(d, filename), '*'*80)
    if filename in skip:
      print 'skipping'
      continue
    if os.path.exists(filename.replace('.pdb', '.updated.pdb')):
      run(filename.replace('.pdb', '.updated.pdb'))
    else:
      run(os.path.join(d, filename))

def run(pdb_filename, model_completion=True):
  print "run",pdb_filename
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

  if 0: # moved adding hydrogens to complete_pdb_hierarchy
    if pdb_filename.endswith('.updated.pdb'):
      pdb_filename_h = pdb_filename
    else:
      pdb_filename_h = pdb_filename.replace('.pdb', '.updated.pdb')
    if not os.path.exists(pdb_filename_h):
      pdb_filename_h = os.path.basename(pdb_filename_h)
    print 'pdb_filename_h',pdb_filename_h
    if not os.path.exists(pdb_filename_h):
      pdb_filename = run_ready_set(pdb_filename)
    else: pdb_filename = pdb_filename_h

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

  # should use cctbx
  hierarchy = remove_alt_loc(ppf.all_chain_proxies.pdb_hierarchy)
  if 0:
    output = hierarchy_utils.write_hierarchy(
      pdb_filename, # uses to get output filename
      ppf.all_chain_proxies.pdb_inp,
      hierarchy,
      "temp0")
    ppf = hierarchy_utils.get_processed_pdb(pdb_filename=output)
  else:
    raw_records = hierarchy_utils.get_raw_records(
      ppf.all_chain_proxies.pdb_inp,
      hierarchy,
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
  )
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
  print "total_charge",total_charge
  ## after no error getting total charge, write the completed pdb file
  hierarchy_utils.write_hierarchy(pdb_filename, # uses to get output filename
                                  ppf.all_chain_proxies.pdb_inp,
                                  ppf.all_chain_proxies.pdb_hierarchy,
                                  fname)

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
  print '''
  model/cluster completion'
    args : %(args)s
    kwds : %(kwds)s
  ''' % locals()
  run(*tuple(args), **kwds)

#  args = sys.argv[1:]
#  del sys.argv[1:]
#  run(*tuple(args))
