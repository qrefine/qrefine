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

def run(model,
        model_completion=True,
        skip_validation=False,
        append_to_end_of_model=False,
        neutron_option=None,
        hydrogen_atom_occupancies=0.,
        use_reduce=True,
        remove_selection=None,
        ):
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

  MC = completion.complete_pdb_hierarchy(
    hierarchy                   = model.get_hierarchy(),
    geometry_restraints_manager = model.get_restraints_manager(),
    crystal_symmetry            = model.crystal_symmetry(),
    use_capping_hydrogens       = use_capping_hydrogens,
    append_to_end_of_model      = append_to_end_of_model,
    use_reduce                  = use_reduce
  )

  # Idealize H as riding
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  h = MC.get_hierarchy()
  asc = h.atom_selection_cache()
  sel = asc.selection("element H or element D")
  model = mmtbx.model.manager(
    model_input               = None,
    pdb_hierarchy             = MC.get_hierarchy(),
    crystal_symmetry          = MC.crystal_symmetry(),
    log                       = null_out())
  model.process(make_restraints=True, grm_normalization=True,
    pdb_interpretation_params = params)
  model.idealize_h_riding()
  hierarchy=model.get_hierarchy()

  if neutron_option=='all_h':
    pass
  elif neutron_option=='all_d':
    from mmtbx.ligands.ready_set_utils import perdeuterate_model_ligands
    perdeuterate_model_ligands(hierarchy)
  else:
    # set for "hd_and_h"
    from mmtbx.hydrogens.neutron_utils import neutron_exchange_hydrogens
    neh_kwds = {"exchange_sites_only" : True,
                "perdeuterate" : False,
                # "cifs" : results.get("cifs", None),
                }
    if neutron_option=="all_hd":
      neh_kwds["exchange_sites_only"] = False
    elif neutron_option=="hd_and_d":
      neh_kwds["perdeuterate"]=True
    hierarchy = neutron_exchange_hydrogens(hierarchy,
                                           **neh_kwds)
    model=mmtbx.model.manager(
    model_input               = None,
    pdb_hierarchy             = hierarchy,
    crystal_symmetry          = model.crystal_symmetry(),
    log                       = null_out())
    h=model.get_hierarchy()
    asc = h.atom_selection_cache()
    sel = asc.selection("element H or element D")

  model.set_occupancies(hydrogen_atom_occupancies, selection=sel)
  return model

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
