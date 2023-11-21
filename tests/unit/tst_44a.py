from __future__ import division
from __future__ import print_function

import os
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import mmtbx.model
from qrefine import qr, refine

from qrefine.fragment import fragment_extracts
from qrefine import  cluster_restraints
from qrefine.utils import hierarchy_utils
from libtbx.utils import null_out
from qrefine.tests.unit import run_tests
from qrefine.command_line import granalyse

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def get_model(file_name):
  file_name = os.path.join(qr_unit_tests,"data_files",file_name)
  pdb_inp = iotbx.pdb.input(file_name)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params)
  return model

def get_restraints_manager(expansion, file_name):
  model = get_model(file_name=file_name)
  params = qr.get_default_params()
  params.restraints="qm"
  params.expansion = expansion
  if(not expansion):
    params.cluster.clustering=True
  else:
    params.cluster.clustering=False
  params.quantum.nproc=1
  params.parallel.nproc=1
  params.quantum.method='PM7'
  params.quantum.basis=''
  params.cluster.two_buffers=False
  params.cluster.maxnum_residues_in_cluster=2

  result = refine.create_restraints_manager(params=params, model=model), \
         model.get_sites_cart()
  return result

def run(prefix, verbose=False):
  """
  Exercise expansion=False / expansion=True
  """
  path = qr_unit_tests+"/data_files/"
  files = ["m2_complete_box_large.pdb",
           "p212121.pdb", 
           "p1_box_small.pdb", 
           "p1_box_large.pdb"
           ]
  for f in files:
    fn = path + f
    if(verbose): print(fn)
    ph = iotbx.pdb.input(fn).construct_hierarchy()
    if list(ph.altloc_indices()) != ['']: continue
    #
    rm1, sites_cart = get_restraints_manager(expansion=False, file_name=fn)
    t1, g1 = rm1.target_and_gradients(sites_cart = sites_cart)
    #
    rm2, sites_cart = get_restraints_manager(expansion=True, file_name=fn)
    t2, g2 = rm2.target_and_gradients(sites_cart = sites_cart)
    #
    if(verbose): 
      atoms = ph.atoms()
      ds = flex.sqrt((g1 - g2).dot())
      for d, g, gg, dist, a in zip((g1-g2), g1, g2, ds, atoms):
        print(["%8.4f"%i for i in d], \
              ["%8.4f"%i for i in g], \
              ["%8.4f"%i for i in gg], "%8.4f"%dist, a.quote())
    #
    rs = flex.double()
    for a, b in zip(g1.as_double(), g2.as_double()):
      r = abs(abs(a)-abs(b))/(abs(a)+abs(b))*2.*100
      rs.append(r)
    if(verbose): print(f, "min/max/mean:", rs.min_max_mean().as_tuple())
    
    d = max(granalyse.get_grad_wdelta(ref=g1, g=g2))
    if(verbose): print(d, "get_grad_wdelta")
    assert d < 10. # XXX Is this reasonable?

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
