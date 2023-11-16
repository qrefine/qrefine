from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import mmtbx.model
from qrefine import qr, refine
from qrefine.fragment import fragment_extracts
from qrefine import cluster_restraints
from qrefine.tests.unit import run_tests
from qrefine.utils import hierarchy_utils
from libtbx.utils import null_out

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def get_model(file_name):
  file_name = os.path.join(qr_unit_tests,"data_files",file_name)
  pdb_inp = iotbx.pdb.input(file_name)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  model.process(make_restraints = True)
  return model

def get_restraints_manager(expansion, file_name, altloc_method):
  model = get_model(file_name=file_name)
  params = qr.get_default_params()
  params.restraints="cctbx"
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
  params.cluster.altloc_method=altloc_method
  params.cluster.bond_with_altloc = False # False to agree with Min's code
  #params.cluster.maxnum_residues_in_cluster=6
  params.cluster.maxnum_residues_in_cluster=999

  result = refine.create_restraints_manager(params=params, model=model), \
         model.get_sites_cart()
  return result

def run(prefix, verbose=False):
  """
  Exercise expansion=False / expansion=True
  
  XXX TEST FAILS (expansion=False / expansion=True)
  
  """
  for altloc_method in ["average", "subtract"]:
    path = qr_unit_tests+"/data_files/"
    files = ["gly2_1.pdb", "altlocs2.pdb", "altlocs.pdb", "gly2_2.pdb"]
    for f in files:
      print(f)
      fn = path + f
      ph = iotbx.pdb.input(fn).construct_hierarchy()
      #
      if(verbose): print("expansion=False ")
      rm1, sites_cart = get_restraints_manager(
        expansion=False, file_name=fn, altloc_method=altloc_method)
      t1, g1 = rm1.target_and_gradients(sites_cart = sites_cart)
      #
      if(verbose): print("expansion=True ")
      rm2, sites_cart = get_restraints_manager(
        expansion=True, file_name=fn, altloc_method=altloc_method)
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
      result = flex.max(flex.abs(rs))
      assert result < 1.e-6, result

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=True)
