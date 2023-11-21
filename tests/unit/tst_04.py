from __future__ import division
from __future__ import print_function

import os
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import mmtbx.model
from qrefine import qr, refine
from qrefine.utils import hierarchy_utils
from libtbx.utils import null_out
from qrefine.tests.unit import run_tests

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

def get_restraints_manager(expansion, clustering, file_name):
  model = get_model(file_name=file_name)
  params = qr.get_default_params()
  params.restraints="cctbx"
  params.expansion = expansion
  params.cluster.clustering=clustering
  return refine.create_restraints_manager(params=params, model=model), \
         model.get_sites_cart()

def run(prefix):
  """
  Exercise expansion=False / expansion=True
  """
  clustering = True
  done = False
  path = qr_unit_tests+"/data_files/"
  for fn in os.listdir(path):
    if(not fn=="1yjp.pdb"): continue
    if not fn.endswith(".pdb"): continue
    fn = path + fn
    ph = iotbx.pdb.input(fn).construct_hierarchy()
    if list(ph.altloc_indices()) != ['']: continue
    print(fn)
    #
    rm1, sites_cart = get_restraints_manager(
      expansion=False, clustering=clustering, file_name=fn)
    t1, g1 = rm1.target_and_gradients(sites_cart = sites_cart)
    #
    rm2, sites_cart = get_restraints_manager(
      expansion=True, clustering=clustering, file_name=fn)
    t2, g2 = rm2.target_and_gradients(sites_cart = sites_cart)
    #
    diff = flex.abs(g1.as_double()-g2.as_double())
    print(diff.min_max_mean().as_tuple())
    #
    assert flex.max(diff) < 1.e-6, flex.max(diff)
    #
    done = True
  assert done

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
