from __future__ import division
from __future__ import print_function

import os
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import mmtbx.model
from qrefine import qr

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def get_model(file_name):
  file_name = os.path.join(qr_unit_tests,"data_files",file_name)
  pdb_inp = iotbx.pdb.input(file_name)
  model = qr.process_model_file(
    pdb_file_name = file_name,
    cif_objects = None,
    crystal_symmetry=pdb_inp.crystal_symmetry())
  return model

def get_restraints_manager(expansion, file_name):
  model = get_model(file_name=file_name)
  params = qr.get_master_phil().extract()
  params.restraints="cctbx"
  params.expansion = expansion
  return qr.create_restraints_manager(params=params, model=model), \
         model.model.get_sites_cart()

def run():
  path = qr_unit_tests+"/data_files/"
  for fn in os.listdir(path):
    if not fn.endswith(".pdb"): continue
    fn = path + fn
    ph = iotbx.pdb.input(fn).construct_hierarchy()
    if list(ph.altloc_indices()) != ['']: continue
    print(fn)
    #
    rm1, sites_cart = get_restraints_manager(expansion=False, file_name=fn)
    t1, g1 = rm1.target_and_gradients(sites_cart = sites_cart)
    #
    rm2, sites_cart = get_restraints_manager(expansion=True, file_name=fn)
    t2, g2 = rm2.target_and_gradients(sites_cart = sites_cart)
    #
    diff = flex.abs(g1.as_double()-g2.as_double())
    print(diff.min_max_mean().as_tuple())
    assert flex.max(diff) < 1.e-6

if(__name__ == "__main__"):
    run()
