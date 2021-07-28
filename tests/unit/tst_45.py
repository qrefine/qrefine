from __future__ import division

import os
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import mmtbx.model
from qrefine import qr
from qrefine.fragment import fragment_extracts
from qrefine import cluster_restraints
import run_tests

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
  if(not expansion):
    params.cluster.clustering=True
  else:
    params.cluster.clustering=False
  params.quantum.nproc=1
  params.parallel.nproc=1
  params.quantum.method='PM7'
  params.quantum.basis=''
  params.cluster.two_buffers=False
  params.cluster.altloc_method="subtract"
  #params.cluster.maxnum_residues_in_cluster=6
  params.cluster.maxnum_residues_in_cluster=999

  result = qr.create_restraints_manager(params=params, model=model), \
         model.model.get_sites_cart()
  return result

def run(prefix, verbose=False):
  path = qr_unit_tests+"/data_files/"
  files = ["gly2_1.pdb", "altlocs2.pdb", "altlocs.pdb", "gly2_2.pdb"]
  for f in files:
    fn = path + f
    ph = iotbx.pdb.input(fn).construct_hierarchy()
    #
    if(verbose): print "expansion=False "
    rm1, sites_cart = get_restraints_manager(expansion=False, file_name=fn)
    t1, g1 = rm1.target_and_gradients(sites_cart = sites_cart)
    #
    if(verbose): print "expansion=True "
    rm2, sites_cart = get_restraints_manager(expansion=True, file_name=fn)
    t2, g2 = rm2.target_and_gradients(sites_cart = sites_cart)
    #
    if(verbose):
      atoms = ph.atoms()
      ds = flex.sqrt((g1 - g2).dot())
      for d, g, gg, dist, a in zip((g1-g2), g1, g2, ds, atoms):
        print ["%8.4f"%i for i in d], \
              ["%8.4f"%i for i in g], \
              ["%8.4f"%i for i in gg], "%8.4f"%dist, a.quote()
    #
    rs = flex.double()
    for a, b in zip(g1.as_double(), g2.as_double()):
      r = abs(abs(a)-abs(b))/(abs(a)+abs(b))*2.*100
      rs.append(r)
    assert flex.max(flex.abs(rs)) < 1.e-6

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
