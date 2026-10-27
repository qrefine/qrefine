from __future__ import division
from __future__ import print_function

import os
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
import mmtbx.model
from qrefine import qr, refine

from qrefine import  cluster_restraints, restraints
from qrefine.utils import hierarchy_utils
from libtbx.utils import null_out
from qrefine.tests.unit import run_tests
from qrefine.command_line import granalyse

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def get_model(file_name):
  file_name = os.path.join(qr_unit_tests,"data_files",file_name)
  pdb_inp = iotbx.pdb.input(file_name)
  return restraints.model_from_hierarchy(
    pdb_hierarchy    = pdb_inp.construct_hierarchy(),
    crystal_symmetry = pdb_inp.crystal_symmetry())

def get_restraints_manager(expansion, file_name, two_buffers=False,
                           distance_buffer=False, distance_buffer_radius=5.0):
  model = get_model(file_name=file_name)
  mpar = mmtbx.model.manager.get_default_pdb_interpretation_params()
  mpar.pdb_interpretation.const_shrink_donor_acceptor=0.6
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = mpar)
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
  params.cluster.two_buffers=two_buffers
  params.cluster.distance_buffer=distance_buffer
  params.cluster.distance_buffer_radius=distance_buffer_radius
  params.cluster.maxnum_residues_in_cluster=2

  result = refine.create_restraints_manager(params=params, model=model,
         altlocs_present=False), model.get_sites_cart()
  return result

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py",""),
        two_buffers=False, distance_buffer=False, distance_buffer_radius=5.0):
  """
  Exercise expansion=False / expansion=True
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  print("two_buffers=%s distance_buffer=%s distance_buffer_radius=%s" % (
    two_buffers, distance_buffer, distance_buffer_radius))
  #
  path = qr_unit_tests+"/data_files/"
  files = ["m2_complete_box_large.pdb",
           "p212121.pdb",
           "p1_box_large.pdb",
           "p1_box_small.pdb",
           ]
  for f in files:
    print(f, "*"*30)
    fn = path + f
    ph = iotbx.pdb.input(fn).construct_hierarchy()
    if list(ph.altloc_indices()) != ['']: continue
    #
    rm1, sites_cart = get_restraints_manager(expansion=False, file_name=fn,
      two_buffers            = two_buffers,
      distance_buffer        = distance_buffer,
      distance_buffer_radius = distance_buffer_radius)
    t1, g1 = rm1.target_and_gradients(sites_cart = sites_cart)
    #
    rm2, sites_cart = get_restraints_manager(expansion=True, file_name=fn)
    t2, g2 = rm2.target_and_gradients(sites_cart = sites_cart)
    #
    if 0:
      atoms = ph.atoms()
      ds = flex.sqrt((g1 - g2).dot())
      for d, g, gg, dist, a in zip((g1-g2), g1, g2, ds, atoms):
        print(["%8.4f"%i for i in d], \
              ["%8.4f"%i for i in g], \
              ["%8.4f"%i for i in gg], "%8.4f"%dist, a.quote())
    #
    # Per-atom absolute differences |g1-g2| vs. gradient magnitudes |g1|
    diffs = flex.sqrt((g1-g2).dot())
    norms = flex.sqrt(g1.dot())
    print(f, "|g1-g2| per atom min/max/mean:", diffs.min_max_mean().as_tuple())
    print(f, "|g1|    per atom min/max/mean:", norms.min_max_mean().as_tuple())

    d = max(granalyse.get_grad_wdelta(ref=g1, g=g2))
    print(d, "get_grad_wdelta")
    assert d < 1.5, d

if(__name__ == "__main__"):
  run(distance_buffer=True, distance_buffer_radius=10.0)
