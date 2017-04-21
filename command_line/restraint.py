from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.restraint
import sys
import time
import os.path
import argparse
import iotbx.pdb
import libtbx.load_env
from scitbx.array_family import flex
from qrefine.core.restraints import from_qm

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

log = sys.stdout

def example():
  print >> log, "No pdb specified, using helix example"
  example_pdb = os.path.join(qrefine_path,"examples/1us0/a87_99_h.pdb")
  run(example_pdb)

def run(pdb_file):
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  sites_cart = ph.atoms().extract_xyz()
  fq = from_qm(
    calculator_name= "mopac",
    pdb_hierarchy=ph,
    crystal_symmetry=cs,
    use_cluster_qm=False)
  energy,gradients = fq.target_and_gradients(sites_cart)
  print >> log,"Energy: ", energy
  print >> log,"Gradients: "
  for gradient in list(gradients):
    print >> log, gradient

if (__name__ == "__main__"):
  print "Restraint for Q|R"
  parser = argparse.ArgumentParser(description='Calculate restraint for Q|R')
  parser.add_argument('--qm', action='store_true',
                      default=False,
                      help='compute the energy and gradient using a QM calculator ')
  parser.add_argument('--cluster', action='store_true',
                      default=False,
                      help='construct a set of clusters, and then calculate the combined gradient')  # nightly build?
  parser.add_argument('--cctbx', action='store_true',
                      default=False,
                      help='''compute the standard cctbx restraint''')
  parser.add_argument('--all', action='store_true',
                      default=False,
                      help='''run the set of restraints for comparison''')
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  if len(args) == 1:
    run(args[0])
  else:
    example()
  print >> log, "Time: %6.4f" % (time.time() - t0)
