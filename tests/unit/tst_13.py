from __future__ import division

import os
import sys
import random
import time
import numpy
import iotbx.pdb
import mmtbx.command_line
import libtbx.load_env
from libtbx.test_utils import approx_equal

from ase.io import write
from ase.io import read as ase_io_read
from ase.optimize.lbfgs import LBFGS
from scitbx.array_family import flex
from qrefine.restraints import from_qm
from qrefine.fragment import fragments
from qrefine.cluster_restraints import from_cluster
from qrefine.clustering import betweenness_centrality_clustering

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests/unit/")

master_params_str ="""
method = *multiprocessing pbs sge lsf threading
.type = choice(multi=False)
nproc = 1
.type = int
qsub_command = None
.type = str
"""

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=master_params_str)

# gradient-only LBFGS without line search
class lbfgs_gradient(object):
  def __init__(self, atoms, restraints):
    self.restraints = restraints
    self.opt = LBFGS(atoms=atoms)

  def step(self):
    pos = self.opt.atoms.get_positions()
    sites_cart = flex.vec3_double(pos)
    e, g = self.restraints.target_and_gradients(sites_cart)
    forces = numpy.array(g) * -1
    self.opt.step(forces)

  def write(self, file):
    write(file, self.opt.atoms)

  def run(self, nstep):
    for i in range(nstep):
      self.step()

# fix random seed in this script
random.seed(0)
flex.set_random_seed(0)

def run(prefix = "tst_13"):
  """
  compare gradients from entire qm and clustered qm.
  """
  pdb_inp = iotbx.pdb.input(os.path.join(qr_unit_tests,"data_files/helix.pdb"))
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()

  ## compare the absolute value of gradients
  g_entire = qm_gradient(cs, ph, clustering=False)
  g_cluster = qm_gradient(cs, ph, clustering = True)
  assert approx_equal(list(g_entire.as_double()),list(g_cluster.as_double()) , g_entire*0.05)
  ## compare the geometry rmsd after 5 steps optimization
  file_entire_qm = "entire_qm.pdb"
  qm_opt(cs, ph,file_entire_qm,cluster=False)
  file_cluster_qm = "cluster_qm.pdb"
  qm_opt(cs, ph,file_cluster_qm,cluster = True)
  sites_cart_entire_qm = iotbx.pdb.input(file_entire_qm).atoms().extract_xyz()
  sites_cart_cluster_qm = iotbx.pdb.input(file_cluster_qm).atoms().extract_xyz()
  rmsd_diff = sites_cart_entire_qm.rms_difference(sites_cart_cluster_qm)
  os.system("rm %s %s"%(file_entire_qm, file_cluster_qm))
  os.system("rm -rf ase")
  assert rmsd_diff<0.02

def qm_gradient(cs, ph, clustering=False):
  fq = from_qm(
    charge_embedding=False,
    pdb_hierarchy=ph,
    qm_engine_name="mopac",
    crystal_symmetry=cs,
    clustering=clustering,
    maxnum_residues_in_cluster=8)

  fm = fragments(
   working_folder             = os.path.split("./ase/tmp_ase.pdb")[0]+ "/",
   clustering_method          = betweenness_centrality_clustering,
   maxnum_residues_in_cluster = 8,
   charge_embedding           = False,
   pdb_hierarchy              = ph,
   qm_engine_name             = "mopac",
   crystal_symmetry           = cs)

  #pass qm to cluster
  cluster = from_cluster(
   restraints_manager = fq,
   fragment_manager =fm,
   parallel_params = get_master_phil().extract())

  sites_cart = ph.atoms().extract_xyz()
  e,g = cluster.target_and_gradients(sites_cart)
  return g

def qm_opt(cs, ph, file, cluster=False):
  fq = from_qm(
    charge_embedding=False,
    pdb_hierarchy=ph,
    qm_engine_name="mopac",
    crystal_symmetry=cs,
    clustering=cluster,
    maxnum_residues_in_cluster=8)



  sys = ase_io_read(os.path.join(qr_unit_tests,"data_files/helix.pdb"))
  opt = lbfgs_gradient(sys, fq)
  opt.run(5)
  opt.write(file)

if(__name__ == "__main__"):
  t0 = time.time()
  log = sys.stdout
  run()
  print "Time: %6.2f"%(time.time()-t0)
  print "OK"
