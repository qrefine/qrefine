from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import random
import time
import numpy as np
import iotbx.pdb
import mmtbx.command_line
import libtbx.load_env
from libtbx.test_utils import approx_equal
from qrefine.tests.unit import run_tests

from ase.io import write
from ase.io import read as ase_io_read
from ase.optimize.lbfgs import LBFGS
from scitbx.array_family import flex
from qrefine.restraints import from_qm
from qrefine.fragment import fragments
from qrefine.cluster_restraints import from_cluster
from qrefine.clustering import betweenness_centrality_clustering

'''
tests the cluster gradient. Modified to use two_buffers Feb2019.
'''

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

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

def approx_equal2(v1,v2,tol):
  # approx_equal between two vectors where the tolerance is multiplied
  # with the vector v2
  happy_bob=True 
  v1=np.array(v1)
  v2=np.array(v2)
  vtol=abs(v1*tol)
  for i in range(v1.shape[0]):
    error=abs(v1[i]-v2[i])
    if( error>=vtol[i] ):
      print('oh dear!',i,error,vtol[i])
      happy_bob=False
  return happy_bob

# gradient-only LBFGS without line search
class lbfgs_gradient(object):
  def __init__(self, atoms, restraints):
    self.restraints = restraints
    self.opt = LBFGS(atoms=atoms)

  def step(self):
    pos = self.opt.atoms.get_positions()
    sites_cart = flex.vec3_double(pos)
    e, g = self.restraints.target_and_gradients(sites_cart)
    forces = np.array(g) * -1
    self.opt.step(forces)

  def write(self, file):
    write(file, self.opt.atoms)

  def run(self, nstep):
    for i in range(nstep):
      self.step()

# fix random seed in this script
random.seed(0)
flex.set_random_seed(0)

def run(prefix):
  """
  compare gradients from entire qm and clustered qm.
  """
  pdb_inp = iotbx.pdb.input(os.path.join(qr_unit_tests,"data_files","helix.pdb"))
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  restraints_entire = generate_restraints(cs, ph, clustering=False)
  ## compare the absolute value of gradients
  g_entire = qm_gradient(ph, restraints_entire)
  restraints_cluster = generate_restraints(cs, ph, clustering=True)
  g_cluster = qm_gradient(ph, restraints_cluster)
  assert approx_equal2(list(g_entire.as_double()),list(g_cluster.as_double()),0.05)
  #old & wrong, but kept to know the original intention: assert approx_equal(list(g_entire.as_double()),list(g_cluster.as_double()) , g_entire*0.05)
  ## compare the geometry rmsd after 5 steps optimization
  file_entire_qm = "entire_qm.pdb"
  qm_opt(restraints_entire,file_entire_qm)
  file_cluster_qm = "cluster_qm.pdb"
  qm_opt(restraints_cluster,file_cluster_qm)
  sites_cart_entire_qm = iotbx.pdb.input(file_entire_qm).atoms().extract_xyz()
  sites_cart_cluster_qm = iotbx.pdb.input(file_cluster_qm).atoms().extract_xyz()
  rmsd_diff = sites_cart_entire_qm.rms_difference(sites_cart_cluster_qm)
  #os.system("rm %s %s"%(file_entire_qm, file_cluster_qm))
  #os.system("rm -rf ase")
  assert rmsd_diff<0.02

def generate_restraints(cs, ph, clustering=False):
  fq = from_qm(
    pdb_hierarchy=ph,
    qm_engine_name="mopac",
    crystal_symmetry=cs,
    clustering=clustering)
  if(clustering):
    fm = fragments(
     working_folder             = os.path.split("./ase/tmp_ase.pdb")[0]+ "/",
     clustering_method          = betweenness_centrality_clustering,
     maxnum_residues_in_cluster = 8,
     charge_embedding           = False,
     two_buffers                = True,
     pdb_hierarchy              = ph,
     qm_engine_name             = "mopac",
     fast_interaction           = True,
     crystal_symmetry           = cs)
    restraints = from_cluster(
     restraints_manager = fq,
     fragment_manager =fm,
     parallel_params = get_master_phil().extract())
  else:
    restraints = fq
  return restraints

def qm_gradient(ph, restraints):
  sites_cart = ph.atoms().extract_xyz()
  e,g = restraints.target_and_gradients(sites_cart)
  return g

def qm_opt(restraints, file):
  sys = ase_io_read(os.path.join(qr_unit_tests,"data_files/helix.pdb"))
  opt = lbfgs_gradient(sys, restraints)
  opt.run(5)
  print("AFTER RUN")
  opt.write(file)

if(__name__ == "__main__"):
  """
  If this test hangs then MOPAC needs to be updated. Run MOPAC command for 
  update instructions.
  """
  disable = False
  # if(os.environ.get("MOPAC_COMMAND") is None): disable = True
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=disable)
