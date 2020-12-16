from __future__ import division
from __future__ import print_function
import iotbx.pdb
from qrefine.restraints import from_qm

class Restraints(object):

  def create(self):
    print(" creating ")
    self.pdb_inp = iotbx.pdb.input(self.pdb)
    self.ph = self.pdb_inp.construct_hierarchy()
    self.cs = self.pdb_inp.crystal_symmetry()
    self.sites_cart = self.ph.atoms().extract_xyz()
    self.manager = from_qm(
             use_cluster_qm             = True,
             pdb_hierarchy              = self.ph,
             crystal_symmetry           = self.cs,
             maxnum_residues_in_cluster = int(self.maxnum_residues_in_cluster)
             )

  def process(self,pdb):
    energy, gradients = self.manager.target_and_gradients(self.sites_cart)
    return Result(self.pdb_code,energy, gradients)

  def __call__(self,pdb):
    print("calling", pdb)
    self.pdb= pdb
    self.create(pdb)
    return self.process(pdb)


class Result(object):
  def __init__(self, pdb_code, energy, gradients):
    self.pdb_code = pdb_code
    self.create()
    self.energy = energy
    self.gradients = gradients
