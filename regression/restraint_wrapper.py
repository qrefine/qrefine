from __future__ import division
import os
import sys
import time
import os.path
import iotbx.pdb
from qrefine.core.restraints import from_qm

class Restraints(object):
  def __init__(self,manager):
    self.manager = manager

  def create(pdb):
    self.pdb =pdb
    self.pdb_inp = iotbx.pdb.input(self.pdb)
    self.ph = self.pdb_inp.construct_hierarchy()
    self.cs = self.pdb_inp.crystal_symmetry()
    self.sites_cart = self.ph.atoms().extract_xyz()
    self.manager = manager

  def process(self,pdb):
    energy, gradients = self.manager.target_and_gradients(self.sites_cart)
    return Result(self.pdb_code,energy, gradients)

  def __call__(self,pdb):
    self.create(pdb)
    return self.process(pdb)

