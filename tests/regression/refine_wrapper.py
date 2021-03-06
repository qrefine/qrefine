from __future__ import division
import iotbx.pdb

class Driver(object):

  def create(self,pdb):
      self.pdb = pdb
      self.pdb_inp = iotbx.pdb.input(self.pdb)
      self.ph = self.pdb_inp.construct_hierarchy()
      self.cs = self.pdb_inp.crystal_symmetry()
      self.sites_cart = self.ph.atoms().extract_xyz()

  def refine(self,pdb):
    self.driver(self.sites_cart)
    return Result(self.pdb_code,
                  self.rsmd_inital,
                  self.rsmd_final,
                  self.r_start,
                  self.r_work)

  def __call__(self,pdb):
    self.create()
    return self.refine()

class Result(object):
  def __init__(self,
                 pdb_code,
                 rsmd_inital,
                 rsmd_final,
                 r_start,
                 r_work):
    self.pdb_code    = pdb_code
    self.rsmd_inital = rsmd_inital
    self.rsmd_final  = rsmd_final
    self.r_start     = r_start
    self.r_work      = r_work
