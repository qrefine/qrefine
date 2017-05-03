from __future__ import division
import os
import sys
import time
import os.path
import iotbx.pdb
import libtbx.load_env
from libtbx.easy_mp import parallel_map
from qrefine.core.restraints import from_qm
from pymongo import MongoClient
from abc import ABCMeta, abstractmethod

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")


qsub_command = 'qsub  -N reg_test_cluster -m ae -q fat  -l nodes=1:ppn=32'

class test_base:
  __metaclass__ = ABCMeta
  
  def __init__(self):
    self.db = MongoClient('localhost', 27017).pyoink


  @abstractmethod
  def check_assertions(self):
    pass

  @abstractmethod
  def insert(self):
    pass

  def pdbs(self,qr_reg_data):
    pdbs=[]
    for pdb_file in os.listdir(qr_reg_data):
      pdbs.append(os.path.join(qr_reg_data,pdb_file))
    return pdbs

  def run(self):
    test_results =parallel_map(
      func=self.func,
      iterable=self.pdbs,
      method='pbs',
      preserve_exception_message=True,
      processes=len(self.pdbs),
      qsub_command=qsub_command,
      use_manager=True)
    for test_result in test_results:
      check_assertions(test_results)

