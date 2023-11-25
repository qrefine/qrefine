from __future__ import division
import os
import os.path
import libtbx.load_env
from libtbx.easy_mp import parallel_map
from abc import ABCMeta, abstractmethod

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")


qr_reg_data_finalise = os.path.join(qrefine_path, "tests/regression/datasets/finalise")
qr_reg_data_cluster  = os.path.join(qrefine_path, "tests/regression/datasets/cluster")



class test_base:
  __metaclass__ = ABCMeta

  def __init__(self):
    pass

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

  #TODO compare with  /home/xuyanting/test_prime/elbow.py
  def run(self):
    print("running regression tests")
    test_results = []

    test_results.append(os.system(f"iotbx.fetch_pdb 3dtj"))
    test_results.append(os.system(f"iotbx.fetch_pdb 1f8t"))


    # Finalise
    for filename in os.listdir(qr_reg_data_finalise):
      test_results.append(os.system(f"qr.finalise {qr_reg_data_finalise}/{filename}"))


    # Clustering
    for filename in os.listdir(qr_reg_data_cluster):
      test_results.append(os.system(f"qr.cluster {qr_reg_data_cluster}/{filename}"))
      test_results.append(os.system(f"qr.fragment {qr_reg_data_cluster}/{filename}"))




    for test_result in test_results:
      self.check_assertions(test_results)
