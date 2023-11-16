from __future__ import division
from __future__ import absolute_import

import time
import os
import os.path
import libtbx.load_env
from libtbx.test_utils import approx_equal
import iotbx.pdb
from mmtbx.pair_interaction import pair_interaction
from qrefine.tests.unit import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def normalize(x):
  r = []
  for x_ in x:
    x_ = list(x_)
    x_.sort()
    r.append(x_)
  r.sort()
  return r

def run(prefix):
  """
  Exercise buffer region of cluster.
  """
  #
  pdb_inp = iotbx.pdb.input(file_name= os.path.join(qr_unit_tests,"data_files","2lvr.pdb"))
  ph = pdb_inp.construct_hierarchy()
  # Clusters
  bc_clusters = [[1, 2],
                  [3, 4, 5, 13, 14, 15, 16, 17, 18, 19, 20, 21],
                  [6, 7, 8, 9, 10, 11, 12, 22, 23, 26, 31],
                  [24, 25, 27, 28, 29, 30]]

  # Framgents = cluster + buffer
  bc_qms = [[1, 2, 3, 4], [1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17,
                           18, 19, 20, 21, 22, 23, 24, 25],
            [3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31],
            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]]


  check_buffer(ph.deep_copy(), bc_clusters, bc_qms)

def check_buffer(ph, clusters, qms):
  #
  qms_calculated = []
  for i in range(len(clusters)):
    fragment = pair_interaction.run(ph.deep_copy(), core=clusters[i])[2]
    qms_calculated.append(fragment)
  qms = normalize(qms)
  qms_calculated = normalize(qms_calculated)
  assert approx_equal(qms, qms_calculated)

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
