from __future__ import division
from __future__ import absolute_import

import time, os
from libtbx.test_utils import approx_equal
import qrefine.clustering as clustering
from qrefine.tests.unit import run_tests

def run(prefix):
  """
  Exercise interaction graph clustering.
  """
  interaction_list = [[24, 27], [19, 22], [5, 11], [18, 22], [25, 26], [23, 26], [2, 3], [5, 7], [9, 11], [27, 29],
                   [5, 10], [9, 10], [18, 21], [11, 12], [24, 25], [5, 12], [6, 10], [6, 7], [20, 21], [10, 11],
                   [21, 25], [5, 6], [21, 24], [17, 21], [6, 11], [12, 13], [5, 13], [17, 20], [4, 5], [3, 12],
                   [29, 30], [8, 9], [20, 24], [11, 13], [4, 13], [5, 19], [28, 29], [6, 9], [6, 19], [6, 8],
                   [13, 19], [7, 19], [6, 13], [13, 18], [15, 19], [7, 8], [4, 19], [16, 19], [6, 22], [15, 18],
                   [15, 16], [9, 31], [6, 31], [26, 27], [13, 14], [11, 22], [18, 19], [4, 15], [17, 18], [14, 15],
                   [16, 17], [22, 31], [8, 26], [4, 16], [8, 23], [15, 17], [6, 23], [7, 23], [13, 15], [22, 23],
                   [2, 4], [13, 22], [23, 24], [19, 23], [10, 12], [19, 20], [24, 30], [21, 22], [9, 26], [23, 27],
                   [4, 14], [1, 2], [16, 20], [26, 31], [25, 28], [27, 28], [22, 26], [24, 28], [20, 23], [17, 19],
                   [27, 30], [16, 18], [20, 22], [1, 3], [6, 26], [28, 30], [3, 13], [3, 5], [3, 4], [22, 25],
                   [3, 14], [9, 22]]
  gn_clusters = [[1, 2],
                  [3, 4, 5, 14, 13],
                  [6, 7, 8, 9, 22, 23, 26, 31],
                  [10, 11, 12],
                  [15, 16, 17, 18, 19, 20, 21],
                  [24, 25, 27, 28, 29, 30]]
  bc_clusters = [[1, 2],
                  [3, 4, 5, 13, 14, 15, 16, 17, 18, 19, 20, 21],
                  [6, 7, 8, 9, 10, 11, 12, 22, 23, 26, 31],
                  [24, 25, 27, 28, 29, 30]]
  for e in interaction_list:
    e.sort()
  interaction_list.sort()
  cc=clustering.betweenness_centrality_clustering(interaction_list)
  assert approx_equal(cc.get_clusters(), bc_clusters)

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
