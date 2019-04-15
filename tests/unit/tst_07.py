from __future__ import division

import os
import time
import iotbx.pdb
import libtbx.load_env
from libtbx.test_utils import approx_equal
from qrefine.utils import yoink_utils
from qrefine.plugin.yoink.pyoink import PYoink
import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix):
  """
  Exercise interaction graph construction.
  """
  pdb_inp = iotbx.pdb.input(file_name= os.path.join(qr_unit_tests,"data_files","2lvr.pdb"))
  ph = pdb_inp.construct_hierarchy()
  yoink_utils.write_yoink_infiles("cluster.xml",
                                  "qmmm.xml",
                                  ph,
                                  os.path.join(qrefine,"plugin","yoink","dat"))
  pyoink=PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"),
                os.path.join(qrefine,"plugin","yoink","dat"),
                "cluster.xml")
  interaction_list, weight = pyoink.get_interactions_list()
  #print interaction_list
  expected_list = [[24, 27], [19, 22], [5, 11], [18, 22], [25, 26], [23, 26], [2, 3], [5, 7], [9, 11], [27, 29],
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
  for e1, e2 in zip(expected_list, interaction_list):
    e1.sort()
    e2.sort()
  expected_list.sort()
  interaction_list.sort()
  assert approx_equal(expected_list, interaction_list)

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
