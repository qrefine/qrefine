from __future__ import division

import time
import os
import os.path
import libtbx.load_env
from libtbx.test_utils import approx_equal
from qrefine.plugin.yoink.pyoink import PYoink
import iotbx.pdb
from qrefine.utils import yoink_utils
import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix = "tst_09"):
  """
  Exercise buffer region of cluster.
  """
  run_tests.assert_folder_is_empty(prefix=prefix)
  pdb_inp = iotbx.pdb.input(file_name= os.path.join(qr_unit_tests,"data_files","2lvr.pdb"))
  ph = pdb_inp.construct_hierarchy()
  yoink_utils.write_yoink_infiles("cluster.xml",
                                  "qmmm.xml",
                                  ph,
                                  os.path.join(qrefine,"plugin","yoink","dat"))
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
  bc_qms = [[1, 2, 3, 4], [1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            [3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31],
            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]]
  gn_qms = [[1, 2, 3, 4], [1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 18, 19, 22],
         [5, 6, 7, 8, 9, 10, 11, 13, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31],
         [3, 5, 6, 9, 10, 11, 12, 13, 22],
         [4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
         [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]]
  pyoink = PYoink(os.path.join(qrefine,"plugin","yoink","Yoink-0.0.1.jar"),
                               os.path.join(qrefine,"plugin","yoink","dat"),
                               "qmmm.xml")
  check_buffer(bc_clusters, bc_qms, pyoink)

def check_buffer(clusters, qms, pyoink):
  qms_calculated = []
  for i in range(len(clusters)):
    pyoink.update(clusters[i])
    qm_atoms, qm_molecules = pyoink.get_qm_indices()
    qms_calculated.append(list(qm_molecules))
  assert approx_equal(qms, qms_calculated)

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_09"
  run(prefix)
  print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
  run_tests.clean_up(prefix)
