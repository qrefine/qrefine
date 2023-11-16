from __future__ import division
from __future__ import absolute_import

import os
import time
import iotbx.pdb
import libtbx.load_env
from libtbx.test_utils import approx_equal
from mmtbx.pair_interaction import pair_interaction
from qrefine.tests.unit import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def get_hierarchy():
  pdb_inp = iotbx.pdb.input(file_name= os.path.join(
    qr_unit_tests,"data_files","2lvr.pdb"))
  return pdb_inp.construct_hierarchy()

def run(prefix):
  """
  Exercise interaction graph construction.
  """
  ph = get_hierarchy()
  interaction_list_cpp = pair_interaction.run(ph)
  expected_list_cpp = [
    (9, 26), (24, 30), (6, 9), (17, 20), (1, 3), (18, 19), (23, 26), (6, 7),
    (4, 19), (24, 27), (6, 10), (11, 22), (25, 26), (7, 19), (27, 28), (15, 18),
    (5, 11), (29, 30), (4, 16), (6, 23), (16, 19), (6, 26), (17, 18), (22, 25),
    (3, 12), (4, 15), (16, 18), (6, 13), (19, 23), (15, 16), (21, 24), (22, 23),
    (22, 26), (8, 9), (17, 21), (20, 21), (24, 28), (6, 11), (13, 19), (18, 21),
    (23, 24), (3, 5), (5, 10), (7, 8), (5, 7), (24, 25), (16, 20), (13, 22),
    (18, 22), (28, 29), (5, 13), (19, 22), (15, 19), (16, 17), (11, 12), (4, 13),
    (9, 11), (11, 13), (20, 22), (13, 15), (2, 3), (8, 26), (6, 8), (20, 24),
    (6, 31), (26, 31), (21, 22), (13, 18), (27, 30), (23, 27), (3, 4), (2, 4),
    (10, 11), (8, 23), (5, 6), (9, 22), (5, 19), (22, 31), (3, 14), (27, 29),
    (1, 2), (28, 30), (5, 12), (10, 12), (4, 5), (7, 10), (6, 22), (7, 23),
    (17, 19), (22, 24), (3, 13), (4, 14), (9, 10), (9, 31), (19, 20), (25, 28),
    (15, 17), (6, 19), (21, 25), (20, 23), (26, 27), (13, 14), (12, 13), (14, 15)]

  result = list(set(expected_list_cpp).symmetric_difference(interaction_list_cpp))
  assert len(result) < 3, result # Allow some differences due to numerics

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
