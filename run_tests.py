from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list = [
  "$D/tests/unit/tst_00.py",
  "$D/tests/unit/tst_01.py",
  "$D/tests/unit/tst_02.py",
  "$D/tests/unit/tst_03.py",
  "$D/tests/unit/tst_04.py",
  "$D/tests/unit/tst_05.py",
  "$D/tests/unit/tst_06.py",
  "$D/tests/unit/tst_07.py",
  "$D/tests/unit/tst_08.py",
  "$D/tests/unit/tst_09.py",
  "$D/tests/unit/tst_10.py",
  "$D/tests/unit/tst_11.py",
  "$D/tests/unit/tst_12.py",
  "$D/tests/unit/tst_13.py",
  "$D/tests/unit/tst_14.py",
  "$D/tests/unit/tst_15.py",
  "$D/tests/unit/tst_16.py",
  "$D/tests/unit/tst_17.py",
  "$D/tests/unit/tst_18.py",
  "$D/tests/unit/tst_19.py",
  "$D/tests/unit/tst_20.py",
  "$D/tests/unit/tst_21.py",
  "$D/tests/unit/tst_22.py",
  "$D/tests/unit/tst_23.py",
  "$D/tests/unit/tst_24.py",
  "$D/tests/unit/tst_25.py",
  "$D/tests/unit/tst_26.py",
  ]

def run():
  build_dir = libtbx.env.under_build("qrefine")
  dist_dir = libtbx.env.dist_path("qrefine")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
