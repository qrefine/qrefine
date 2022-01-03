from __future__ import division
from __future__ import absolute_import

import os
import time

from qrefine.tests.unit import run_tests
import iotbx.pdb
import libtbx.load_env
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix):
  """
  Exercise refine_sites=False.
  """
  xrs_good,xrs_poor,f_obs,r_free_flags = run_tests.setup_helix_example()
  run_tests.run_cmd(prefix,args = ["restraints=cctbx","refine_sites=False"])
  xrs_start = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests,"data_files","m00_poor.pdb")).xray_structure_simple()
  xrs_refined = iotbx.pdb.input(
    file_name = os.path.join(prefix,"m00_poor_refined.pdb")).xray_structure_simple()
  d = flex.sqrt((xrs_start.sites_cart() - xrs_refined.sites_cart()).dot())
  assert approx_equal(d.min_max_mean().as_tuple(), [0,0,0])

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
  