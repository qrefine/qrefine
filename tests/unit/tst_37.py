from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import libtbx.load_env
from scitbx.array_family import flex
from libtbx import easy_pickle
from libtbx.test_utils import approx_equal
from qrefine.tests.unit import run_qrefine

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Exercise gradients match: using clustering vs not using clustering.
  Non-P1 case (P212121), cctbx restraints.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  common = [
    os.path.join(qr_unit_tests_data,"2olx.pdb"),
    "mode=opt",
    "restraints=cctbx",
    "debug=True",
    "number_of_micro_cycles=1",
    "max_iterations_refine=1",
    "parallel.nproc=1"]
  for clustering in ["False", "True"]:
    r_qr = run_qrefine(prefix = prefix,
      args = common+["clustering=%s"%clustering])
  #
  # Gradients at the starting coordinates match. Later gradients are not
  # comparable: the two optimizations diverge after the first step.
  #
  g1 = easy_pickle.load("debug_optimization_clusteringFalse.pkl").gradients[0]
  g2 = easy_pickle.load("debug_optimization_clusteringTrue.pkl").gradients[0]
  assert g1.size() == g2.size()
  diff = flex.sqrt((g1-g2).dot())
  assert approx_equal(flex.max(diff), 0, 1.e-6)

if(__name__ == '__main__'):
  run()
