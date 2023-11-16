from __future__ import division
from __future__ import absolute_import

import os
import time
from qrefine.tests.unit import run_tests
import libtbx.load_env
from mmtbx import monomer_library
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation
import iotbx.pdb
from libtbx.utils import null_out
import mmtbx.restraints
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from qrefine import restraints

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix):
  """
  CCTBX gradients in presence of altlocs
  
  XXX TEST FAILS. CRITICAL TO FIX. The test exercises the fundamental logic of
  XXX handling gradients in presence of alternative conformations.
  
  """
  for file_name in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
                    "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:
    file_name = os.path.join(qr_unit_tests,"data_files",file_name)

    pi = iotbx.pdb.input(file_name = file_name)
    ph = pi.construct_hierarchy()
    cs = pi.crystal_symmetry()
    # Gradients from CCTBX (directly)
    g1 = restraints.get_cctbx_gradients(ph = ph, cs = cs).gradients
    # Gradients from CCTBX (via decomposing hierarchy into concormers A, B, blanc)
    g2 = restraints.from_cctbx_altlocs(ph = ph, cs = cs)
    # Make sure they match!
    assert approx_equal(g1, g2)

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
