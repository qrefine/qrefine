from __future__ import division
from __future__ import absolute_import

import os
from qrefine.tests.unit import run_tests
import libtbx.load_env
import iotbx.pdb
import mmtbx.restraints
from libtbx.test_utils import approx_equal
from qrefine import restraints
from scitbx.array_family import flex

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix):
  """
  CCTBX gradients in presence of altlocs
  """
  method = "subtract"
  for file_name in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
                    "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:
    file_name = os.path.join(qr_unit_tests,"data_files",file_name)
    print(file_name)
    pi = iotbx.pdb.input(file_name = file_name)

    ph = pi.construct_hierarchy()
    cs = pi.crystal_symmetry()
    # Gradients from CCTBX (directly)
    g1 = restraints.get_cctbx_gradients(ph = ph, cs = cs).gradients
    # Gradients from CCTBX (via decomposing hierarchy into concormers A, B, blanc)
    g2 = restraints.from_cctbx_altlocs(ph = ph, cs = cs, option=1, method=method)
    g3 = restraints.from_cctbx_altlocs(ph = ph, cs = cs, option=2, method=method)

    _, gX = restraints.from_altlocs2(ph = ph, cs = cs, method=method
      ).target_and_gradients(sites_cart = ph.atoms().extract_xyz())
    assert approx_equal(g1, gX)

    # Make sure they match!
    assert approx_equal(g1, g2)
    if method == "subtract": assert approx_equal(g2, g3)

    # Individual differences
    #g1 = g1.as_double()
    #g2 = g2.as_double()
    #diff = flex.abs(g1-g2)
    #for d, atom in zip(diff, ph.atoms()):
    #  print(atom.format_atom_record(), d)

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
