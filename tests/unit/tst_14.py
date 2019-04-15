from __future__ import division

import time, os
import iotbx.pdb
from qrefine import super_cell
from libtbx.test_utils import approx_equal
import run_tests

pdb_str = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      9  O   HOH Z 333       5.000   5.000   5.000  1.00137.30           O
TER
"""

pdb_str_super_sphere_answer = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
SCALE1      0.100000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.100000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.100000        0.00000
ATOM      1  O   HOH Z 333       5.000   5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHSS   0       5.000   5.000  15.000  1.00137.30           O
ATOM      1  O   HOHSS   1       5.000  -5.000   5.000  1.00137.30           O
ATOM      1  O   HOHSS   2       5.000   5.000  -5.000  1.00137.30           O
ATOM      1  O   HOHSS   3       5.000  15.000   5.000  1.00137.30           O
ATOM      1  O   HOHSS   4      -5.000   5.000   5.000  1.00137.30           O
ATOM      1  O   HOHSS   5      15.000   5.000   5.000  1.00137.30           O
TER
"""

def run(prefix):
  """
  Exercise supercell.
  """
  of = open("%s.pdb"%prefix,"w")
  print >> of, pdb_str
  of.close()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph = pdb_inp.construct_hierarchy()
  sites_cart_start = ph.atoms().extract_xyz()
  o = super_cell.expand(
    pdb_hierarchy        = ph,
    crystal_symmetry     = pdb_inp.crystal_symmetry(),
    select_within_radius = 11)
  #
  o.write_super_cell_selected_in_sphere(file_name="%s_super_sphere.pdb"%prefix)
  sites_cart = ph.atoms().extract_xyz()
  o.update(sites_cart = sites_cart)
  o.write_super_cell_selected_in_sphere(file_name="%s_super_sphere.pdb"%prefix)
  #
  sites_cart_super_sphere_answer = iotbx.pdb.input(source_info=None,
    lines=pdb_str_super_sphere_answer).atoms().extract_xyz()
  super_sphere_answer = list(sites_cart_super_sphere_answer.as_double())
  super_sphere = list(o.ph_super_sphere.atoms().extract_xyz().as_double())
  super_sphere_answer.sort()
  super_sphere.sort()
  assert approx_equal(super_sphere_answer,super_sphere)

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
