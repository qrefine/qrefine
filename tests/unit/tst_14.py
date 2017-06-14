from __future__ import division

import time
import iotbx.pdb
from qrefine import super_cell
from libtbx.test_utils import approx_equal

pdb_str = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      9  O   HOH Z 333       5.000   5.000   5.000  1.00137.30           O
TER
"""

pdb_str_super_cell_answer = """
ATOM      1  O   HOH Z 333       5.000   5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHG7 333      -5.000  -5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHG6 333      -5.000  -5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHG5 333      -5.000  -5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHG4 333      -5.000   5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHG3 333      -5.000   5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHG2 333      -5.000   5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHG1 333      -5.000  15.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHG0 333      -5.000  15.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHG9 333      -5.000  15.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHG8 333       5.000  -5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGW 333       5.000  -5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGV 333       5.000  -5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHGU 333       5.000   5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGT 333       5.000   5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHGS 333       5.000  15.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGR 333       5.000  15.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGQ 333       5.000  15.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHGP 333      15.000  -5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGZ 333      15.000  -5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGY 333      15.000  -5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHGX 333      15.000   5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGF 333      15.000   5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGE 333      15.000   5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHGD 333      15.000  15.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGC 333      15.000  15.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGB 333      15.000  15.000  15.000  1.00137.30           O
TER
"""

pdb_str_super_sphere_answer = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
SCALE1      0.100000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.100000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.100000        0.00000
ATOM      1  O   HOH Z 333       5.000   5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHG3 333      -5.000   5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGW 333       5.000  -5.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGU 333       5.000   5.000  -5.000  1.00137.30           O
TER
ATOM      1  O   HOHGT 333       5.000   5.000  15.000  1.00137.30           O
TER
ATOM      1  O   HOHGR 333       5.000  15.000   5.000  1.00137.30           O
TER
ATOM      1  O   HOHGF 333      15.000   5.000   5.000  1.00137.30           O
TER
"""

def run(prefix = "tst_14"):
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
  o.write_super_cell(file_name="%s_super_cell.pdb"%prefix)
  o.write_p1(file_name="%s_p1.pdb"%prefix)
  o.write_super_cell_selected_in_sphere(file_name="%s_super_sphere.pdb"%prefix)
  sites_cart = ph.atoms().extract_xyz()
  o.update(sites_cart = sites_cart)
  o.write_super_cell(file_name="%s_super_cell.pdb"%prefix)
  #
  sites_cart_super_sphere_answer = iotbx.pdb.input(source_info=None,
    lines=pdb_str_super_sphere_answer).atoms().extract_xyz()
  assert approx_equal(sites_cart_super_sphere_answer,
    o.ph_super_sphere.atoms().extract_xyz())
  #
  sites_cart_super_cell_answer = iotbx.pdb.input(source_info=None,
    lines=pdb_str_super_cell_answer).atoms().extract_xyz()
  assert approx_equal(sites_cart_super_cell_answer,
    o.ph_super_cell.atoms().extract_xyz())

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_14"
  run(prefix)
  print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
