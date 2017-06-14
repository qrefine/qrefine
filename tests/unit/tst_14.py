from __future__ import division

import time
import iotbx.pdb
from qrefine import super_cell

#pdb_str = """
#CRYST1   10.000   10.000   10.000  70.00  80.00 120.00 P 1
#ATOM      1  N   SER A 149       3.000   3.000   3.000  1.00137.30           N
#ATOM      9  O   HOH A 333       4.000   4.000   4.000  1.00137.30           N
#TER
#ATOM      1  N   SER Z 149       5.000   5.000   5.000  1.00137.30           N
#TER
#"""

pdb_str = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      9  O   HOH Z 333       5.000   5.000   5.000  1.00137.30           O
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

if(__name__ == "__main__"):
  t0 = time.time()
  prefix = "tst_14"
  run(prefix)
  print prefix + ":  OK  " + "Time: %6.2f (s)" % (time.time() - t0)
