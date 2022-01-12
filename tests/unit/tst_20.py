from __future__ import division
from __future__ import absolute_import
import iotbx.pdb
import os
import time
import libtbx.load_env
from qrefine.super_cell import expand
from qrefine.tests.unit import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def run(prefix):
  """
  Make sure expand works for this particular file.
  """
  pdb_name = os.path.join(qr_unit_tests_data,
    "1bdw_ala_refine_001_complete_minimized.pdb_modified.pdb")
  pdb_inp = iotbx.pdb.input(pdb_name)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  super_cell = expand(
    pdb_hierarchy=ph,
    crystal_symmetry=cs)

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
