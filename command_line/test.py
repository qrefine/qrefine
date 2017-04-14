from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.test
import sys
import argparse
from qrefine.tests      import regression_tests
from qrefine.core.tests import run_tests as unit_tests

if __name__ == '__main__':
  print " Testing Quantum Refinement Program"
  parser = argparse.ArgumentParser(description='Test runners for Q|R code')
  parser.add_argument('--unit',       action='store_true',
                                      default=True,
                                      help='run unit tests ')
  parser.add_argument('--regression', action='store_true',
                                      default=False,
                                      help='run regression tests, will take a long time') # nightly build?
  parser.add_argument('--pdb',        action='store_true',
                                      default=False,
                                      help='''run comprehensive tests on entire pdb,
                                              only run on a HPC cluster''')

  args = parser.parse_args()
  if (args.unit)      : unit_tests.run()
  if (args.regression): regression_tests.run(args=sys.argv[1:])
  if (args.pdb)       : regression_tests.run(args=sys.argv[1:])
