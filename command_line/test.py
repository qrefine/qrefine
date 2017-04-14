from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.test
import sys
import argparse
from qrefine.tests      import run_tests
from qrefine.core.tests import run_tests

if __name__ == '__main__':
  print " Test Quantum Refinement Program"
  parser = argparse.ArgumentParser(description='Test runners for Q|R code')
  parser.add_argument('--unit',       action='store_true',
                                      default=True,
                                      help='run unit tests ')
  parser.add_argument('--regression', action='store_true',
                                      default=False,
                                      help='run regression tests, will take a long time')
  parser.add_argument('--pdb',        action='store_true',
                                      default=False,
                                      help='run comprehensive tests on entire pdb')

  args = parser.parse_args()
  if (args.unit)      : run_tests.run()
  if (args.regression): run_tests.run_regression_tests(args=sys.argv[1:])
  if (args.pdb)       : run_tests.run_pdb_tests(args=sys.argv[1:])
  if (not args.unit and not args.regression and not args.pdb):
    success = run_tests.run_all_tests(args=sys.argv[1:])
