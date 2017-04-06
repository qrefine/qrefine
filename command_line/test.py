from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.test

import argparse
from qrefine.core.tests import run_tests

def run_regression():
  import sys
  success = test_runner.run_multiple_tests(args=sys.argv[1:])
  return success

if __name__ == '__main__':
  print " Test Quantum Refinement Program"
  parser = argparse.ArgumentParser(description='Shaz.IR')
  parser.add_argument('--unit', action='store_true', default=True, help='run unit tests    ')
  parser.add_argument('--regression', action='store_true', default=False, help='run regression tests, will take a long time')
  parser.add_argument('--pdb', action='store_true', default=False, help='test optimizer code only ')

  args = parser.parse_args()
  if (args.unit): run_tests.unit()
  if (args.regression): run_tests.regression()
  if (args.pdb): run_tests.minmizer()
  if (not args.unit and not args.regression and not args.pdb):
      run_tests.run()