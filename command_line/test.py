from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.test
import sys
import argparse
from qrefine.core.tests import run_tests as unit_tests
from qrefine.regression import regression_tests

log = sys.stdout

if __name__ == '__main__':
  print "Testing Q|R"
  parser = argparse.ArgumentParser(description='Test runners for Q|R code')
  parser.add_argument('--unit',
                      action='store_true',
                      default=False,
                      help='run unit regression ')
  parser.add_argument('--regression',
                      action='store_true',
                      default=False,
                      help='run regression regression, will take a long time')
  parser.add_argument('--pdb',
                      action='store_true',
                      default=False,
                      help='''run comprehensive regression on entire pdb,
                              only run on a HPC cluster''')
  parser.add_argument('--all',
                      action='store_true',
                      default=False,
                      help='''run all tests, only run on a HPC cluster''')
  args = parser.parse_args()
  if (args.unit)      :
    print >> log,"Running unit tests"
    unit_tests.run()
  if (args.regression):
    print >> log,"Running regression tests"
    regression_tests.run()
  if (args.pdb)       :
    print >> log,"Running pdb tests"
    regression_tests.run(args=sys.argv[1:])
  if(args.all):
    print >> log,"Running all regression tests"
    unit_tests.run()
    regression_tests.run()
    regression_tests.run(args=sys.argv[1:])
  if(not args.all and not args.unit and not args.regression and not args.pdb):
    print >> log,"Running unit tests"
    unit_tests.run()
