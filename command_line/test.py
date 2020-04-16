from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.test
# LIBTBX_SET_DISPATCHER_NAME qr.run_tests
import sys
import argparse
import qrefine.tests.unit.run_tests as unit_tests
from qrefine.tests.regression import regression_tests

log = sys.stdout

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Test runners for Q|R code')
  parser.add_argument('--unit',
                      action='store_true',
                      default=False,
                      help='run unit regression ')
  parser.add_argument('--reg',
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
  parser.add_argument('--nproc',
                      action='store',
                      default=1,
                      help='nprocs',
                      type=int)
  parser.add_argument('--only',
                      action='store',
                      default=None,
                      help='only this test',
                      type=int)
  parser.add_argument('--non_mopac_only',
                      action='store_true',
                      default=False,
                      help="only tests that don't use MOPAC")
  args = parser.parse_args()
  if (args.unit)      :
    print >> log,"Running Q|R unit tests"
    unit_tests.run(nproc=args.nproc)
  if (args.reg):
    print >> log,"Running Q|R regression tests"
    regression_tests.run(nproc=args.nproc)
  if (args.pdb)       :
    print >> log,"Running Q|R pdb tests"
    # swith to acceptance_tests.py
    regression_tests.run(args=sys.argv[1:])
  if(args.all):
    print >> log,"Running all regression tests"
    unit_tests.run()
    regression_tests.run()
    regression_tests.run(args=sys.argv[1:])
  if(not args.all and not args.unit and not args.reg and not args.pdb):
    print >> log,"Running Q|R unit tests"
    rc = unit_tests.run(nproc=args.nproc,
                        only_i=args.only,
                        non_mopac_only=args.non_mopac_only)
    assert not rc, 'qr.test rc : %s' % rc
