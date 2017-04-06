from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.show
import argparse

def run(args):
   print "running "
   job.time()
   os.resources()
   pbs.cluster()
   qm_ready.valid()

if (__name__ == "__main__"):
  import sys
  parser = argparse.ArgumentParser(description='qr.stats')
  parser.add_argument('--all', action='store_true', default=True, help=' show all information   ')
  parser.add_argument('--time', action='store_true', default=False, help='show timing data')
  parser.add_argument('--resources', action='store_true', default=False, help='show available resources ')
  parser.add_argument('--q', action='store_true', default=False, help='show pbs queue  ')
  parser.add_argument('--qmready', action='store_true', default=False, help='test structure preparation code     ')
  args = parser.parse_args()

  if (args.time): job.time()
  if (args.resources): os.resources()
  if (args.q): pbs.cluster()
  if (args.qmready): qm_ready.valid()
  if (not args.all and not args.regression and not args.minimizer and not args.cluster):
    run(args=sys.argv[1:])
