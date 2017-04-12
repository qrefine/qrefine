from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.show
import argparse

if (__name__ == "__main__"):
  import sys
  parser = argparse.ArgumentParser(description='qr.show')
  parser.add_argument('--refine',    action='store_true',
                                     default=False,
                                     help='show current refinement metrics ')
  parser.add_argument('--time',      action='store_true',
                                     default=False,
                                     help='show timing data')
  parser.add_argument('--resources', action='store_true',
                                     default=False,
                                     help='show available resources ')
  args = parser.parse_args()

  args = sys.argv[1:]
  job_id = args[1]

  #TODO
  if (args.refine)   : job_manager.refine(job_ids)
  if (args.time)     : job_manager.time(job_ids)
  if (args.resources): job_manager.resources(job_ids)
  if (not args.refine and not args.time and not args.resources ):
      job_manager.refine(job_ids)
      job_manager.time(job_ids)
      job_manager.resources(job_ids)