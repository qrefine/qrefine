from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.finalise
import os, sys, time
import argparse
import libtbx.load_env
from libtbx import easy_run

log = sys.stdout

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

def example():
  cmd = "phenix.python " +  \
        os.path.join(qr_path,"finalise.py ") +  \
        os.path.join(qrefine_path,"examples/1us0/a87_99_h.pdb ") + \
        "model_completion=False > 1us0.log"
  cmd = cmd.replace("\n", "")
  print "Running example:", cmd
  easy_run.fully_buffered(command=cmd)

def run(args):
  cmd = "phenix.python " +  \
        os.path.join(qr_path,"finalise.py ") + \
        " ".join(args).replace("\n", "")
  cmd = cmd.replace("\n", "")
  print "Running example:", cmd
  easy_run.fully_buffered(command=cmd)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='Finalise a model before quantum refinement '
  )
  parser.add_argument('--example',
                      action='store_true',
                      default=False,
                      help='run finalise example.')
  known, unknown = parser.parse_known_args()
  t0 = time.time()
  print >> log,"Starting Q|R"
  if(known.example):
    example()
  else:
    run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
