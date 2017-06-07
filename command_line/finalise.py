from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.finalise
import os, sys, time
import argparse
import libtbx.load_env
import qrefine.finalise

log = sys.stdout

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

legend = """\
Finalise a model before quantum refinement
"""


def run(args):
  cmd = "phenix.python " +  \
        os.path.join(qr_path,"finalise.py ") + \
        " ".join(args).replace("\n", "")
  cmd = cmd.replace("\n", "")
  print "Running example:", cmd
  easy_run.fully_buffered(command=cmd)

if __name__ == '__main__':

  t0 = time.time()
  print >> log,"Starting Q|R"
  run(params, log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
