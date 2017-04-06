from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.start

import os
import sys
import time

def start_qr():
  print " Quantum Refinement Program"
  import libtbx.load_env
  top_dir = os.path.dirname(libtbx.env.dist_path("qrefine"))
  qrefine_path = libtbx.env.find_in_repositories("qrefine")
  qrefine_core_path = os.path.join(qrefine_path, "core")
  print os.listdir(qrefine_core_path)

if __name__ == '__main__':
  from qrefine.core import qr
  t0 = time.time()
  log = sys.stdout
  #we want a reference to the running job.
  qr.run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
