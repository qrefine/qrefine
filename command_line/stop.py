from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.stop

import os
import sys
import time
import libtbx.load_env
from qrefine.core import qr

top_dir = os.path.dirname(libtbx.env.dist_path("qrefine"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
qrefine_core_path = os.path.join(qrefine_path, "core")
print os.listdir(qrefine_core_path)

def stop(args, log):
  jobid = args[1]
  print "Are you sure you want to stop job {}".format(jobid)
  # how do we get a reference to all of the running jobs.
  # This could be hundreds of qdels, must automate this.

if __name__ == '__main__':
  t0 = time.time()
  log = sys.stdout
  stop(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
