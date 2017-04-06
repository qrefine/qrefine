from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.stop

import os
import sys
import time

def setup():
  from qrefine.core import qr
  import libtbx.load_env
  top_dir = os.path.dirname(libtbx.env.dist_path("qrefine"))
  qrefine_path = libtbx.env.find_in_repositories("qrefine")
  qrefine_core_path = os.path.join(qrefine_path, "core")
  print os.listdir(qrefine_core_path)

if __name__ == '__main__':
  id=0
  print "Are you sure you want to stop job {}".format(id)
  t0 = time.time()
  log = sys.stdout
  qr.run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
