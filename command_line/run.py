from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.run

from qrefine.core import qr
import libtbx.load_env

import sys
import time


print " Quantum Refinement Program"

top_dir = os.path.dirname(libtbx.env.dist_path("qrefine"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
qrefine_core_path = os.path.join(qrefine_path, "core")





if __name__ == '__main__':
    t0 = time.time()
    log = sys.stdout
    print os.listdir(qrefine_core_path)
    qr.run(args=sys.argv[1:], log=log)
    print >> log, "Time: %6.4f" % (time.time() - t0)
