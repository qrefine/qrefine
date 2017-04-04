from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.run

from qr.core import qr
import sys
import time
print " Quantum Refinement Program"

if __name__ == '__main__':
    t0 = time.time()
    log = sys.stdout
    qr.run(args=sys.argv[1:], log=log)
    print >> log, "Time: %6.4f" % (time.time() - t0)
