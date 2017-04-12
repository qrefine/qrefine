from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.pause

import os
import sys
import time
from libtbx.command_line import easy_qsub, easy_run

# We need some sort of a job manager
from qrefine.core import qr


if __name__ == '__main__':
  id = 0
  print "Are you sure you want to pause job {}".format(id)
  t0 = time.time()
  log = sys.stdout
  #qr.run(args=sys.argv[1:], log=log)
  cmds= []
  # example job ids, just for getting the code to work
  sub_job_ids =['305.mu01','306.mu01']
  for sub_job_id in sub_job_ids:
    cmds.append(" qhold" +  sub_job_id)
  easy_run.call(cmds)
  print >> log, "Time: %6.4f" % (time.time() - t0)