from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.restraint
import sys
import time
import os.path
import argparse
import iotbx.pdb
import libtbx.load_env
import mmtbx.command_line
import qrefine.qr as qr
from qrefine.restraints import from_qm, from_cctbx

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

log = sys.stdout

legend = """\
Compute energy and gradient for a system
"""

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=qr.master_params_str)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, get_master_phil().show()

def run(args, log):
  print "args are" , args
  print_legend_and_usage(log)
  cmdline = mmtbx.command_line.load_model_and_data(
        args=args,
        master_phil=get_master_phil(),
        create_fmodel=False,
        out=log)
  qr.run(cmdline=cmdline, log=log)

if (__name__ == "__main__"):
  print "Restraint for Q|R"
  t0 = time.time()
  print >> log, "Starting Q|R"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
