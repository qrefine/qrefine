from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.refine
import os
import sys
import time
import argparse
import libtbx.load_env
from libtbx import  easy_run
from libtbx.command_line import easy_qsub
import mmtbx.command_line
from qrefine import qr

phenix_source = os.path.dirname(libtbx.env.dist_path("phenix"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
example_path = os.path.join(qrefine_path,"examples")

log = sys.stdout

legend = """
Refine a model using restraints from Quantum Chemistry
"""

master_params_str ="""
import scope qrefine.qr.master_params_str
"""

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=master_params_str)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, "                               phenix.polder"
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, get_master_phil().show()

def run_cmd(cmd):
  # we want a reference to the running job.
  if cmd.find("qsub_command") is not -1:
    easy_qsub.run(
      phenix_source=phenix_source,
      where="working_dir",
      commands=cmd
    )
  else:
    # we want a reference to the running job.
    result = easy_run.fully_buffered(cmd).raise_if_errors()

def example():
  cmd = "phenix.python " +  \
      os.path.join(qrefine_core_path,"qr.py ") +  \
      os.path.join(qrefine_path,"examples/1us0/data.mtz ") + \
      os.path.join(qrefine_path,"examples/1us0/a87_99_h.pdb ") + \
      "output_folder_name = 1us0 > 1us0.log"
  cmd = cmd.replace("\n", "")
  print "Running example:", cmd
  run_cmd(cmd)

def run(args, log):
  print "running refine "
  print_legend_and_usage(log)
  cmdline = mmtbx.command_line.load_model_and_data(
       args          = args,
       master_phil   = get_master_phil(),
       create_fmodel = False,
       out           = log)
  params = cmdline.params
  qr.run(params.refine,log)

if __name__ == '__main__':
  t0 = time.time()
  print >> log,"Starting Q|R"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
