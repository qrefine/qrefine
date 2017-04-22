from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.refine
import os
import sys
import time
import argparse
import libtbx.load_env
from libtbx import  easy_run
from libtbx.command_line import easy_qsub

phenix_source = os.path.dirname(libtbx.env.dist_path("phenix"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
qrefine_core_path = os.path.join(qrefine_path, "core")
example_path = os.path.join(qrefine_path,"examples")

log = sys.stdout

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
  cmd = " phenix.python " +  \
        os.path.join(qrefine_core_path,"qr.py ") + \
        " ".join(args).replace("\n", "")
  print "Running example:", cmd
  easy_run.call(cmd)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Refine a model using restraints from Quantum Chemistry')
  parser.add_argument('--example', action='store_true', default=False, help='run refinement example.')
  known, unknown = parser.parse_known_args()
  t0 = time.time()
  print >> log,"Starting Q|R"
  if(known.example):
    example()
  else:
    run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
