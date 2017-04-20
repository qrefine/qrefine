from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.refine
import os
import sys
import time
import argparse
import libtbx.load_env
import libtbx
from libtbx import easy_run
from libtbx import command_line
from libtbx import  easy_run
from libtbx.command_line import easy_qsub
from qrefine.core import qr

phenix_source = os.path.dirname(libtbx.env.dist_path("phenix"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
qrefine_core_path = os.path.join(qrefine_path, "core")
example_path = os.path.join(qrefine_path,"examples")

log = sys.stdout

def help():
  """  
  Options and keywords in qrefine:
       - qm_calculator
       - macro_cycles
       - micro_cycles
       - max_bond_rmsd
       - refine_sites
       - refine_adp
       - cluster_qm
       - charge_embedding
       - cluster
   """

def stop(args, log):
  """how do we get a reference to all of the running jobs.
  This could be hundreds of qdels,we must automate this."""
  jobid = args[1]
  print "Are you sure you want to stop job {}".format(jobid)

def pause():
  id = 0
  print "Are you sure you want to pause job {}".format(id)
  # example job ids, just for getting the code to work
  sub_job_ids =['305.mu01','306.mu01']
  cmds= []
  for sub_job_id in sub_job_ids:
    cmds.append(" qhold" +  sub_job_id)
  easy_run.call(cmds)

def run_cmd(cmd):
  # we want a reference to the running job.
  if cmd.find("qsub_command") is not -1:
    easy_qsub.run(
      phenix_source=phenix_source,
      where="working_dir",
      commands=cmd)
  else:
    # we want a reference to the running job.
    result = easy_run.fully_buffered(cmd).raise_if_errors()

def example():
  cmd = """
          phenix.python """+  \
      os.path.join(qrefine_core_path,"qr.py ") +  \
      os.path.join(qrefine_path,"examples/1us0/data.mtz ") + \
      os.path.join(qrefine_path,"examples/1us0/a87_99_h.pdb ") + \
       """output_folder_name = refine
          > refine.log"""
  run_cmd(cmd)

def run(args, log):
  cmd = " phenix.python " +  \
      os.path.join(qrefine_core_path,"qr.py ") +  \
      os.path.join(qrefine_path,"examples/1us0/data.mtz ") + \
      os.path.join(qrefine_path,"examples/1us0/a87_99_h.pdb ") + \
      """max_iterations = 90
      maxnum_residues_in_cluster = 5
      max_atoms = 10000
      number_of_micro_cycles = 20
      qm_calculator = terachem
      mode = refine
      stpmax = 0.9
      restraints = qm
      gradient_only = true
      line_search = true
      shake_sites = false
      use_cluster_qm = true
      qsub_command = "qsub"
      shared_disk = false
      output_folder_name = cluster_glr
      restraints_weight_scale = 32
      > cluster_glr.log"""
  print cmd
  easy_run.call(cmd)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='qrefine.example')
  parser.add_argument('--example', action='store_true', default=False, help='run refinement example.     ')
  t0 = time.time()
  log = sys.stdout
  print "Starting Q|R"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
