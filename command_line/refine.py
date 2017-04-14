from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.refine
import os
import sys
import time
import argparse
import libtbx.load_env
from libtbx.command_line import easy_qsub, easy_run
from qrefine.core import qr

phenix_source = os.path.dirname(libtbx.env.dist_path("phenix"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
qrefine_core_path = os.path.join(qrefine_path, "core")

def help():
  """ Commands in qrefine:
      - start (including restart, check{dry run, and 1SCF, qmready} )
      - pause
      - stop
      - ≈show
      - help
      - test
      - example
  """
  """ Options and keywords in qrefine:
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
  jobid = args[1]
  print "Are you sure you want to stop job {}".format(jobid)
  # how do we get a reference to all of the running jobs.
  # This could be hundreds of qdels, must automate this.


def pause():
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
          phenix.python qr.py 1uso.mtz 1uso.pdb
          max_iterations = 90
          max_atoms = 10000
          number_of_micro_cycles = 20
          qm_calculator = terachem
          mode = refine
          stpmax = 0.9
          restraints = qm
          gradient_only = true
          line_search = true
          shake_sites = false
          output_folder_name = refine
          restraints_weight_scale = 32
          > refine.log"""
  run_cmd(cmd)


def run(args, log):
  cmd = """
      phenix.python qr.py input.mtz input.pdb
      max_iterations = 90
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

if __name__ == '__main__':
  print "refine"
  parser = argparse.ArgumentParser(description='qrefine.example')
  parser.add_argument('--example', action='store_true', default=False, help='run refinement example.     ')
  t0 = time.time()
  log = sys.stdout
  print "Starting Q|R"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
