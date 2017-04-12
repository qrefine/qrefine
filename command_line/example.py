from __future__ import division
import argparse
import start
# LIBTBX_SET_DISPATCHER_NAME qr.example

def run_refine():
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
  start.run_cmd(cmd)

def run_cluster():
  cmd = """
          phenix.python qr.py 3dtj.mtz 3dtj.pdb
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
  start.run_cmd(cmd)


if __name__ == '__main__':
  print " Run Q|R Examples"
  parser = argparse.ArgumentParser(description='qrefine.example')
  parser.add_argument('--refine', action='store_true', default=False, help='run refinement example.     ')
  parser.add_argument('--cluster', action='store_true', default=False, help='run clustering example')

  args = parser.parse_args()
  if (args.refine): run_refine()
  if (args.cluster): run_cluster()
  if (not args.unit and not args.regression and not args.pdb):
    run_refine()
    run_cluster()
