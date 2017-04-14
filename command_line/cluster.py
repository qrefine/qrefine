from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.cluster
import sys
import time
import os.path
import iotbx.pdb
import argparse
try:
  from jpype import startJVM
except ImportError, e:
  raise Sorry(str(e))
from restraints import from_qm
import completion

qr_path = os.environ["QR_REPO_PARENT"]

log = sys.stdout

def example():
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
  run_cmd(cmd)


def run(pdb_file, maxnum_residues_in_cluster=15):
  #  completion.run(pdb_file)
  #  pdb_file = os.path.basename(pdb_file)[:-4] + "_complete.pdb"
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  sites_cart = ph.atoms().extract_xyz()
  fq = from_qm(
    pdb_hierarchy=ph,
    crystal_symmetry=cs,
    use_cluster_qm=True,
    maxnum_residues_in_cluster=int(maxnum_residues_in_cluster),
    yoink_jar_path=qr_path + "/./qr-core/plugin/yoink/Yoink-0.0.1.jar",
    yoink_dat_path=qr_path + "/./qr-core/plugin/yoink/dat")
  print >> log, "molecular indices in clusters:(the molecular index starts from 1)", fq.fragments.clusters
  print >> log, "pdb hierarchy for each fragment (cluster+buffer)", fq.fragments.qm_pdb_hierarchies

if (__name__ == "__main__"):

  if (__name__ == "__main__"):
    parser = argparse.ArgumentParser(description='qr.help')
    parser.add_argument('--commands', action='store_true',
                        default=False,
                        help='display the set of available commands ')
    parser.add_argument('--keywords', action='store_true',
                        default=False,
                        help='display the set of available keywords ')
    parser.add_argument('--settings', action='store_true',
                        default=False,
                        help='display the set of default settings   ')
    args = parser.parse_args()

    if (args.commands):   print commands.__doc__
    if (args.keywords):   print keywords.__doc__
    if (args.settings):   print settings.__doc__
    if (not args.commands and not args.keywords and not args.settings):
      print run.__doc__
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
  print >> log, "Time: %6.4f" % (time.time() - t0)

