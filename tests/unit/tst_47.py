from __future__ import absolute_import, division, print_function
import os
import time
import libtbx.load_env
from libtbx import easy_run
import numpy as np

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine, "tests", "unit", "data_files")

def run_command(cmd, log_name):
  cmd = "%s > %s 2>&1" % (cmd, log_name)
  print(cmd)
  rc = easy_run.call(cmd)
  text = open(log_name).read()
  assert rc == 0, text
  assert text.find("Sorry") == -1 and text.find("Traceback") == -1, text
  return text

def exercise_cluster(pdb):
  """
  qr.cluster: residues of the model partitioned into clusters.
  """
  text = run_command("qr.cluster %s maxnum_residues_in_cluster=5" % pdb,
    "cluster.log")
  n_residues = 13 # helix.pdb
  n_clusters = int([l for l in text.splitlines()
                    if l.startswith("# clusters:")][0].split()[-1])
  sizes = [l for l in text.splitlines()
           if l.startswith("# residues in clusters:")][0]
  sizes = [int(x) for x in sizes.split("[")[1].split("]")[0].split(",")]
  assert len(sizes) == n_clusters and n_clusters > 1, (sizes, n_clusters)
  assert sum(sizes) == n_residues and max(sizes) <= 5, sizes

def exercise_charges(pdb):
  """
  qr.charges: total formal charge of a complete model (0 for this peptide).
  """
  text = run_command("qr.charges %s verbose=True" % pdb, "charges.log")
  charge = [l for l in text.splitlines() if l.startswith("Charge:")][0]
  assert charge.split()[-1] == "0", charge

def exercise_restraint(pdb):
  """
  qr.restraint: one gradient evaluation and minimizer step (cctbx here).
  """
  text = run_command("qr.restraint %s restraints=cctbx" % pdb, "restraint.log")
  assert text.find("Convergence reached") > -1 or \
         text.find("cycle   0:") > -1, text
  assert os.path.isfile(os.path.join("pdb", "helix_refined.pdb"))

def exercise_gtest(pdb):
  """
  qr.gtest: gradients of the clustered restraints for two cluster sizes,
  saved as <g_mode>-<size>.npy, same size as the model.
  """
  import iotbx.pdb
  run_command("qr.gtest %s restraints=cctbx g_scan='3 6'" % pdb, "gtest.log")
  n_atoms = iotbx.pdb.input(file_name=pdb).atoms().size()
  for fn in ["1-3.npy", "1-6.npy"]:
    assert os.path.isfile(fn), fn
    g = np.load(fn)
    assert g.shape == (3*n_atoms,), (g.shape, n_atoms)
    assert np.linalg.norm(g) > 0

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  The auxiliary command lines run: qr.cluster, qr.charges, qr.restraint and
  qr.gtest (cctbx restraints) on helix.pdb.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb = os.path.join(qr_unit_tests_data, "helix.pdb")
  for exercise in [exercise_cluster, exercise_charges, exercise_restraint,
                   exercise_gtest]:
    folder = exercise.__name__
    os.makedirs(folder, exist_ok=True)
    os.chdir(folder)
    try:
      exercise(pdb)
    finally:
      os.chdir("..")
    print(folder, "OK")

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
