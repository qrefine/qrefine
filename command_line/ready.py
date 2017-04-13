from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.ready

import sys
import os
import time
import os.path
import iotbx.pdb

try:
  from jpype import startJVM
except ImportError, e:
  raise Sorry(str(e))
from restraints import from_qm
import completion

qr_path = os.environ["QR_REPO_PARENT"]

log = sys.stdout

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
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
  print >> log, "Time: %6.4f" % (time.time() - t0)

