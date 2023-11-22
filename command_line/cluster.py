from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.cluster
import sys
import time
import os.path
import libtbx
import iotbx.pdb
from qrefine.fragment import fragments

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_yoink_path =os.path.join(qrefine_path, "plugin","yoink","yoink")

log = sys.stdout

legend = """\
Cluster a system into many small pieces
"""

def run(pdb_file, log,  maxnum_residues_in_cluster=4):
  print("max number of residues in each cluster:\n", maxnum_residues_in_cluster, file=log)
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  fq = fragments(
    pdb_hierarchy=ph,
    crystal_symmetry=cs,
    maxnum_residues_in_cluster=int(maxnum_residues_in_cluster),
    qm_run=False)# not run qm_calculation, just the clustering result
  print("Residue indices for each cluster:\n", fq.clusters, file=log)
  print('# clusters  : ',len(fq.clusters), file=log)


if (__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args[0], log)
  print("Time: %6.4f" % (time.time() - t0), file=log)
