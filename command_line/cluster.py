from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.cluster
import sys
import time
import os.path
import libtbx
import iotbx.pdb
from qrefine.restraints import from_qm

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_yoink_path =os.path.join(qrefine_path, "plugin","yoink","yoink")

log = sys.stdout

legend = """\
Cluster a system into many small pieces
"""

def run(pdb_file, maxnum_residues_in_cluster=15):
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  fq = from_qm(
    pdb_hierarchy=ph,
    crystal_symmetry=cs,
    use_cluster_qm=True,
    maxnum_residues_in_cluster=int(maxnum_residues_in_cluster))
  chunks = []
  chunk_sizes = []
  #This should be in chunk
  for chunk in fq.fragments.qm_pdb_hierarchies:
    res_in_chunk = []
    atom_tot_per_residue = 0
    for chain in chunk.only_model().chains():
      for residue_group in chain.residue_groups():
        res_in_chunk.append(residue_group.resid())
        for atom_group in residue_group.atom_groups():
          atom_tot_per_residue += atom_group.atoms_size()
          chunk_sizes.append(atom_tot_per_residue)
    chunks.append(res_in_chunk)
  print >> log, "Residue indices in each chunk:", fq.fragments.clusters
  print >> log, "pdb hierarchy for each fragment (cluster+buffer)", chunks

if (__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args[0], log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
