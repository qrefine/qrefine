from __future__ import absolute_import, division, print_function
import os
import time
import libtbx.load_env
import iotbx.pdb
from mmtbx.pair_interaction import pair_interaction
from qrefine.clustering import betweenness_centrality_clustering
from qrefine.fragment import bonded_residue_pairs

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine, "tests", "unit", "data_files")

def is_connected(residues, interaction_list):
  residues = set(residues)
  neighbors = dict((r, set()) for r in residues)
  for i, j in interaction_list:
    if(i in residues and j in residues):
      neighbors[i].add(j)
      neighbors[j].add(i)
  seen = set()
  stack = [min(residues)]
  while(stack):
    r = stack.pop()
    if(r in seen): continue
    seen.add(r)
    stack.extend(neighbors[r] - seen)
  return len(seen) == len(residues)

def check_clusters(clusters, interaction_list, n_residues, maxnum):
  """
  The clusters partition the residues 1..n_residues, none has more than
  maxnum residues and each is connected in the interaction graph.
  """
  flat = sorted(r for cluster in clusters for r in cluster)
  assert flat == list(range(1, n_residues+1)), flat
  sizes = [len(c) for c in clusters]
  assert max(sizes) <= maxnum, (sizes, maxnum)
  for cluster in clusters:
    assert is_connected(cluster, interaction_list), cluster

def exercise_clique():
  """
  A clique (all residues interact with all): betweenness centrality cannot
  split it, the clusters come from the fallback split. Three times the
  maximum size: one split in half is not enough.
  """
  n, maxnum = 13, 4
  interaction_list = [(i, j) for i in range(1, n+1) for j in range(i+1, n+1)]
  clusters = betweenness_centrality_clustering(interaction_list,
    maxnum_residues_in_cluster=maxnum, bcc_threshold=9, size=n).get_clusters()
  check_clusters(clusters, interaction_list, n, maxnum)

def exercise_1yjp():
  """
  1yjp (7 residues and 7 waters): at maxnum 4 betweenness centrality leaves
  a component of 6 residues (3-6 and the waters 13, 14) that the fallback
  split has to cut, and a cut by residue index (3-5 | 6, 13, 14) leaves
  residue 6 without an interaction partner in its cluster. With the bonded
  pairs the cut is the one that cuts no peptide bond: 3-6 | 13, 14.
  """
  maxnum = 4
  ph = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,
    "1yjp.pdb")).construct_hierarchy()
  n = len(list(ph.residue_groups()))
  interaction_list = pair_interaction.run(ph.deep_copy())
  clusters = betweenness_centrality_clustering(interaction_list,
    maxnum_residues_in_cluster=maxnum, bcc_threshold=9, size=n).get_clusters()
  check_clusters(clusters, interaction_list, n, maxnum)
  bonded_pairs = bonded_residue_pairs(ph)
  assert (3, 4) in bonded_pairs and (5, 6) in bonded_pairs, bonded_pairs
  clusters = betweenness_centrality_clustering(interaction_list,
    maxnum_residues_in_cluster=maxnum, bcc_threshold=9, size=n,
    bonded_pairs=bonded_pairs).get_clusters()
  check_clusters(clusters, interaction_list, n, maxnum)
  assert [3, 4, 5, 6] in clusters and [13, 14] in clusters, clusters

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Clusters larger than maxnum_residues_in_cluster that betweenness centrality
  clustering cannot split are split into connected clusters of at most
  maxnum_residues_in_cluster residues, cutting as few covalent bonds as
  possible.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_clique()
  exercise_1yjp()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
