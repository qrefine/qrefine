from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from qrefine import qr
import itertools
from libtbx.program_template import ProgramTemplate
from qrefine.fragment import fragments

def connected_pieces(residues, neighbors):
  """
  The connected components of the residues (a list) in the graph given by
  neighbors (residue -> set of residues), each sorted, in order of their
  smallest residue.
  """
  left = set(residues)
  pieces = []
  while(left):
    stack = [min(left)]
    piece = set()
    while(stack):
      r = stack.pop()
      if(r in piece): continue
      piece.add(r)
      stack.extend((neighbors[r] & left) - piece)
    pieces.append(sorted(piece))
    left -= piece
  return pieces

BOND_WEIGHT = 10 # a covalent bond versus one interaction, ordering only

def bisect(residues, neighbors, bonded):
  """
  Cut the connected cluster residues (a list) in two parts, each of at least
  a third of the residues, cutting the fewest covalent bonds (bonded:
  residue -> set of residues it is bonded to) and, among those cuts, the
  fewest interactions (edges of neighbors) per pair of residues separated
  (the ratio cut: a cut through a dense block is balanced, a cut of a chain
  is at its weakest link); ties go to the more balanced cut. The candidate
  cuts are the prefixes of the residues ordered along the Fiedler vector
  (the eigenvector of the second smallest eigenvalue of the graph Laplacian,
  bonds weighing BOND_WEIGHT), which orders a chain along the chain and a
  compact blob along its longest extent. Residue indices break ties in the
  order, and the sign of the vector (arbitrary from the eigensolver) is
  fixed by its first non-zero component.
  """
  n = len(residues)
  index = dict((r, i) for i, r in enumerate(residues))
  laplacian = np.zeros((n, n))
  for r in residues:
    for s in neighbors[r] | bonded[r]:
      if(s in index):
        w = BOND_WEIGHT if s in bonded[r] else 1
        laplacian[index[r], index[s]] -= w
        laplacian[index[r], index[r]] += w
  fiedler = np.linalg.eigh(laplacian)[1][:, 1]
  first = fiedler[np.abs(fiedler) > 1.e-9][0]
  if(first < 0): fiedler = -fiedler
  order = sorted(range(n), key=lambda i: (fiedler[i], residues[i]))
  position = dict((residues[i], p) for p, i in enumerate(order))
  def crossing(pairs, k): # pairs of positions separated by a cut after k
    return sum(1 for p, q in pairs if p < k <= q)
  edges = [(position[r], position[s]) for r in residues for s in neighbors[r]
           if s in index and position[r] < position[s]]
  bonds = [(position[r], position[s]) for r in residues for s in bonded[r]
           if s in index and position[r] < position[s]]
  smallest_part = max(1, (n+2)//3)
  best = None # (bonds cut, ratio cut, imbalance, size of the first part)
  for k in range(smallest_part, n-smallest_part+1):
    candidate = (crossing(bonds, k), crossing(edges, k)/(k*(n-k)),
                 abs(2*k-n), k)
    if(best is None or candidate < best): best = candidate
  k = best[3]
  return [sorted(residues[i] for i in order[:k]),
          sorted(residues[i] for i in order[k:])]

def split_cluster(residues, neighbors, bonded, maxnum):
  """
  Split the cluster residues (a list) into connected clusters of at most
  maxnum residues: a cluster that is not connected falls apart into its
  pieces, a connected one too large is bisected, and the parts are split
  again as needed.
  """
  assert maxnum >= 1
  pieces = connected_pieces(residues, neighbors)
  if(len(pieces) == 1 and len(residues) <= maxnum): return [sorted(residues)]
  if(len(pieces) == 1): pieces = bisect(residues, neighbors, bonded)
  result = []
  for piece in pieces:
    result.extend(split_cluster(piece, neighbors, bonded, maxnum))
  return result

class betweenness_centrality_clustering(object):
  """
  Clusters (lists of 1-based residue indices) of a system of size residues
  from its interaction_list (pairs of residue indices): betweenness
  centrality clustering of the interaction graph, with the threshold
  lowered from bcc_threshold (down to 4) until no cluster has more than
  maxnum_residues_in_cluster residues; clusters still too large (a dense
  block of residues has no edge of high centrality) are split by
  split_cluster, which cuts as few of the bonded_pairs (covalently bonded
  residue pairs, optional) as it can.
  """
  def __init__(self, interaction_list,maxnum_residues_in_cluster=20,bcc_threshold=9,size=None,
               bonded_pairs=None):
    self.interaction_list = interaction_list
    self.bonded_pairs = bonded_pairs or []
    self.maxnum_residues_in_cluster=maxnum_residues_in_cluster
    self.bcc_threshold=bcc_threshold
    self.size = size

  def get_clusters(self):
    """
    The clusters (lists of 1-based residue indices), in order of their
    smallest residue, each of at most maxnum_residues_in_cluster residues.
    """
    from mmtbx.pair_interaction import ext
    from scitbx.array_family import flex
    if self.size is None:
      merged = list(itertools.chain.from_iterable( self.interaction_list))
      self.size = max(merged)
    edges = flex.int()
    for i, j in self.interaction_list:
      edges.append(i - 1)
      edges.append(j - 1)
    residues = list(range(1, self.size + 1))
    threshold = self.bcc_threshold
    while True:
      # Girvan-Newman (mmtbx.pair_interaction.ext.girvan_newman, the boost
      # betweenness_centrality_clustering with component-local recompute)
      # continues on the edges the previous threshold left
      edges = ext.girvan_newman(n_vertices=self.size, edges=edges,
        threshold=threshold)
      remaining = dict((r, set()) for r in residues)
      for k in range(edges.size() // 2):
        a, b = edges[2*k] + 1, edges[2*k+1] + 1
        remaining[a].add(b)
        remaining[b].add(a)
      components = connected_pieces(residues, remaining)
      if max(len(c) for c in components) <= self.maxnum_residues_in_cluster:
        break
      if threshold <= 4: break
      threshold -= 1
    neighbors = dict((r, set()) for r in residues)
    for i, j in self.interaction_list:
      neighbors[i].add(j)
      neighbors[j].add(i)
    bonded = dict((r, set()) for r in residues)
    for i, j in self.bonded_pairs:
      bonded[i].add(j)
      bonded[j].add(i)
    final_components = []
    for component in components:
      final_components.extend(split_cluster(
        residues  = component,
        neighbors = neighbors,
        bonded    = bonded,
        maxnum    = self.maxnum_residues_in_cluster))
    return final_components

class Program(ProgramTemplate):

  description = """
  qr.cluster Cluster a system into many small pieces


  Example:
  qr.cluster model.pdb  [<param_name>=<param_value>] ...
  """


  datatypes = ['model', 'phil', ]

  master_phil_str = qr.master_phil_str

  def validate(self):
    print('Validate inputs:', file=self.logger)
    self.data_manager.has_models(
      expected_n=1,
      exact_count=True,
      raise_sorry=True)


  def run(self):
    log = self.logger
    self.header("Clustering start")
    print("max number of residues in each cluster: ", self.params.cluster.maxnum_residues_in_cluster, file=log)
    print("bcc threshold value: ", self.params.cluster.bcc_threshold, file=log)
    ph = self.data_manager.get_model().get_hierarchy()
    cs = self.data_manager.get_model().crystal_symmetry()
    fq = fragments(
      pdb_hierarchy=ph,
      crystal_symmetry=cs,
      maxnum_residues_in_cluster=self.params.cluster.maxnum_residues_in_cluster,
      bcc_threshold = self.params.cluster.bcc_threshold,
      select_within_radius = self.params.cluster.select_within_radius,
      clusters_only = True)
    self.clusters = fq.clusters
    print('# clusters: ',len(fq.clusters), file=log)
    res_sizes = ','.join(str(len(cluster)) for cluster in fq.clusters)
    print("# residues in clusters: [ %s ]" % res_sizes, file=log)

  def get_results(self):
    return self.clusters

