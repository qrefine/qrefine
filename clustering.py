from __future__ import division
from boost_adaptbx import graph
from boost_adaptbx.graph import clustering_algorithm
from graph import connected_component_algorithm as cca
from libtbx.utils import Sorry

class girvan_nweman_clustering(object):

  def __init__(self, interaction_list,maxnum_residues_in_cluster=15,size = None):
        self.interaction_list = interaction_list
        self.size = size
        self.maxnum_residues_in_cluster = maxnum_residues_in_cluster

  def get_clusters(self):
   self.clusters = self.first_clustering()
   self.clusters = self.second_clustring(self.clusters)
   return self.clusters

  def sublistExists(self,big_list, small_list):
    return set(small_list).issubset(set(big_list))

  def first_clustering(self):
    try:
      import networkx as nx
      import plugin.cmty as gn
    except ImportError, e:
      raise Sorry(str(e))
    Gx = nx.Graph()
    if self.size is not None:
      for item in range(self.size):
        Gx.add_node(item+1)
    weight = [1.0]*len(self.interaction_list)
    es_triple = []
    for i in range(len(weight)):
      temp = (self.interaction_list[i][0],self.interaction_list[i][1],weight[i])
      es_triple.append(temp)
    Gx.add_weighted_edges_from( es_triple)
    clusters =gn.run_graph(Gx)
    return clusters

  def sub_clustering(self,cluster,final):
    try:
      import networkx as nx
      import plugin.cmty as gn
    except ImportError, e:
      raise Sorry(str(e))
    if len(cluster) < self.maxnum_residues_in_cluster:
      final.append(cluster)
    if len(cluster) >= self.maxnum_residues_in_cluster:
      es_triple = []
      for interaction in self.interaction_list:
        if self.sublistExists(cluster,interaction):
          temp = (interaction[0],interaction[1],1.0)
          es_triple.append(temp)
      Gx = nx.Graph()
      Gx.add_weighted_edges_from( es_triple)
      subclusters = gn.run_graph(Gx,5)
      subcluster_is_list = True
      for subcluster in subclusters:
        subcluster_is_list = isinstance(subcluster, list)
      if subcluster_is_list == True:
        for subcluster in subclusters:
          self.sub_clustering(subcluster,final)
      else:
         final.append(subclusters)

  def second_clustring(self,clusters) :
    final_clusters = []
    for cluster in clusters:
       self.sub_clustering(cluster,final_clusters)
    return final_clusters

class betweenness_centrality_clustering(object):
  def __init__(self, interaction_list,maxnum_residues_in_cluster=20,size=None):
    self.interaction_list = interaction_list
    self.g = graph.adjacency_list(
      graph_type="undirected",
      vertex_type="vector",
      edge_type="set")
    self.maxnum_residues_in_cluster=maxnum_residues_in_cluster
    self.size = size

  def get_clusters(self):
    self.build_graph()
    threshold = 9
    clustering = True
    while (clustering and threshold >=4):
      edge_centrality_map = clustering_algorithm.\
        betweenness_centrality_clustering(graph=self.g, threshold=threshold)
      components = cca.connected_components(graph=self.g)
      components_size = [ len(component) for component in components]
      #print "max(components_size): ",max(components_size)
      if max(components_size) <= self.maxnum_residues_in_cluster:
        clustering = False
      else:
        threshold -= 1
    final_components = []
    for component in components:
      component[:] = [x + 1 for x in component]
      if len(component) > self.maxnum_residues_in_cluster:
        component.sort()
        index = int(len(component)/2)
        final_components.append(component[:index])
        final_components.append(component[index:])
      else:
        final_components.append(component)
    return final_components

  def build_graph(self):
    import itertools
    vertices = []
    if self.size is None:
      merged = list(itertools.chain.from_iterable( self.interaction_list))
      self.size = max(merged)
    for i in range(self.size):
      vertices.append(self.g.add_vertex())
    for pair in self.interaction_list:
      self.g.add_edge(vertex1 = vertices[pair[0]-1],
        vertex2 = vertices[pair[1]-1], weight = 1)
