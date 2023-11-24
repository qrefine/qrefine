from __future__ import division
from __future__ import absolute_import
import sys
import time
from qrefine import qr
import iotbx.pdb
from boost_adaptbx import graph
from boost_adaptbx.graph import clustering_algorithm
from boost_adaptbx.graph import connected_component_algorithm as cca
from libtbx.utils import Sorry
from libtbx.program_template import ProgramTemplate
from qrefine.fragment import fragments



class girvan_nweman_clustering(object):

  def __init__(self, interaction_list,maxnum_residues_in_cluster=15, 
                     size = None):
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
    except ImportError as e:
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
    except ImportError as e:
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
  def __init__(self, interaction_list,maxnum_residues_in_cluster=20,bcc_threshold=9,size=None):
    self.interaction_list = interaction_list
    self.g = graph.adjacency_list(
      graph_type="undirected",
      vertex_type="vector",
      edge_type="set")
    self.maxnum_residues_in_cluster=maxnum_residues_in_cluster
    self.bcc_threshold=bcc_threshold
    self.size = size

  def get_clusters(self):
    self.build_graph()
    threshold = self.bcc_threshold
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
    print "max(components_size): ",max(components_size)    
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

class Program(ProgramTemplate):

  description = """
  qr.cluster Cluster a system into many small pieces


  Example:
  qr.cluster model.pdb  [<param_name>=<param_value>] ...
  """

  local_phil ="""
  maxnum_residues_in_cluster = 25
      .type = int
      .help = maximum number of residues in a cluster
  bcc_threshold = 9
      .type = int
      .help = threshold value for bcc
  """


  datatypes = ['model', 'phil', ]

  master_phil_str = qr.master_phil_str + local_phil

  def validate(self):
    print('Validate inputs:', file=self.logger)
    self.data_manager.has_models(
      expected_n=1,
      exact_count=True,
      raise_sorry=True)


  def run(self):
    self.header("Refinement start")
    print("max number of residues in each cluster:\n", self.params.maxnum_residues_in_cluster, file=log)
    print("bcc threshold value:\n", self.params.bcc_threshold, file=log)
    ph = self.data_manager.get_model().get_hierarchy()
    cs = self.data_manager.get_model().crystal_symmetry()
    fq = fragments(
      pdb_hierarchy=ph,
      crystal_symmetry=cs,
      maxnum_residues_in_cluster=self.params.maxnum_residues_in_cluster,
      bcc_threshold = self.params.bcc_threshold,
      clusters_only = True)
    print("Residue indices for each cluster:\n", fq.clusters, file=log)
    print('# clusters  : ',len(fq.clusters), file=log)

