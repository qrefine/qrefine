from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from qrefine import qr
from boost_adaptbx import graph
from boost_adaptbx.graph import clustering_algorithm
from boost_adaptbx.graph import connected_component_algorithm as cca
from libtbx.program_template import ProgramTemplate
from qrefine.fragment import fragments

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

