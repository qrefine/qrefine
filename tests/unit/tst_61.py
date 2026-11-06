from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
import iotbx.pdb
from mmtbx.pair_interaction import pair_interaction
from qrefine.clustering import betweenness_centrality_clustering

qr_unit_tests_data = os.path.join(libtbx.env.find_in_repositories("qrefine"),
  "tests", "unit", "data_files")

def hierarchy(file_name):
  return iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,
    file_name)).construct_hierarchy()

def as_sets(result):
  return tuple(sorted(x) for x in result)

def exercise_prepared(file_name, maxnum):
  """
  run_prepared(prepare(ph), cluster) gives what run(ph.deep_copy(), cluster)
  gives, for every cluster of the model; prepare leaves ph unchanged.
  """
  ph = hierarchy(file_name)
  xyz_before = ph.atoms().extract_xyz().deep_copy()
  n = len(list(ph.residue_groups()))
  interaction_list = pair_interaction.run(ph.deep_copy())
  clusters = betweenness_centrality_clustering(interaction_list,
    maxnum_residues_in_cluster=maxnum, bcc_threshold=9, size=n).get_clusters()
  assert len(clusters) > 1, clusters
  prepared = pair_interaction.prepare(ph)
  assert ph.atoms().extract_xyz().rms_difference(xyz_before) == 0
  assert prepared.n_residues == n
  for cluster in clusters:
    expected = as_sets(pair_interaction.run(ph.deep_copy(), cluster))
    got = as_sets(pair_interaction.run_prepared(prepared, cluster))
    assert got == expected, (file_name, cluster, got, expected)
  # the whole-model interaction list from the prepared data too
  assert sorted(pair_interaction.run_prepared(prepared, None)) == \
    sorted(interaction_list)
  assert ph.atoms().extract_xyz().rms_difference(xyz_before) == 0

def cluster_route(file_name, maxnum):
  """The cctbx cluster-route restraints manager of a model and the model."""
  from qrefine.tests.unit.tst_11 import get_model
  from qrefine import qr, refine
  model = get_model(file_name)
  params = qr.get_default_params()
  params.restraints = "cctbx"
  params.expansion = False
  params.cluster.clustering = True
  params.altlocs.handling = "native"
  params.cluster.maxnum_residues_in_cluster = maxnum
  params.parallel.nproc = 1
  return refine.create_restraints_manager(params=params, model=model), model

def exercise_grm_cache():
  """
  The cctbx cluster route selects the restraints of every fragment once:
  the second evaluation selects nothing and gives the same gradients; the
  super sphere is not written at every evaluation.
  """
  import mmtbx.restraints
  from scitbx.array_family import flex
  from qrefine.cluster_restraints import from_cluster
  rm, model = cluster_route("2lvr.pdb", maxnum=5)
  assert isinstance(rm, from_cluster)
  super_sphere_file = os.path.join("ase", "tmp_ase.pdb")
  if(os.path.exists(super_sphere_file)): os.remove(super_sphere_file)
  n_fragments = len(rm.fragment_manager.fragment_selections)
  assert n_fragments > 3, n_fragments
  sites = model.get_sites_cart()
  calls = []
  orig = mmtbx.restraints.manager.select
  def counting(self, *args, **kw):
    calls.append(1)
    return orig(self, *args, **kw)
  mmtbx.restraints.manager.select = counting
  try:
    g1 = rm.target_and_gradients(sites_cart=sites)[1]
    n1 = len(calls)
    g2 = rm.target_and_gradients(sites_cart=sites)[1]
    n2 = len(calls)
  finally:
    mmtbx.restraints.manager.select = orig
  assert n1 == n_fragments, (n1, n_fragments)
  assert n2 == n1, (n1, n2)
  assert flex.max(flex.abs(g1.as_double() - g2.as_double())) == 0
  assert not os.path.exists(super_sphere_file)

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_prepared("1yjp.pdb", maxnum=4)
  exercise_prepared("2olx.pdb", maxnum=3)
  exercise_grm_cache()

if(__name__ == "__main__"):
  run()
  print("OK")
