from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex
from qrefine import conformers, refine, restraints, cluster_restraints
from qrefine.fragment import fragments
from qrefine.tests.unit import tst_46
from qrefine.tests.unit.tst_64 import three_label_model, hetero_model
from qrefine.tests.unit.tst_11 import qr_unit_tests

data_files = os.path.join(qr_unit_tests, "data_files")

altloc_files = ["altlocs.pdb", "altlocs2.pdb", "h_altconf.pdb",
  "h_altconf_2.pdb", "h_altconf_2_complete.pdb",
  "h_altconf_2_complete_modified.pdb", "h_altconf_complete.pdb"]

def own_box_model(file_name):
  """
  Processed model of data_files/file_name in the box of the file.
  """
  pdb_inp = iotbx.pdb.input(file_name=os.path.join(data_files, file_name))
  return tst_46.process(mmtbx.model.manager(model_input=pdb_inp,
    log=null_out()))

def contact_model(pdb_hierarchy, contact, name):
  """
  tst_46.get_model for a hierarchy (written as <name>.pdb first): processed
  model in a P1 box with the closest model-image contact at most contact A.
  """
  local = name + ".pdb"
  pdb_hierarchy.write_pdb_file(file_name=local)
  lines = open(local).read().splitlines()
  lines = [l for l in lines if l.startswith(("ATOM", "HETATM", "TER"))]
  pdb_inp = iotbx.pdb.input(source_info=None, lines=lines)
  xyz = pdb_inp.atoms().extract_xyz()
  extent = [b-a for a, b in zip(xyz.min(), xyz.max())]
  for gap in [25.0] + [6.0 - 0.5*i for i in range(19)]:
    box = [e+gap for e in extent]
    d = tst_46.closest_image_distance(xyz, box)
    if(d <= contact): break
  assert d <= contact and d > 1.0, (gap, d)
  cryst1 = "CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1" % tuple(box)
  file_name = "%s_contact%.1f.pdb" % (name, contact)
  with open(file_name, "w") as fo:
    fo.write("\n".join([cryst1] + lines) + "\n")
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  return tst_46.process(mmtbx.model.manager(model_input=pdb_inp,
    log=null_out()))

def clustered(model, handling, maxnum):
  """
  cctbx clustered restraints manager of model, altlocs handling native or
  conformers (fragments with conformer variants); no files.
  """
  params = tst_46.get_cctbx_params(expansion=False, clustering=True,
    handling=handling, maxnum=maxnum)
  return refine.create_restraints_manager(params, model, working_folder=None,
    log=null_out())

def expansion(model, handling):
  params = tst_46.get_cctbx_params(expansion=True, clustering=False,
    handling=handling, maxnum=999)
  return refine.create_restraints_manager(params, model, working_folder=None,
    log=null_out())

def gradients(rm, model):
  return rm.target_and_gradients(model.get_sites_cart())[1]

def max_diff(g1, g2):
  return flex.max((g1 - g2).norms())

def check_multiplicity(rm, model):
  m = conformers.multiplicity(model.get_hierarchy())
  assert rm.multiplicity.size() == m.size()
  assert rm.multiplicity.all_eq(m), list(rm.multiplicity)

def exercise_one_cluster_limit():
  """
  maxnum=999 (one cluster): the clustered conformers route equals
  from_conformers over expansion, with images in the super sphere.
  """
  models = [
    ("h_altconf_complete", tst_46.get_model("h_altconf_complete.pdb",
      contact=2.0)),
    ("three_label", contact_model(three_label_model().get_hierarchy(),
      contact=2.0, name="three_label"))]
  for name, model in models:
    rm_c = clustered(model, "conformers", maxnum=999)
    rm_e = expansion(model, "conformers")
    assert isinstance(rm_e, conformers.from_conformers)
    g_c = gradients(rm_c, model)
    g_e = gradients(rm_e, model)
    d = max_diff(g_c, g_e)
    print("%s one cluster: max |clustered - expansion| = %.2e" % (name, d))
    assert d < 1.e-6, d
    n_conf = conformers.n_conformers(model.get_hierarchy())
    assert rm_c.fragment_manager.n_fragments == 1
    assert rm_c.n_jobs == n_conf, (rm_c.n_jobs, n_conf)
    check_multiplicity(rm_c, model)
  # Target: without images in the super sphere (with images, the buffer of
  # a variant is that of the altloc fragment, a superset, and buffer-buffer
  # terms enter the fragment energy)
  models = [
    ("h_altconf_complete", tst_46.get_model("h_altconf_complete.pdb",
      contact=1000.)),
    ("three_label", tst_46.process(three_label_model()))]
  for name, model in models:
    rm_c = clustered(model, "conformers", maxnum=999)
    assert rm_c.fragment_manager.pdb_hierarchy_super.atoms_size() == \
      model.get_number_of_atoms()
    check_target_one_cluster(rm_c, model, name)

def conformer_mean_target(model):
  """
  Independent reference for the clustered conformers target with one
  cluster: the mean over conformers of the native clustered target of each
  single-conformer model (one fragment: the cctbx energy of that
  conformer's super sphere).
  """
  targets = []
  for c in conformers.split(model):
    m = tst_46.process(c.model)
    rm = clustered(m, "native", maxnum=999)
    assert rm.fragment_manager.n_fragments == 1
    targets.append(rm.target_and_gradients(m.get_sites_cart())[0])
  return sum(targets) / len(targets), targets

def check_target_one_cluster(rm_c, model, name):
  t_c = rm_c.target_and_gradients(model.get_sites_cart())[0]
  t_ref, targets = conformer_mean_target(model)
  assert len(set(targets)) == len(targets) # conformers differ
  rel = abs(t_c - t_ref) / abs(t_ref)
  print("%s one cluster: target %.6f, conformer mean %.6f (rel %.1e)" % (
    name, t_c, t_ref, rel))
  assert rel < 1.e-9, (t_c, t_ref)

def exercise_fragments_without_altlocs():
  """
  tst_64.hetero_model (blank peptide, far water with conformers A, B),
  maxnum=1: fragments without altloc atoms (one variant, weight
  n_conformers; fragment 0 is one of them) and with them; every
  fragment's weights sum to n_conformers, and the target is the mean over
  conformers of the sum of fragment energies: sum over jobs of
  w*e / n_conformers.
  """
  model = tst_46.process(hetero_model())
  rm = clustered(model, "conformers", maxnum=1)
  fm = rm.fragment_manager
  t = rm.target_and_gradients(model.get_sites_cart())[0]
  n_conf = conformers.n_conformers(model.get_hierarchy())
  assert n_conf == 2
  single = [i for i in range(fm.n_fragments)
    if fm.variant_fragment_index.count(i) == 1]
  multi = [i for i in range(fm.n_fragments)
    if fm.variant_fragment_index.count(i) > 1]
  assert 0 in single and len(multi) > 0
  energies = []
  sites_super = fm.pdb_hierarchy_super.atoms().extract_xyz()
  for j, selection in enumerate(fm.fragment_selections):
    energies.append(rm.restraints_manager([selection, sites_super, j])[0])
  for i in range(fm.n_fragments):
    jobs = [j for j, k in enumerate(fm.variant_fragment_index) if k == i]
    assert sum([fm.variant_weights[j] for j in jobs]) == n_conf
  expected = sum([w*e for w, e in zip(fm.variant_weights, energies)])
  expected /= n_conf
  print("hetero maxnum=1: %d fragments (%d without altlocs, fragment 0 "
    "among them), target %.6f, expected %.6f" % (fm.n_fragments,
    len(single), t, expected))
  assert abs(t - expected) <= 1.e-9 * abs(expected), (t, expected)
  # n_conformers 1 (weights of one variant) would be off by a factor 2
  assert abs(t*n_conf - expected) > 1.e-3 * abs(expected)

def blank_and_altloc(ph):
  ci = ph.get_conformer_indices().conformer_indices
  blank = ci == 0
  return blank, ~blank

def exercise_relation_per_atom():
  """
  Native clustered (a) vs conformers clustered (c), same clustering.
  One cluster: altloc atoms equal, blank atoms m*g_c == g_a +
  (m-1)*g_blank_only (tst_31's relation; g_blank_only from expansion on the
  blank atoms). Real clusters (maxnum=3): altloc atoms equal.
  """
  for file_name in ["h_altconf_complete.pdb", "h_altconf_2_complete.pdb"]:
    model = own_box_model(file_name)
    ph = model.get_hierarchy()
    blank, altloc = blank_and_altloc(ph)
    m = conformers.multiplicity(ph)
    assert m.select(altloc).all_eq(1.0)
    assert blank.count(True) > 0 and altloc.count(True) > 0
    bsel = blank.iselection()
    model_blank = tst_46.process(mmtbx.model.manager(model_input=None,
      pdb_hierarchy=ph.select(bsel), crystal_symmetry=model.crystal_symmetry(),
      log=null_out()))
    g_blank_only = gradients(expansion(model_blank, "native"), model_blank)
    for maxnum in [999, 3]:
      rm_a = clustered(model, "native", maxnum)
      rm_c = clustered(model, "conformers", maxnum)
      g_a = gradients(rm_a, model)
      g_c = gradients(rm_c, model)
      check_multiplicity(rm_c, model)
      assert rm_a.multiplicity.all_eq(1.0)
      assert rm_a.n_jobs == rm_a.fragment_manager.n_fragments
      assert rm_c.n_jobs > rm_c.fragment_manager.n_fragments
      d_alt = max_diff(g_c.select(altloc), g_a.select(altloc))
      print("%s maxnum=%d: %d fragments, %d jobs; altloc atoms "
        "max |c - a| = %.2e" % (file_name, maxnum,
        rm_c.fragment_manager.n_fragments, rm_c.n_jobs, d_alt))
      assert d_alt < 1.e-9, d_alt
      if(maxnum != 999): continue
      mb = m.select(bsel)
      lhs = g_c.select(bsel) * mb
      rhs = g_a.select(bsel) + g_blank_only * (mb - 1)
      d_rel = max_diff(lhs, rhs)
      print("%s one cluster: blank atoms max |m*g_c - (g_a + "
        "(m-1)*g_blank_only)| = %.2e" % (file_name, d_rel))
      assert d_rel < 1.e-6, d_rel
      # the relation is not trivial: blank atoms bonded to altlocs differ
      assert max_diff(g_c.select(bsel), g_a.select(bsel)) > 1.e-3

def exercise_clustering_vs_expansion():
  """
  Altloc data files in their own boxes, maxnum=3, conformers handling:
  clustering vs expansion (from_conformers) within the frozen-buffer
  tolerance of single-conformer clustering.
  """
  n_done = 0
  for file_name in altloc_files:
    model = own_box_model(file_name)
    ph = model.get_hierarchy()
    try:
      conformers.check_completed(ph)
      conformers.check_all_links(ph)
    except Sorry as e:
      print("%s: skipped (%s)" % (file_name, str(e)[:60]))
      continue
    rm_c = clustered(model, "conformers", maxnum=3)
    rm_e = expansion(model, "conformers")
    g_c = gradients(rm_c, model)
    g_e = gradients(rm_e, model)
    d = max_diff(g_c, g_e)
    print("%s: %d fragments, %d jobs, max |clustered - expansion| = %.2e" % (
      file_name, rm_c.fragment_manager.n_fragments, rm_c.n_jobs, d))
    assert d < 1.e-4, d
    check_multiplicity(rm_c, model)
    n_done += 1
  assert n_done >= 4, n_done

def unweighted_sum(rm, model):
  """
  The combination before conformer variants: plain sum over the jobs of
  the job gradients, buffer atoms zeroed per job.
  """
  fm = rm.fragment_manager
  fm.update_xyz(model.get_sites_cart())
  rm.restraints_manager.fragment_extracts = fm.get_fragment_extracts(
    with_super_sphere_grm=True)
  sites_super = fm.pdb_hierarchy_super.atoms().extract_xyz()
  n = model.get_number_of_atoms()
  target = 0
  result = flex.vec3_double(n)
  for index, selection in enumerate(fm.fragment_selections):
    t, g = rm.restraints_manager([selection, sites_super, index])
    g_i = flex.vec3_double(n).set_selected(selection, g)
    g_i = g_i.set_selected(fm.buffer_selections[index], [0,0,0])
    result += g_i
    target += t
  return target, result

def direct(model, selections, maxnum):
  """
  cctbx clustered restraints manager built from fragments directly with the
  given conformer_selections (None or []: create_restraints_manager never
  passes [], it dispatches a model without conformers as None).
  """
  params = tst_46.get_cctbx_params(expansion=False, clustering=True,
    handling="conformers", maxnum=maxnum)
  fm = fragments(
    cif_objects                = model.get_restraint_objects(),
    working_folder             = None,
    make_working_folder        = False,
    clustering_method          = "bcc",
    maxnum_residues_in_cluster = maxnum,
    pdb_hierarchy              = model.get_hierarchy(),
    crystal_symmetry           = model.crystal_symmetry(),
    qm_engine_name             = "mopac",
    cap_fragments              = False,
    conformer_selections       = selections)
  return cluster_restraints.from_cluster(
    restraints_manager = restraints.from_cctbx(restraints_manager=None),
    fragment_manager   = fm,
    parallel_params    = params.parallel)

def exercise_single_conformer_bit_identical():
  """
  No altlocs: conformer_selections None and [] give the same gradients as
  the plain sum, to the bit; one job per fragment, multiplicity 1.
  """
  model = own_box_model("helix.pdb")
  ph = model.get_hierarchy()
  assert conformers.conformer_selections(ph) == []
  for maxnum in [1, 999]:
    rm_n = direct(model, None, maxnum)
    rm_c = direct(model, conformers.conformer_selections(ph), maxnum)
    t_n, g_n = rm_n.target_and_gradients(model.get_sites_cart())
    t_c, g_c = rm_c.target_and_gradients(model.get_sites_cart())
    t_r, g_r = unweighted_sum(rm_n, model)
    assert flex.max((g_n - g_c).norms()) == 0
    assert flex.max((g_n - g_r).norms()) == 0
    assert list(g_n) == list(g_r)
    assert t_n == t_c == t_r
    for rm in [rm_n, rm_c]:
      assert rm.n_jobs == rm.fragment_manager.n_fragments
      assert rm.multiplicity.all_eq(1.0)
    print("helix maxnum=%d: %d fragments, bit identical" % (maxnum,
      rm_n.n_jobs))

def exercise_single_conformer_only_water():
  """
  tst_31.model_with_b_only_water (HOH 100 in conformer B only) in a small
  box (images in the super sphere): maxnum 999 (two clusters, HOH 101
  does not interact) equals from_conformers over expansion (both route
  through conformers._selections), and with
  small clusters the multiplicity is the membership of the rule (1 on the
  water). The expansion result itself is checked against cctbx on
  single-conformer models in tst_31.
  """
  from qrefine.tests.unit import tst_31
  model = contact_model(tst_31.model_with_b_only_water().get_hierarchy(),
    contact=2.0, name="b_only_water")
  ph = model.get_hierarchy()
  water = flex.size_t([a.i_seq for a in ph.atoms()
    if a.parent().parent().resseq_as_int() == 100])
  assert conformers.multiplicity(ph).select(water).all_eq(1.0)
  rm_c = clustered(model, "conformers", maxnum=999)
  rm_e = expansion(model, "conformers")
  assert isinstance(rm_e, conformers.from_conformers)
  assert [c.n_absent for c in rm_e.conformers] == [1, 0]
  d = max_diff(gradients(rm_c, model), gradients(rm_e, model))
  print("B-only water maxnum 999: max |clustered - expansion| = %.2e" % d)
  assert d < 1.e-6, d
  print("B-only water maxnum 999: %d fragments, %d jobs" % (
    rm_c.fragment_manager.n_fragments, rm_c.n_jobs))
  check_multiplicity(rm_c, model)
  rm_c = clustered(model, "conformers", maxnum=1)
  g = gradients(rm_c, model)
  check_multiplicity(rm_c, model)
  fm = rm_c.fragment_manager
  print("B-only water maxnum 1: %d fragments, %d jobs" % (fm.n_fragments,
    rm_c.n_jobs))
  assert rm_c.n_jobs < 2*fm.n_fragments

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  """
  Weighted combination of conformer variant gradients in the clustered
  route (cluster_restraints.from_cluster), cctbx restraints.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_single_conformer_bit_identical()
  exercise_one_cluster_limit()
  exercise_fragments_without_altlocs()
  exercise_relation_per_atom()
  exercise_clustering_vs_expansion()
  exercise_single_conformer_only_water()

if(__name__ == "__main__"):
  run()
  print("OK")
