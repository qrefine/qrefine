from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import math
import time
from scitbx.array_family import flex
from qrefine.tests.unit import run_fmodel
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.tst_06 import write_inputs, r_work_of
from qrefine.tests.unit.tst_06 import check_refined
from qrefine.tests.unit.tst_06 import refine_fast
from qrefine.tests.unit.tst_06 import get_model
from qrefine.tests.unit import tst_25
from qrefine.fragment import fragments
from qrefine import conformers, refine
from libtbx.test_utils import approx_equal

clustering_args = ["restraints=cctbx", "clustering=True",
  "altlocs.handling=conformers", "maxnum_residues_in_cluster=3",
  "parallel.nproc=1"]

def rmsd(model_1, model_2):
  """
  rms coordinate difference of two models with the same atom order.
  """
  h1 = model_1.get_hierarchy()
  h2 = model_2.get_hierarchy()
  assert [a.id_str() for a in h1.atoms()] == [a.id_str() for a in h2.atoms()]
  d = model_1.get_sites_cart() - model_2.get_sites_cart()
  return math.sqrt(flex.mean(d.dot()))

def fragment_counts(pdb_file, maxnum):
  """
  Number of fragments and of conformer variant jobs that clustering (cctbx
  restraints, bcc) makes of the model in pdb_file.
  """
  model = get_model(pdb_file)
  ph = model.get_hierarchy()
  fm = fragments(
    cif_objects                = model.get_restraint_objects(),
    working_folder             = None,
    make_working_folder        = False,
    clustering_method          = "bcc",
    maxnum_residues_in_cluster = maxnum,
    pdb_hierarchy              = ph,
    crystal_symmetry           = model.crystal_symmetry(),
    qm_engine_name             = "mopac",
    cap_fragments              = False,
    conformer_selections       = conformers.conformer_selections(ph))
  return fm.n_fragments, len(fm.variant_fragment_index)

def exercise_routes(pdb_poor, mtz, r_fm, r_start, prefix):
  """
  The three routes on the altloc toy model: (a) native, whole model;
  (b) conformers, whole model; (c) conformers with clustering (fragments
  evaluated per conformer variant). (b) and (c) evaluate the same per
  conformer restraints, (a) is cctbx's own altloc-aware treatment.
  """
  n_frag, n_jobs = fragment_counts(pdb_poor, maxnum=3)
  print("maxnum_residues_in_cluster=3: %d fragments, %d jobs" % (
    n_frag, n_jobs))
  assert n_frag > 1, n_frag
  assert n_jobs > n_frag, (n_frag, n_jobs)
  r_a = run_qrefine_in(folder="a_native", prefix=prefix+"_a",
    args=[pdb_poor, mtz, "restraints=cctbx", "clustering=False",
      "altlocs.handling=native"] + refine_fast)
  check_refined(r_qr=r_a, pdb_poor=pdb_poor, r_fm=r_fm, r_start=r_start)
  r_b = run_qrefine_in(folder="b_conformers", prefix=prefix+"_b",
    args=[pdb_poor, mtz, "restraints=cctbx", "clustering=False",
      "altlocs.handling=conformers"] + refine_fast)
  check_refined(r_qr=r_b, pdb_poor=pdb_poor, r_fm=r_fm, r_start=r_start)
  r_c = run_qrefine_in(folder="c_clustering", prefix=prefix+"_c",
    args=[pdb_poor, mtz] + clustering_args + refine_fast)
  check_refined(r_qr=r_c, pdb_poor=pdb_poor, r_fm=r_fm, r_start=r_start)
  assert r_c.std_out_text.find(
    "Clustering with 2 conformers: per-fragment variants") > -1
  assert r_b.std_out_text.find("per-fragment variants") == -1
  assert r_c.std_out_text.find("Job complete") > -1
  # cctbx fragments write no job folders; no per-conformer scratch folders
  for r in [r_b, r_c]:
    assert not os.path.exists(os.path.join(r.folder, "ase", "A"))
    assert not os.path.exists(os.path.join(r.folder, "ase", "B"))
  # (b) and (c) evaluate the same per-conformer cctbx restraints (measured
  # 0 at PDB precision); (a) treats the altlocs natively (measured 0.10 A).
  rmsd_ab = rmsd(r_a.model, r_b.model)
  rmsd_ac = rmsd(r_a.model, r_c.model)
  rmsd_bc = rmsd(r_b.model, r_c.model)
  print("rmsd(a,b): %8.5f rmsd(a,c): %8.5f rmsd(b,c): %8.5f" % (
    rmsd_ab, rmsd_ac, rmsd_bc))
  assert rmsd_bc < 0.01, rmsd_bc
  assert rmsd_ab < 0.2, rmsd_ab
  assert rmsd_ac < 0.2, rmsd_ac

def exercise_reclustering(pdb_poor, mtz, r_fm, r_start, prefix):
  """
  Route (c) with re-clustering forced at every check.
  """
  r = run_qrefine_in(folder="reclustering", prefix=prefix+"_reclustering",
    args=[pdb_poor, mtz] + clustering_args + [
      "cluster.re_calculate_rmsd_tolerance=0",
      "number_of_weight_search_cycles=1", "number_of_refine_cycles=2"])
  n = r.std_out_text.count("need to redo clustering")
  print("re-clustering events:", n)
  assert n >= 2, n
  assert r.std_out_text.count("per-fragment variants") >= 1
  check_refined(r_qr=r, pdb_poor=pdb_poor, r_fm=r_fm, r_start=r_start)

def exercise_hd_exchange():
  """
  H/D exchange sites in A/B atom groups next to heavy-atom altlocs A/B
  (tst_25 model), neutron data, route (c): hd_mapper and per-fragment
  conformer variants together.
  """
  labels = ("A","B")
  ph, cs = tst_25.read_hierarchy()
  model_h = tst_25.get_model(ph, cs)
  model_hd, _ = tst_25.split_hd(model_h, n_residues=3, labels=labels)
  pdb = os.path.abspath("model_hd.pdb")
  with open(pdb, "w") as fo:
    fo.write(model_hd.model_as_pdb())
  r_fm = run_fmodel(args=[pdb, "high_resolution=2.0",
    "scattering_table=neutron"], prefix="model_hd")
  r = run_qrefine_in(folder="hd_exchange", prefix="hd_exchange",
    args=[pdb, os.path.abspath(r_fm.mtz), "scattering_table=neutron"] +
      clustering_args + refine_fast)
  assert r.std_out_text.find("Using hd_mapper") > -1
  assert r.std_out_text.find("per-fragment variants") > -1
  assert r.std_out_text.find("Best Rwork, Rfree (after refinement)") > -1
  assert r.std_out_text.find("Job complete") > -1
  assert r.model is not None
  h = r.model.get_hierarchy()
  assert sorted(h.altloc_indices().keys()) == ["", "A", "B"]
  assert h.atoms_size() == model_hd.size()
  # H and its D copy stay together (hd_mapper)
  pairs = refine.exchangeable_hd_pairs(h)
  assert len(pairs) == 3, pairs
  xyz = h.atoms().extract_xyz()
  for pair in pairs:
    assert approx_equal(xyz[pair[0][0]], xyz[pair[1][0]])

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine, cctbx restraints, model with alternative conformations: routes
  (a) native, (b) conformers and (c) conformers with clustering agree;
  forced re-clustering; H/D exchange with heavy-atom altlocs.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb_good, pdb_poor = write_inputs(prefix)
  r_fm = run_fmodel(prefix=prefix, args=[pdb_good, "high_res=2"])
  mtz = os.path.abspath(r_fm.mtz)
  r_start = r_work_of(pdb_poor, r_fm)
  exercise_routes(pdb_poor=pdb_poor, mtz=mtz, r_fm=r_fm, r_start=r_start,
    prefix=prefix)
  exercise_reclustering(pdb_poor=pdb_poor, mtz=mtz, r_fm=r_fm,
    r_start=r_start, prefix=prefix)
  exercise_hd_exchange()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
