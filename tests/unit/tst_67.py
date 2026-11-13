from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import time
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx.array_family import flex
from qrefine import refine, conformers, restraints, cluster_restraints
from qrefine.tests.unit import tst_46, altloc_scenarios, run_qrefine_in
from qrefine.tests.unit.tst_31 import duplicate
from qrefine.tests.unit.tst_65 import contact_model, gradients, max_diff

def toy_hierarchy():
  """
  tst_00 peptide (GLY1 ASN2 ASN3 GLN4 GLN5 ASN6 TYR7, H, OXT) without
  waters.
  """
  ph = altloc_scenarios.base_model().get_hierarchy()
  sel = ph.atom_selection_cache().selection("not water")
  ph = ph.select(sel, copy_atoms=True) # edited in place later
  ph.atoms().reset_i_seq()
  ph.atoms_reset_serial()
  return ph

def single_toy():
  """
  The peptide in a P1 box with no image within the super sphere (contact
  30 A), processed.
  """
  return contact_model(toy_hierarchy(), contact=30., name="single_toy")

def super_toy():
  """
  The peptide with the side chain of ASN3 in conformers A, B (chi1 0 and
  -60, occupancies 0.5), in a P1 box with no image within the super
  sphere, processed; written to super_toy.pdb. The HB atoms stay blank, so
  a chi1 rotation of 120 puts CG of B on HB3 (MOPAC gradient 3.7e4
  kcal/mol/A); -60 is between them.
  """
  ph = toy_hierarchy()
  m = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph, log=null_out())
  altloc_scenarios.side_chain_conformers(
    rg           = altloc_scenarios.residue(m, 3),
    labels       = ["A", "B"],
    occupancies  = [0.5, 0.5],
    chi1_degrees = [0, -60])
  model = contact_model(m.get_hierarchy(), contact=30., name="super_toy_box")
  assert conformers.conformer_labels(model.get_hierarchy()) == ["A", "B"]
  altloc_scenarios.write(model, "super_toy.pdb")
  return model

def mopac_params(expansion, clustering, maxnum):
  params = tst_46.get_mopac_params(expansion=expansion,
    clustering=clustering)
  params.cluster.maxnum_residues_in_cluster = maxnum
  params.parallel.nproc = 1
  params.altlocs.handling = "conformers"
  return params

def qm_managers(rm):
  """
  Every from_qm the manager holds (from_conformers over expansion, or
  from_cluster).
  """
  if(isinstance(rm, conformers.from_conformers)):
    result = [c.restraints_manager.restraints_manager for c in rm.conformers]
  else:
    assert isinstance(rm, cluster_restraints.from_cluster)
    result = [rm.restraints_manager]
  for r in result: assert isinstance(r, restraints.from_qm)
  return result

def vacuum(rm):
  """
  MOPAC without COSMO: with EPS=78.4 its gradients are not consistent with
  its energies (tst_45), and routes are compared in vacuum only.
  """
  for r in qm_managers(rm): r.qm_engine.set(EPS=None)
  return rm

def manager(model, expansion, clustering, maxnum):
  return vacuum(refine.create_restraints_manager(
    params = mopac_params(expansion, clustering, maxnum),
    model  = model,
    log    = null_out()))

def touching(fm, residue_selection):
  """
  Number of fragments whose full super selection has an atom of
  residue_selection (no images: super atoms are the model atoms).
  """
  assert fm.pdb_hierarchy_super.atoms_size() == residue_selection.size()
  return [i for i, s in enumerate(fm.fragment_full_super_selections)
    if (s & residue_selection).count(True) > 0]

def residue_selection(model, resseq):
  return flex.bool([a.parent().parent().resseq_as_int() == resseq
    for a in model.get_hierarchy().atoms()])

def exercise_one_cluster_limit(model):
  """
  maxnum=999 (one cluster): the clustered conformers route evaluates the
  same capped sphere per conformer as from_conformers over expansion, so
  the gradients agree to the engine's noise; one job per conformer.
  """
  rm_c = manager(model, expansion=False, clustering=True, maxnum=999)
  rm_e = manager(model, expansion=True, clustering=False, maxnum=999)
  assert isinstance(rm_e, conformers.from_conformers)
  g_c = gradients(rm_c, model)
  g_e = gradients(rm_e, model)
  d = max_diff(g_c, g_e)
  print("one cluster: %d fragment(s), %d jobs, max |clustered - expansion| "
    "= %.2e kcal/mol/A (max |g| %.1f)" % (rm_c.fragment_manager.n_fragments,
    rm_c.n_jobs, d, flex.max(g_e.norms())))
  assert rm_c.fragment_manager.n_fragments == 1
  assert rm_c.n_jobs == 2, rm_c.n_jobs
  assert flex.max(g_e.norms()) > 1 # not vacuous
  assert d < 1.e-4, d

def exercise_duplicate():
  """
  The single-conformer peptide clustered (maxnum=1) vs the same with ASN3
  duplicated into identical conformers A, B: same gradients; one job per
  fragment plus one per fragment touching ASN3.
  """
  single = single_toy()
  rm_s = manager(single, expansion=False, clustering=True, maxnum=1)
  g_s = gradients(rm_s, single)
  assert rm_s.n_jobs == rm_s.fragment_manager.n_fragments
  dup, original = duplicate(model=single, labels=["A", "B"],
    residue_indices=[2])
  dup = tst_46.process(dup)
  assert conformers.conformer_labels(dup.get_hierarchy()) == ["A", "B"]
  rm_d = manager(dup, expansion=False, clustering=True, maxnum=1)
  g_d = gradients(rm_d, dup)
  fm = rm_d.fragment_manager
  n_touching = len(touching(fm, residue_selection(dup, 3)))
  d = max_diff(g_d, g_s.select(original))
  print("duplicate: %d fragments (single %d), %d touching ASN3, %d jobs, "
    "max |g_dup - g_single| = %.2e kcal/mol/A" % (fm.n_fragments,
    rm_s.fragment_manager.n_fragments, n_touching, rm_d.n_jobs, d))
  assert fm.n_fragments == rm_s.fragment_manager.n_fragments
  assert 0 < n_touching < fm.n_fragments
  assert rm_d.n_jobs == fm.n_fragments + n_touching
  assert flex.max(g_s.norms()) > 1 # not vacuous
  assert d < 1.e-4, d

def rmsd(m1, m2):
  x1, x2 = m1.get_sites_cart(), m2.get_sites_cart()
  assert x1.size() == x2.size()
  return x1.rms_difference(x2)

def exercise_opt(model):
  """
  qr.refine mode=opt with MOPAC on super_toy.pdb, whole model (b) and
  clustered (c, maxnum=1): the same final model; job folders ase/<i>_A and
  ase/<i>_B for the fragments touching ASN3, ase/<i> for the others.
  Two optimization cycles (number_of_micro_cycles; max_iterations_refine
  alone does not bound mode=opt). COSMO (the qr.refine default) in both.
  """
  common = ["super_toy.pdb", "mode=opt", "restraints=qm",
    "quantum.engine_name=mopac", "max_iterations_refine=3",
    "number_of_micro_cycles=2", "parallel.nproc=1"]
  r_b = run_qrefine_in(folder="opt_b", prefix="opt_b",
    args=common + ["clustering=False"])
  r_c = run_qrefine_in(folder="opt_c", prefix="opt_c",
    args=common + ["clustering=True", "maxnum_residues_in_cluster=1"])
  # expected job folders from the fragments of the same clustering (no
  # engine call, no files)
  params = mopac_params(expansion=False, clustering=True, maxnum=1)
  fm = refine.create_restraints_manager(params=params, model=model,
    working_folder=None, log=null_out()).fragment_manager
  t = touching(fm, residue_selection(model, 3))
  expected = set()
  for i in range(fm.n_fragments):
    if(i in t): expected |= set(["%d_A" % i, "%d_B" % i])
    else:       expected.add("%d" % i)
  ase = os.path.join(r_c.folder, "ase")
  found = set([n for n in os.listdir(ase)
    if os.path.isdir(os.path.join(ase, n))])
  for m in [r_b.model, r_c.model]:
    assert m is not None
    assert m.get_number_of_atoms() == model.get_number_of_atoms()
  d = rmsd(r_b.model, r_c.model)
  d0 = rmsd(r_b.model, model)
  print("opt: rmsd(whole, clustered) = %.4f A (whole moved %.4f A); "
    "%d fragments, %d touching ASN3, job folders %s" % (d, d0,
    fm.n_fragments, len(t), " ".join(sorted(found))))
  assert found == expected, (sorted(found), sorted(expected))
  assert d < 0.02, d

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  """
  MOPAC (PM7) on a peptide with one residue in two side-chain conformers
  (super toy): the clustered conformers route (per-fragment variants)
  against the whole-model and duplicate-conformer references, and qr.refine
  mode=opt end to end.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  model = super_toy()
  for exercise, args in [
      (exercise_one_cluster_limit, [model]),
      (exercise_duplicate,         []),
      (exercise_opt,               [model])]:
    t0 = time.time()
    exercise(*args)
    print("%s: %.1f s" % (exercise.__name__, time.time()-t0))

if(__name__ == "__main__"):
  run()
  print("OK")
