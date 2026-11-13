from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import iotbx.pdb
from scitbx.array_family import flex
from qrefine import conformers
from qrefine.tests.unit import tst_46
from qrefine.tests.unit.tst_31 import duplicate
from qrefine.tests.unit.tst_65 import clustered, gradients, max_diff, \
  contact_model, data_files

def single_model():
  """
  h_altconf_complete without altlocs, in a box with no images within the
  super sphere (contact 30 A).
  """
  pdb_inp = iotbx.pdb.input(file_name=os.path.join(data_files,
    "h_altconf_complete.pdb"))
  ph = pdb_inp.construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  return contact_model(ph, contact=30., name="single")

def cluster_sets(fm, original=None):
  """
  Clusters as sets of (original) atom indices.
  """
  result = []
  for s in fm.fragment_cluster_selections:
    idx = list(s.iselection())
    if(original is not None): idx = [original[i] for i in idx]
    result.append(frozenset(idx))
  return sorted(result, key=lambda x: sorted(x))

def exercise(labels, residue=2):
  """
  Duplicate-conformer identity through the clustered route: one residue of
  a single-conformer model copied into identical conformers (labels; then
  with other occupancies, which are not used) gives the single-conformer
  gradients, the same clusters, and one more job per extra label for every
  fragment that touches the residue.
  """
  model = single_model()
  rm_s = clustered(model, "conformers", maxnum=3)
  g_single = gradients(rm_s, model)
  fm_s = rm_s.fragment_manager
  n_frag_single = fm_s.n_fragments
  assert rm_s.n_jobs == n_frag_single
  # residue atoms in the single model
  rg = list(model.get_hierarchy().residue_groups())[residue]
  res_atoms = set(rg.atoms().extract_i_seq())
  sets_single = cluster_sets(fm_s)
  # same model with duplicated residue, then with altered occupancies
  for altered in [False, True]:
    dup_model, original = duplicate(model=model, labels=labels,
      residue_indices=[residue])
    if(altered):
      ph = dup_model.get_hierarchy()
      for ag in ph.atom_groups():
        if(ag.altloc == labels[0]): ag.atoms().set_occ(
          flex.double(ag.atoms_size(), 0.9))
        elif(ag.altloc == labels[1]): ag.atoms().set_occ(
          flex.double(ag.atoms_size(), 0.1))
        elif(ag.altloc != ""): ag.atoms().set_occ(
          flex.double(ag.atoms_size(), 0.0))
    dup_model = tst_46.process(dup_model)
    assert conformers.conformer_labels(dup_model.get_hierarchy()) == labels
    rm = clustered(dup_model, "conformers", maxnum=3)
    g = gradients(rm, dup_model)
    fm = rm.fragment_manager
    d = max_diff(g, g_single.select(original))
    # Atoms of the residue in the new model (all conformers), no images in
    # the super sphere: super atoms are the model atoms
    assert fm.pdb_hierarchy_super.atoms_size() == \
      dup_model.get_number_of_atoms()
    res_new = flex.bool([original[i] in res_atoms
      for i in range(original.size())])
    n_touching = 0
    for s in fm.fragment_full_super_selections:
      if((s & res_new).count(True) > 0): n_touching += 1
    n_expected = fm.n_fragments + (len(labels)-1)*n_touching
    print("labels %s occ %s: n_fragments %d (single %d), n_touching %d, "
      "n_jobs %d (expected %d), max |g - g_single| %.2e" % (
      "".join(labels), "0.9/0.1" if altered else "1/n", fm.n_fragments,
      n_frag_single, n_touching, rm.n_jobs, n_expected, d))
    assert d < 1.e-9, d
    assert fm.n_fragments == n_frag_single
    assert cluster_sets(fm, original) == sets_single
    assert n_touching > 0
    assert rm.n_jobs == n_expected, (rm.n_jobs, n_expected)
    # Number of conformers an atom belongs to: 1 for the copies of the
    # residue (each in its own conformer), len(labels) for blank atoms
    m = flex.double([1. if r else float(len(labels)) for r in res_new])
    assert rm.multiplicity.all_eq(m)
    assert rm.multiplicity.all_eq(conformers.multiplicity(
      dup_model.get_hierarchy()))

def run(prefix="qrefine_"+os.path.basename(__file__).replace(".py", "")):
  """
  Duplicate-conformer identity (cctbx) of the clustered conformers route
  with two and three labels.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  for labels in [["A", "B"], ["A", "B", "C"]]:
    exercise(labels)

if (__name__ == "__main__"):
  run()
  print("OK")
