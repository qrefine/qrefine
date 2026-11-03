from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from qrefine import conformers
from qrefine.tests.unit import run_qrefine_in, run_fmodel
from qrefine.tests.unit.tst_31 import gradients_from_conformers, \
  cctbx_calculator
from qrefine.tests.unit.altloc_scenarios import base_model, residue, \
  side_chain_conformers, write, shaken, r_work, labels_and_occ

def build(angles_4, labels_4=("A","B","C"), occupancies_4=(0.4,0.35,0.25)):
  """
  Residues 3 and 5 with A/B chi1 conformers, residue 4 with labels_4.
  """
  m = base_model()
  side_chain_conformers(residue(m,3), ["A","B"], [0.5,0.5], [0,120])
  side_chain_conformers(residue(m,5), ["A","B"], [0.5,0.5], [0,-120])
  side_chain_conformers(residue(m,4), list(labels_4), list(occupancies_4),
    angles_4)
  return m

def keys(model):
  """
  (resseq, name, altloc) per atom, in model order.
  """
  result = []
  for a in model.get_hierarchy().atoms():
    ag = a.parent()
    result.append((ag.parent().resseq_as_int(), a.name.strip(), ag.altloc))
  return result

def exercise_gradients():
  """
  C side chain placed exactly on A: conformer C is a repeat of conformer A.
  """
  m_id = build(angles_4=[0,120,0])
  m_ab = build(angles_4=[0,120], labels_4=["A","B"], occupancies_4=[0.5,0.5])
  g_abc = gradients_from_conformers(m_id)
  g_ab  = gradients_from_conformers(m_ab)
  part_a = conformers.split(m_id)[0]
  assert part_a.label == "A"
  g_a = flex.vec3_double(m_id.get_hierarchy().atoms_size(), (0,0,0))
  g_a = g_a.set_selected(part_a.iselection, cctbx_calculator(
    ph=part_a.model.get_hierarchy(),
    crystal_symmetry=part_a.model.crystal_symmetry()))
  k_id = keys(m_id)
  i_ab = dict([(k, i) for i, k in enumerate(keys(m_ab))])
  n_labelled, n_blank, d_blank = 0, 0, 0
  for i, k in enumerate(k_id):
    if(k[2] == "C"): continue
    j = i_ab[k]
    if(k[2] == ""):
      n_blank += 1
      d_blank = max(d_blank, max([abs(x-y) for x, y in zip(g_abc[i], g_ab[j])]))
      assert approx_equal(
        [3*x for x in g_abc[i]],
        [2*x+y for x, y in zip(g_ab[j], g_a[i])], 1.e-6)
    else:
      n_labelled += 1
      assert approx_equal(g_abc[i], g_ab[j], 1.e-6)
  assert n_labelled > 0 and n_blank > 0
  # The test is not vacuous: blank gradients differ between A/B/C and A/B.
  assert d_blank > 1.e-3, d_blank

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  1US0 pattern (ARG 40, SER 282): residues with A/B side chains in a model
  where one residue also has a C side chain. Residues 3 and 5 (ASN, GLN) get
  A/B chi1 conformers (0.5/0.5), residue 4 (GLN) gets A/B/C
  (0.4/0.35/0.25). The model is complete in its own conformers, so
  completion changes nothing and reduce2 would not touch it; qr.refine
  builds conformer C from conformer A of residues 3 and 5 (first own label).

  Part 1, in process with cctbx restraints: with the C side chain placed
  exactly on A, every labelled atom gets the gradient of the A/B-only model
  and blank atoms obey 3 g_ABC = 2 g_AB + g_A (g_A: gradient of the
  single-conformer A model).

  Part 2: qr.refine restraints=cctbx on a shaken copy against fmodel data
  from the model: it runs, prints the fallback line, keeps every label and
  occupancy and lowers R-work.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  m = build(angles_4=[0,120,-120])
  assert not m.complete_altlocs(check_only=True).changed
  assert conformers.conformer_labels(m.get_hierarchy()) == ["A","B","C"]
  # Part 1
  exercise_gradients()
  # Part 2
  pdb_good = "%s_good.pdb"%prefix
  write(m, pdb_good)
  m_poor = shaken(m)
  pdb_poor = "%s_poor.pdb"%prefix
  write(m_poor, pdb_poor)
  r_fm = run_fmodel(prefix = prefix, args=[pdb_good, "high_res=1.2"])
  r_qr = run_qrefine_in(prefix = prefix, args=[pdb_poor, r_fm.mtz,
    "restraints=cctbx", "clustering=False", "stpmax=0.3"])
  line = "2 residue groups without label C take their first own label"
  assert r_qr.std_out_text.find(line) > -1
  lo_in, lo_out = labels_and_occ(m), labels_and_occ(r_qr.model)
  assert sorted(lo_out.keys()) == sorted(lo_in.keys())
  for k, occ in lo_in.items():
    assert approx_equal(lo_out[k], occ, 1.e-3), k
  r_start = r_work(r_fm.f_obs, r_fm.flags, m_poor)
  r_final = r_work(r_fm.f_obs, r_fm.flags, r_qr.model)
  print("R-work: shaken %.4f refined %.4f" % (r_start, r_final))
  assert r_final < 0.5 * r_start, (r_start, r_final)

if(__name__ == "__main__"):
  run()
  print("OK")
