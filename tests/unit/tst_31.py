from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import libtbx.load_env
import iotbx.pdb
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from scitbx.array_family import flex
from qrefine import restraints
from qrefine.tests.unit import run_tests

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine,"tests","unit","data_files")

def gradients_from_conformers(ph, crystal_symmetry, calculator):
  """
  Gradients for a model with altlocs from a calculator that only accepts
  single-conformer models (e.g. QM).

  Conformer X = blank atoms + atoms with altloc X. Each conformer is passed
  to calculator(ph, crystal_symmetry) as a complete molecule; it must return
  gradients in the order of ph.atoms(). Per atom i:

    g_i = (1/m_i) * sum over conformers c containing i of g_c,i

  where m_i is the number of conformers atom i belongs to: altloc atoms get
  the gradient of their own conformer, blank atoms the mean over conformers.
  Only gradients are defined (no energy). Occupancies are not used.
  """
  ci = ph.get_conformer_indices().conformer_indices
  n_conformers = flex.max(ci)
  if(n_conformers == 0):
    return calculator(ph, crystal_symmetry)
  # Each conformer must be a complete molecule: every residue with altlocs
  # must have all altloc labels present in the model.
  all_labels = set([a.parent().altloc for a in ph.atoms()]) - set([""])
  for rg in ph.residue_groups():
    labels = set([ag.altloc for ag in rg.atom_groups()]) - set([""])
    if(len(labels) > 0 and labels != all_labels):
      raise Sorry("Residue %s has altlocs %s but model has %s: conformers "
        "would be incomplete." % (rg.id_str(), sorted(labels),
        sorted(all_labels)))
  blank = ci == 0
  g = flex.vec3_double(ci.size(), [0,0,0])
  m = flex.double(ci.size(), 0)
  for c in range(1, n_conformers+1):
    isel = ((ci == c) | blank).iselection()
    g_c = calculator(ph.select(isel), crystal_symmetry)
    assert g_c.size() == isel.size()
    g = g.set_selected(isel, g.select(isel) + g_c)
    m = m.set_selected(isel, m.select(isel) + 1)
  assert m.all_gt(0)
  return g * (1/m)

def cctbx_calculator(ph, crystal_symmetry):
  """
  Stand-in for QM: builds its own restraints from the hierarchy it is given.
  """
  return restraints.get_cctbx_gradients(ph = ph, cs = crystal_symmetry).gradients

def exercise_vs_cctbx():
  """
  cctbx counts each restraint once (blank-blank once, plus every blank-X
  term), which cannot be recovered from whole-conformer gradients. So:
    - altloc atoms: must equal cctbx exactly;
    - blank atoms: m*g = g_cctbx + (m-1)*g_blank_only, which also means
      equal to cctbx for blank atoms that do not interact with altlocs.
  g_blank_only is available from cctbx (not from QM) and is used only here.
  """
  for file_name in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
                    "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:
    pi = iotbx.pdb.input(file_name = os.path.join(qr_unit_tests_data,file_name))
    ph = pi.construct_hierarchy()
    cs = pi.crystal_symmetry()
    g_cctbx = cctbx_calculator(ph, cs)
    g = gradients_from_conformers(
      ph = ph, crystal_symmetry = cs, calculator = cctbx_calculator)
    ci = ph.get_conformer_indices().conformer_indices
    n_conformers = flex.max(ci)
    assert n_conformers == 2
    blank = ci == 0
    altloc = ~blank
    assert approx_equal(g.select(altloc), g_cctbx.select(altloc))
    if(blank.count(True) == 0): continue
    bsel = blank.iselection()
    g_blank_only = cctbx_calculator(ph.select(bsel), cs)
    m = n_conformers
    assert approx_equal(g.select(bsel)*m,
      g_cctbx.select(bsel) + g_blank_only*(m-1))
    # Blank atoms bonded to altloc atoms really differ from cctbx
    assert not approx_equal(g.select(bsel), g_cctbx.select(bsel), out=None)

def exercise_no_altlocs():
  """
  Without altlocs the calculator is called once on the whole model.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  ph = pi.construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  cs = pi.crystal_symmetry()
  g = gradients_from_conformers(
    ph = ph, crystal_symmetry = cs, calculator = cctbx_calculator)
  assert approx_equal(g, cctbx_calculator(ph, cs))

def exercise_duplicate():
  """
  Duplicate a residue into identical conformers A and B: every atom must get
  the gradient of the original single-conformer model.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  ph = pi.construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  cs = pi.crystal_symmetry()
  g_single = cctbx_calculator(ph, cs)
  # Map each atom of the duplicated model to its original atom
  for a in ph.atoms(): a.tmp = a.i_seq
  rg = list(ph.residue_groups())[2]
  ag_a = rg.only_atom_group()
  ag_b = ag_a.detached_copy()
  for a_b, a_a in zip(ag_b.atoms(), ag_a.atoms()): # detached_copy resets tmp
    a_b.tmp = a_a.tmp
  ag_a.altloc = "A"
  ag_b.altloc = "B"
  rg.append_atom_group(ag_b)
  ph.atoms_reset_serial()
  original = flex.size_t([a.tmp for a in ph.atoms()])
  g = gradients_from_conformers(
    ph = ph, crystal_symmetry = cs, calculator = cctbx_calculator)
  assert flex.max(ph.get_conformer_indices().conformer_indices) == 2
  assert approx_equal(g, g_single.select(original))

def exercise_incomplete_conformers():
  """
  Residue with A/B/C in a model where other residues have only A/B: conformer
  C would lack atoms, which QM cannot handle.
  """
  pi = iotbx.pdb.input(
    file_name = os.path.join(qr_unit_tests_data,"h_altconf_complete.pdb"))
  ph = pi.construct_hierarchy()
  cs = pi.crystal_symmetry()
  rgs = [rg for rg in ph.residue_groups() if rg.atom_groups_size() > 1]
  assert len(rgs) > 1
  ag_c = rgs[0].atom_groups()[-1].detached_copy()
  ag_c.altloc = "C"
  rgs[0].append_atom_group(ag_c)
  try:
    gradients_from_conformers(
      ph = ph, crystal_symmetry = cs, calculator = cctbx_calculator)
  except Sorry as e:
    assert str(e).find("conformers would be incomplete") > -1
  else: raise AssertionError("Expected Sorry for incomplete conformers")

def run(prefix):
  """
  Combine gradients from single-conformer calculations (the QM route) for
  models with altlocs. cctbx is used as a stand-in for the QM calculator.
  """
  exercise_vs_cctbx()
  exercise_no_altlocs()
  exercise_duplicate()
  exercise_incomplete_conformers()

if(__name__ == '__main__'):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
