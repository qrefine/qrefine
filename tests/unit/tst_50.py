from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from scitbx import matrix
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from qrefine.tests.unit import run_qrefine_in, run_fmodel
from qrefine.tests.unit.altloc_scenarios import base_model, residue, \
  side_chain_conformers, label_atoms, write, shaken, r_work, labels_and_occ, \
  assert_scattering_neutral, _rotated

REDUCE2 = os.path.join(os.environ["LIBTBX_BUILD"], "bin", "mmtbx.reduce2")

def build(single_source=True):
  """
  Residue 3: CG onward A/B (chi1 0/120), CB with HB2, HB3 A only (0.46).
  Residue 4: CB onward A/B (chi2 rotation, CB-B equals CB-A), N, CA, C, O, H,
  HA A only (0.52). Without single_source the same atoms stay blank, with the
  occupancies 0.46 / 0.52 (the hand-completed model).
  """
  m = base_model()
  rg3, rg4 = residue(m, 3), residue(m, 4)
  side_chain_conformers(rg3, ["A","B"], [0.46,0.54], [0,120])
  side_chain_conformers(rg4, ["A","B"], [0.52,0.48], [0,0])
  cb = matrix.col([a for a in rg4.atoms() if a.name.strip()=="CB"][0].xyz)
  cg = matrix.col([a for a in rg4.atoms() if a.name.strip()=="CG"][0].xyz)
  for ag in rg4.atom_groups():
    if(ag.altloc == "B"):
      for a in ag.atoms():
        a.xyz = _rotated(a.xyz, cb, cg - cb, -120)
  names_3, names_4 = ["CB","HB2","HB3"], ["N","CA","C","O","H","HA"]
  if(single_source):
    label_atoms(rg3, names_3, "A", 0.46)
    label_atoms(rg4, names_4, "A", 0.52)
  else:
    for rg, names, occ in [(rg3, names_3, 0.46), (rg4, names_4, 0.52)]:
      for a in rg.atoms():
        if(a.name.strip() in names): a.occ = occ
  return m

def coordinates(model):
  """
  {(resseq, name, altloc): xyz} of the non-hydrogen atoms.
  """
  result = {}
  for a in model.get_hierarchy().atoms():
    if(a.element.strip().upper() in ("H","D")): continue
    ag = a.parent()
    result[(ag.parent().resseq_as_int(), a.name.strip(), ag.altloc)] = a.xyz
  return result

def compare_occupancies(lo_1, lo_2, tolerance):
  assert sorted(lo_1.keys()) == sorted(lo_2.keys())
  for k, occ in lo_1.items():
    assert approx_equal(lo_2[k], occ, tolerance), k

def hydrogens(model, resseq):
  return sorted([(a.parent().altloc, a.name.strip())
    for a in model.get_hierarchy().atoms()
    if a.parent().parent().resseq_as_int() == resseq and
       a.element.strip().upper() in ("H","D")])

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  1PWM pattern (LYS 85, ASN 129): atoms that exist under one label only
  while the rest of the residue has two. Residue 3 (ASN): CB with HB2, HB3
  under A only (0.46), CG, OD1, ND2 and their H under A (0.46) and B (0.54,
  chi1 + 120). Residue 4 (GLN): N, CA, C, O, H, HA under A only (0.52), CB
  onward under A (0.52) and B (0.48, chi2 rotation, so CB-B equals CB-A).
  Completion shares the single-source atoms (blank, occupancies kept) because
  each is bonded to the other conformer.

  Asserted: complete_altlocs shares exactly CB, HB2, HB3 of residue 3 and N,
  CA, C, O, H, HA of residue 4, all from A, and copies nothing; the number of
  atoms and the structure factors are unchanged; the result equals the model
  made by hand (same atoms, labels, occupancies, coordinates);
  mmtbx.reduce2 qmready=true writes the same heavy atoms, reports the shared
  atoms in its log and keeps HB2, HB3 of residue 3 and H, HA of residue 4;
  qr.refine on the uncompleted model stops with a message pointing to
  reduce2; on the completed shaken model qr.refine keeps labels and
  occupancies and lowers R-work to below 0.6 of the start value.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  m_in, m_hand = build(), build(single_source=False)
  m_done = build()
  r = m_done.complete_altlocs()
  got = sorted([(" ".join(rid.split()), name.strip(), src.strip())
    for rid, name, src in r.shared])
  expected = sorted(
    [("ASN A 3", n, "A") for n in ["CB","HB2","HB3"]] +
    [("GLN A 4", n, "A") for n in ["N","CA","C","O","H","HA"]])
  assert len(r.copied) == 0 and len(r.truncated) == 0 and r.changed
  assert got == expected, got
  assert m_done.get_hierarchy().atoms_size() == \
    m_in.get_hierarchy().atoms_size()
  assert_scattering_neutral(m_in, m_done, copies_expected=False)
  # equals the hand-made model
  assert m_done.get_hierarchy().atoms_size() == \
    m_hand.get_hierarchy().atoms_size()
  lo_done, lo_hand = labels_and_occ(m_done), labels_and_occ(m_hand)
  compare_occupancies(lo_hand, lo_done, 1.e-6)
  c_done, c_hand = coordinates(m_done), coordinates(m_hand)
  d = max([max([abs(x-y) for x, y in zip(c_done[k], c_hand[k])])
    for k in c_hand.keys()])
  assert d < 1.e-6, d
  # reduce2
  in_pdb = "%s_in.pdb"%prefix
  out_pdb = "%s_reduce2.pdb"%prefix
  log_file = "%s_reduce2.log"%prefix
  write(m_in, in_pdb)
  easy_run.fully_buffered("%s %s qmready=True output.filename=%s > %s 2>&1"%(
    REDUCE2, in_pdb, out_pdb, log_file)).raise_if_errors()
  log = open(log_file).read()
  heading = ("Altloc completion (atoms borrowed between conformers of a "
    "residue):")
  assert heading in log, log
  lines = [" ".join(l.split()) for l in log.splitlines()]
  for name in ["CB", "N"]:
    rid = "ASN A 3" if name == "CB" else "GLN A 4"
    line = "%s %s: altloc A -> shared (needed by B)"%(rid, name)
    assert line in lines, line
  m_r2 = mmtbx.model.manager(
    model_input=iotbx.pdb.input(file_name=out_pdb), log=null_out())
  compare_occupancies(lo_done, labels_and_occ(m_r2), 1.e-3)
  h3, h4 = hydrogens(m_r2, 3), hydrogens(m_r2, 4)
  for name in ["HB2", "HB3"]:
    assert [x for x in h3 if x[1] == name], (name, h3)
  for name in ["H", "HA"]:
    assert [x for x in h4 if x[1] == name], (name, h4)
  # qr.refine on the uncompleted model
  r_fm = run_fmodel(prefix = prefix, args=[in_pdb, "high_res=1.2"])
  r_sorry = run_qrefine_in(prefix = prefix+"_sorry", folder="sorry",
    args=[in_pdb, r_fm.mtz, "restraints=cctbx", "clustering=False"],
    sorry_expected=True)
  assert "mmtbx.reduce2 qmready=true" in r_sorry.std_out_text
  # refinement of the completed shaken model
  m_poor = shaken(m_done)
  poor = "%s_poor.pdb"%prefix
  write(m_poor, poor)
  r_qr = run_qrefine_in(prefix = prefix, args=[poor, r_fm.mtz,
    "restraints=cctbx", "clustering=False", "stpmax=0.3"])
  compare_occupancies(lo_done, labels_and_occ(r_qr.model), 1.e-3)
  r_start = r_work(r_fm.f_obs, r_fm.flags, m_poor)
  r_final = r_work(r_fm.f_obs, r_fm.flags, r_qr.model)
  print("R-work: shaken %.4f refined %.4f" % (r_start, r_final))
  assert r_final < 0.6 * r_start, (r_start, r_final)

if(__name__ == "__main__"):
  run()
  print("OK")
