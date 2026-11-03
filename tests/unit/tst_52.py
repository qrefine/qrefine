from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx.test_utils import approx_equal
from scitbx import matrix
from qrefine import conformers
from qrefine.tests.unit import run_qrefine_in, run_fmodel
from qrefine.tests.unit.altloc_scenarios import base_model, residue, \
  side_chain_conformers, backbone_conformers, label_atoms, write, shaken, \
  r_work, labels_and_occ, assert_scattering_neutral

BACKBONE_NAMES = ["N","CA","C","O","H","HA"]
CB_NAMES       = ["CB","HB2","HB3"]
SIDE_CHAIN     = ["CG","OD1","ND2","HD21","HD22"]
BOND_LIMIT     = 1.92 # C-C: covalent radii + 0.4 A

def _atom(rg, name, altloc):
  found = [a for ag in rg.atom_groups() if ag.altloc == altloc
    for a in ag.atoms() if a.name.strip() == name]
  assert len(found) == 1, (name, altloc)
  return found[0]

def dist(x, y):
  return (matrix.col(x.xyz) - matrix.col(y.xyz)).length()

def build():
  """
  (model, dists). Waters are left out (as in tst_51).
  """
  m = base_model()
  for chain in m.get_hierarchy().only_model().chains():
    for rg in list(chain.residue_groups()):
      if(rg.atom_groups()[0].resname == "HOH"): chain.remove_residue_group(rg)
  m.get_hierarchy().atoms().reset_i_seq()
  m.get_hierarchy().atoms_reset_serial()
  backbone_conformers(m, 3, "B", 0.18, 90)
  rg = residue(m, 3)
  side_chain_conformers(rg, ["A","C"], [0.41, 0.41], [0, 120])
  label_atoms(rg, CB_NAMES, "A", 0.82)
  for a in rg.atoms():
    l, n = a.parent().altloc, a.name.strip()
    if(n in BACKBONE_NAMES): a.occ = {"A":0.82, "B":0.18}[l]
    elif(n in CB_NAMES):     a.occ = 0.82
    else:                    a.occ = 0.41
  m.get_hierarchy().atoms().reset_i_seq()
  m.get_hierarchy().atoms_reset_serial()
  dists = {"CA-B..CB-A": dist(_atom(rg,"CA","B"), _atom(rg,"CB","A")),
           "CB-A..CG-C": dist(_atom(rg,"CB","A"), _atom(rg,"CG","C")),
           "CB-A..CG-A": dist(_atom(rg,"CB","A"), _atom(rg,"CG","A"))}
  return m, dists

def sums(model):
  s = {}
  for a in residue(model, 3).atoms():
    if(a.element.strip().upper() in ("H","D")): continue
    s[a.name.strip()] = s.get(a.name.strip(), 0) + a.occ
  return s

def check_sums(model):
  for k, v in sums(model).items():
    expected = 1.0 if k in BACKBONE_NAMES else 0.82
    assert approx_equal(v, expected, 1.e-6), (k, v, expected)

def atoms_size(m):
  return m.get_hierarchy().atoms_size()

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  1R6J pattern (ASP 195): backbone A/B, CB under A only, side chain beyond CB
  under A/C. Residue 3 (ASN): N, CA, C, O, H, HA under A (0.82) and B (0.18,
  rotated 90 degrees about C(2)...N(4) so that CA-B is more than 1.92 A from
  CB-A); CB, HB2, HB3 under A only (0.82); CG, OD1, ND2 with H under A (0.41)
  and C (0.41, chi1 + 120). Conformer B cannot borrow CB (CB-A is not bonded
  to anything of B) and is truncated at CA; conformer C borrows CB-A
  (anchored by CG-C) and then the A backbone (CA-B is not bonded to CB), each
  a single candidate, copied with 0.41/0.41.

  Asserted: the precondition, on the built model, that CA-B..CB-A is more
  than 1.92 A and CB-A..CG-C less than 1.92 A; per-name occupancy sums before
  and after are 1 for the backbone names and 0.82 for CB and the side chain;
  complete_altlocs reports exactly one truncated entry (residue 3, label B,
  CB, CG, HB2, HB3, HD21, HD22, ND2, OD1; the reason names CB, "A copy", and
  the distance to the nearest atom of B, with its name, computed from the
  model), copies exactly C, CA, CB,
  H, HA, HB2, HB3, N, O to C, all from A and nothing to other targets, shares
  nothing; occupancies of CB become A 0.41, C 0.41 and of the backbone A 0.41,
  B 0.18, C 0.41; the sums are unchanged; 9 atoms are added and the
  structure factors are unchanged (assert_scattering_neutral);
  conformers.check_completed raises nothing (a truncated conformer is not a
  change); three conformers, none needs the fallback, and conformer B has
  exactly 8 atoms fewer than A and C (which are equal); qr.refine with cctbx
  restraints on the shaken completed model succeeds with the truncated
  conformer, keeps labels and occupancies and lowers R-work to below 0.6 of
  the start value. A QM run would need reduce2 capping=True for the
  truncated B conformer.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  m_in, dists = build()
  print("distances:", dists)
  assert dists["CA-B..CB-A"] > BOND_LIMIT, dists
  assert dists["CB-A..CG-C"] < BOND_LIMIT, dists
  check_sums(m_in)
  in_pdb = "%s_in.pdb"%prefix
  write(m_in, in_pdb)
  m_done, _ = build()
  r = m_done.complete_altlocs()
  print("copied:", sorted(r.copied))
  print("truncated:", r.truncated)
  assert len(r.truncated) == 1, r.truncated
  rid, label, names, reason = r.truncated[0]
  assert " ".join(rid.split()) == "ASN A 3" and label == "B", r.truncated
  assert sorted(names) == ["CB","CG","HB2","HB3","HD21","HD22","ND2","OD1"], \
    names
  # the nearest atom of B to CB-A (hydrogens included) is named in the reason
  rg = residue(m_in, 3)
  cb_a = _atom(rg, "CB", "A")
  nearest = min([(dist(a, cb_a), a.name.strip()) for ag in rg.atom_groups()
    if ag.altloc == "B" for a in ag.atoms()])
  print("nearest atom of B to CB-A:", nearest)
  for s in ["CB", "A copy", "%.1f A from %s B"%nearest]:
    assert s in reason, (s, reason)
  got = sorted([n.strip() for _, n, s, t, _ in r.copied
    if t.strip() == "C" and s.strip() == "A"])
  assert got == ["C","CA","CB","H","HA","HB2","HB3","N","O"], got
  assert len(r.copied) == 9, r.copied
  assert len(r.shared) == 0 and len(r.inconsistent) == 0 and r.changed
  assert atoms_size(m_done) == atoms_size(m_in) + 9
  assert_scattering_neutral(m_in, m_done)
  lo = labels_and_occ(m_done)
  expected_occ = {("CB","A"):0.41, ("CB","C"):0.41}
  for n in ["N","CA","C","O"]:
    for l, o in [("A",0.41), ("B",0.18), ("C",0.41)]: expected_occ[(n,l)] = o
  for (n, l), o in expected_occ.items():
    assert approx_equal(lo[(3, n, l)], o, 1.e-6), (n, l, lo[(3, n, l)])
  assert (3, "CB", "B") not in lo
  check_sums(m_done)
  ph_done = m_done.get_hierarchy()
  conformers.check_completed(ph_done)
  sel = conformers.conformer_selections(ph_done)
  assert len(sel) == 3, len(sel)
  fb = [s.n_fallback for s in sel]
  assert fb == [0, 0, 0], fb
  sizes = [len(p.model.get_hierarchy().atoms())
    for p in conformers.split(m_done)]
  print("conformer sizes:", sizes)
  assert len(sizes) == 3 and sizes[0] == sizes[2] and \
    sizes[0] - sizes[1] == 8, sizes
  # refinement
  r_fm = run_fmodel(prefix = prefix, args=[in_pdb, "high_res=1.2"])
  m_poor = shaken(m_done)
  poor = "%s_poor.pdb"%prefix
  write(m_poor, poor)
  r_qr = run_qrefine_in(prefix = prefix, args=[poor, r_fm.mtz,
    "restraints=cctbx", "clustering=False", "stpmax=0.3",
    "number_of_weight_search_cycles=1", "data_weight=10"])
  lo_qr = labels_and_occ(r_qr.model)
  assert sorted(lo_qr.keys()) == sorted(lo.keys())
  for k, o in lo.items(): assert approx_equal(lo_qr[k], o, 1.e-3), k
  r_start = r_work(r_fm.f_obs, r_fm.flags, m_poor)
  r_final = r_work(r_fm.f_obs, r_fm.flags, r_qr.model)
  print("R-work: shaken %.4f refined %.4f" % (r_start, r_final))
  assert r_final < 0.6 * r_start, (r_start, r_final)

if(__name__ == "__main__"):
  run()
  print("OK")
