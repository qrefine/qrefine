from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry, null_out
from scitbx import matrix
import iotbx.pdb
from iotbx.pdb import altloc_completion
import mmtbx.model
from qrefine import conformers
from qrefine.tests.unit import run_qrefine_in, run_fmodel
from qrefine.tests.unit.altloc_scenarios import base_model, residue, \
  side_chain_conformers, backbone_conformers, write, shaken, r_work, \
  labels_and_occ, assert_scattering_neutral, _rotated

REDUCE2 = os.path.join(os.environ["LIBTBX_BUILD"], "bin", "mmtbx.reduce2")

BACKBONE_CB = ["N","CA","C","O","H","HA","CB","HB2","HB3"]
SIDE_CHAIN  = ["CG","OD1","ND2","HD21","HD22"]
ROTATION = 20

def _atom(rg, name, altloc):
  found = [a for ag in rg.atom_groups() if ag.altloc == altloc
    for a in ag.atoms() if a.name.strip() == name]
  assert len(found) == 1, (name, altloc)
  return found[0]

def _named(rg, name):
  found = [a for a in rg.atoms() if a.name.strip() == name]
  assert len(found) == 1, name
  return found[0]

def dist(x, y):
  return (matrix.col(x.xyz) - matrix.col(y.xyz)).length()

def build(occ_d=0.20):
  """
  (model, d_degrees, distances). The D side chain angle is searched: CG-D must
  be bonded to both CB-A and CB-B, farthest from CG-A, CG-B and CG-C.
  """
  m = base_model()
  # D side chain is a symmetry-mate clash with water 8 (reduce2 qmready=true
  # refuses bad clashes), so the waters are left out.
  for chain in m.get_hierarchy().only_model().chains():
    for rg in list(chain.residue_groups()):
      if(rg.atom_groups()[0].resname == "HOH"): chain.remove_residue_group(rg)
  m.get_hierarchy().atoms().reset_i_seq()
  m.get_hierarchy().atoms_reset_serial()
  backbone_conformers(m, 3, "B", 0.19, ROTATION, names=BACKBONE_CB)
  rg = residue(m, 3)
  side_chain_conformers(rg, ["A","B","C","D"], [0.43,0.10,0.18,occ_d],
    [0,0,120,0])
  ca, cb = [matrix.col(_atom(rg, n, "A").xyz) for n in ["CA","CB"]]
  c_prev = matrix.col(_named(residue(m, 2), "C").xyz)
  n_next = matrix.col(_named(residue(m, 4), "N").xyz)
  for a in rg.atoms():
    if(a.parent().altloc == "B" and a.name.strip() in SIDE_CHAIN):
      a.xyz = _rotated(a.xyz, c_prev, n_next - c_prev, ROTATION)
  cb_a, cb_b = _atom(rg, "CB", "A"), _atom(rg, "CB", "B")
  cg = {l: _atom(rg, "CG", l) for l in "ABCD"}
  best = None
  for d in range(0, 360, 10):
    trial = cg["A"].detached_copy()
    trial.xyz = _rotated(cg["A"].xyz, cb, cb - ca, d)
    if(not (altloc_completion.bonded(trial, cb_a) and
            altloc_completion.bonded(trial, cb_b))): continue
    far = min([dist(trial, cg[l]) for l in "ABC"])
    if(best is None or far > best[0]): best = (far, d)
  assert best is not None
  d = best[1]
  for a in rg.atoms():
    if(a.parent().altloc == "D" and a.name.strip() in SIDE_CHAIN):
      a.xyz = _rotated(_atom(rg, a.name.strip(), "A").xyz, cb, cb - ca, d)
  # explicit occupancies per name set
  for a in rg.atoms():
    l, n = a.parent().altloc, a.name.strip()
    if(n in BACKBONE_CB): a.occ = {"A":0.81, "B":0.19}[l]
    else: a.occ = {"A":0.43, "B":0.10, "C":0.18, "D":occ_d}[l]
  m.get_hierarchy().atoms().reset_i_seq()
  m.get_hierarchy().atoms_reset_serial()
  od1 = {l: _atom(rg, "OD1", l) for l in "AD"}
  dists = {"CG-D..CB-A": dist(cg["D"], cb_a),
           "CG-D..CB-B": dist(cg["D"], cb_b),
           "CG-C..CB-A": dist(cg["C"], cb_a),
           "CG-C..CB-B": dist(cg["C"], cb_b),
           "CG-D..OD1-D": dist(cg["D"], od1["D"]),
           "CG-A..OD1-A": dist(cg["A"], od1["A"]),
           "CG-D..nearest CG": min([dist(cg["D"], cg[l]) for l in "ABC"])}
  return m, d, dists

def sums(model):
  s = {}
  for a in residue(model, 3).atoms():
    if(a.element.strip().upper() in ("H","D")): continue
    s[a.name.strip()] = s.get(a.name.strip(), 0) + a.occ
  return s

def check_sums(model, side_chain_sum):
  """
  Backbone and CB names sum to 1; CG, OD1, ND2 to 0.43 + 0.10 + 0.18 + 0.20.
  """
  for k, v in sums(model).items():
    expected = side_chain_sum if k in SIDE_CHAIN else 1.0
    assert approx_equal(v, expected, 1.e-6), (k, v, expected)

def atoms_size(m):
  return m.get_hierarchy().atoms_size()

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  1R6J pattern (SER 205): a side chain with four conformers on a backbone
  with two. Residue 3 (ASN): N, CA, C, O, H, HA, CB, HB2, HB3 under A (0.81)
  and B (0.19, rotated 20 degrees about C(2)...N(4)); CG, OD1, ND2 with H
  under A (0.43, original), B (0.10, rotated with the B backbone), C (0.18,
  chi1 + 120 on A) and D (0.20, chi1 rotation about the CA-CB axis of A,
  searched in steps of 10 degrees so that CG-D is bonded to both CB-A and
  CB-B and is farthest from CG-A/B/C; all D side-chain atoms are the A ones
  rotated by that angle).
  Completion copies the backbone and CB (nine names) from A into C and D:
  the bond test leaves only CB-A for C, the occupancy budget picks CB-A over
  CB-B for D (B has 0.19 for 0.20 + 0.10).

  Asserted: the precondition, on the built model, that CG-D is bonded
  (< 1.92 A) to both CB-A and CB-B and that CG-D..OD1-D equals CG-A..OD1-A;
  per-name occupancy sums before and after are 1 for backbone and CB and
  0.91 for CG, OD1, ND2 (0.43 + 0.10 + 0.18 + 0.20); complete_altlocs
  copies exactly the nine names to C and to D, all from A, and shares and
  truncates nothing, no inconsistency; occupancies of CB and of every heavy
  backbone atom become A 0.43, B 0.19, C 0.18, D 0.20 and the sums are
  unchanged; 18 atoms are added and the structure factors are unchanged
  (assert_scattering_neutral); conformers.split gives four parts of equal
  atom count and no conformer needs the fallback; mmtbx.reduce2 qmready=true
  on the input writes the same heavy atoms and occupancies and logs "copied
  A -> C" and "copied A -> D" lines naming the nine names of residue 3;
  qr.refine on the shaken completed model keeps labels and occupancies and
  lowers R-work to below 0.6 of the start value. Ambiguous variant (D side
  chain 0.45): complete_altlocs(check_only=True) raises Sorry naming
  conformer D, ASN A 3, CB and both candidates altloc A and altloc B, and
  leaves the hierarchy unchanged.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  m_in, d, dists = build()
  print("D chi1 angle: %d; distances: %s"%(d, dists))
  assert dists["CG-D..CB-A"] < 1.92 and dists["CG-D..CB-B"] < 1.92, dists
  assert approx_equal(dists["CG-D..OD1-D"], dists["CG-A..OD1-A"], 1.e-6), dists
  check_sums(m_in, side_chain_sum = 0.91)
  in_pdb = "%s_in.pdb"%prefix
  write(m_in, in_pdb)
  m_done, _, _ = build()
  r = m_done.complete_altlocs()
  print("copied:", sorted(r.copied))
  got = sorted([(" ".join(rid.split()), n.strip(), s.strip(), t.strip())
    for rid, n, s, t, o in r.copied])
  expected = sorted([("ASN A 3", n, "A", t) for n in BACKBONE_CB
    for t in ["C","D"]])
  assert got == expected, got
  assert len(r.shared) == 0 and len(r.truncated) == 0 and \
    len(r.inconsistent) == 0 and r.changed
  assert atoms_size(m_done) == atoms_size(m_in) + 18
  assert_scattering_neutral(m_in, m_done)
  lo = labels_and_occ(m_done)
  for n in BACKBONE_CB:
    if(n.startswith("H")): continue
    for l, o in [("A",0.43), ("B",0.19), ("C",0.18), ("D",0.20)]:
      assert approx_equal(lo[(3, n, l)], o, 1.e-6), (n, l, lo[(3, n, l)])
  check_sums(m_done, side_chain_sum = 0.91)
  parts = conformers.split(m_done)
  assert len(parts) == 4, len(parts)
  sizes = [len(p.model.get_hierarchy().atoms()) for p in parts]
  assert len(set(sizes)) == 1, sizes
  fb = [s.n_fallback for s in
    conformers.conformer_selections(m_done.get_hierarchy())]
  assert fb == [0, 0, 0, 0], fb
  # reduce2
  out_pdb = "%s_reduce2.pdb"%prefix
  log_file = "%s_reduce2.log"%prefix
  easy_run.fully_buffered("%s %s qmready=True output.filename=%s > %s 2>&1"%(
    REDUCE2, in_pdb, out_pdb, log_file)).raise_if_errors()
  log = open(log_file).read()
  heading = \
    "Altloc completion (atoms borrowed between conformers of a residue):"
  assert heading in log, log
  lines = [" ".join(l.split()) for l in log.splitlines()]
  for t in ["C", "D"]:
    assert [l for l in lines if l.startswith("ASN A 3 ") and
      l.endswith(": copied A -> %s"%t) and
      sorted(l.split(": ")[0][len("ASN A 3 "):].split(", ")) ==
        sorted(BACKBONE_CB)], (t, log)
  print("\n".join([l for l in log.splitlines() if "copied" in l or
    heading in l]))
  m_r2 = mmtbx.model.manager(
    model_input=iotbx.pdb.input(file_name=out_pdb), log=null_out())
  lo_r2 = labels_and_occ(m_r2)
  assert sorted(lo_r2.keys()) == sorted(lo.keys())
  for k, o in lo.items(): assert approx_equal(lo_r2[k], o, 1.e-3), k
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
  # ambiguous variant
  m_amb, _, _ = build(occ_d=0.45)
  before = m_amb.get_hierarchy().as_pdb_string()
  try:
    m_amb.complete_altlocs(check_only=True)
  except Sorry as e:
    text = str(e)
  else:
    raise AssertionError("Sorry expected")
  print(text)
  for s in ["Cannot complete conformer D", " ASN A   3 ", " CB ", "altloc A",
            "altloc B"]:
    assert s in text, (s, text)
  assert m_amb.get_hierarchy().as_pdb_string() == before

if(__name__ == "__main__"):
  run()
  print("OK")
