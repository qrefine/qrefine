# Altloc Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `mmtbx.reduce2 qmready=true` and `qr.refine` accept the partial altlocs of real PDB entries by completing residues from their own conformers (borrowing, scattering-neutral) and by letting qr.refine's conformer models fall back to a residue's first label.

**Architecture:** One pure-hierarchy module in iotbx implements the borrowing rule and edits the hierarchy; mmtbx.model wraps it and reduce2 qmready calls it before its atom-complete check, which stops reporting chain-inherited conformers. qrefine's `conformers.py` builds conformer models with the first-label fallback, counts multiplicity from membership, and qr.refine only checks (Sorry pointing to reduce2), never edits.

**Tech Stack:** cctbx (iotbx.pdb.hierarchy, mmtbx.model, mmtbx.programs.reduce2), qrefine (conformers.py, qr.py), libtbx tests. Dispatchers from `<root>/build/bin/` (`libtbx.python`, `mmtbx.reduce2`, `qr.refine`, `mmtbx.python -m mmtbx.command_line.fmodel`).

**Spec:** `qrefine/docs/superpowers/specs/2026-11-02-altloc-completion-design.md` (approved, D1-D6 settled). Scratch evidence: `<root>/rt/run_1793639952_reduce2/` (classify.log, parent_choice.py, reduce2_list_clashes.py, the three entries and the survey).

## Global Constraints

- cctbx conventions: 2-space indent, 80 columns, imports inside functions where feasible, `Sorry` for user errors, print to `log`, no `import *`, numpydoc docstrings on new public functions. iotbx must not import mmtbx (dependency direction).
- Constants (spec D1): `BOND_TOLERANCE = 0.4` (same value as `mmtbx.hydrogens.capping.GUESS_BOND_TOLERANCE`), `OCCUPANCY_TOLERANCE = 0.05`. Not parameters.
- Bonded: `d < r_a + r_b + 0.4` with covalent radii from `cctbx.eltbx.covalent_radii` (fallback 0.77 as `altloc_reconciliation._covalent_radius`).
- Scattering neutrality (spec 5): no edit changes a site's total occupancy; copies keep xyz, B, uij, segid, hetero.
- No PHIL parameter for completion (D2). An ambiguous borrow is a Sorry (D3). Parent for a residue without label X is its first own label; a broken labelled link is a Sorry (D5).
- Every qrefine test writes only inside its own `qrefine_tst_XX` folder (keep the `os.makedirs` / `os.chdir` lines of `run()`). Tests run from an empty scratch directory outside the checkout.
- Nothing is committed by the implementer: leave changes uncommitted and report the diff; the user commits.
- Run tests with the build's dispatchers, e.g. `<root>/build/bin/libtbx.python <file>`.

## Review Focus

1. A residue group with exactly one own label (1PWM ASP 125, A only): the completion must do nothing and `missing_atoms_by_residue` must not report conformer B for it. Test in Task 1 (`exercise_single_label`) and Task 2 (`exercise_inherited_conformer`).
2. A borrowed heavy atom whose H are under the source label (HB2/HB3 of a CB-A): the H must be borrowed with it, never left under A alone. Test in Task 1 (`exercise_hydrogens_follow`).
3. Hydrogen-only labels (neutron H/D sites, label with no non-H atom): not conformer labels, never a borrowing target, never counted in the fallback. Test in Task 1 (`exercise_hydrogen_only_label`) and Task 4 (existing `exercise_split` H-only part kept).
4. Waters and ligands with partial labels (a water O under A and B, a ligand with one atom under C): the same rule, no amino-acid assumptions. Test in Task 1 (`exercise_water_and_ligand`).
5. Microheterogeneity (PRO under A, SER under B/C, shared N as 1EJG 22): SER never borrows PRO atoms. Test in Task 1 (`exercise_microheterogeneity`).

---

### Task 1: `iotbx.pdb.altloc_completion` and `hierarchy.root.complete_altlocs`

**Files:**
- Create: `cctbx_project/iotbx/pdb/altloc_completion.py`
- Modify: `cctbx_project/iotbx/pdb/hierarchy.py` (add method next to `reconcile_altlocs`, line ~2669)
- Create: `cctbx_project/iotbx/pdb/tst_altloc_completion.py`
- Modify: `cctbx_project/iotbx/run_tests.py` (register after `tst_hierarchy_reconcile_altlocs.py`, line 127)

**Interfaces:**
- Consumes: `iotbx.pdb.altloc_reconciliation._element`, `_covalent_radius` (import them; do not copy).
- Produces:
  - `BOND_TOLERANCE = 0.4`, `OCCUPANCY_TOLERANCE = 0.05`
  - `bonded(atom_a, atom_b) -> bool`
  - `conformer_labels(residue_group) -> list[str]` sorted labels with at least one non-H atom
  - `label_occupancy(residue_group, label) -> float` smallest occupancy among non-H atoms of that label (spec occ_R)
  - `run(hierarchy, check_only=False, log=None) -> libtbx.group_args` with fields `shared` (list of `(rg_id_str, atom_name, source_label)`), `copied` (list of `(rg_id_str, atom_name, source_label, target_label, occupancy)`), `truncated` (list of `(rg_id_str, label, [atom_names], reason)`), `inconsistent` (list of `rg_id_str`), `changed` (bool). Raises `Sorry` on an ambiguous borrow. With `check_only=True` the hierarchy is untouched and `changed` says whether it would be.
  - `report_lines(result) -> list[str]`: the lines of spec section 6, exactly: `"%s %s: altloc %s -> shared (needed by %s)"`, `"%s %s: copied %s -> %s"` (names joined by ", ", one line per (source, target) pair), `"%s %s: %s copy is %.1f A from %s, not borrowed; truncated in %s"`. Empty list when nothing happened.
  - `hierarchy.root.complete_altlocs(self, check_only=False, log=None)` returning `run(...)`'s result; when `log` is given it prints `report_lines`.

- [ ] **Step 1: Write the failing tests** in `tst_altloc_completion.py`, building hierarchies from lines with the `atom_line` helper of `tst_hierarchy_reconcile_altlocs.py` (copy the helper; it takes name, altloc, resname, resseq, xyz, element and writes occupancy 0.50; add an `occ` argument). Use a fake residue `XXX` with carbon atoms on a line 1.5 A apart (`C1` at x=0, `C2` at 1.5, `C3` at 3.0, `C4` at 4.5) unless a real residue is needed.

```python
def exercise_single_source_shared():
  # C1 blank; C2 under A only (occ 0.46); C3 under A (0.46) and B (0.54), B at y=1.5
  r = h.complete_altlocs()
  assert r.shared == [(rg_id, "C2", "A")] and r.copied == [] and r.changed
  a = atom_named(h, "C2"); assert a.parent().altloc == "" and approx_equal(a.occ, 0.46)
  assert altloc_completion.report_lines(r) == ["%s C2: altloc A -> shared (needed by B)" % rg_id]

def exercise_multi_source_copies():
  # SER 205 pattern on XXX: C1 blank; C2 under A (0.81) at x=1.5 and B (0.19) at (1.5, 0.7);
  # C3 under A 0.43 (2.9,0,0), B 0.10 (2.9,0.7,0), C 0.18 (2.2,1.1,0), D 0.20 (2.5,1.3,0)
  # distances chosen so that C3-C is bonded to C2-A only and C3-D to both C2-A and C2-B.
  r = h.complete_altlocs()
  assert [c[:4] for c in sorted(r.copied)] == [(rg_id, "C2", "A", "C"), (rg_id, "C2", "A", "D")]
  assert approx_equal([c[4] for c in sorted(r.copied)], [0.18, 0.20])
  assert occupancies_of(h, "C2") == {"A": 0.43, "B": 0.19, "C": 0.18, "D": 0.20}  # approx_equal
  assert r.inconsistent == []

def exercise_ambiguous_sorry():
  # as above with C3-D at occupancy 0.45: no candidate passes the budget test
  try: h.complete_altlocs()
  except Sorry as e:
    s = str(e); assert "XXX" in s and " C2 " in s and "conformer D" in s and "altloc A" in s and "altloc B" in s
  else: raise AssertionError

def exercise_truncated():
  # ASP 195 pattern: C1 under A (0,0,0) and B (0,3,0); C2 under A only at (1.5,0,0)
  r = h.complete_altlocs()
  assert r.truncated == [(rg_id, "B", ["C2"], "C2: A copy is 3.4 A from C1 B")]  # distance to 1 decimal in the reason
  assert r.shared == [] and r.copied == [] and not r.changed

def exercise_hydrogens_follow():
  # single-source C2-A with H21-A, H22-A bonded to it (1.0 A): all three shared
  assert sorted(n for _, n, _ in r.shared) == ["C2", "H21", "H22"]

def exercise_hydrogen_only_label():
  # labels A/B complete; one extra H under label C only: conformer_labels(rg) == ["A", "B"], nothing happens
  assert altloc_completion.conformer_labels(rg) == ["A", "B"] and not r.changed

def exercise_single_label():
  # residue with atoms under A only plus blank atoms: nothing happens
  assert not r.changed and r.truncated == []

def exercise_water_and_ligand():
  # water O under A and B at different positions: nothing; ligand LIG with C1 blank, C2 A-only, C3 A/B: C2 shared
  ...

def exercise_microheterogeneity():
  # N blank; CA, C under A (PRO) and under B (SER); CB under A (PRO) only: SER B does not borrow CB
  assert not r.changed

def exercise_check_only_and_idempotent():
  r1 = h.complete_altlocs(check_only=True); assert r1.changed and h_unchanged
  r2 = h.complete_altlocs(); r3 = h.complete_altlocs()
  assert r2.shared == r1.shared and not r3.changed and r3.shared == [] and r3.copied == []

def exercise_scattering_neutral():
  # multi-source hierarchy with CRYST1 10 10 10 90 90 90 P 1: f_calc amplitudes before and after agree to 1e-6
  xrs0 = pdb_inp.xray_structure_simple(); ...; assert flex.max(flex.abs(f1.data()-f0.data())) < 1e-6 * flex.max(flex.abs(f0.data()))

def exercise_label_occupancy_and_bonded():
  assert approx_equal(altloc_completion.label_occupancy(rg, "A"), 0.43)
  assert altloc_completion.bonded(c_atom_at_0, c_atom_at_1_9) and not altloc_completion.bonded(c_atom_at_0, c_atom_at_1_93)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `<root>/build/bin/libtbx.python cctbx_project/iotbx/pdb/tst_altloc_completion.py`
Expected: `ImportError` / `AttributeError: 'root' object has no attribute 'complete_altlocs'`.

- [ ] **Step 3: Implement `altloc_completion.py`**

Module docstring: the rule of spec section 4 in prose (as `altloc_reconciliation.py` does). Algorithm per residue group R with `conformer_labels(R)` of size >= 2 (size 0 or 1: skip):

```
for X in labels:
  group_X = atom group of R with altloc X (resname rx)
  in_X = atoms of R with altloc "" or X            # by atom name
  names_elsewhere = names of atoms in atom groups with altloc in labels - {X} and resname == rx, not in in_X
  gaps = names_elsewhere; borrowers_of[(atom)] accumulates labels
  repeat until no acceptance:
    for name in sorted(gaps):
      candidates = [(P, atom)] from labels P != X, resname rx, with that name
      passing = [c for c in candidates if any(bonded(c.atom, b) for b in in_X)]
      if not passing: continue
      if len(passing) > 1:
        passing = [c for c in passing if c.atom.occ >= occ_R(X) + occ_R(P) + sum(occ_R(Y) for Y in borrowers_of[c.atom]) - OCCUPANCY_TOLERANCE]
        if len(passing) != 1: raise Sorry(...)        # message names rg.id_str(), name, X, each candidate: "altloc P (d.dd A from <anchor name>, occ o.oo)"
      accept: borrowers_of[atom].append(X); in_X.append(atom); gaps.remove(name)
  truncated.append((rg.id_str(), X, sorted(gaps), reason)) if gaps     # reason: the nearest candidate of the first gap name and its distance to the nearest in_X atom, 1 decimal
edits (skipped when check_only):
  for atom, borrowers in borrowers_of.items():
    P = atom.parent().altloc
    if set(borrowers) == set(labels) - {P}:  -> move atom to the blank atom group of R (resname rx; create it with rg.insert_atom_group(0, ag) if none); occupancy unchanged; shared.append
    else: weights = {Y: occ_R(Y) for Y in borrowers + [P]}; total = atom.occ
          if abs(sum(weights) - total) > OCCUPANCY_TOLERANCE: inconsistent.append(rg.id_str()) (once per rg)
          for Y in borrowers: copy = atom.detached_copy(); copy.set_occ(total*weights[Y]/sum(weights)); group_Y.append_atom(copy); copied.append(...)
          atom.set_occ(total*weights[P]/sum(weights))
  hierarchy.atoms().reset_i_seq() when anything changed
```

`occ_R` is `label_occupancy` evaluated on the input occupancies (before any edit), computed once per residue group. H atoms take part like any atom (the bond test with covalent radius of H 0.31 accepts 1.0 A X-H). Keep `check_only` by collecting edits first and applying them in a second phase.

- [ ] **Step 4: Add `complete_altlocs` to `hierarchy.root`** right after `reconcile_altlocs` (docstring in the same style), delegating to `altloc_completion.run` and printing `report_lines` to `log`.

- [ ] **Step 5: Run the test to verify it passes**

Run: `<root>/build/bin/libtbx.python cctbx_project/iotbx/pdb/tst_altloc_completion.py`
Expected: `OK`.

- [ ] **Step 6: Register and run neighbours**

Add `"$D/pdb/tst_altloc_completion.py"` to `iotbx/run_tests.py`. Run `tst_hierarchy_reconcile_altlocs.py` and `tst_hierarchy.py` to check nothing in the package broke. Run `libtbx.find_clutter` on the new files.

---

### Task 2: `mmtbx.model.manager.complete_altlocs` and the missing-atoms report

**Files:**
- Modify: `cctbx_project/mmtbx/model/model.py` (`_missing_atoms_records` lines ~2440-2500 in the loop building `key`; new method next to `add_terminal_oxt`, line ~4698)
- Modify: `cctbx_project/mmtbx/regression/model/tst_model_atom_completeness.py` (add three exercises, call them from `run()`)

**Interfaces:**
- Consumes: Task 1 `hierarchy.root.complete_altlocs`.
- Produces: `mmtbx.model.manager.complete_altlocs(self, check_only=False, log=None)` returning Task 1's `group_args`; after a change it calls `self._pdb_hierarchy.atoms().reset_serial()` and `self.reset_after_changing_hierarchy()` (as `add_terminal_oxt` does). `missing_atoms_by_residue()` keys: a residue group with own labels is reported only under its own labels; a residue group without labels is reported once under its blank key (`ag.id_str()`), even in a chain with altlocs.

- [ ] **Step 1: Write the failing tests** in `tst_model_atom_completeness.py`, built from `tripeptide` (GLN 1, ASN 2, TYR 3) with the existing `remove` helper and the line-editing idiom of `exercise_alt_confs`:

```python
def exercise_inherited_conformer():
  # TYR OH under A and B (complete), ASN ND2 under A, B and C: no entry for TYR in conformer C
  m = get_model(...)
  assert m.is_atom_complete(), m.missing_atoms_by_residue()
  # the same, with TYR OH under A only: still reported only for B
  assert m.missing_atoms_by_residue() == {"BTYR A   3 ": ["OH"]}

def exercise_plain_residue_in_altloc_chain():
  # GLN 1 without CD, OE1, NE2 (blank residue) while TYR OH is under A and B
  assert m.missing_atoms_by_residue() == {" GLN A   1 ": ["CD", "NE2", "OE1"]}
  assert list(m.get_missing_atoms().keys()) == [" GLN A   1 "]

def exercise_complete_altlocs():
  # TYR CZ under A and B, OH under A only (the incomplete case of exercise_alt_confs)
  assert not m.is_atom_complete()
  r = m.complete_altlocs()
  assert r.shared == [(" TYR A   3 ", "OH", "A")]
  assert m.is_atom_complete()
  assert m.get_xray_structure().scatterers().size() == m.get_hierarchy().atoms_size()
  assert m.get_hierarchy().atoms().extract_serial()[-1].strip() == str(m.get_hierarchy().atoms_size())  # reset_serial applied
```

- [ ] **Step 2: Run to verify failure**

Run: `<root>/build/bin/libtbx.python cctbx_project/mmtbx/regression/model/tst_model_atom_completeness.py`
Expected: the first assertion of `exercise_inherited_conformer` fails (TYR reported under C), then `AttributeError` for `complete_altlocs`.

- [ ] **Step 3: Change `_missing_atoms_records`**: after `ag = expected[0].parent()`, compute `own = set(g.altloc for g in ag.parent().atom_groups() if g.altloc.strip())`; `if mm.conf_altloc and own and mm.conf_altloc not in own: continue`; when `own` is empty use the blank key (`key = legacy_key = id_str`) regardless of `mm.conf_altloc`. Docstring of `missing_atoms_by_residue`: replace "Residues are checked per conformer." with the new rule (own conformers only; plain residues once).

- [ ] **Step 4: Add `complete_altlocs`** to `manager` with a numpydoc docstring (what it does, the group_args fields, the Sorry).

- [ ] **Step 5: Run the test to verify it passes**, then the module's neighbours: `mmtbx/regression/model/tst_model_atom_completeness.py`, `mmtbx/regression/tst_reduce2_capping.py`, and `libtbx.run_tests_parallel module=mmtbx nproc=8 run_in_unique_dirs=True` from an empty scratch directory (the key change can affect other reports; note any test that changed expectation and why).

---

### Task 3: reduce2 qmready integration

**Files:**
- Modify: `cctbx_project/mmtbx/programs/reduce2.py` (`run()`, the `qmready` block after `add_terminal_oxt()`, line ~2381)
- Create: `cctbx_project/mmtbx/regression/tst_reduce2_altloc_completion.py`
- Modify: `cctbx_project/mmtbx/run_tests.py` (register next to `tst_reduce2_capping.py`, line 499)

**Interfaces:**
- Consumes: Task 2 `model.complete_altlocs(log=...)`, Task 1 `report_lines`.
- Produces: reduce2 log block headed `Altloc completion (atoms borrowed between conformers of a residue):` followed by `report_lines` indented two spaces, printed only when the list is not empty; the written model is the completed one.

- [ ] **Step 1: Write the failing test** using `run`, `check_sorry` from `tst_reduce2_capping.py` (import them) and models derived from `tst_model_atom_completeness.tripeptide` (add `CRYST1` as `tst_reduce2_capping` models have one):

```python
def test_single_source_completed():
  # TYR CZ A/B, OH A only: before this change qmready stops with "Model must be atom-complete"
  run(model, "tst_r2ac_1", ["qmready=True"], fn_out=out)
  h = iotbx.pdb.input(out).construct_hierarchy()
  oh = [a for a in h.atoms() if a.name.strip() == "OH"]; assert len(oh) == 1 and oh[0].parent().altloc == ""
  hh = [a for a in h.atoms() if a.name.strip() == "HH"]; assert len(hh) == 1 and hh[0].parent().altloc == ""
  log text contains "Altloc completion" and " OH: altloc A -> shared (needed by B)"

def test_multi_source_copies():
  # ASN 2 with CB under A (0.81) and B (0.19, shifted 0.7 A), CG/OD1/ND2 under A 0.43, B 0.10, C 0.18, D 0.20
  # (C and D: chi1 rotations of A by +120 and -120 degrees about CA-CB)
  output has CB under A/B/C/D with occupancies 0.43/0.19/0.18/0.20 (backbone stays blank), with reduce2's HB2/HB3 under each of the four labels

def test_ambiguous_sorry():
  check_sorry(model_with_D_at_0_45, "tst_r2ac_3", ["qmready=True"], "Altloc completion", "ASN", "CB", "altloc A", "altloc B")

def test_truncated_still_needs_capping():
  # ASP 195 pattern: CB under A only, 3 A from CA-B: qmready=True gives the atom-complete Sorry, capping=True runs
  check_sorry(model, "tst_r2ac_4", ["qmready=True"], "Model must be atom-complete")
  run(model, "tst_r2ac_4b", ["qmready=True", "capping=True"], fn_out=out)

def test_constants_agree():
  from mmtbx.hydrogens import capping; from iotbx.pdb import altloc_completion
  assert altloc_completion.BOND_TOLERANCE == capping.GUESS_BOND_TOLERANCE
```

- [ ] **Step 2: Run to verify failure** (`test_single_source_completed` stops with the atom-complete Sorry).

- [ ] **Step 3: Implement**: in `run()` after `self.model.add_terminal_oxt()`:

```python
r = self.model.complete_altlocs()
lines = altloc_completion.report_lines(r)
if lines:
  print("Altloc completion (atoms borrowed between conformers of a residue):", file=self.logger)
  for l in lines: print("  " + l, file=self.logger)
```

The `Sorry` of an ambiguous borrow propagates as is (raised before anything is written, like the atom-complete Sorry). Update the `qmready` PHIL help: add "complete residues whose conformers lack atoms their other conformers have (altloc completion)".

- [ ] **Step 4: Run the test to verify it passes**; register it in `mmtbx/run_tests.py`; run `tst_reduce2_capping.py` and `tst_reduce2*.py` neighbours; `libtbx.find_clutter` on changed files.

---

### Task 4: qrefine `conformers.py`: fallback, membership multiplicity, guard, check

**Files:**
- Modify: `qrefine/conformers.py` (`multiplicity`, `_check_complete`, `split`, `from_conformers.__init__` print; module docstring)
- Modify: `qrefine/tests/unit/tst_31.py` (`exercise_incomplete_conformers` replaced; `exercise_split` extended)

**Interfaces:**
- Consumes: Task 1 `altloc_completion.bonded`, `altloc_completion.run(check_only=True)`.
- Produces:
  - `conformer_selections(pdb_hierarchy) -> list[group_args(label, iselection, n_fallback)]`: per conformer label X (in `conformer_labels` order) the ascending iselection of blank atoms, X atoms, hydrogen-only-label atoms, and, for every residue group with own labels lacking X, the atoms of its first own label; `n_fallback` is the number of such residue groups.
  - `multiplicity(pdb_hierarchy) -> flex.double`: per atom the number of selections containing it (never 0: no labels gives 1 everywhere), replacing the label rule.
  - `split(model)` unchanged signature and return; built on `conformer_selections`; runs the link guard; no `_check_complete`.
  - `check_links(pdb_hierarchy, label, parents)` raising `Sorry("conformer %s: residue %s takes its first label %s, but its %s is not bonded to %s of residue %s (%.2f A)")`; `parents` maps residue group memory_id to the chosen label.
  - `check_completed(pdb_hierarchy)`: raises `Sorry("Altloc completion needed: %d atoms to share and %d to copy (%s): run mmtbx.reduce2 qmready=true on the model first")` listing up to five `residue atom` pairs, when `altloc_completion.run(ph, check_only=True)` has `shared` or `copied`.
  - `from_conformers` prints, after its per-conformer line, `"  %d residue groups without label %s take their first own label"` when `n_fallback > 0`.

- [ ] **Step 1: Write the failing tests** in `tst_31.py`:

```python
def exercise_fallback():
  # h_altconf_complete.pdb (TYR 94 and GLY 95 with A/B) plus a C atom group on TYR 94 only
  # (copy of its B group, altloc C, occupancies A 0.4, B 0.35, C 0.25)
  sel = conformers.conformer_selections(ph)
  assert [s.label for s in sel] == ["A", "B", "C"] and [s.n_fallback for s in sel] == [0, 0, 1]
  gly95_a = i_seqs of GLY 95 altloc A atoms
  assert all(i in set(sel[2].iselection) for i in gly95_a)           # conformer C takes GLY 95 A
  m = conformers.multiplicity(ph)
  assert m.select(gly95_a).all_eq(2.0) and m.select(gly95_b).all_eq(1.0)
  assert m.select(tyr94_any_label).all_eq(1.0) and m.select(blank).all_eq(3.0)
  parts = conformers.split(model); assert len(parts) == 3
  assert parts[2].model.get_hierarchy().atoms_size() == sel[2].iselection.size()
  # old Sorry is gone
  conformers.split(model)   # no exception

def exercise_link_guard():
  # tst_41.pdb_str_swapped-style dipeptide: residue 1 N..O under A/B/C (C = A shifted 4 A along z,
  # occupancy 0.2), residue 2 under A/B with B shifted the same 4 A: conformer C of residue 1 bonds to
  # residue 2 B, but the first label is A -> Sorry naming both residues and C/N
  try: conformers.split(model)
  except Sorry as e: assert "conformer C" in str(e) and "residue" in str(e) and "not bonded" in str(e)
  else: raise AssertionError
  # the same with residue 2's labels swapped (A shifted): no Sorry

def exercise_check_completed():
  # TYR 94 of h_altconf_complete.pdb with OH removed from B: Sorry naming reduce2
  try: conformers.check_completed(ph)
  except Sorry as e: assert "mmtbx.reduce2 qmready=true" in str(e) and "OH" in str(e)
  # complete model: no exception
```

Replace `exercise_incomplete_conformers` with `exercise_fallback`; keep every other exercise. In `exercise_split`, the existing assertion `m.select(ci == 0).all_eq(2.0) and m.select(ci != 0).all_eq(1.0)` stays valid.

- [ ] **Step 2: Run to verify failure**: `<root>/build/bin/libtbx.python qrefine/tests/unit/tst_31.py` (AttributeError `conformer_selections`).

- [ ] **Step 3: Implement** in `conformers.py`. `conformer_selections`: one pass over residue groups computing `own = conformer_labels of the rg` (reuse `altloc_completion.conformer_labels`) and `parent = own[0]`; then for each label X and each atom: include if `altloc` is blank or hydrogen-only or `== X`, or (`X not in own` and `altloc == parent`). Build `parents_by_label` for the guard. `check_links(ph, X, parents)`: for consecutive residue groups of a chain (i, i+1) where one lacks X (parent P) and the other has X: take the lacking one's link atom under P (C if it is residue i, N if i+1) and the other's link atom under X; if both exist (labelled, not blank), require `altloc_completion.bonded`. Call it inside `split` for each label. Update the module docstring (the paragraph "Every residue group that has conformer labels must have all of the model's conformer labels" becomes the fallback rule and the guard; multiplicity paragraph: membership).

- [ ] **Step 4: Run tst_31, then the altloc tests**: `tst_05`, `tst_06`, `tst_11`, `tst_25`, `tst_31`, `tst_40`-`tst_46` with `libtbx.python`, each from the scratch directory. Expected: all `OK`.

---

### Task 5: qr.refine check-only load, scenario helpers, `tst_49` (inherited conformer)

**Files:**
- Modify: `qrefine/qr.py` (after the `Reconciling altlocs` block, before `params = ...get_default_pdb_interpretation_params()`, line ~460)
- Create: `qrefine/tests/unit/altloc_scenarios.py`
- Create: `qrefine/tests/unit/tst_49.py`
- Modify: `qrefine/run_tests.py` (append `tst_49`)

**Interfaces:**
- Consumes: Task 4 `conformers.check_completed`, `conformers.n_conformers`; `tst_00.pdb_str_good`; `tst_31.build_cctbx`, `gradients_from_conformers`; `qrefine.tests.unit.run_fmodel`, `run_qrefine_in`.
- Produces `altloc_scenarios.py`:
  - `base_model() -> mmtbx.model.manager` from `tst_00.pdb_str_good` (CRYST1 21.937 4.866 23.477 90 107.08 90 P 1 21 1; GLY 1 ASN 2 ASN 3 GLN 4 GLN 5 ASN 6 TYR 7 with H, 7 waters).
  - `residue(model, resseq) -> residue_group`
  - `side_chain_conformers(rg, labels, occupancies, chi1_degrees)`: for the atoms beyond CB (CG onward and their H, HD/HE/HG names), create one atom group per label with the atoms rotated about the CA-CB axis by the given angle (0 keeps the original), occupancies as given; atoms not beyond CB stay blank. `labels[0]` reuses the existing atoms (relabelled).
  - `label_atoms(rg, names, label, occupancy)`: move the named blank atoms into the atom group of `label` (created if absent) and set their occupancy (the "A only" patterns).
  - `backbone_conformers(model, resseq, label, occupancy, degrees)`: adds a copy of N, CA, C, O, H, HA of residue `resseq` under `label` rotated by `degrees` about the axis C(resseq-1) -> N(resseq+1) (both peptide bond lengths preserved), relabels the originals to the first existing label (or "A") and sets occupancies.
  - `write(model, file_name)`, `shaken(model, mean_distance=0.3) -> model` (`xray_structure.shake_sites_in_place`, seed fixed), `r_work(f_obs, flags, model) -> float` (as tst_00's fmodel block: `mmtbx.f_model.manager` + `update_all_scales()`, so bulk solvent and scaling are included), `f_calc_rel_diff(model_1, model_2, d_min=1.2) -> float`, `labels_and_occ(model) -> dict[(resseq, name, altloc)] = occ` for comparing output with input.
  - `assert_scattering_neutral(model_in, model_done, copies_expected=True)`: the one check of spec section 5 used by every completion test. Asserts `f_calc_rel_diff(model_in, model_done) < 1e-6` (amplitudes of `xray_structure.structure_factors(d_min=1.2)` on identical index sets) and, when `copies_expected`, `model_done.get_hierarchy().atoms_size() > model_in.get_hierarchy().atoms_size()` so the check is not vacuous. Fcalc being unchanged is the whole claim: R-factors follow from it.
- Produces in `qr.py`: the check-only Sorry at load.

- [ ] **Step 1: Write the failing test** `tst_49.py` with this docstring (adapt wording, keep the content):

> 1US0 pattern (ARG 40, SER 282): residues with A/B side chains in a model where one residue also has a C side chain. Residues 3 and 5 (ASN, GLN) get A/B chi1 conformers (0.5/0.5), residue 4 (GLN) gets A/B/C (0.4/0.35/0.25). The model is complete in its own conformers, so completion changes nothing and reduce2 would not touch it; qr.refine builds conformer C from conformer A of residues 3 and 5 (first own label). Part 1, in process with cctbx: with the C side chain placed exactly on A, every labelled atom gets the gradient of the A/B-only model and blank atoms obey 3 g_ABC = 2 g_AB + g_A. Part 2: qr.refine restraints=cctbx on a shaken copy against fmodel data from the model: runs, prints the fallback line, keeps every label and occupancy, lowers R-work.

```python
def run(prefix=...):
  os.makedirs(prefix, exist_ok=True); os.chdir(prefix)
  m = base_model(); side_chain_conformers(residue(m,3), ["A","B"], [0.5,0.5], [0,120])
  side_chain_conformers(residue(m,5), ["A","B"], [0.5,0.5], [0,-120])
  side_chain_conformers(residue(m,4), ["A","B","C"], [0.4,0.35,0.25], [0,120,-120])
  assert not m.complete_altlocs(check_only=True).changed
  # part 1
  m_id = same with residue 4 angles [0,120,0]; m_ab = same residues 3,5 and residue 4 ["A","B"],[0.5,0.5],[0,120]
  g_abc = gradients_from_conformers(m_id); g_ab = gradients_from_conformers(m_ab); g_a = cctbx gradient of conformers.split(m_id)[0].model mapped back
  assert approx_equal(g_abc.select(labelled_in_both), g_ab.select(labelled_in_both), 1e-6)   # same atoms by (resseq, name, altloc)
  assert approx_equal(3*g_abc.select(blank), 2*g_ab.select(blank) + g_a.select(blank), 1e-6)
  # part 2
  write(m, good); r_fm = run_fmodel(prefix=prefix, args=[good, "high_res=1.2"])
  write(shaken(m), poor)
  r_qr = run_qrefine_in(prefix=prefix, args=[poor, r_fm.mtz, "restraints=cctbx", "clustering=False", "stpmax=0.3"])
  assert "2 residue groups without label C take their first own label" in r_qr.std_out_text
  assert labels_and_occ(r_qr.model) == labels_and_occ(m)
  assert r_work(r_fm.f_obs, r_fm.flags, r_qr.model) < 0.5 * r_work(r_fm.f_obs, r_fm.flags, shaken model)
```

- [ ] **Step 2: Run to verify failure** (`ImportError: altloc_scenarios`).

- [ ] **Step 3: Implement `altloc_scenarios.py`** (rotation with `scitbx.matrix.col` and `rotate_around_origin` about the normalised axis through the first axis atom; occupancies set on every atom of the group; `atoms().reset_i_seq()` and `reset_serial()` after each edit; rebuild the model with `mmtbx.model.manager(model_input=None, pdb_hierarchy=..., crystal_symmetry=...)` so the X-ray structure matches). Then in `qr.py` add, after the reconcile block:

```python
if(conformers.n_conformers(self.model.get_hierarchy()) > 1):
  conformers.check_completed(self.model.get_hierarchy())
```

- [ ] **Step 4: Run tst_49** (`libtbx.python qrefine/tests/unit/tst_49.py` from the scratch dir). Expected `OK`, under 60 s. Tune `stpmax` / shake distance only if the R assertion is not met for a reason you can state. Register `tst_49` in `qrefine/run_tests.py`.

---

### Task 6: `tst_50` single-source atoms (1PWM)

**Files:**
- Create: `qrefine/tests/unit/tst_50.py`
- Modify: `qrefine/run_tests.py`

**Interfaces:** Task 5 helpers; Task 2 `complete_altlocs`; Task 3 reduce2.

- [ ] **Step 1: Write the test** with this docstring:

> 1PWM pattern (LYS 85, ASN 129): atoms that exist under one label only while the rest of the residue has two. Residue 3 (ASN): CB with HB2, HB3 under A only (0.46), CG, OD1, ND2 and their H under A (0.46) and B (0.54, chi1 + 120). Residue 4 (GLN): N, CA, C, O, H, HA under A only (0.52), CB onward under A (0.52) and B (0.48, chi2 rotation, so CB-B equals CB-A). Completion shares the single-source atoms (blank, occupancies kept) because each is bonded to the other conformer: structure factors are unchanged, the result equals the hand-made model, mmtbx.reduce2 qmready=true writes it, qr.refine on the uncompleted model stops and points to reduce2, and on the completed model refines.

Assertions: `r.shared` names exactly CB, HB2, HB3 of residue 3 and N, CA, C, O, H, HA of residue 4 with source "A"; `assert_scattering_neutral(m_in, m_done, copies_expected=False)` (sharing adds no atoms); `labels_and_occ(m_done) == labels_and_occ(m_hand)` and max coordinate difference < 1e-6, where `m_hand` is built with `label_atoms` never called (the atoms simply left blank with occupancies 0.46 / 0.52); `mmtbx.reduce2 <in> qmready=True` via `easy_run.fully_buffered(...).raise_if_errors()` writes a model whose `labels_and_occ` of heavy atoms equals `m_done`'s; `run_qrefine_in(args=[in, mtz, ...], sorry_expected=True)` and `"mmtbx.reduce2 qmready=true" in std_out_text`; refinement of the shaken completed model lowers R-work as in tst_49.

- [ ] **Step 2: Run to verify failure, Step 3: adjust helpers only if needed, Step 4: run to pass, register.**

---

### Task 7: `tst_51` multi-source atoms and the ambiguous Sorry (1R6J SER 205)

**Files:**
- Create: `qrefine/tests/unit/tst_51.py`
- Modify: `qrefine/run_tests.py`

- [ ] **Step 1: Write the test** with this docstring:

> 1R6J pattern (SER 205): a side chain with four conformers on a backbone with two. Residue 3 (ASN): backbone and CB under A (0.81) and B (0.19, backbone rotated 20 degrees about C(2)...N(4), CB following CA-B); CG, OD1, ND2 with H under A (0.43, original), B (0.10, rotated with the B backbone), C (0.18, chi1 + 120 on A) and D (0.20, chi1 angle searched in steps of 10 degrees so that CG-D is bonded to both CB-A and CB-B and farthest from CG-A/B/C). Completion gives C and D copies of CB and of the backbone from A, chosen by the bond test (CB-A only for C) and by the occupancy budget (CB-A over CB-B for D: B has 0.19 for 0.20 + 0.10). Occupancies: CB A 0.43, B 0.19, C 0.18, D 0.20; structure factors unchanged although 10 heavy atoms and their H were added (`assert_scattering_neutral(m_in, m_done)`); reduce2 qmready writes the same; split gives four complete conformer models; refinement lowers R-work. With the D side chain at 0.45 no candidate passes the budget: Sorry naming the residue, CB and both candidates.

Assertions as the docstring states, plus `[len(p.model.get_hierarchy().atoms()) for p in conformers.split(m_done)]` all equal (every conformer model has the same atom count) and `[s.n_fallback for s in conformers.conformer_selections(ph)] == [0, 0, 0, 0]`. The ambiguous variant uses `run_qrefine_in(..., sorry_expected=True)` only after `m.complete_altlocs(check_only=True)` raised `Sorry` in process with "ASN" and "CB" and "altloc A" and "altloc B" in the message.

- [ ] **Steps 2-4 as Task 6.**

---

### Task 8: `tst_52` truncated in one conformer (1R6J ASP 195)

**Files:**
- Create: `qrefine/tests/unit/tst_52.py`
- Modify: `qrefine/run_tests.py`

- [ ] **Step 1: Write the test** with this docstring:

> 1R6J pattern (ASP 195): backbone A/B, CB under A only, side chain beyond CB under A/C. Residue 3 (ASN): backbone under A (0.82) and B (0.18, rotated 90 degrees about C(2)...N(4) so CA-B is more than 2.5 A from CB-A); CB, HB2, HB3 under A only (0.82); CG, OD1, ND2 with H under A (0.41) and C (0.41, chi1 + 120). Conformer B cannot borrow CB (not bonded to anything of B) and is truncated at CA; conformer C borrows CB-A (anchored by CG-C) and then the A backbone (CA-B is not bonded to CB), each a single candidate, copied with 0.41/0.41. Structure factors unchanged; qr.refine with cctbx restraints tolerates the truncated conformer and lowers R-work (a QM run would need reduce2 capping=True).

Assertions: `r.truncated == [(id, "B", [CB, CG, HB2, HB3, HD21, HD22, ND2, OD1 sorted], reason)]` with "not borrowed" and "CA B" in `reason`; `sorted(n for _, n, _, t, _ in r.copied if t == "C")` equals `["C", "CA", "CB", "H", "HA", "HB2", "HB3", "N", "O"]`; occupancies of CB: A 0.41, C 0.41; backbone A 0.41, B 0.18, C 0.41; `assert_scattering_neutral(m_in, m_done)`; refinement lowers R-work.

- [ ] **Steps 2-4 as Task 6.**

---

### Task 9: Acceptance on the three entries and the full suites

**Files:** none changed. Working directory: a fresh `<root>/rt/run_$(date +%s)_acceptance`.

- [ ] **Step 1: reduce2** on `1us0.pdb`, `1pwm.pdb`, `1r6j.pdb` (copy from `rt/run_1793639952_reduce2/`) using the clash-listing driver `rt/run_1793639952_reduce2/reduce2_list_clashes.py` (the water clash gate is separate work) with `qmready=true`, plus `capping=true` for 1US0 and 1R6J. Expected: 1PWM log lists four shared atoms; 1R6J log lists copies for SER 205 (C, D), HIS 202 (C), ASP 195 (C) and ASP 195 truncated in B; 1US0 lists nothing under Altloc completion and only HIS 312 / GLU 313 as missing. Record the three logs.

- [ ] **Step 2: scattering neutrality on the real entries**: for 1PWM and 1R6J, in process, `f_calc_rel_diff` between the deposited model and the model after `complete_altlocs()` alone (before H are added), at d_min 1.0 A, must be below 1e-6; report the two numbers. Then **qr.refine** `restraints=cctbx clustering=False number_of_macro_cycles=1` on each reduce2 output with its own data. Expected: no Sorry about completion; the log shows the fallback counts (1US0: 78 residue groups for C); R-work after one cycle not above the start by more than 0.005.

- [ ] **Step 3: Full suites**: `libtbx.run_tests_parallel module=iotbx nproc=8 run_in_unique_dirs=True`, the same for `mmtbx`, and `qr.test` (or `libtbx.python qrefine/run_tests.py`) from empty scratch directories. Expected: all pass; list any failure with its output.

- [ ] **Step 4: Report**: the diff summary (`git -C cctbx_project status --short`, `git -C qrefine status --short`), the acceptance numbers, and `libtbx.find_clutter` output for every changed file. Leave everything uncommitted.
