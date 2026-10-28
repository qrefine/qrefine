# Reduce2-based fragment capping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** qr.refine caps its QM fragments with `mmtbx.hydrogens.capping` (shared with `mmtbx.reduce2 capping=True`), `finalise.py`/`completion.py` are deleted, and `qr.finalise` runs `mmtbx.reduce2 qmready=True`.

**Architecture:** `mmtbx/hydrogens/capping.py` gains `cut_bonds` (one H per bond a fragment cut severed), `cut_bonds_by_distance` and `cap_model` (the find / charge / cap sequence reduce2 runs inline today); reduce2 calls `cap_model`. qrefine gets `fragment_capping.capper`, which processes a fragment once, caches the dangling-atom records and the charge, and re-places the caps on new coordinates without re-interpretation; `fragment.py` and `restraints.py` consume the capped hierarchy through a `fragment_selection` mask instead of index slicing.

**Tech Stack:** cctbx/mmtbx (Python 3.11, `iotbx.pdb.hierarchy`, `mmtbx.model.manager`, `mmtbx.ligands.electrons`), qrefine test harness (`qrefine.python tests/unit/tst_NN.py`).

**Spec:** `docs/superpowers/specs/2026-10-27-reduce2-capping-design.md` (this repo). Read it first.

## Global Constraints

- Two git repos: `modules/cctbx_project` (mmtbx) and `modules/qrefine`. **Do not commit**: the user commits. Each task ends by listing the changed files instead.
- Python dispatchers (full paths): `/Users/pafonine/Desktop/all/phenix/build/bin/libtbx.python` for mmtbx tests, `/Users/pafonine/Desktop/all/phenix/build/bin/qrefine.python` for qrefine tests. Run tests from a scratch directory under `/private/tmp/claude-502/-Users-pafonine-Desktop-all-phenix/56c077e3-afa4-4664-b497-da86370e2815/scratchpad/`, never inside the repos.
- reduce2's behaviour and output text do not change: `mmtbx/regression/tst_reduce2_capping.py`, `tst_reduce2_qmready_charge.py`, `tst_reduce2_qmready_clash.py` must pass unchanged. The Sorry text `capping=True: no dictionary for %s; set capping_guess_for_unknowns=True to cap them by element valence` is asserted by a test and stays verbatim.
- Neutral link-atom caps at cuts (NH2 on N, one H on C, HG on SG); a cut through a phosphodiester is a Sorry (hypervalent P). Never guess protonation.
- No capped hierarchy is ever passed to pdb_interpretation / `charges_class`; charges come from `electrons.run(model=, extra_sigma_bonds=)` before capping.
- Cap names come from `capping.cap_names`: `H2` on a backbone N, `H` on a PRO N, `HG` on SG, `HOXT` on a cut carbonyl C, `H` + missing atom name on truncations.
- Follow the existing style of each file (2-space indent, `from __future__` lines, docstrings as in capping.py).

## Review Focus

1. A fragment atom bonded to a metal that lies outside the fragment (CYS SG with Zn outside): no cap, a coordination bond satisfies no valence. Test: `test_cut_bonds_by_distance_metal` (Task 1).
2. A fragment with alternate conformations: `electrons` miscounts silently. Expected: Sorry from the capper. Test: `test_altloc_sorry` (Task 3).
3. A master not prepared with `reduce2 qmready=True` (a backbone N with only CA and H and no cut bond): the charge would silently be off by one. Expected: Sorry telling the user to run reduce2 qmready=True. Test: `test_not_qmready_sorry` (Task 3).
4. A cut through a nucleic-acid phosphodiester: the P side is hypervalent. Expected: Sorry naming the P; the O3' side alone gets `HO3'`. Test: `test_nucleic_cut` (Task 3).
5. `capper.cap()` called at every gradient evaluation: caps must not accumulate and the cached fragment must not change. Test: `test_cap_twice` (Task 3).

---

### Task 1: `cut_bonds` in `find_dangling_atoms` and `cut_bonds_by_distance` (mmtbx)

**Files:**
- Modify: `modules/cctbx_project/mmtbx/hydrogens/capping.py` (module docstring; `find_dangling_atoms` at line 405; new function after `extra_sigma_bonds`, line 589)
- Test: `modules/cctbx_project/mmtbx/hydrogens/tst_capping.py` (add tests before the `# ---` builder section near line 765; register them in the file's `run()`/`__main__` the way the existing tests are)

**Interfaces:**
- Produces: `find_dangling_atoms(model, guess_for_unknowns=False, log=None, cut_bonds=None)`; `cut_bonds` is `{i_seq: n}` in `model`.
- Produces: `cut_bonds_by_distance(hierarchy, selection)` returning `{i_seq: n}` where `i_seq` indexes the atoms of `hierarchy.select(selection)` in selection order (what `hierarchy.select(selection).atoms().reset_i_seq()` gives); `selection` is a `flex.bool` over `hierarchy.atoms()`.

- [ ] **Step 1: Write the failing tests** in `tst_capping.py`

```python
def i_seq_of(model, resseq, name):
  for a in model.get_hierarchy().atoms():
    if (a.parent().parent().resseq.strip() == str(resseq) and
        a.name.strip() == name):
      return a.i_seq
  raise AssertionError((resseq, name))

def test_cut_bonds_backbone_n():
  '''ALA 2 of pdb_ala_ala cut from ALA 1: N 2 (CA, H) lost its peptide C.
  cut_bonds gives it one H, named H2 (its free dictionary H); without
  cut_bonds it is not examined (no heavy dictionary partner missing).'''
  model = model_with_h(pdb_ala_ala)
  model = model.select(model.selection('resseq 2'))
  n = i_seq_of(model, 2, 'N')
  assert ('ALA', '2', '', 'N') not in records_by_name(capping.find_dangling_atoms(model))
  result = capping.find_dangling_atoms(model, cut_bonds={n: 1})
  r = records_by_name(result)[('ALA', '2', '', 'N')]
  assert (r.n_fill, r.coordination, r.missing_names, r.h_names) == (1, 3, [], ['H2']), r
  assert capping.extra_sigma_bonds(result)[n] == 1
  capping.cap_dangling_atoms(model, result)
  assert h_names_on(model, 2, 'N') == ['H', 'H2'], h_names_on(model, 2, 'N')

def test_cut_bonds_keeps_dictionary_record():
  '''ALA 1 cut from ALA 2: C 1 misses OXT (leaving, nothing replaces it):
  n_fill 1 already; a cut bond on it changes nothing and adds no record.'''
  model = model_with_h(pdb_ala_ala)
  model = model.select(model.selection('resseq 1'))
  c = i_seq_of(model, 1, 'C')
  plain = records_by_name(capping.find_dangling_atoms(model))[('ALA', '1', '', 'C')]
  result = capping.find_dangling_atoms(model, cut_bonds={c: 1})
  assert len([r for r in result.records if r.i_seq == c]) == 1
  assert records_by_name(result)[('ALA', '1', '', 'C')].n_fill == plain.n_fill == 1

def test_cut_bonds_ignored_on_h_and_metal():
  model = model_with_h(pdb_cys_zn)
  zn = [a.i_seq for a in model.get_hierarchy().atoms() if a.element.strip().upper() == 'ZN'][0]
  h = [a.i_seq for a in model.get_hierarchy().atoms() if a.element_is_hydrogen()][0]
  result = capping.find_dangling_atoms(model, cut_bonds={zn: 1, h: 1})
  assert not [r for r in result.records if r.i_seq in (zn, h)], records_by_name(result)

def test_cut_bonds_by_distance_disulfide():
  '''CYS 232 of pdb_disulfide selected alone: its SG lost SG 128.'''
  model = model_with_h(pdb_disulfide)
  h = model.get_hierarchy()
  sel = h.atom_selection_cache().selection('resseq 232')
  cut = capping.cut_bonds_by_distance(h, sel)
  frag = h.select(sel); frag.atoms().reset_i_seq()
  sg = [a.i_seq for a in frag.atoms() if a.name.strip() == 'SG'][0]
  assert cut == {sg: 1}, cut
  # the whole model: nothing is cut
  assert capping.cut_bonds_by_distance(h, flex.bool(h.atoms().size(), True)) == {}

def test_cut_bonds_by_distance_metal():
  '''pdb_cys_zn without its ZN: the coordination bond is no cut bond.'''
  model = model_with_h(pdb_cys_zn)
  h = model.get_hierarchy()
  sel = ~h.atom_selection_cache().selection('element ZN')
  assert capping.cut_bonds_by_distance(h, sel) == {}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd <scratch> && /Users/pafonine/Desktop/all/phenix/build/bin/libtbx.python /Users/pafonine/Desktop/all/phenix/modules/cctbx_project/mmtbx/hydrogens/tst_capping.py`
Expected: FAIL with `TypeError: find_dangling_atoms() got an unexpected keyword argument 'cut_bonds'`.

- [ ] **Step 3: Implement `cut_bonds` in `find_dangling_atoms`**

In the per-atom loop of `find_dangling_atoms`, after the dictionary pass over a residue's atoms, add records for the atoms of `cut_bonds` that are in that atom group, are not H, not metals, and got no record in the pass: `n_fill = n`, neighbours `_conformer_neighbours(atom, table.get(atom.i_seq, ()), atoms, is_metal)`, `missing_names = []`, `h_names` = the dictionary H partners of the atom absent from the conformer (as in the pass), `source = 'dictionary'`, via `_record`. Atoms of residues without a dictionary are handled only in guess mode (`_guess_atom_group` already counts every bonded neighbour; a cut partner is absent, so add `n` to its `n_fill` there too). Add rule 9 to the module docstring: "Cut bonds (cut_bonds, {i_seq: n}): bonds of the model's atoms to atoms not in the model, a fragment cut out of a larger model. An atom without a record gets one with n_fill n and its free dictionary H names (H2 on a backbone N, HG on a cut SG, HO3' on a cut O3'); an atom with a record keeps it, its n_fill already counts the absent partner. H and metals are ignored."

- [ ] **Step 4: Implement `cut_bonds_by_distance(hierarchy, selection)`**

Candidate pairs: `hierarchy.extract_xray_structure(crystal_symmetry=cs)` where `cs` is `uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(sites_cart=hierarchy.atoms().extract_xyz(), buffer_layer=5)` as a P1 `crystal.symmetry`, then `xrs.pair_asu_table(distance_cutoff=GUESS_BOND_CUTOFF).extract_pair_sym_table()`. Keep pairs with exactly one atom selected, both non-H and non-metal (`atom_property().is_metal`), and distance below `_covalent_radius(e1) + _covalent_radius(e2) + GUESS_BOND_TOLERANCE` (skip a pair whose element has no radius). Count per selected atom and map its hierarchy i_seq to its position in `selection.iselection()`. Docstring states the key convention.

- [ ] **Step 5: Run the tests to verify they pass**

Run: same command as Step 2. Expected: `OK` (the file's final print), no assertion errors.

- [ ] **Step 6: Report** the changed files (`capping.py`, `tst_capping.py`); no commit.

---

### Task 2: `cap_model` and reduce2 switched to it (mmtbx)

**Files:**
- Modify: `modules/cctbx_project/mmtbx/hydrogens/capping.py` (`cap_dangling_atoms` line 923; new `cap_model` after it)
- Modify: `modules/cctbx_project/mmtbx/programs/reduce2.py` (`_CalculateTotalCharge` 1912-1953, `_CapDanglingAtoms` 1955-2006, the qmready block in `run()` 2553-2582)
- Test: `modules/cctbx_project/mmtbx/hydrogens/tst_capping.py`

**Interfaces:**
- Consumes: Task 1.
- Produces: `cap_dangling_atoms(model, result, use_neutron_distances=True, mon_lib_srv=None)` (`mon_lib_srv` overrides `model.get_mon_lib_srv()` for the dictionary lookups of the report; the server of the processed model the records came from).
- Produces: `cap_model(model, guess_for_unknowns=False, cut_bonds=None, use_neutron_distances=True, calculate_charge=True, dangling=None, mon_lib_srv=None, log=None)` returning `group_args(dangling, capped, total_charge, formal_charges, warnings)`: `dangling` the `find_dangling_atoms` result used (computed unless given), `capped` the `cap_dangling_atoms` result, `total_charge` an int or None, `formal_charges` the `electrons` distribution or None, `warnings` the open-shell lines `'open-shell cap on <id_str>'`. With `dangling` given, `find_dangling_atoms` is not called (the model may be unprocessed); `calculate_charge=True` requires a processed model (`electrons.run(model=...)`). Sorry on `dangling.unknown_residues` when `guess_for_unknowns` is False, with the verbatim text in Global Constraints.

- [ ] **Step 1: Write the failing tests**

```python
def test_cap_model_matches_two_step():
  '''cap_model on pdb_arg (ARG 3 at CB): same caps and charge as
  find + extra_sigma_bonds + electrons + cap_dangling_atoms by hand.'''
  a, b = model_with_h(pdb_arg), model_with_h(pdb_arg)
  dangling = capping.find_dangling_atoms(a)
  charge = electrons.run(model=a, extra_sigma_bonds=capping.extra_sigma_bonds(dangling))
  capping.cap_dangling_atoms(a, dangling)
  r = capping.cap_model(b)
  assert r.total_charge == charge, (r.total_charge, charge)
  assert r.capped.n_added == 3 and r.warnings == [], (r.capped.n_added, r.warnings)
  assert [x.name for x in a.get_hierarchy().atoms()] == [x.name for x in b.get_hierarchy().atoms()]
  assert a.get_hierarchy().atoms().extract_xyz().max_distance(
    b.get_hierarchy().atoms().extract_xyz()) < 1.e-6

def test_cap_model_cached_dangling_unprocessed():
  '''The records of a processed model placed on an unprocessed manager of a
  copy of its hierarchy (no restraints): identical caps, no charge.'''
  model = model_with_h(pdb_arg)
  dangling = capping.find_dangling_atoms(model)
  copy = mmtbx.model.manager(model_input=None,
    pdb_hierarchy=model.get_hierarchy().deep_copy(),
    crystal_symmetry=model.crystal_symmetry(), log=null_out())
  r = capping.cap_model(copy, dangling=dangling, calculate_charge=False,
    mon_lib_srv=model.get_mon_lib_srv())
  ref = capping.cap_dangling_atoms(model, dangling)
  assert r.total_charge is None and r.formal_charges is None
  assert r.capped.lines == ref.lines, (r.capped.lines, ref.lines)
  a, b = model.get_hierarchy().atoms(), copy.get_hierarchy().atoms()
  assert list(a.extract_name()) == list(b.extract_name())
  assert a.extract_xyz().max_distance(b.extract_xyz()) < 1.e-6

def test_cap_model_unknown_residue_sorry():
  try:
    capping.cap_model(model_with_h(pdb_unl))
  except Sorry as e:
    assert 'capping=True: no dictionary for' in str(e), str(e)
  else:
    raise AssertionError('no Sorry')
```

- [ ] **Step 2: Run to verify they fail**

Run: the Task 1 command. Expected: `AttributeError: module 'mmtbx.hydrogens.capping' has no attribute 'cap_model'`.

- [ ] **Step 3: Implement `mon_lib_srv` on `cap_dangling_atoms` and `cap_model` in `capping.py`**

`cap_model`: `dangling = dangling or find_dangling_atoms(model, guess_for_unknowns, log, cut_bonds)`; Sorry on unknown residues when not guessing; `total_charge, formal_charges, warnings = None, None, []`; if `calculate_charge`: `extra = extra_sigma_bonds(dangling)`; `formal_charges = electrons.run(model=model, extra_sigma_bonds=extra, return_formal_charges=True)`; `total_charge = formal_charges.get_total_charge()`; `warnings` = `'open-shell cap on %s' % atoms[i].id_str()` for each `i` in `sorted(extra)` with `formal_charges.get(i, 0)` truthy (the loop now in reduce2's `_CalculateTotalCharge`). Then `capped = cap_dangling_atoms(model, dangling, use_neutron_distances, mon_lib_srv)`. Document in the module docstring under "Building".

- [ ] **Step 4: Run to verify they pass**

Run: the Task 1 command. Expected: OK.

- [ ] **Step 5: Switch reduce2 to `cap_model`**

In `run()` replace the block from `dangling = None` to the end of `if self.params.capping: outString += ...` with:

```python
      if self.params.capping:
        r = capping.cap_model(qm_model,
          guess_for_unknowns    = self.params.capping_guess_for_unknowns,
          use_neutron_distances = self.params.use_neutron_distances,
          calculate_charge      = self.params.approach == 'add',
          log                   = self.logger)
        if self.params.approach == 'add':
          self._ApplyTotalCharge(r.total_charge)
        outString += self._CapDanglingAtoms(qm_model, r)
      elif self.params.approach == 'add':
        self._CalculateTotalCharge(qm_model)
```

`_CalculateTotalCharge(self, model)`: computes `electrons.run(model=model, return_formal_charges=True).get_total_charge()` and calls `_ApplyTotalCharge`. New `_ApplyTotalCharge(self, total_charge)`: the `make_sub_header('Calculating total charge')`, `self.total_charge = ...`, the print and the output-file renaming that `_CalculateTotalCharge` does today. `_CapDanglingAtoms(self, model, r)`: the output manager and the report as today, reading `r.dangling.records`, `r.dangling.skipped`, `r.capped.lines`, `r.capped.n_atoms/n_residues/n_added`, `list(r.capped.warnings) + list(r.warnings)`. Headers print in the same order as before (charge, then capping).

- [ ] **Step 6: Run the reduce2 tests**

Run (from scratch): `libtbx.python .../mmtbx/regression/tst_reduce2_capping.py && libtbx.python .../mmtbx/regression/tst_reduce2_qmready_charge.py && libtbx.python .../mmtbx/regression/tst_reduce2_qmready_clash.py` and the Task 1 command. Expected: each ends with OK, no assertion errors.

- [ ] **Step 7: Report** changed files (`capping.py`, `tst_capping.py`, `programs/reduce2.py`); no commit.

---

### Task 3: `qrefine/fragment_capping.py` and `tests/unit/tst_38.py`

**Files:**
- Create: `modules/qrefine/fragment_capping.py`
- Create: `modules/qrefine/tests/unit/tst_38.py`
- Modify: `modules/qrefine/utils/hierarchy_utils.py` (add `model_from_hierarchy`), `modules/qrefine/restraints.py:30-47` (remove the definition, `from qrefine.utils.hierarchy_utils import model_from_hierarchy`), `modules/qrefine/run_tests.py` (add `"$D/tests/unit/tst_38.py"`)

**Interfaces:**
- Consumes: Task 1 and Task 2 functions.
- Produces:

```python
class capper(object):
  def __init__(self, parent_hierarchy, selection, crystal_symmetry, cif_objects=None, log=None)
  # attributes: fragment_hierarchy (parent.select(selection), i_seq reset, never changed),
  #   crystal_symmetry, charge (int), n_caps (int), lines (list of str), dangling
  def cap(self, sites_cart=None):
    # -> group_args(pdb_hierarchy, fragment_selection)
    # pdb_hierarchy: a capped deep copy of fragment_hierarchy at sites_cart (None: setup sites),
    #   caps inside their residues; i_seq and serial reset
    # fragment_selection: flex.bool over pdb_hierarchy.atoms(), True for the fragment atoms,
    #   in fragment order: pdb_hierarchy.select(fragment_selection) is the fragment at sites_cart
```
  `hierarchy_utils.model_from_hierarchy(pdb_hierarchy, crystal_symmetry, cif_objects=None, log=None)`: the function now in restraints.py, unchanged.

- [ ] **Step 1: Write the failing test `tests/unit/tst_38.py`** following the structure of `tst_35.py` (`run(prefix=...)`, `os.makedirs(prefix)`, `os.chdir(prefix)`, `if __name__ == '__main__': run()`), with these cases as functions called from `run()`:

```python
def fragment_of(file_name, sel_str):
  pi = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data, file_name))
  ph = pi.construct_hierarchy()
  return ph, ph.atom_selection_cache().selection(sel_str), pi.crystal_symmetry()

def cap_names(r):
  return sorted(a.name.strip() for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection))

def check_fragment_atoms(c, r, ph, sel):
  frag = ph.select(sel)
  kept = r.pdb_hierarchy.atoms().select(r.fragment_selection)
  assert kept.size() == frag.atoms().size()
  assert list(kept.extract_name()) == list(frag.atoms().extract_name())
  assert kept.extract_xyz().max_distance(frag.atoms().extract_xyz()) < 1.e-6

def test_peptide_cut():
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  r = c.cap()
  assert (c.n_caps, c.charge) == (2, 0), (c.n_caps, c.charge)
  assert cap_names(r) == ['H2', 'HOXT'], cap_names(r)
  check_fragment_atoms(c, r, ph, sel)
  for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection):
    parent = {'H2': 'N', 'HOXT': 'C'}[a.name.strip()]
    p = [x for x in a.parent().atoms() if x.name.strip() == parent][0]
    assert 0.95 < a.distance(p) < 1.15, (a.name, a.distance(p))

def test_truncated_side_chain():
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb",
    "resseq 3:5 and not (resseq 4 and (name NE2 or name HE21 or name HE22))")
  c = fragment_capping.capper(ph, sel, cs)
  assert (c.n_caps, c.charge) == (3, 0), (c.n_caps, c.charge)
  assert cap_names(c.cap()) == ['H2', 'HNE2', 'HOXT']

def test_pro_cut():
  ph, sel, cs = fragment_of("2lvr.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  r = c.cap()
  assert (c.n_caps, c.charge) == (2, 0), (c.n_caps, c.charge)
  assert cap_names(r) == ['H', 'HOXT'], cap_names(r)
  h = [a for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection) if a.name.strip() == 'H'][0]
  assert h.parent().resname == 'PRO'

def test_disulfide_cut():
  from mmtbx.hydrogens.tst_capping import pdb_disulfide, model_with_h
  model = model_with_h(pdb_disulfide)
  ph = model.get_hierarchy()
  sel = ph.atom_selection_cache().selection("resseq 232")
  c = fragment_capping.capper(ph, sel, model.crystal_symmetry())
  assert (c.n_caps, c.charge) == (1, 0), (c.n_caps, c.charge)
  assert cap_names(c.cap()) == ['HG']

def test_whole_model_no_caps():
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "all")
  c = fragment_capping.capper(ph, sel, cs)
  assert c.n_caps == 0
  assert c.charge == charges_class(pdb_hierarchy=ph.deep_copy(), crystal_symmetry=cs).get_total_charge()
  r = c.cap()
  assert r.fragment_selection.all_eq(True)

def test_cap_moves_with_sites():
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  r0 = c.cap()
  shift = (1.0, 2.0, 3.0)
  r1 = c.cap(sites_cart=c.fragment_hierarchy.atoms().extract_xyz() + shift)
  assert cap_names(r0) == cap_names(r1)
  d = r1.pdb_hierarchy.atoms().extract_xyz() - r0.pdb_hierarchy.atoms().extract_xyz()
  assert (d - shift).max_distance(flex.vec3_double(d.size(), (0,0,0))) < 1.e-6  # rigid shift, caps too

def test_cap_twice():
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  n = c.fragment_hierarchy.atoms().size()
  a, b = c.cap(), c.cap()
  assert a.pdb_hierarchy.atoms().size() == b.pdb_hierarchy.atoms().size() == n + 2
  assert c.fragment_hierarchy.atoms().size() == n

def test_altloc_sorry():
  from mmtbx.hydrogens.tst_capping import pdb_ser_altloc, model_with_h
  model = model_with_h(pdb_ser_altloc)
  ph = model.get_hierarchy()
  try: fragment_capping.capper(ph, flex.bool(ph.atoms().size(), True), model.crystal_symmetry())
  except Sorry as e: assert 'alternate conformations' in str(e), str(e)
  else: raise AssertionError('no Sorry')

def test_not_qmready_sorry():
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "not (resseq 1 and (name H1 or name H2 or name H3))")
  try: fragment_capping.capper(ph.select(sel), flex.bool(sel.count(True), True), cs)
  except Sorry as e: assert 'qmready' in str(e), str(e)
  else: raise AssertionError('no Sorry')

def test_nucleic_cut():
  from mmtbx.hydrogens.tst_capping import pdb_dna, model_with_h
  model = model_with_h(pdb_dna)
  ph = model.get_hierarchy()
  resseqs = [rg.resseq.strip() for rg in ph.residue_groups()]
  first, last = resseqs[0], resseqs[-1]
  # O3' side alone: HO3'
  c = fragment_capping.capper(ph, ph.atom_selection_cache().selection("resseq %s" % first), model.crystal_symmetry())
  assert "HO3'" in cap_names(c.cap())
  # P side: hypervalent, Sorry naming P
  try: fragment_capping.capper(ph, ph.atom_selection_cache().selection("resseq %s" % last), model.crystal_symmetry())
  except Sorry as e: assert ' P ' in str(e) or "P " in str(e), str(e)
  else: raise AssertionError('no Sorry')
```
Imports: `os, iotbx.pdb, libtbx.load_env, flex, Sorry, from qrefine import fragment_capping, from qrefine.charges import charges_class`; `qr_unit_tests_data` as in tst_35. If `pdb_dna`'s first residue carries a 5' phosphate instead of a 5'-OH, pick the first residue whose P is absent; if no such residue exists, use `pdb_rna3` from `mmtbx/regression/tst_reduce2_capping.py` instead and note it in the test docstring.

- [ ] **Step 2: Run to verify it fails**

Run: `cd <scratch> && qrefine.python /Users/pafonine/Desktop/all/phenix/modules/qrefine/tests/unit/tst_38.py`
Expected: `ImportError: cannot import name 'fragment_capping'`.

- [ ] **Step 3: Move `model_from_hierarchy`** from `restraints.py` to `utils/hierarchy_utils.py` verbatim; `restraints.py` imports it from there (nothing else changes in restraints.py yet).

- [ ] **Step 4: Implement `fragment_capping.capper`**

`__init__`: `self.fragment_hierarchy = parent_hierarchy.select(selection)` + `atoms().reset_i_seq()`; Sorry `'Capping: the fragment has alternate conformations; QM fragments must be single-conformer'` if `fragment_hierarchy.altloc_indices()` has a non-blank key; `cut = capping.cut_bonds_by_distance(parent_hierarchy, selection)`; `self._model = model_from_hierarchy(self.fragment_hierarchy.deep_copy(), crystal_symmetry, cif_objects, log)`; `self.dangling = capping.find_dangling_atoms(self._model, cut_bonds=cut, log=log)`; Sorry `'Capping: no dictionary for %s in a QM fragment; supply its restraint CIF to qr.refine'` on `unknown_residues`; the qmready check: for every amino-acid residue (`iotbx.pdb.common_residue_names_get_class` in `common_amino_acid`, `modified_amino_acid`, `d_amino_acid`) backbone `N` with `len(dangling.bond_table.get(i_seq, ())) + cut.get(i_seq, 0) < 3`, Sorry `'Capping: %s has an unsatisfied valence that is not a fragment cut; the model is not QM-ready, run mmtbx.reduce2 qmready=True first'`; `r = capping.cap_model(self._model, dangling=self.dangling, log=log)`: this caps `self._model`'s hierarchy in place, which is fine because `self._model` is only needed for the records, the charge and its server; `self.charge = r.total_charge`, `self.n_caps = r.capped.n_added`, `self.lines = r.capped.lines`, `self._mon_lib_srv = self._model.get_mon_lib_srv()`. `self.fragment_hierarchy` is a separate copy and stays uncapped.

`cap(sites_cart=None)`: `h = self.fragment_hierarchy.deep_copy()`; if `sites_cart` is not None: `assert sites_cart.size() == h.atoms().size()`, `h.atoms().set_xyz(sites_cart)`; `m = mmtbx.model.manager(model_input=None, pdb_hierarchy=h, crystal_symmetry=self.crystal_symmetry, log=null_out())`; `capping.cap_model(m, dangling=self.dangling, calculate_charge=False, mon_lib_srv=self._mon_lib_srv)`; before capping set `atom.tmp = 1` on every atom of `h` (caps are created by `capping._new_h` with `tmp` 0); after capping `fragment_selection = h.atoms().extract_tmp_as_size_t() == 1` (a `flex.bool`). Return `group_args(pdb_hierarchy=h, fragment_selection=...)`.

- [ ] **Step 5: Run to verify it passes**

Run: Step 2 command. Expected: the test prints OK for every case and exits 0. If `test_pro_cut` or `test_truncated_side_chain` disagree only on a count, investigate with the capper's `lines` before touching the expected values (they come from probes on the same files).

- [ ] **Step 6: Register** `"$D/tests/unit/tst_38.py"` in `modules/qrefine/run_tests.py` `tst_list`.

- [ ] **Step 7: Report** changed/created files; no commit.

---

### Task 4: `fragment.py` and `restraints.from_qm` use the capper

**Files:**
- Modify: `modules/qrefine/fragment.py` (imports 13-16; `__init__` 72-139 to store `cif_objects`; `get_fragment_hierarchies_and_charges` 387-435; `get_fragment_extracts` 436-461; `get_qm_file_name_and_pdb_hierarchy` 462-485; `write_cluster_and_fragments_pdbs` 563-605)
- Modify: `modules/qrefine/restraints.py` (`from_qm.target_and_gradients` 558-626)

**Interfaces:**
- Consumes: `fragment_capping.capper` (Task 3).
- Produces: `fragment_extracts.fragment_cappers` (list, one `capper` per fragment, same index as `fragment_super_selections`); `get_qm_file_name_and_pdb_hierarchy(fragment_extracts, index)` returns `(qm_pdb_file, pdb_hierarchy, fragment_selection)` (was a 2-tuple).

- [ ] **Step 1: `fragment.py`**

`from . import completion` goes; `from . import fragment_capping` comes. `__init__`: `self.cif_objects = cif_objects`. `get_fragment_hierarchies_and_charges`: `self.fragment_cappers = []`; per fragment `c = fragment_capping.capper(self.pdb_hierarchy_super, fragment_super_selection, self.expansion.cs_box, cif_objects=self.cif_objects)`; `self.fragment_cappers.append(c)`; `self.fragment_capped_initial.append(c.cap().pdb_hierarchy)`; `fragment_charge = c.charge`; delete the `charges_class(...)` call and `raw_records`. `get_fragment_extracts`: add `fragment_cappers = self.fragment_cappers`. `get_qm_file_name_and_pdb_hierarchy`: after writing `qm_pdb_file`, `r = fragment_extracts.fragment_cappers[index].cap(sites_cart=fragment_hierarchy.atoms().extract_xyz())`; write `r.pdb_hierarchy` to `complete_qm_pdb_file`; `return os.path.abspath(complete_qm_pdb_file), r.pdb_hierarchy, r.fragment_selection`. `write_cluster_and_fragments_pdbs`: `capped_hierarchy = F.fragment_cappers[index].cap(sites_cart=index_frag.atoms().extract_xyz()).pdb_hierarchy`; drop the `expansion_file` local. `command_line/fragmentation.py` ignores the return value: no change.

- [ ] **Step 2: `restraints.from_qm.target_and_gradients`**

Clustering branch: `qm_pdb_file, ph, fragment_selection = get_qm_file_name_and_pdb_hierarchy(...)`. Non-clustering branch: `fragment_selection = flex.bool(self.system_size, True)` (rename the existing `selection =` line). After the engine call: `gradients = flex.vec3_double(ase_gradients).select(fragment_selection)`, then the existing `gradients_scale` multiplication (assert `gradients.size() == len(gradients_scale)`). Remove the `selection` parameter's use for slicing (keep the signature).

- [ ] **Step 3: Run the clustering QM tests**

Run from scratch, each: `qrefine.python .../tests/unit/tst_NN.py` for NN in 08, 09, 10, 12, 35, 36 (mopac is installed as `phenix.mopac`). Expected: each exits 0; `tst_35` prints `get_grad_wdelta (cluster vs no cluster): <d>` with `d < 10`. Record `d` in the report (the spec expects a small change from the new cap geometry).

- [ ] **Step 4: Run `tst_38`** again (unchanged expectation) and report changed files; no commit.

---

### Task 5: `restraints.from_expansion` uses the capper; `update(charge=)`

**Files:**
- Modify: `modules/qrefine/restraints.py` (`restraints.update` 68-107; `h_diff_sel` 109-124 deleted; `from_expansion` 126-239)

**Interfaces:**
- Consumes: `fragment_capping.capper`; `expand.SCM.super_cell_hierarchy`, `expand.SCM.keep_selection` (cctbx `super_cell.manager`), `expand.ph_super_sphere`.
- Produces: `restraints.update(self, pdb_hierarchy, crystal_symmetry, charge=None)`: `from_qm(charge=self.params.quantum.charge if self.params.quantum.charge is not None else charge, ...)`.

- [ ] **Step 1: `from_expansion.__init__`**

QM branch: `self.capper = fragment_capping.capper(self.expansion.SCM.super_cell_hierarchy, self.expansion.SCM.keep_selection, self.expansion.cs_box, cif_objects=restraints_source.cif_objects)`; `assert self.capper.fragment_hierarchy.is_similar_hierarchy(self.expansion.ph_super_sphere)` (same atoms in the same order; `is_similar_hierarchy` ignores coordinates); `r = self.capper.cap(sites_cart=self.expansion.ph_super_sphere.atoms().extract_xyz())` (always the sphere's current coordinates: the super cell is in the crystal frame while the initial sphere may be boxed, and `SCM.update` re-selects the sphere from the super cell); `self.pdb_hierarchy_super_completed = r.pdb_hierarchy`; `n = self.model.get_hierarchy().atoms().size()`; `self.selection = flex.bool(r.pdb_hierarchy.atoms().size(), r.fragment_selection.iselection()[:n])`; keep the master assertion as `assert self.pdb_hierarchy_super_completed.select(self.selection).is_similar_hierarchy(self.model.get_hierarchy())`. Non-QM branch: `self.capper = None`, completed = sphere deep copy, `self.selection` = first `n` atoms as today. Delete the `h_diff_sel` call and function, and the `selection = ~master_selection` block. `self.restraints_source.update(..., charge=self.capper.charge if self.capper else None)`.

- [ ] **Step 2: `from_expansion._update`**

After `self.expansion.update(...)`: `xyz = self.expansion.ph_super_sphere.atoms().extract_xyz()`; if `self.capper`: `xyz = self.capper.cap(sites_cart=xyz).pdb_hierarchy.atoms().extract_xyz()`; `assert xyz.size() == self.pdb_hierarchy_super_completed.atoms().size()`; `self.pdb_hierarchy_super_completed.atoms().set_xyz(xyz)`.

- [ ] **Step 3: Run the expansion tests**

Run from scratch: `tst_04, tst_24, tst_25, tst_26, tst_29, tst_32, tst_33, tst_34` (all exercise expansion; those with `restraints=qm` use mopac). Expected: each exits 0.

- [ ] **Step 4: Report** changed files; no commit.

---

### Task 6: `qr.finalise` wrapper; `tst_20` and `tst_18` rewritten

**Files:**
- Modify: `modules/qrefine/command_line/finalise.py` (replace whole file)
- Modify: `modules/qrefine/tests/unit/tst_20.py:560-610`, `modules/qrefine/tests/unit/tst_18.py:327-358`

**Interfaces:**
- Produces: `qr.finalise <model> [reduce2 options]` == `mmtbx.reduce2 <model> [options] qmready=True`.

- [ ] **Step 1: Rewrite `tst_20.run`**: write `test_remove_hg_from_zn_cys.pdb` as now; `rc = easy_run.go('qr.finalise %s' % fn)`; `assert rc.return_code == 0, rc.stdout_lines + rc.stderr_lines`; `outs = glob.glob('test_remove_hg_from_zn_cysH_charge*.pdb'); assert len(outs) == 1, outs`; read it; keep the six `' HG  CYS A ...' not in lines` assertions and `' HE  ARG A  10' in lines`; assert `'qr.finalise now runs mmtbx.reduce2 qmready=True' in '\n'.join(rc.stdout_lines)`. Delete the disabled `qr.charges` block. Docstring: "Exercise qr.finalise (runs mmtbx.reduce2 qmready=True)".

- [ ] **Step 2: Run `tst_20`** to see it fail: `ModuleNotFoundError`/`Sorry` from the old wrapper (expected failure).

- [ ] **Step 3: Replace `command_line/finalise.py`**

```python
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME qr.finalise
import sys
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

note = """\
qr.finalise now runs mmtbx.reduce2 qmready=True. The former qr.finalise
options (action, calculate_charge, reduce, options.neutron,
hydrogen_atom_occupancies, append_to_end_of_model, skip_validation,
stop_if_poor_model) are gone: the total charge is always calculated and put
into the output file name (<prefix>H_charge<q>.pdb). All other arguments are
passed to mmtbx.reduce2 (for example capping=True); run mmtbx.reduce2
directly for its full option list.
"""

def run(args):
  print(note)
  return run_program(program_class=reduce2.Program,
                     args=list(args) + ['qmready=True'])

if __name__ == '__main__':
  run(sys.argv[1:])
```

- [ ] **Step 4: Run `tst_20`**. Expected: exit 0. If reduce2 stops with a Sorry on this input (it already carries H), strip the H from `pdb_lines` in the test before writing the file (keep heavy atoms and Zn) and rerun; record the reason in the test docstring.

- [ ] **Step 5: Rewrite the capping half of `tst_18`** (lines 342-355): remove the `completion.py` command and the `mmtbx.pdb_interpretation` check; instead

```python
  pi = iotbx.pdb.input(file_name='expansion.pdb')
  ph = pi.construct_hierarchy()
  sel = ph.atom_selection_cache().selection("chain A and resseq 89:91")
  c = fragment_capping.capper(ph, sel, pi.crystal_symmetry())
  r = c.cap()
  names = sorted(a.name.strip() for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection))
  assert names == ['H2', 'HOXT'], names
  assert c.charge == 0, c.charge
```
Also remove the writing of `test_alt_loc_original.pdb` if nothing else in the test uses it. Run `tst_18`: exit 0.

- [ ] **Step 6: Report** changed files; no commit.

---

### Task 7: Removals, `use_reduce`, dead helpers, docs

**Files:**
- Delete: `modules/qrefine/finalise.py`, `modules/qrefine/completion.py`, `modules/qrefine/tests/unit/batch_run_finalise.py`, `modules/qrefine/tests/unit/broken_to_fix/tst_48.py`, `modules/qrefine/tests/unit/broken_to_fix/tst_49.py`, `modules/qrefine/tests/unit/babel_pdbs/` (directory), `modules/qrefine/tests/regression/test_reg_01_finalise.py`, `modules/qrefine/tests/regression/datasets/finalise/` (directory)
- Modify: `modules/qrefine/tests/regression/regression_tests.py` (drop the import and entry), `modules/qrefine/tests/regression/test_reg_00_base.py` (drop `qr_reg_data_finalise` and the finalise loop), `modules/qrefine/qr.py:33,347,381`, `modules/qrefine/utils/hierarchy_utils.py`, `modules/qrefine/charges.py` (`c_capping`, lines 290-292 and 593-595), `modules/qrefine/rst_files/qr/finalise.txt`, `modules/qrefine/rst_files/qr.txt:60-70`, `modules/qrefine/rst_files/qr/tests.txt`, `modules/qrefine/CHANGES.md`

- [ ] **Step 1: Delete the files and directories listed**, then verify nothing references them:

Run: `cd modules/qrefine && grep -rn "completion\b\|finalise\|batch_run_finalise\|babel_pdbs\|qr_reg_data_finalise" --include=*.py . | grep -v "command_line/finalise.py\|tests/unit/tst_20.py"`
Expected: no output (the wrapper and tst_20 are the only legitimate mentions).

- [ ] **Step 2: Remove `use_reduce`** from `qr.py` (PHIL line 33, `self.params.use_reduce=False` at 347, the print at 381). Run: `grep -rn use_reduce modules/qrefine --include=*.py --include=*.txt`. Expected: no output.

- [ ] **Step 3: Dead `hierarchy_utils` functions.** For each of `write_hierarchy`, `add_hydrogens_using_reduce`, `add_hydrogens_using_ReadySet`, `smart_add_atoms`, `get_pdb_interpretation_params`, `remove_atom_from_chain`, `is_n_terminal_residue`: run `grep -rn "<name>" modules/qrefine --include=*.py | grep -v utils/hierarchy_utils.py`; delete the function if the grep is empty, keep it otherwise and list it in the report. Remove imports that become unused in `hierarchy_utils.py`.

- [ ] **Step 4: `charges.py`**: in both `c_capping` definitions make the check `_terminal(atom_names, [' HC ']) or _terminal(atom_names, ['HOXT'])`.

- [ ] **Step 5: Docs.** `rst_files/qr/finalise.txt`: title "Preparing a QM-ready model", text: `qr.finalise model.pdb` runs `mmtbx.reduce2 model.pdb qmready=True` (adds H including water, OXT and NH3+ at chain ends and breaks, neutron X-H distances, checks for clashes, calculates the total charge, writes `modelH_charge<q>.pdb`); extra arguments are passed to reduce2 (`capping=True` caps truncated side chains for QM, but such a model cannot be refined); replace `{{phil:qrefine.command_line.finalise}}` with `{{phil:mmtbx.programs.reduce2}}`. `rst_files/qr.txt` step 2: output `1yjp_completeH_charge<q>.pdb`... use the real name: `1yjpH_charge<q>.pdb`, and step 3 accordingly. `rst_files/qr/tests.txt`: update the `tst_18`/`tst_20` rows if present and add `tst_38 | Exercise fragment capping (fragment_capping.capper) |`. `CHANGES.md` under `[Unreleased]`: `### Changed` "- qr.finalise runs mmtbx.reduce2 qmready=True; QM fragment capping uses mmtbx.hydrogens.capping (fragment_capping.capper); caps sit inside residues, cut carbonyl cap is named HOXT"; `### Removed` "- finalise.py, completion.py, the use_reduce parameter". Leave `doc/*.html` alone (generated) and say so in the report.

- [ ] **Step 6: Import check.** Run: `qrefine.python -c "import qrefine.fragment, qrefine.restraints, qrefine.charges, qrefine.qr, qrefine.utils.hierarchy_utils, qrefine.command_line.finalise"`. Expected: no output.

- [ ] **Step 7: Report** deleted and changed files; no commit.

---

### Task 8: Full verification

- [ ] **Step 1: mmtbx tests.** From scratch: `libtbx.python .../mmtbx/hydrogens/tst_capping.py`, `tst_reduce2_capping.py`, `tst_reduce2_qmready_charge.py`, `tst_reduce2_qmready_clash.py`. Expected: all OK.

- [ ] **Step 2: qrefine suite.** From a fresh scratch directory: `qrefine.python -c "from qrefine.tests.unit import run_tests; run_tests.run(nproc=4)"` (the file's `__main__` block is broken; `run()` is the entry). Expected: every `qrefine_tst_NN:  OK` line for tst_00 to tst_38, no `FAILED`. Record total time.

- [ ] **Step 3: Report** the test output summary (per-test OK lines, tst_35's wdelta), the complete list of changed/added/deleted files in both repos, and anything left out, for the user to review and commit.
