# Capping H Recognition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Models capped by `mmtbx.reduce2 qmready=true capping=true` load in every program that builds restraints, their caps are identifiable by chemistry alone, selectable as `capping_h`, carry occupancy 0, and stay out of geometry statistics.

**Architecture:** One function, `capping.find_capping_h`, implements the definition of a cap (an H the dictionary cannot attach, on a parent with a dangling valence that accommodates it). pdb_interpretation calls it after proxies are built and adds bond and angle proxies with a new origin id `capping` plus the generic H nonbonded type. The `capping_h` selection keyword and the occupancy guards call the same function. reduce2 places caps at occupancy 0 and with names that never collide with the dictionary.

**Tech Stack:** cctbx_project (iotbx.pdb.atom_selection, mmtbx.hydrogens.capping, mmtbx.monomer_library.pdb_interpretation, cctbx.geometry_restraints), qrefine tests, phenix_regression.

**Spec:** `qrefine/docs/superpowers/specs/2026-11-04-capping-h-recognition-design.md`

## Global Constraints

- Repositories: `modules/cctbx_project` (code and most tests), `modules/qrefine` (one test), `modules/phenix_regression` (one test). The user commits; leave changes uncommitted. Never run `rm -rf`; scratch goes in a fresh `rt/run_$(date +%s)` folder.
- Python to run anything: `/Users/pafonine/Desktop/all/phenix/build/bin/libtbx.python`.
- Definition of a cap: spec section 3, conditions 1-4, exactly. Counts and chemistry, never coordinates.
- Bonding distance: sum of covalent radii + `GUESS_BOND_TOLERANCE` (0.4 A), the existing constant.
- Cap restraints: bond ideal `capping.xh_length(element, coordination, use_neutron_distances)`, sigma 0.02 A; angles by `capping.steric_number` (109.5, 120, 180 degrees), sigma 3.0 degrees; nonbonded type `H`; origin id `capping`; no dihedral, planarity or chirality proxies.
- Caps get occupancy 0; B, segid, hetero flag stay the parent's.
- Selection keyword name: `capping_h`.
- Speed: work proportional to the number of candidate residues; one pass over H atoms with set lookups; no pair search over the whole model.
- Every qrefine test writes only into its own folder `qrefine_tst_XX` (keep the `os.makedirs`/`os.chdir` lines in `run()`).
- Existing tests must keep passing: `mmtbx/hydrogens/tst_capping.py`, `mmtbx/regression/tst_reduce2_capping.py`, `iotbx/pdb/tst_atom_selection.py`, qrefine `tst_38.py`.

## Review Focus

1. A capped model with two MODEL records: caps of MODEL 2 must find their parents in MODEL 2 only. Test `test_two_models` in Task 3.
2. A cap on a blank-altloc parent whose residue has altloc atom groups (the `pdb_ser_altloc` fixture): the cap is in the blank group, candidates per parent are counted once per conformer, not twice. Test `test_altloc_parent` in Task 3.
3. A model with caps but processed with `stop_for_unknowns=False`: caps must still get proxies, not be silently skipped (reduce2 itself builds such a manager). Test `test_process_no_stop` in Task 4.
4. A hierarchy whose i_seqs are stale (not reset) given to `find_capping_h`: selection indices must follow `hierarchy.atoms()` order. Test `test_stale_i_seq` in Task 3.
5. `phenix.pdbtools remove=capping_h` on a model without any H: the keyword must return an empty selection quickly, not raise. Test `test_keyword_no_h` in Task 5.

---

### Task 1: Origin id `capping` and the test file

**Files:**
- Modify: `cctbx_project/cctbx/geometry_restraints/auto_linking_types.py:150-213`
- Create: `cctbx_project/mmtbx/hydrogens/tst_capping_h.py`
- Modify: `cctbx_project/mmtbx/run_tests.py:197` (register after `tst_capping.py`)

**Interfaces:**
- Produces: `linking_class().get_origin_id('capping')` (bonds and angles, internals `[0, 1]`); `capping` in `not_covalent` so `iterate_covalent()` skips it.

- [ ] **Step 1: Write the failing test**

`tst_capping_h.py` starts with the imports of `tst_capping.py` plus `from mmtbx.hydrogens import tst_capping as fixtures` and ends with `run()` calling each test and printing OK, like `tst_capping.py`.

```python
def test_origin_id():
  from cctbx.geometry_restraints.linking_class import linking_class
  from cctbx.geometry_restraints.auto_linking_types import iterate_covalent
  ids = linking_class()
  i = ids.get_origin_id('capping')
  assert i > 0
  assert ids.data['capping'][3] == ['Capping H', 'Capping H']
  assert 'capping' not in list(iterate_covalent())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd modules/cctbx_project && libtbx.python mmtbx/hydrogens/tst_capping_h.py`
Expected: AssertionError from `get_origin_id` ("linking origin id not found").

- [ ] **Step 3: Add the entry**

In the main `for link_info in [...]` list of `auto_linking_types.py`, after `'reference hydrogen bonds'` and before the closing `]:`, add the spec section 6 entry verbatim (short description `'capping'`, headers `['Capping H', 'Capping H']`, internals `[0, 1]`). Add `'capping'` to `not_covalent` under a comment `# necessary: QM caps, excluded from geometry statistics`.

- [ ] **Step 4: Run test to verify it passes**

Run: the command of Step 2. Expected: OK.

- [ ] **Step 5: Register and run the neighbours**

Add `"$D/hydrogens/tst_capping_h.py",` after line 197 of `mmtbx/run_tests.py`. Run `libtbx.python cctbx/geometry_restraints/tst_ext.py` and `libtbx.python mmtbx/regression/tst_reduce2_capping.py`. Expected: both OK (ids are looked up by name, so the shift of the standard cif link ids is harmless; if a test compares .geo text with hard-coded ids, report it instead of editing the expectation).

---

### Task 2: capping.py groundwork: shared server, valence budget helper, occupancy 0, dictionary-wide names

**Files:**
- Modify: `cctbx_project/mmtbx/hydrogens/capping.py` (`find_dangling_atoms` 413-531, `cap_names` 846-880, `_new_h` 932-946, module docstring Limitations 164-180)
- Modify: `cctbx_project/mmtbx/hydrogens/tst_capping.py:1140` (occupancy assertion)
- Test: `cctbx_project/mmtbx/hydrogens/tst_capping_h.py`

**Interfaces:**
- Produces:
  - `shared_mon_lib_srv() -> mmtbx.monomer_library.server.server`, module-level cache, built on first call.
  - `_budget(atom, own, order, present, neighbours, rg) -> group_args(n_fill, missing, h_names)`: the body of the per-atom loop of `find_dangling_atoms` from `missing = ...` through `h_names = ...` (lines 473-500), unchanged in behaviour; `find_dangling_atoms` calls it.
  - `_dictionary_names(resname, mon_lib_srv, cache) -> None | group_args(names: set, h_partner: {h_name: heavy_name})` from `_dictionary_atoms` and its bonds (an H's partner is its one bonded heavy atom).
  - `_new_h` sets occupancy 0.
  - `cap_names(atom_group, h_names, missing_names, n_fill, parent=None, dictionary_names=())`: candidates in `dictionary_names` are skipped unless they are in `h_names`; `cap_dangling_atoms` passes the residue's dictionary names.

- [ ] **Step 1: Write the failing tests**

```python
def test_shared_server():
  s1 = capping.shared_mon_lib_srv()
  assert s1 is capping.shared_mon_lib_srv()

def test_cap_occupancy_zero():
  model, result = fixtures.cap(fixtures.pdb_arg)
  cb = fixtures.atom_named(model, 3, 'CB')
  for name in fixtures.h_names_on(model, 3, 'CB'):
    h = fixtures.atom_named(model, 3, name)
    assert h.occ == 0.0, (name, h.occ)
    assert (h.b, h.segid, h.hetero) == (cb.b, cb.segid, cb.hetero)
  assert cb.occ == 1.0

def test_cap_name_avoids_dictionary():
  # DC 4 of pdb_dna without O3': the cap on C3' must not be HO3' (the
  # dictionary H of the absent O3'); the fallback is HC1
  model, result = fixtures.cap(fixtures.without(fixtures.pdb_dna, 4, ["O3'"]))
  names = fixtures.h_names_on(model, 4, "C3'")
  assert "HO3'" not in names, names
  assert 'HC1' in names, names

def test_cap_name_keeps_free_dictionary_h():
  model, result = fixtures.cap(fixtures.pdb_arg)
  assert fixtures.h_names_on(model, 3, 'CB') == ['HB2', 'HB3', 'HCG']

def test_budget_refactor_unchanged():
  # find_dangling_atoms before and after the refactor: same records
  r = fixtures.find(fixtures.pdb_lys)
  by = fixtures.records_by_name(r)
  assert list(by) == [('LYS', '2', '', 'CD')]
  assert (by[('LYS', '2', '', 'CD')].n_fill,
          by[('LYS', '2', '', 'CD')].missing_names) == (2, ['CE'])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `libtbx.python mmtbx/hydrogens/tst_capping_h.py`
Expected: `test_shared_server` AttributeError; `test_cap_occupancy_zero` assertion `h.occ == 0.0`; `test_cap_name_avoids_dictionary` assertion on HO3'.

- [ ] **Step 3: Implement**

In `capping.py`: `shared_mon_lib_srv` with a module global `_SHARED_MON_LIB_SRV = None`; `_dictionary_names`; extract `_budget`; `_new_h` with `atom.set_occ(new_occ=0.0)`; `cap_names` with the `dictionary_names` parameter; in `cap_dangling_atoms` build the residue's dictionary names once per resname (cache dict) and pass them. Update the docstrings of `cap_names` (rule: "never a dictionary name of the residue unless it is a free H name of the parent") and the module docstring: Names of the caps rule, Building paragraph (occupancy 0), Limitations (replace the first bullet with: processed by pdb_interpretation through `find_capping_h`, see its docstring; add the old-file HG case of spec section 8).

- [ ] **Step 4: Fix the existing assertion and run**

`tst_capping.py:1140`: expect `h.occ == 0.0` and the other three fields equal to CB's. Run `libtbx.python mmtbx/hydrogens/tst_capping.py`, `libtbx.python mmtbx/hydrogens/tst_capping_h.py`, `libtbx.python mmtbx/regression/tst_reduce2_capping.py` (its line 543 area mentions occupancies of altloc caps; update any `occ` expectation for caps to 0). Expected: OK, OK, OK.

---

### Task 3: `find_capping_h`

**Files:**
- Modify: `cctbx_project/mmtbx/hydrogens/capping.py` (new public function after `extra_sigma_bonds`, plus two private helpers)
- Test: `cctbx_project/mmtbx/hydrogens/tst_capping_h.py`

**Interfaces:**
- Consumes: `_dictionary_names`, `_dictionary_partners`, `_conformer_names`, `_present_names`, `_conformer_neighbours`, `_budget`, `steric_number`, `reduce_hydrogen._compatible`, `reduce_hydrogen._covalent_radius`, `atom_property().is_metal`.
- Produces:

```python
def find_capping_h(hierarchy, mon_lib_srv, bond_table=None, log=None):
  '''
  -> group_args(selection=flex.bool over hierarchy.atoms() (i_seq reset
     on a working copy, never on the caller's hierarchy),
     records=[group_args(i_seq, parent_i_seq, coordination, steric)],
     rejected=[(atom label, reason)])
  '''
def _neighbours_by_distance(atom, atoms, sites, is_metal):
  '''Bonded non-metal neighbours of atom by covalent radii + 0.4 A among
  its residue group and the previous and next residue groups of its chain,
  one per conformer (_conformer_neighbours).'''
def capping_h_selection(hierarchy, mon_lib_srv=None) -> flex.bool:
  '''find_capping_h(...).selection with shared_mon_lib_srv() when None.'''
```

Algorithm of `find_capping_h` (the one place that decides):
1. `atoms = hierarchy.atoms()`; indices are positions in this array (build `i_seq` lookup via `atoms.reset_i_seq()` on a `hierarchy.deep_copy()` if any `atoms.extract_i_seq()` is not `0..n-1`, else use the hierarchy as is).
2. Per atom group with a dictionary (`_dictionary_names` not None): candidates = H/D atoms whose stripped name is not in `names`, or is in `h_partner` and its partner is not in `_present_names(atom, _conformer_names(rg))`.
3. Parent of a candidate: the closest non-H, non-metal atom of the same residue group that `_compatible` with it and within bonding distance, else `rejected` with `'no heavy atom of its residue within bonding distance'`.
4. Group candidates by `(parent i_seq, candidate altloc)`. For each group: `own = partners.get(parent name)`; if None, reject every candidate with `'parent <name> is not in the dictionary'`. `neighbours` from `_conformer_neighbours(parent, bond_table[parent], ...)` when `bond_table` is given, else `_neighbours_by_distance`; drop the group's candidates from `neighbours`. `b = _budget(...)`. Reject with `'parent <name> misses no heavy atom'` when `b.missing` is empty, with `'%d H on <name>, its valence allows %d'` when `len(group) > b.n_fill`. Otherwise every candidate is a cap: `coordination = len(neighbours) + len(group)`, `steric = steric_number(element, coordination)`.
5. `log`: one line per rejected candidate when log is not None.

- [ ] **Step 1: Write the failing tests**

```python
def capped_hierarchy(pdb_str):
  model, result = fixtures.cap(pdb_str)
  return model, result, model.get_hierarchy()

def added_i_seqs(hierarchy, result):
  '''
  i_seqs of the caps cap_dangling_atoms placed: result.added is
  [(atom_key of the parent, [cap names])] (capping.py:1088); the caps are
  the atoms of the parent's atom group with those names.
  '''
  want = set()
  for key, names in result.added:
    model_id, chain_id, resseq, icode, resname, altloc, parent_name = key
    for a in hierarchy.atoms():
      ag = a.parent(); rg = ag.parent()
      if ((rg.parent().id.strip(), rg.resseq.strip(), ag.altloc.strip())
          == (chain_id, resseq, altloc) and a.name.strip() in names):
        want.add(a.i_seq)
  return want

def test_find_matches_added():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  added = added_i_seqs(h, result)
  assert set(r.selection.iselection()) == added, (r.selection.iselection(), added)
  assert not r.rejected, r.rejected
  # HB2 and HB3 on the truncated CB are dictionary atoms: not caps
  names = set(h.atoms()[i].name.strip() for i in r.selection.iselection())
  assert names == {'HCG'}, names

def test_find_lys_two_caps():
  model, result, h = capped_hierarchy(fixtures.pdb_lys)
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  rec = [x for x in r.records if h.atoms()[x.parent_i_seq].name.strip() == 'CD']
  assert len(rec) == 2 and all(x.coordination == 4 and x.steric == 4 for x in rec)

def test_find_old_ho3_name():
  # a file written before Task 2: the cap on C3' named HO3' (dictionary H of
  # the absent O3') is recognized by condition 2
  model, result, h = capped_hierarchy(fixtures.without(fixtures.pdb_dna, 4, ["O3'"]))
  cap = [a for a in h.atoms() if a.name.strip() == 'HC1'][0]
  cap.set_name(new_name="HO3'")
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection[cap.i_seq]

def test_find_rejects_wrong_name_complete_residue():
  model = fixtures.model_with_h(fixtures.pdb_ala_ala)
  h = model.get_hierarchy()
  a = fixtures.atom_named(model, 2, 'HB1'); a.set_name(new_name=' HX1')
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection.count(True) == 0
  assert r.rejected and 'misses no heavy atom' in r.rejected[0][1], r.rejected

def test_find_rejects_excess():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  cb = fixtures.atom_named(model, 3, 'CB')
  extra = fixtures.atom_named(model, 3, 'HCG').detached_copy()
  extra.set_name(new_name=' HX9'); extra.set_xyz(new_xyz=(cb.xyz[0]+1.0, cb.xyz[1], cb.xyz[2]))
  cb.parent().append_atom(extra); h.atoms().reset_i_seq()
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert r.selection.count(True) == 0
  assert any('valence allows' in reason for label, reason in r.rejected), r.rejected

def test_find_with_bond_table_equals_without():
  model, result, h = capped_hierarchy(fixtures.pdb_disulfide)
  model.process(make_restraints=True)   # Task 4 makes this pass; until then skip the table branch
  a = capping.find_capping_h(h, model.get_mon_lib_srv())
  b = capping.find_capping_h(h, model.get_mon_lib_srv(), bond_table=capping.bond_table(model))
  assert list(a.selection) == list(b.selection)

def test_two_models():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  two = iotbx.pdb.input(lines=(model.model_as_pdb().replace('END', '') +
    'ENDMDL\n' + model.model_as_pdb()).split('\n'), source_info=None).construct_hierarchy()
  r = capping.find_capping_h(two, model.get_mon_lib_srv())
  assert r.selection.count(True) == 2 and not r.rejected

def test_altloc_parent():
  model, result, h = capped_hierarchy(fixtures.pdb_ser_altloc)
  r = capping.find_capping_h(h, model.get_mon_lib_srv())
  assert set(r.selection.iselection()) == added_i_seqs(h, result)
  assert not r.rejected, r.rejected

def test_stale_i_seq():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  expected = list(capping.find_capping_h(h, model.get_mon_lib_srv()).selection)
  for a in h.atoms(): a.set_i_seq(a.i_seq + 100)
  assert list(capping.find_capping_h(h, model.get_mon_lib_srv()).selection) == expected

def test_speed():
  import time
  # 20 000 atoms without caps: 40 copies of pdb_lys_full chained as separate chains
  h = big_hierarchy_without_caps()   # helper: concatenate fixtures.pdb_lys_full with chain ids
  srv = capping.shared_mon_lib_srv()
  t0 = time.time(); r = capping.find_capping_h(h, srv); dt = time.time() - t0
  assert r.selection.count(True) == 0
  assert dt < 2.0, dt
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `libtbx.python mmtbx/hydrogens/tst_capping_h.py`
Expected: AttributeError `find_capping_h`.

- [ ] **Step 3: Implement `find_capping_h`, `_neighbours_by_distance`, `capping_h_selection`**

Docstring of `find_capping_h` states the definition (spec section 3) and the permissive neighbour counting without a bond table (spec section 4).

- [ ] **Step 4: Run tests to verify they pass**

Run: the command of Step 2 (with `test_find_with_bond_table_equals_without` left out of `run()` until Task 4). Expected: OK. Also `libtbx.python mmtbx/hydrogens/tst_capping.py` OK.

---

### Task 4: pdb_interpretation restrains caps

**Files:**
- Modify: `cctbx_project/mmtbx/hydrogens/capping.py` (new `add_capping_h_restraints`)
- Modify: `cctbx_project/mmtbx/monomer_library/pdb_interpretation.py:6225-6227` (in `process.__init__`, before `# attempt to fix pH problems`)
- Test: `cctbx_project/mmtbx/hydrogens/tst_capping_h.py`

**Interfaces:**
- Consumes: `find_capping_h`, `xh_length`, `steric_number`, `linking_class().get_origin_id('capping')`, registries `geometry_proxy_registries.bond_simple` / `.angle` (`process(source_info=, proxy=)`, `__init__.py:179, 264`), `nonbonded_energy_type_registry.assign_directly(i_seq, symbol)`.
- Produces:

```python
class _capping_source_info(object):
  def labels(self): return ['capping H']
  def n_expected_atoms(self): return 0

def add_capping_h_restraints(hierarchy, mon_lib_srv, geometry_proxy_registries,
                             nonbonded_energy_type_registry,
                             use_neutron_distances, log=None):
  '''
  find_capping_h with the bond table of geometry_proxy_registries.bond_simple
  (simple proxies only), then per cap: the bond proxy, angle proxies to
  every bonded neighbour of the parent (other caps included), nonbonded
  type 'H'. Returns group_args(i_seqs=flex.size_t of the caps,
  rejected=[...]). Logs 'Capping H: N recognized' and the rejected lines.
  '''
```

In `process.__init__`, right after `self._xray_structure = None` and before the pH block:

```python
    from mmtbx.hydrogens import capping
    acp = self.all_chain_proxies
    self.capping_h = capping.add_capping_h_restraints(
      hierarchy                      = acp.pdb_hierarchy,
      mon_lib_srv                    = mon_lib_srv,
      geometry_proxy_registries      = acp.geometry_proxy_registries,
      nonbonded_energy_type_registry = acp.nonbonded_energy_type_registry,
      use_neutron_distances          = acp.params.use_neutron_distances,
      log                            = log)
```

Run it only when `acp.nonbonded_energy_type_registry.n_unknown_type_symbols()` or the hierarchy has an H without a bond proxy (cheap check: collect i_seqs of H atoms absent from `bond_simple.table`); otherwise `self.capping_h = None`. This keeps models without caps at zero extra cost.

- [ ] **Step 1: Write the failing tests**

```python
def processed(pdb_str, **kw):
  model, result, h = capped_hierarchy(pdb_str)
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out(), **kw)
  m.process(make_restraints=True)
  return m, result

def cap_proxies(m):
  from cctbx.geometry_restraints.linking_class import linking_class
  oid = linking_class().get_origin_id('capping')
  grm = m.get_restraints_manager().geometry
  bonds = [p for p in grm.get_all_bond_proxies(sites_cart=m.get_sites_cart())[0]
           if p.origin_id == oid]
  angles = [p for p in grm.get_all_angle_proxies() if p.origin_id == oid]
  return bonds, angles

def test_process_capped_model():
  m, result = processed(fixtures.pdb_arg)      # default stop_for_unknowns=True
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1
  assert len(angles) == 3                      # HCG-CB-CA, HCG-CB-HB2, HCG-CB-HB3
  reg = m.get_restraints_manager().geometry   # no unknown nonbonded types:
  assert m.selection('capping_h').count(True) == 1   # Task 5; until then use capping.capping_h_selection
  assert abs(bonds[0].distance_ideal - capping.xh_length('C', 4, True)) < 1e-6

def test_process_lys_angles():
  m, result = processed(fixtures.pdb_lys)
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 2
  # each cap: CG, HD2, HD3, the other cap -> 4 angles, 8 in all, 109.5
  assert len(angles) == 8 and all(abs(a.angle_ideal - 109.5) < 1e-6 for a in angles)

def test_statistics_exclude_caps():
  m, result = processed(fixtures.pdb_arg)
  uncapped = fixtures.model_with_h(fixtures.pdb_arg)
  s_c = m.geometry_statistics().result()
  s_u = uncapped.geometry_statistics().result()
  assert s_c.bond.n == s_u.bond.n and s_c.angle.n == s_u.angle.n
  assert abs(s_c.bond.mean - s_u.bond.mean) < 1e-6
  geo = m.restraints_as_geo()
  assert 'Capping H' in geo

def test_idempotent():
  m, result = processed(fixtures.pdb_arg)
  r = capping.find_dangling_atoms(m)
  assert not r.records, [x.key for x in r.records]

def test_still_stops_wrong_name():
  model = fixtures.model_with_h(fixtures.pdb_ala_ala)
  fixtures.atom_named(model, 2, 'HB1').set_name(new_name=' HX1')
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out())
  fixtures.assert_sorry(lambda: m.process(make_restraints=True), 'unknown')

def test_process_no_stop():
  m, result = processed(fixtures.pdb_arg, stop_for_unknowns=False)
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1

def test_neutron_and_deuterium():
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  cap = fixtures.atom_named(model, 3, 'HCG'); cap.set_element('D'); cap.set_name(new_name=' DCG')
  text = model.model_as_pdb()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(lines=text.split('\n'),
    source_info=None), log=null_out())
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.use_neutron_distances = True
  m.process(make_restraints=True, pdb_interpretation_params=p)
  bonds, angles = cap_proxies(m)
  assert len(bonds) == 1 and abs(bonds[0].distance_ideal - capping.xh_length('C', 4, True)) < 1e-6

def test_riding_minimization():
  m, result = processed(fixtures.pdb_lys)
  m.setup_riding_h_manager()
  m.idealize_h_riding()
  from libtbx import easy_run
  open('capped_lys.pdb', 'w').write(m.model_as_pdb())
  assert easy_run.call('mmtbx.geometry_minimization capped_lys.pdb '
    'output_file_name_prefix=capped_lys_min > capped_lys.log') == 0
  out = iotbx.pdb.input(file_name='capped_lys_min.pdb').construct_hierarchy()
  # caps still bonded to CD at the X-H length
  cd = [a for a in out.atoms() if a.name.strip() == 'CD'][0]
  for a in out.atoms():
    if a.name.strip() in ('HCE', 'HC1', 'HC2') or a.name.strip().startswith('HC'):
      d = (matrix.col(a.xyz) - matrix.col(cd.xyz)).length()
      assert abs(d - capping.xh_length('C', 4, True)) < 0.05, (a.name, d)
```

Use the `restraints_as_geo` and `geometry_statistics` method names as they exist on `mmtbx.model.manager` (`grep -n "def restraints_as_geo\|def geometry_statistics" mmtbx/model/model.py`); adjust the attribute path to `bond.n` / `bond.mean` per `mmtbx/model/statistics.py:89-101`. The riding test writes into the test's working directory; run the test from a scratch folder.

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd rt/run_$(date +%s) && libtbx.python ../../modules/cctbx_project/mmtbx/hydrogens/tst_capping_h.py`
Expected: `test_process_capped_model` raises Sorry "unknown nonbonded energy type".

- [ ] **Step 3: Implement `add_capping_h_restraints` and the hook**

Bond proxy: `geometry_restraints.bond_simple_proxy(i_seqs=(parent, cap), distance_ideal=xh, weight=1/0.02**2, origin_id=oid)`. Angle proxy: `geometry_restraints.angle_proxy(i_seqs=(x, parent, cap), angle_ideal=ideal, weight=1/3.0**2, origin_id=oid)` with ideal from `steric_number(parent element, record.coordination)`: 4 -> 109.5, 3 -> 120.0, 2 -> 180.0. Neighbours of the parent for the angles: all i_seqs bonded to it in the registry table, caps of the same parent included, each pair once. Register with `process(source_info=_capping_source_info(), proxy=...)`. Nonbonded type: `assign_directly(i_seq=cap, symbol='H')`.

- [ ] **Step 4: Run tests to verify they pass**

Run: the command of Step 2, with `test_find_with_bond_table_equals_without` restored in `run()`. Expected: OK. Then `libtbx.python mmtbx/regression/tst_reduce2_capping.py`, `libtbx.python mmtbx/hydrogens/tst_capping.py`, `libtbx.python mmtbx/monomer_library/tst_pdb_interpretation.py` (if present) and `libtbx.python mmtbx/conformation_dependent_library/tst_pH_mechanism.py`. Expected: all OK.

- [ ] **Step 5: Docstrings**

`capping.py` module docstring: a new paragraph "Reading a capped model" describing `find_capping_h` and `add_capping_h_restraints` in two sentences, pointing to the spec sections 3 and 5. `mmtbx/programs/reduce2.py` `capping` PHIL help (line ~155): drop the words saying the output cannot be processed, say caps have occupancy 0 and are selectable as `capping_h`. `_CapDanglingAtoms` comment (line 2102-2106): same correction.

---

### Task 5: Selection keyword `capping_h`

**Files:**
- Modify: `cctbx_project/iotbx/pdb/atom_selection.py` (`__slots__` 173-202, `__init__` 206-224, `sel_*` 555-575, parser branch ~762)
- Test: `cctbx_project/mmtbx/hydrogens/tst_capping_h.py`

**Interfaces:**
- Consumes: `mmtbx.hydrogens.capping.capping_h_selection(hierarchy)`.
- Produces: keyword `capping_h` in every selection string; `cache.get_capping_h()` / `cache.sel_capping_h()`.

- [ ] **Step 1: Write the failing tests**

```python
def test_keyword_selection():
  m, result = processed(fixtures.pdb_arg)
  expected = capping.capping_h_selection(m.get_hierarchy())
  assert list(m.selection('capping_h')) == list(expected)
  asc = m.get_hierarchy().atom_selection_cache()
  assert list(asc.selection('capping_h')) == list(expected)
  assert asc.selection('not capping_h').count(False) == 1
  assert asc.selection('capping_h and resseq 3').count(True) == 1

def test_keyword_no_h():
  h = iotbx.pdb.input(lines=fixtures.pdb_arg.split('\n'), source_info=None).construct_hierarchy()
  assert h.atom_selection_cache().selection('capping_h').count(True) == 0

def test_pdbtools_remove():
  from libtbx import easy_run
  model, result, h = capped_hierarchy(fixtures.pdb_arg)
  uncapped = fixtures.model_with_h(fixtures.pdb_arg)
  open('capped_arg.pdb', 'w').write(model.model_as_pdb())
  assert easy_run.call('phenix.pdbtools capped_arg.pdb remove=capping_h '
                       'output.file_name=stripped.pdb > pdbtools.log') == 0
  out = iotbx.pdb.input(file_name='stripped.pdb').construct_hierarchy()
  names = [(a.parent().parent().resseq.strip(), a.name.strip()) for a in out.atoms()]
  want = [(a.parent().parent().resseq.strip(), a.name.strip()) for a in uncapped.get_hierarchy().atoms()]
  assert names == want, (names, want)
  assert easy_run.call('phenix.pdbtools capped_arg.pdb remove="not capping_h" '
                       'output.file_name=only_caps.pdb > pdbtools2.log') == 0
  only = iotbx.pdb.input(file_name='only_caps.pdb').construct_hierarchy()
  assert [a.name.strip() for a in only.atoms()] == ['HCG']
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `libtbx.python .../tst_capping_h.py` from a scratch folder.
Expected: `AtomSelectionError` on `capping_h` (unknown keyword).

- [ ] **Step 3: Implement**

In `atom_selection.py`: add `'capping_h'` to `__slots__`, `self.capping_h = None` in `__init__`, `get_capping_h()` that does `from mmtbx.hydrogens.capping import capping_h_selection` inside a try/except ImportError raising `Sorry('capping_h selection needs mmtbx')`, stores `capping_h_selection(self.root).iselection()` and returns `[self.capping_h]`; `sel_capping_h()` like `sel_water`; a parser branch `elif (lword == "capping_h"): result_stack.append(self.sel_capping_h())` next to `single_atom_residue`. Update the class docstring's keyword list.

- [ ] **Step 4: Run tests to verify they pass**

Run: the Step 2 command, then `libtbx.python iotbx/pdb/tst_atom_selection.py` and `libtbx.python iotbx/regression/tst_atom_selections_10k.py`. Expected: OK.

---

### Task 6: Occupancy guards

**Files:**
- Modify: `cctbx_project/mmtbx/model/model.py:3591-3660` (`reset_occupancy_for_hydrogens_simple`, `reset_occupancies_for_hydrogens`)
- Modify: `cctbx_project/mmtbx/refinement/occupancies.py:315-365` (`occupancy_selections`)
- Test: `cctbx_project/mmtbx/hydrogens/tst_capping_h.py`

**Interfaces:**
- Consumes: `model.selection('capping_h')`.
- Produces: both resets leave caps untouched; `occupancy_selections` returns groups without cap i_seqs.

- [ ] **Step 1: Write the failing test**

```python
def test_occupancy_survives():
  from mmtbx.refinement import occupancies
  m, result = processed(fixtures.pdb_ser_altloc)      # caps inside altloc groups
  caps = m.selection('capping_h')
  assert caps.count(True) > 0
  m.reset_occupancy_for_hydrogens_simple()
  assert all(o == 0.0 for o in m.get_occ().select(caps))
  m.reset_occupancies_for_hydrogens()
  assert all(o == 0.0 for o in m.get_occ().select(caps))
  groups = occupancies.occupancy_selections(model=m)
  flat = set(i for g in groups for gg in g for i in gg)
  assert not (flat & set(caps.iselection())), flat & set(caps.iselection())
  # structure factors: caps add nothing
  uncapped = fixtures.model_with_h(fixtures.pdb_ser_altloc)
  for mm in (m, uncapped): mm.setup_scattering_dictionaries(scattering_table='n_gaussian')
  fc_c = m.get_xray_structure().structure_factors(d_min=2.0).f_calc()
  fc_u = uncapped.get_xray_structure().structure_factors(d_min=2.0).f_calc()
  assert flex.max(flex.abs(fc_c.data() - fc_u.common_set(fc_c).data())) < 1e-6
```

If `reset_occupancies_for_hydrogens` needs `refinement_flags`, set `m.refinement_flags = None`-safe path or skip that call with a comment naming the attribute it needs; do not leave the assertion out.

- [ ] **Step 2: Run test to verify it fails**

Expected: after `reset_occupancy_for_hydrogens_simple` a cap has the parent's occupancy.

- [ ] **Step 3: Implement**

In both resets: `caps = self.selection('capping_h')` once; `continue` for any `i_h` with `caps[i_h]`. In `occupancy_selections`, after the `result.extend(pogl)` block and before "add partial occupancies": `caps = model.selection('capping_h').iselection()`; if non-empty, `result = mmtbx.utils.remove_selections(selection=result, other=caps, size=model.get_number_of_atoms())` and exclude them from `sel` (`sel = sel & ~model.selection('capping_h')`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `tst_capping_h.py`; `libtbx.python mmtbx/refinement/tst_occupancy_selections.py`; `libtbx.python mmtbx/regression/tst_reduce2_altloc_completion.py`; `libtbx.python mmtbx/regression/tst_reduce2_neutron_output.py`. Expected: OK.

---

### Task 7: qrefine end-to-end test

**Files:**
- Create: `qrefine/tests/unit/tst_53.py`
- Modify: `qrefine/run_tests.py:58` (add `"$D/tests/unit/tst_53.py",`)
- Modify: `qrefine/CHANGES.md` (one line, same style as the latest entry)

**Interfaces:**
- Consumes: `qrefine.tests.unit.run_qrefine_in(args, prefix, folder=None)` (attributes `pdb`, `model`, `std_out_text`), `run_fmodel(args, prefix)` (attribute `mtz`); `mmtbx.reduce2` command line as in `tst_23.py:55`; `phenix.pdbtools remove=capping_h`.

- [ ] **Step 1: Write the test**

```python
pdb_str = <the ARG-truncated fixture: tst_capping.pdb_arg text, inline>

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine on a model capped by mmtbx.reduce2 qmready=true capping=true:
  loads, refines with cctbx restraints and with mopac, keeps the caps at
  occupancy 0; phenix.pdbtools remove=capping_h strips exactly the caps.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  open("in.pdb", "w").write(pdb_str)
  assert easy_run.call("mmtbx.reduce2 qmready=true capping=true in.pdb > reduce2.log") == 0
  capped = glob.glob("inH_charge*.pdb")[0]          # charge in the name, tst_23
  charge = int(capped.split("charge")[1].split(".pdb")[0])
  r = run_fmodel(args=[capped, "high_res=1.5"], prefix="data")
  caps_before = capping_atoms(capped)                # (resseq, name) of capping_h, via model.selection
  assert caps_before == [("3", "HCG")], caps_before
  for engine_args in (["restraints=cctbx"], ["restraints=qm", "quantum.engine_name=mopac"]):
    out = run_qrefine_in(args=[capped, r.mtz, "number_of_macro_cycles=1",
      "max_iterations_refine=5"] + engine_args, prefix=prefix + "_" + engine_args[0][-5:])
    assert capping_atoms(out.pdb) == caps_before
    occ = [a.occ for a in out.model.get_hierarchy().atoms() if a.name.strip() == "HCG"]
    assert occ == [0.0], occ
    if engine_args[0] == "restraints=qm":
      assert ("charge: %d" % charge) in out.std_out_text or str(charge) in out.std_out_text
  assert easy_run.call("phenix.pdbtools %s remove=capping_h output.file_name=stripped.pdb > pdbtools.log" % out.pdb) == 0
  names = [a.name.strip() for a in iotbx.pdb.input(file_name="stripped.pdb").construct_hierarchy().atoms()]
  assert "HCG" not in names and "HB2" in names and "HB3" in names
```

`capping_atoms(file_name)`: `mmtbx.model.manager` on the file, `selection('capping_h')`, return sorted `(resseq, name)`. Parameter names for cycles and the mopac engine: copy them from `tst_38.py` / `qr.get_default_params()` (`grep -n "number_of_macro_cycles\|engine_name" qrefine/qr.py`). Where qr.refine prints the total charge, find the exact phrase in `qrefine/refine.py` or `charges.py` (`grep -n "charge" qrefine/refine.py`) and assert on it.

- [ ] **Step 2: Run it**

Run: `cd rt/run_$(date +%s) && libtbx.python ../../modules/qrefine/tests/unit/tst_53.py`
Expected: OK with Tasks 1-6 in place (before them, qr.refine stops on unknown atoms: that is the failing state).

- [ ] **Step 3: Register and run the suite**

Add to `run_tests.py`; run `qr.test` or `libtbx.python qrefine/run_tests.py`. Expected: every test that passed before still passes (`tst_38.py` in particular).

---

### Task 8: phenix.refine regression test

**Files:**
- Create: `phenix_regression/refinement/tst_173.py` (numbered after the last `tst_NNN.py`; the directory is discovered, no list to edit)

**Interfaces:**
- Consumes: `phenix_regression.refinement.run_phenix_refine(args, prefix)` (attributes `pdb`, `log`, `geo`), `run_fmodel(args, prefix)` (attribute `mtz`), `get_r_factors(file_name)` (`__init__.py:10`).

- [ ] **Step 1: Write the test**

```python
def run(prefix = os.path.basename(__file__).replace(".py","")):
  """
  phenix.refine on a model capped by mmtbx.reduce2 qmready=true capping=true:
  runs with the default strategy, keeps the caps at occupancy 0, geometry
  statistics and R factors at cycle 0 equal the uncapped model's, and
  phenix.pdbtools remove=capping_h strips the caps.
  """
  open("%s_in.pdb" % prefix, "w").write(pdb_str)   # pdb_arg fixture, inline
  assert easy_run.call("mmtbx.reduce2 qmready=true capping=true %s_in.pdb > %s_reduce2.log" % (prefix, prefix)) == 0
  capped = glob.glob("%s_inH_charge*.pdb" % prefix)[0]
  r = run_fmodel(args=[capped, "high_res=1.5", "type=real", "r_free=0.1"], prefix=prefix)
  args = [capped, r.mtz, "main.number_of_mac=1", "strategy=individual_sites+individual_adp"]
  c = run_phenix_refine(args=args, prefix=prefix + "_capped")
  # uncapped: strip, refine with the same args
  assert easy_run.call("phenix.pdbtools %s remove=capping_h output.file_name=%s_stripped.pdb > %s_pdbtools.log" % (capped, prefix, prefix)) == 0
  u = run_phenix_refine(args=["%s_stripped.pdb" % prefix, r.mtz, "main.number_of_mac=1",
                              "strategy=individual_sites+individual_adp"], prefix=prefix + "_uncapped")
  rc, ru = get_r_factors(c.log), get_r_factors(u.log)
  assert approx_equal(rc.r_work_start, ru.r_work_start, 1e-4)   # attribute names per get_r_factors
  hc = iotbx.pdb.input(file_name=c.pdb).construct_hierarchy()
  caps = [a for a in hc.atoms() if a.name.strip() == "HCG"]
  assert len(caps) == 1 and caps[0].occ == 0.0
  assert "Capping H" in open(c.geo).read()
  # bond rmsd at start equal (log line of cycle 0)
  assert bond_rmsd_start(c.log) == bond_rmsd_start(u.log)
```

`bond_rmsd_start(log)`: parse the first "bonds" rmsd line of the refinement log (look at `tst_000.py`-style tests that read `get_r_factors` for the log format, `__init__.py:10-38`). Use the attribute names `get_r_factors` defines.

- [ ] **Step 2: Run it**

Run: `cd rt/run_$(date +%s) && libtbx.python ../../modules/phenix_regression/refinement/tst_173.py`
Expected: OK. If phenix.refine's riding H manager refuses a cap and prints a warning, the test still passes (free H); record the warning text in the test docstring.

---

### Task 9: Other Phenix programs on a capped model

**Files:**
- Create: `phenix_regression/refinement/tst_174.py` (same discovery rule as Task 8)

**Interfaces:**
- Consumes: the capped model and fmodel data made as in Task 8 (duplicate the ten lines, do not import from `tst_173.py`), `phenix_regression.real_space_refine.run_real_space_refine(args, prefix)` (attributes `pdb`, `std_out`; see `real_space_refine/__init__.py:12-33` and `tst_00.py` for how the map is passed), `easy_run.call` for the command lines.

- [ ] **Step 1: Write the test**

One `run()` with five blocks, each on the capped model `capped` and, where a comparison is needed, on `stripped.pdb` made with `phenix.pdbtools remove=capping_h`. Helper `caps_of(file_name)` returns `[(name, occ)]` of `selection('capping_h')` atoms of an `mmtbx.model.manager` on the file; helper `cap_bonded(file_name)` asserts each cap is within 0.05 A of `capping.xh_length('C', 4, True)` from its nearest heavy atom in the residue.

```python
  # 1. phenix.geometry_minimization
  assert not easy_run.call("phenix.geometry_minimization %s output_file_name_prefix=%s_gm macro_cycles=1 max_iterations=20 > %s_gm.log" % (capped, prefix, prefix))
  assert caps_of("%s_gm.pdb" % prefix) == [("HCG", 0.0)]
  cap_bonded("%s_gm.pdb" % prefix)
  assert "Capping H" in open("%s_gm.geo" % prefix).read()
  # 2. phenix.model_statistics: bond and angle rmsd equal to the stripped model's
  for tag, f in (("capped", capped), ("stripped", "stripped.pdb")):
    assert not easy_run.call("phenix.model_statistics %s use_hydrogens=true > %s_%s.model_statistics" % (f, prefix, tag))
  assert bond_angle("%s_capped.model_statistics" % prefix) == bond_angle("%s_stripped.model_statistics" % prefix)
  # 3. phenix.real_space_refine on the fmodel map coefficients
  rs = run_real_space_refine(args=[capped, r.mtz, "resolution=1.5", "macro_cycles=1"], prefix=prefix + "_rsr")
  assert caps_of(rs.pdb) == [("HCG", 0.0)]
  cap_bonded(rs.pdb)
  # 4. phenix.dynamics
  assert not easy_run.call("phenix.dynamics %s number_of_steps=50 output_file_name_prefix=%s_dyn > %s_dyn.log" % (capped, prefix, prefix))
  assert caps_of("%s_dyn.pdb" % prefix) == [("HCG", 0.0)]
  cap_bonded("%s_dyn.pdb" % prefix)
  # 5. phenix.maps: map coefficients identical with and without caps
  for tag, f in (("capped", capped), ("stripped", "stripped.pdb")):
    assert not easy_run.call("phenix.maps %s %s output.prefix=%s_%s > %s_%s_maps.log" % (f, r.mtz, prefix, tag, prefix, tag))
  mc = mtz.object(file_name="%s_capped_map_coeffs.mtz" % prefix).as_miller_arrays()
  ms = mtz.object(file_name="%s_stripped_map_coeffs.mtz" % prefix).as_miller_arrays()
  for a, b in zip(mc, ms):
    assert a.info().label_string() == b.info().label_string()
    assert flex.max(flex.abs(a.data() - b.data())) < 1e-6, a.info().label_string()
```

`bond_angle(file_name)` parses the "  Bond      :" and "  Angle     :" lines as `refinement/tst_147.py:44-48` does. Output file names of phenix.maps and phenix.dynamics: take them from the programs' logs on the first run and fix the names in the test (`phenix.maps` writes `<prefix>_map_coeffs.mtz` by default; confirm). Parameter names (`number_of_steps`, `macro_cycles`, `resolution`) are checked against each program's `--show-defaults` before use; replace any that does not exist with the right one and leave a comment naming the program.

- [ ] **Step 2: Run it**

Run: `cd rt/run_$(date +%s) && libtbx.python ../../modules/phenix_regression/refinement/tst_174.py`
Expected: OK. A program that refuses the capped model is a defect of Tasks 4-6, not of the test: report it with the program's message.

---

## Self-review notes

- Spec coverage: section 3 -> Task 3; 4 -> Task 3; 5 -> Task 4; 6 -> Task 1; 7 -> Task 5; 8 -> Task 2; 9 -> Tasks 2 and 6; 10 -> Tasks 7, 8 and 9; 11 tests 1-11 -> Tasks 2-6 (1: `test_find_matches_added` + `test_cap_occupancy_zero`; 2: `test_process_capped_model`; 3: `test_statistics_exclude_caps`; 4: `test_idempotent`; 5: Task 5; 6: `test_cap_name_avoids_dictionary` + `test_find_old_ho3_name`; 7: `test_still_stops_wrong_name` + `test_find_rejects_excess`; 8: `test_riding_minimization`; 9: `test_neutron_and_deuterium`; 10: `test_speed`; 11: `test_occupancy_survives`); 12 -> task order; 13 decisions carried as constants above.
- Names used across tasks: `find_capping_h`, `capping_h_selection`, `shared_mon_lib_srv`, `add_capping_h_restraints`, `_budget`, `_dictionary_names`, `_neighbours_by_distance`, `_capping_source_info`; keyword `capping_h`; origin id `'capping'`.
