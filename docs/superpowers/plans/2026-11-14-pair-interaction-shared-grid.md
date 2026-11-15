# pair_interaction shared grid Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the per-view interaction search of qr.refine's conformer
clustering route cost seconds per extra view instead of ~57 s, with results
identical by construction to independent per-view searches on one shared
grid.

**Architecture:** `points_and_pairs` (C++) learns a region of differing
atoms and reports per pair whether it was found at a point far from that
region; `mmtbx.pair_interaction` gains a shared grid and one function that
runs the reference view in full and the other views on near points only,
assembling each view's exact result; qrefine's `fragment.py` prepares its
conformer views on the super sphere's grid and calls that function.

**Tech Stack:** cctbx (Boost.Python extension `mmtbx_pair_interaction_ext`,
`scitbx.array_family.flex`), Python 3 under `libtbx.python`.

**Spec:**
`modules/qrefine/docs/superpowers/specs/2026-11-14-pair-interaction-shared-grid-design.md`
(read it first; section 2 is the exactness argument every task relies on).

## Global Constraints

- Root is `/Users/pafonine/Desktop/all/phenix`; run everything through
  `<root>/build/bin/` dispatchers (`libtbx.python`, `qr.refine`,
  `libtbx.scons`, `libtbx.run_tests_parallel`). Never source setpaths, never
  plain `python`.
- Rebuild the extension after any C++ change:
  `cd <root>/build && ./bin/libtbx.scons -j8 lib/mmtbx_pair_interaction_ext.so`
  (exit 0, then re-run the tests that use it).
- Without `region_sites` and without a `grid`, `points_and_pairs`,
  `get_interactions`, `run` and `run_prepared` are bit-identical to today:
  mmtbx `tst_02` to `tst_07` pass unchanged, the non-altloc qrefine route is
  untouched.
- Memory: no new dense matrix in `points_and_pairs`; the dedup matrix stays
  one byte per residue pair.
- Occupancies are never used anywhere in this work.
- No new PHIL parameters, no new log lines.
- Edited and new lines are at most 80 columns.
- mmtbx `tst_08` writes no files and runs under 60 s; qrefine tests write
  only into their own `qrefine_tst_XX/` folder.
- Manual runs go into fresh scratch directories under
  `/private/tmp/claude-502/-Users-pafonine-Desktop-all-phenix/6cfbe5ff-cf2b-4471-8e70-dedb1b8b13e1/scratchpad/run_$(date +%s)_<tag>`;
  never `rm -rf`.
- Leave both trees (`cctbx_project`, `qrefine`) uncommitted at the end of
  every task; the user commits. "Commit" steps below mean: report the diff
  (`git -C <repo> diff --stat`) and stop.

> Amended 2026-11-14 (ledger, Task 4 ruling): the grid is per core
> (`views_grid`, the union box of the views' sub-selections), `grid_of` and
> `prepare(grid)` do not exist, `run_prepared` takes `grid=None`; tst_64
> and tst_67 pin the duplicate-conformer identity. Tasks 2-4 below are the
> original text; tst_08 and the code follow the amended spec 3.2.

## Review Focus

1. An atom of a view outside the shared grid (caller passed a grid of the
   wrong hierarchy): `_grid_search` must fail with an assertion naming the
   grid, not silently clamp the atom into the last cell. Test in Task 2.
2. Views prepared on different grids passed to `run_prepared_views`:
   assertion, not a wrong union. Test in Task 3.
3. A core with no residue in any view, or in one view only: `None` for the
   views without it, exact results for the others. Test in Task 3.
4. A single view passed to `run_prepared_views`: identical to
   `run_prepared`. Test in Task 3.
5. `near_only=True` with an empty region: nothing is searched and the
   result is empty (no fallback to a full search). Test in Task 1.

---

### Task 1: near/far points and flags in `points_and_pairs`, Python face

**Files:**
- Modify: `cctbx_project/mmtbx/pair_interaction/interaction.h`
  (`points_and_pairs`, lines 442-558, and its comment block)
- Modify: `cctbx_project/mmtbx/pair_interaction/interaction_ext.cpp`
  (the `points_and_pairs` def, lines 115-135)
- Modify: `cctbx_project/mmtbx/pair_interaction/pair_interaction.py`
  (`get_interactions`, lines 219-266)
- Create: `cctbx_project/mmtbx/pair_interaction/tst_08.py`
- Modify: `cctbx_project/mmtbx/run_tests.py` (add
  `"$D/pair_interaction/tst_08.py"` after tst_07)

**Interfaces:**
- Produces (C++, `interaction.h`):
  `points_and_pairs(ngrid, step_size, xyz, xyz_min, atom_in_residue,
  element_flags, wfc_obj, silva_type, residue_in_core,
  region_sites=af::shared<vec3<double> >(), near_only=false)` returning
  `af::shared<vec3<int> >` of `(i, j, flag)`, `i<j` 1-based residue
  indices, `flag` 0 (found at a far point) or 1 (found only at near
  points); a pair may appear twice (first with 1, later with 0).
- Produces (Python): `_grid_search(ph, atom_in_residue, step_size,
  silva_type, core, element_wfc_dict, grid=None, region_sites=None,
  near_only=False) -> dict {(i, j): flag}` with the minimum flag per pair;
  `grid` is `group_args(xyz_min, xyz_max)` in Bohr (origin fixed to
  `grid.xyz_min`, `ngrid[k] = floor((xyz_max[k]-xyz_min[k])/step)+1`), else
  the box of `ph`'s atoms as today; asserts every atom lies within
  `[xyz_min, xyz_max]` when a grid is given. `get_interactions` keeps its
  signature and return (`list(set(pairs))`) and calls `_grid_search`.
- Produces (test fixture, `tst_08.py`): `altloc_fixture()` returning
  `group_args(ph, views)`: `ph` is the 2lvr hierarchy of tst_04 (`from
  mmtbx.pair_interaction.tst_04 import pdb_str`) with the side-chain atoms
  (all atoms except N, CA, C, O, CB) of TYR 4, GLN 12 and ARG 21 duplicated
  as altloc A (original coordinates) and B (shifted by +0.7 A along x), and
  every atom of HIS 26 given altloc B (present in view B only); `views` is
  the list for labels A, B of `group_args(ph=parent.select(sel),
  iselection=sel.iselection(), residues=[1-based parent residue group
  index of each view residue group, ascending])` where `sel` selects atoms
  with altloc blank or the label. Build the fixture from PDB text (edit the
  lines of `pdb_str`: altloc column 17, x column 31-38), not through the
  hierarchy API, so atom order is the parent order.

- [ ] **Step 1: Write the failing test `exercise_flags` in `tst_08.py`**

```python
def exercise_flags():
  f = altloc_fixture()
  p = pair_interaction.prepare(f.views[0].ph)        # view A, no grid
  air = list(p.atom_in_residue)
  wfc = p.element_wfc_dict
  step = 0.3*pair_interaction.A2B
  full = pair_interaction._grid_search(p.ph, air, step, 'sedd', None, wfc)
  assert len(full) > 0 and set(full.values()) == set([0])
  ref = pair_interaction.get_interactions(p.ph, air, silva_type='sedd',
    element_wfc_dict=wfc)
  assert sorted(full) == sorted(ref)                 # bit identity
  # region: the side chain of TYR 4 (parent residue 4) in view A, Bohr
  region = p.xyz.select(p.residue_atoms[3])
  flagged = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    wfc, region_sites=region)
  assert sorted(flagged) == sorted(full)             # same pairs
  only_near = set([k for k, v in flagged.items() if v == 1])
  assert len(only_near) > 0
  near = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    wfc, region_sites=region, near_only=True)
  assert set(near.values()) == set([1])
  assert only_near <= set(near) < set(full)          # strict subset
  # with a core, as run_prepared uses it
  core = [3, 4, 5]
  c_full = pair_interaction._grid_search(p.ph, air, step, 'sedd', core, wfc)
  c_flag = pair_interaction._grid_search(p.ph, air, step, 'sedd', core, wfc,
    region_sites=region)
  assert sorted(c_full) == sorted(c_flag)
  # Review Focus 5: near_only with nothing near searches nothing
  empty = pair_interaction._grid_search(p.ph, air, step, 'sedd', None, wfc,
    region_sites=flex.vec3_double(), near_only=True)
  assert empty == {}
```

`run()` of tst_08 calls `exercise_flags()` and prints `OK`; the test is
executed as `libtbx.python tst_08.py` from a scratch directory.

- [ ] **Step 2: Run it, expect failure**

Run: `cd <scratch> && <root>/build/bin/libtbx.python <root>/modules/cctbx_project/mmtbx/pair_interaction/tst_08.py`
Expected: `AttributeError: ... has no attribute '_grid_search'`.

- [ ] **Step 3: Implement the C++**

In `points_and_pairs`: add the two parameters; build
`cell_list region(region_sites.const_ref(), xyz_min, cutoff)` when
`region_sites.size() > 0`; `near(point)` is true when a region site among
`region.candidates(point)` has `dist_sq < cutoff*cutoff`; with
`near_only`, `continue` on far points before the atom candidate loop.
Replace the bool `pairs` matrix by `af::versa<signed char, af::c_grid<2> >`
with states 0 unseen, 1 found near, 2 found far, and apply the rules of
spec 3.1 (skip on 2; skip on 1 when the point is near; on an interaction
set the state to 2 for a far point, 1 for a near point, and push
`(ia1, ia2, point_is_near ? 1 : 0)`). Update the comment block above the
function (grid, region, flag, duplicate emission). In `interaction_ext.cpp`
extend the function-pointer cast and add
`arg("region_sites")=af::shared<vec3<double> >(), arg("near_only")=false`.

- [ ] **Step 4: Rebuild**

Run: `cd <root>/build && ./bin/libtbx.scons -j8 lib/mmtbx_pair_interaction_ext.so`
Expected: exit 0.

- [ ] **Step 5: Implement `_grid_search` and refactor `get_interactions`**

Move the body of `get_interactions` into `_grid_search` with the signature
above; `get_interactions(ph, atom_in_residue, step_size=0.5*A2B,
silva_type='dori', core=None, element_wfc_dict=None)` keeps its sedd step
override and returns `list(set(_grid_search(...).keys()))`. The flag of a
pair emitted twice is the minimum. `region_sites` is passed to the
extension as `flex.vec3_double` (empty when None).

- [ ] **Step 6: Run tst_08 and the existing tests**

Run: tst_08 as in Step 2; then
`for t in 02 03 04 06 07; do <root>/build/bin/libtbx.python <root>/modules/cctbx_project/mmtbx/pair_interaction/tst_$t.py; done`
Expected: `OK` from tst_08; every existing test exits 0 (tst_06 prints
nothing on success, check `$?`).

- [ ] **Step 7: Register tst_08 in `mmtbx/run_tests.py`; report the diff**

---

### Task 2: shared grid in `prepare` and `run_prepared`

**Files:**
- Modify: `cctbx_project/mmtbx/pair_interaction/pair_interaction.py`
  (`prepare` lines 115-160, `run_prepared` lines 180-217, new `grid_of`)
- Modify: `cctbx_project/mmtbx/pair_interaction/tst_08.py`

**Interfaces:**
- Consumes: `_grid_search(..., grid=...)` (Task 1).
- Produces: `grid_of(ph) -> group_args(xyz_min, xyz_max)`, the bounding
  box of `ph.atoms().extract_xyz() * A2B` (tuples of 3 floats, Bohr);
  `prepare(ph, grid=None)` storing `grid` as `prepared.grid`;
  `run_prepared(prepared, core=None)` passing `prepared.grid` to every
  `_grid_search` call (both passes, and the whole-hierarchy call when
  `core` is None).

- [ ] **Step 1: Write the failing test `exercise_grid`**

```python
def exercise_grid():
  f = altloc_fixture()
  g = pair_interaction.grid_of(f.ph)
  p = pair_interaction.prepare(f.views[0].ph, grid=g)
  assert p.grid is g
  assert pair_interaction.prepare(f.views[0].ph).grid is None
  air = list(p.atom_in_residue)
  step = 0.3*pair_interaction.A2B
  # a grid equal to the atoms' own box samples the same points
  own = group_args(xyz_min=p.xyz.min(), xyz_max=p.xyz.max())
  d0 = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    p.element_wfc_dict)
  d1 = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    p.element_wfc_dict, grid=own)
  assert d0 == d1
  # run_prepared with a grid runs and returns the three lists
  c, a, m = pair_interaction.run_prepared(p, [3, 4, 5])
  assert len(c) > 0 and set(c) <= set(a) and set([3, 4, 5]) <= set(m)
  # Review Focus 1: an atom outside the grid is an error
  small = group_args(xyz_min=g.xyz_min,
    xyz_max=tuple([v-5.0 for v in g.xyz_max]))
  try:
    pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
      p.element_wfc_dict, grid=small)
  except AssertionError as e:
    assert str(e).find("grid") > -1, str(e)
  else:
    raise AssertionError("atom outside the grid was accepted")
```

- [ ] **Step 2: Run tst_08, expect `AttributeError: ... 'grid_of'`**

- [ ] **Step 3: Implement `grid_of`, `prepare(ph, grid=None)`,
  `run_prepared` grid pass-through**

- [ ] **Step 4: Run tst_08 (OK) and tst_02/03/04/06 (exit 0: no grid, no
  change); report the diff**

---

### Task 3: `run_prepared_views`

**Files:**
- Modify: `cctbx_project/mmtbx/pair_interaction/pair_interaction.py`
- Modify: `cctbx_project/mmtbx/pair_interaction/tst_08.py`

**Interfaces:**
- Consumes: `_grid_search` (Task 1), `prepare(ph, grid)`, `grid_of`
  (Task 2); the `prepared` fields `ph, xyz, atom_in_residue,
  residue_atoms, n_residues, element_wfc_dict, grid`.
- Produces: `run_prepared_views(views, core) -> list`: `views` is a list
  of `group_args(prepared, iselection, residues)` (see spec 3.2; all
  `prepared.grid` must be the same object or equal boxes, asserted), `core`
  a list of 1-based parent residue indices; one entry per view: `None`
  when the view has no residue of `core`, else `(core_atoms,
  fragment_atoms, fragment_residues)` as sorted lists in parent terms
  (1-based parent atom positions `iselection[i-1]+1`; parent residue
  indices `residues[k-1]`), equal to `run_prepared(view.prepared,
  core_in_view)` mapped to the parent. A private helper
  `_core_selection(prepared, core_v)` returning `(core_atoms, selection,
  sub_ph, atom_in_residue_selected)` is shared with `run_prepared` so the
  15 Bohr sub-selection cannot diverge between the two paths.

- [ ] **Step 1: Write the failing test `exercise_views`**

```python
def mapped_run_prepared(view, core):
  """run_prepared of one view, mapped to the parent (the reference)."""
  in_view = dict([(r, k+1) for k, r in enumerate(view.residues)])
  core_v = [in_view[r] for r in core if r in in_view]
  if(len(core_v) == 0): return None
  c, a, m = pair_interaction.run_prepared(view.prepared, core_v)
  isel = view.iselection
  return (sorted([isel[i-1]+1 for i in c]),
          sorted([isel[i-1]+1 for i in a]),
          sorted([view.residues[i-1] for i in m]))

def exercise_views():
  f = altloc_fixture()
  g = pair_interaction.grid_of(f.ph)
  views = [group_args(prepared=pair_interaction.prepare(v.ph, grid=g),
    iselection=v.iselection, residues=v.residues) for v in f.views]
  cores = [[1, 2], [3, 4, 5], [13, 14, 15], [24, 25, 26, 27]]
  for core in cores:
    got = pair_interaction.run_prepared_views(views, core)
    expected = [mapped_run_prepared(v, core) for v in views]
    assert got == expected, (core, got, expected)
  # the altloc cores do distinguish the views (non-vacuous)
  for core in [[3, 4, 5], [24, 25, 26, 27]]:
    a, b = pair_interaction.run_prepared_views(views, core)
    assert a[1] != b[1], core
  # Review Focus 4: one view equals run_prepared
  assert pair_interaction.run_prepared_views(views[:1], [3, 4, 5]) == \
    [mapped_run_prepared(views[0], [3, 4, 5])]
  # Review Focus 3: a core absent from a view, or from all views
  got = pair_interaction.run_prepared_views(views, [26])   # HIS 26: B only
  assert got[0] is None and got[1] == mapped_run_prepared(views[1], [26])
  assert pair_interaction.run_prepared_views(views, [999]) == [None, None]
  # Review Focus 2: views on different grids
  other = group_args(prepared=pair_interaction.prepare(f.views[1].ph),
    iselection=f.views[1].iselection, residues=f.views[1].residues)
  try:
    pair_interaction.run_prepared_views([views[0], other], [3, 4, 5])
  except AssertionError:
    pass
  else:
    raise AssertionError("views on different grids were accepted")
```

- [ ] **Step 2: Run tst_08, expect `AttributeError: ... 'run_prepared_views'`**

- [ ] **Step 3: Implement `_core_selection` and `run_prepared_views`**

Follow spec 3.2 steps 1-5 exactly. Decisions the spec leaves open:
the reference view is the first view with a core residue; region sites
are the coordinates (Bohr, from the view holding the atom) of the union
over participating views of `(S_v ^ S_ref)` where `S_v = set(
view.iselection.select(selection_v))`; parent-pair sets are Python sets
of `(min, max)` parent residue tuples; a view's pass-1 pairs are
`{p for p in far_ref if p touches core_v} | near_v`; `new_core_v` and
`interaction_mols_v` follow `run_prepared` literally (pass-2 membership is
tested against the original `core_v`, not `new_core_v`). The reference
view's pass-2 call uses `core = sorted(union of new_core_v mapped into the
reference, absent residues dropped)`, and its own pairs are its flag-0 and
flag-1 pairs both filtered by `new_core_ref`. Return sorted lists.

- [ ] **Step 4: Run tst_08 (OK, under 60 s: print the wall time) and
  tst_02/03/04/06/07 (exit 0); report the diff**

---

### Task 4: qrefine uses the shared grid

**Files:**
- Modify: `qrefine/fragment.py` (`_prepared_sphere_views` lines 400-417,
  `_run_prepared` lines 419-440)

**Interfaces:**
- Consumes: `pair_interaction.grid_of`, `prepare(ph, grid)`,
  `run_prepared_views(views, core)` (Tasks 2-3).
- Produces: `_prepared_sphere_views()` returns a list of
  `group_args(prepared, iselection, residues)`; without conformer
  selections one entry with `iselection=None, residues=None`.
  `_run_prepared(prepared, residues)` keeps its return
  `(core atoms, fragment atoms, fragment residues)` (sorted lists, 1-based
  positions in the super sphere); callers in `get_fragments` are unchanged.

- [ ] **Step 1: Edit `_prepared_sphere_views`**

With conformer selections: `grid = pair_interaction.grid_of(
self.pdb_hierarchy_super)` once, `pair_interaction.prepare(h, grid=grid)`
per view; keep the existing `assert len(list(h.residue_groups())) ==
len(residues)`. Update the docstring (shared grid of the full altloc
super sphere, why).

- [ ] **Step 2: Edit `_run_prepared`**

Non-view entry: `return pair_interaction.run_prepared(entry.prepared,
residues)`. Views: `results = pair_interaction.run_prepared_views(
prepared, residues)`; union of `r[0]` and `r[2]` over the non-None
results; fragment atoms as today (whole residues of the union via
`_sphere_residue_atoms`). Drop the per-view core mapping and the
`in_view` dict. Update the docstring.

- [ ] **Step 3: Run the conformer-route tests**

Run from `<scratch>/run_<ts>_t4`:
`for t in 65 68 69; do <root>/build/bin/libtbx.python <root>/modules/qrefine/tests/unit/tst_$t.py; done`
Expected: `OK` from each (tst_68 about 80 s, tst_69 about 40 s).

- [ ] **Step 4: Run the full qrefine suite**

Run: `mkdir <scratch>/run_<ts>_suite/d && cd <scratch>/run_<ts>_suite/d && <root>/build/bin/libtbx.run_tests_parallel module=qrefine nproc=8 > ../suite.log 2>&1; tail -5 ../suite.log`
Expected: every test passes (75 at the time of writing). Report the diff.

---

### Task 5: 1PWM timing and documentation

**Files:**
- Modify: `qrefine/docs/superpowers/specs/2026-11-14-pair-interaction-shared-grid-design.md`
  (section 5 Results)
- Modify: `qrefine/docs/superpowers/specs/2026-11-07-altloc-clustering-variants-design.md`
  (section 7, "Known limitations and follow-ups", item (1): mark done,
  one sentence with the numbers)
- Modify: `qrefine/CHANGES.md` (one entry under `## [Unreleased]` /
  `### Changed`)

- [ ] **Step 1: Time 1PWM after the change**

In `<scratch>/run_<ts>_1pwm`:
```
/usr/bin/time -l <root>/build/bin/qr.refine \
  <root>/rt3/1PWMH_charge3.pdb <root>/rt3/1PWM.mtz restraints=cctbx \
  clustering=True altlocs.handling=conformers \
  number_of_weight_search_cycles=1 number_of_refine_cycles=1 \
  parallel.nproc=1 > run.log 2> time.log
```
Record: wall (`real`), max RSS, the final `Best Rwork, Rfree (after
refinement)` line, the fragment and job counts from the log. Reference
before the change (2026-11-12, same command): 474 s wall, 2.72 GB,
R-free 0.1350, 190 fragments, 373 jobs; setup 246 s.

- [ ] **Step 2: Profile the setup**

Same scratch directory, the same arguments on
`<root>/modules/qrefine/command_line/qr.py` under
`libtbx.python -m cProfile -o setup.prof`; then
`libtbx.python -c "import pstats; pstats.Stats('setup.prof').sort_stats('cumulative').print_stats(25)"`.
Record the cumulative time of `create_restraints_manager`,
`get_fragments`, `run_prepared_views`, `_grid_search` and the tottime of
`points_and_pairs` (reported under `_grid_search`/`get_interactions`), and
compare with the Task 10 profile of the earlier project (198 s setup,
`get_interactions` 108 s in 766 calls).

- [ ] **Step 3: Write the results**

Spec section 5: a table before/after (wall, setup, grid search, max RSS,
R-free, fragments, jobs), one sentence on whether the fragments changed
(the grid moved) and the cProfile numbers. 2026-11-07 spec item (1): done
2026-11-14, numbers, pointer to this spec. CHANGES.md entry: the
conformer clustering route searches interactions on one shared grid with
the extra views restricted to the altloc neighbourhoods; identical union
by construction; numbers.

- [ ] **Step 4: Report the diffs of both repositories and stop**
