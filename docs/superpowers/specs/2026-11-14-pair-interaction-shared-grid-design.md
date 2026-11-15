# pair_interaction on a shared grid: exact near-region search for views

Date: 2026-11-14. Status: design approved in conversation, spec for review.
Repositories: cctbx_project (`mmtbx/pair_interaction`), qrefine
(`fragment.py`). Follow-up of the altloc clustering project
(`2026-11-07-altloc-clustering-variants-design.md`, section 7, "Known
limitations and follow-ups", item (1)).

## 1. Problem

With `altlocs.handling=conformers` and `clustering=True`, qr.refine runs the
interaction search of `mmtbx.pair_interaction` once per single-conformer
view of the super sphere, for every cluster and both passes (sedd, dori).
On 1PWM (two views, 191 clusters) that is 766 calls and 108 s of the 198 s
cctbx setup; each extra view costs about 57 s. The views differ only in the
altloc residue groups (and in residue groups absent from a conformer), so
almost all of the second view's work repeats the first view's.

The repetition cannot be skipped today because `get_interactions` sets the
grid origin to the minimum coordinate of the sub-hierarchy of each call: the
views sample different points, so no point of one view corresponds to a
point of another. An exact skip of whole clusters (run one view when no
altloc atom is within reach of the cluster) was tried in the project and
never applied on 1PWM.

## 2. Why a shared grid makes the second view cheap

In `points_and_pairs` (`interaction.h`) a grid point is assigned to the pair
of residues of its two nearest atoms within the cutoff of 10 Bohr, and the
density at the point is built from the atoms within 10 Bohr (candidates are
the atoms of the point's cell and the 26 neighbouring cells, visited in
ascending atom index). The outcome at a point therefore depends only on

- the point's position, and
- the atoms within 10 Bohr of it: coordinates, elements, residue membership
  and relative order (ties of the integer-truncated distances are resolved
  by order).

Two conformer views are selections of the same parent (the super sphere):
common atoms have the same coordinates and elements in both, their relative
order is the parent order, and their view residues map to the same parent
residues. Hence, on one shared grid, a point with no differing atom within
10 Bohr gives the same parent-residue pair in every view. Only points within
10 Bohr of a differing atom ("near" points) have to be searched in a second
view; the rest ("far" points) are read off the first view.

Two more facts keep the construction exact:

- The core bounding box (`lo`/`hi`, core atoms plus or minus the cutoff)
  only saves work: a pair that touches a core residue comes from a point
  within 10 Bohr of a core atom, so a larger box finds nothing new.
- `run_prepared` restricts each call to the atoms within 15 Bohr of the core
  atoms (`sub_ph`). This truncation affects results at the edge of the box
  (an existing approximation, kept). The per-view sub-selections differ, so
  the set of differing atoms is taken as the symmetric difference of the
  sub-selections, not just the altloc atoms.
- A residue group absent from a view has its atoms among the differing
  atoms, so every pair that involves it comes from a near point.

## 3. Design

### 3.1 C++: `points_and_pairs`

Two new optional arguments, both defaulting to "off":

- `region_sites` (`af::shared<vec3<double> >`, Bohr, same frame as `xyz`):
  the differing atoms. A grid point is near when a region site lies within
  the cutoff of 10 Bohr of it, found through a second `cell_list` over the
  region sites with the same origin and cell size as the atom cell list.
- `near_only` (bool): skip far points before any density work.

The returned pairs keep their type `af::shared<vec3<int> >`; the third
component (always 0 today, dropped by the Python wrapper) becomes a flag:
0 when the pair was found at a far point, 1 when it has so far been found
only at near points. The dedup matrix holds three states per residue pair
(unseen, found near, found far) in one byte, as the bool matrix does today.
Rules per point:

- state far: skip (as today);
- state near and the point is near: skip;
- state near and the point is far: evaluate the density; if it interacts,
  set the state to far and emit the pair again with flag 0 (Python keeps the
  minimum flag per pair);
- unseen: evaluate; if it interacts, set the state by the point and emit.

Without region sites every point is far, no pair is emitted twice and the
result is bit-identical to today, flag 0 included. Memory is unchanged;
the extra cost with region sites is one cell lookup per point plus the
re-evaluation of pairs first found near and then found far.

### 3.2 Python (`pair_interaction.py`): shared grid and the view function

Amended 2026-11-14 during implementation (ledger ruling, Task 4): the grid
is per core, not per super sphere. tst_64 and tst_67 of qrefine pin the
duplicate-conformer identity: a model whose A and B conformers have the
same coordinates must give the single-conformer fragments bit for bit. A
grid fixed to the whole super sphere samples other points than the
single-conformer route's per-call box and breaks that identity. The box of
the union of the per-view sub-selections of a core equals each view's own
box for duplicate conformers (identity kept) and is one grid shared by all
views of that core for real altlocs (all the near/far argument needs).

- `_grid_search(...)`: the body of `get_interactions` with the search
  parameters, returning the pairs and their flags. It takes `grid`
  (`group_args(xyz_min, xyz_max)` in Bohr; origin fixed to `xyz_min`,
  `ngrid` from `xyz_max`; else the box of the atoms it is given, as today),
  `region_sites` and `near_only`. It asserts that every atom lies inside
  the grid. `get_interactions` calls it and returns the pairs as before.
- `prepare(ph)`: unchanged.
- `run_prepared(prepared, core=None, grid=None)`: `grid` is passed to
  every search of the call; default behaviour unchanged.
- `views_grid(views, core)`: the box (Bohr) of the union over the views
  with a core residue of the atoms within 15 Bohr of their core atoms
  (`_core_selection`, shared with `run_prepared`); None when no view has
  a core residue.
- `run_prepared_views(views, core)`: `views` is a list of
  `group_args(prepared, iselection, residues)`: `prepare(view)`, the
  parent atom position of each view atom (`flex.size_t`, ascending) and
  the 1-based parent residue group index of each view residue (list).
  `core` is a list of 1-based parent residue indices. Returns one
  `(core_atoms, fragment_atoms, fragment_residues)` per view, in parent
  terms (1-based parent atom positions, parent residue indices), equal to
  what `run_prepared(view, core in view indices, grid=views_grid(views,
  core))` gives, mapped to the parent; `None` for a view that has no
  residue of the core.

  Algorithm (reference view: the first view with a core residue):

  1. Per participating view: core atoms (core residues present in the
     view), the 15 Bohr sub-selection, `sub_ph` and `atom_in_residue`, as
     in `run_prepared`; the sub-selection as parent atoms. The grid is the
     box of the union of the sub-selections.
  2. Region: the union over views of the symmetric differences between the
     view's parent sub-selection and the reference's; coordinates (Bohr)
     from whichever view holds the atom.
  3. Pass 1 (sedd): the reference view searches all points with flags and
     its own core (every view's core, mapped into the reference, is a
     subset of the reference's core). Other views search near points only
     with their own core. A view's pairs are the reference's far pairs that
     touch the view's core plus its own near pairs, all in parent residue
     indices. `new_core` per view as in `run_prepared`.
  4. Pass 2 (dori): the reference view searches all points with flags and
     core = union over views of `new_core` mapped into the reference
     (residues absent from the reference are dropped: their pairs come from
     near points). Other views search near points only with their own
     `new_core`. A view's pairs are the reference's far pairs that touch the
     view's `new_core` plus its own near pairs (for the reference view: its
     flag-1 pairs, filtered the same way, since it searched with the union
     core).
  5. Per view: `interaction_mols = new_core` plus the residues of pass-2
     pairs that touch the original core (as `run_prepared`); fragment atoms
     are the view's atoms of those residues, mapped to the parent.

  Pairs are kept as sets of parent residue pairs; the flag of a pair emitted
  twice is the minimum (far wins).

### 3.3 qrefine (`fragment.py`)

- `_prepared_sphere_views`: the returned entries become the
  `group_args(prepared, iselection, residues)` the mmtbx function takes
  (no grid: `run_prepared_views` derives it per core).
- `_run_prepared`: with views, calls `run_prepared_views(views, residues)`
  and unites the per-view results as today (core atoms, fragment residues;
  fragment atoms are whole parent residues of the union). The per-view core
  mapping and the skip of a view lacking the core move into mmtbx.
- Unchanged: `get_clusters` (two whole-model `run` calls, 2.6 s), the
  non-altloc route (one view, no grid), `altlocs.handling=native`, PHIL,
  logs. Conformer-route fragments of real-altloc models may change once
  (the grid of a core is the union box of its views, not each view's own
  box); duplicate-conformer models keep the single-conformer fragments bit
  for bit (tst_64, tst_67).

## 4. Testing

- New `mmtbx/pair_interaction/tst_08.py` (registered in
  `mmtbx/run_tests.py`): the 2lvr fixture of tst_04 (31 residues, 492
  atoms) with the side chains of three residues duplicated as A and B with
  a shifted B copy and one residue present in B only (absent from view A).
  Views are `ph.select(blank or label)` with their iselections and parent
  residue lists built in the test. Checks: (a) for four cores, including
  one that contains an altloc residue and one that contains the B-only
  residue, with either view as the reference, `run_prepared_views` equals
  the independent `run_prepared` of each view on `views_grid`, mapped to
  the parent; (b) with region sites, far pairs plus near pairs equal the
  pairs found without region sites, and the `near_only` pairs are a strict
  subset; (c) without region sites the flags are all 0 and the pairs equal
  those of the current code (tst_04's expectations stay as they are);
  (d) duplicate conformers (shift 0, no B-only residue): every view equals
  the single-conformer `run_prepared` on its own box, bit for bit. Under
  60 s.
- Existing `tst_02` to `tst_07` pass unchanged (the no-region, no-grid path
  is bit-identical).
- qrefine: `tst_65`, `tst_68`, `tst_69` (conformer routes, frozen-buffer and
  route-equivalence checks) and the full suite pass.
- 1PWM, cctbx restraints, `parallel.nproc=1`: setup time before and after;
  expected: the grid search drops from 108 s to the first view's share plus
  seconds for the second view. Recorded in section 5 below and in section 7
  of the 2026-11-07 spec.

## 5. Results (2026-11-14)

Implemented as amended in 3.2 (per-core grid). mmtbx `tst_08` passes:
every view equals its independent run on `views_grid` for four cores with
either view as the reference, duplicate conformers equal the
single-conformer run bit for bit, flags and `near_only` behave as 3.1
says. mmtbx `tst_02` to `tst_07` unchanged. qrefine suite 74 of 75: the
failure is `tst_49`, which fails on this tree without this change too
(bisected by reverting `fragment.py` alone; the tree carries unrelated
uncommitted cctbx edits of another session).

1PWM, cctbx restraints, `clustering=True altlocs.handling=conformers`,
one weight-search and one refinement cycle, `parallel.nproc=1`, same
command as the three-route comparison of the 2026-11-07 spec. "Before" is
the tree with `fragment.py` of HEAD swapped in (the mmtbx changes are
inert without it); setup stages timed in process around
`create_fragment_manager`:

| | before | after |
|---|---|---|
| qr.refine wall / user (s) | 474 (2026-11-12) | 457 / 450 |
| fragment manager (s) | 184.5 | 177.7 |
| of it: grid search (`_grid_search`, s) | 113.8 | 112.5 |
| of it: Girvan-Newman + whole-model searches (`get_clusters`, s) | 48.2 | 48.2 |
| fragments / jobs | 190 / 373 | 190 / 373 |
| fragment atoms (super sphere), sum / median | 66188 / 293 | 65972 / 291 |
| R-free after refinement | 0.1350 | 0.1347 |
| max RSS (GB) | - | 1.69 |

cProfile of the whole qr.refine run (after): `create_restraints_manager`
180.6 s, `get_fragments` 110.4 s, `run_prepared_views` 108.7 s in 153
calls, of which `_grid_search` (the C++ search) 105.9 s in 614 calls and
the Python of the view union (`_core_selection`, sets, mapping) 2.6 s;
`get_clusters` 50.4 s, Girvan-Newman 45.3 s of it.

The saving is about 1 % of the grid search, far from the "seconds per
extra view" the 2026-11-07 spec predicted. The prediction assumed the
altloc neighbourhoods are a small part of a cluster's box. On 1PWM they
are not: 168 of 912 residue groups carry altlocs (3140 of 8452 atoms), 93 %
of the atoms and 72 % of the occupied grid points lie within 10 Bohr of an
altloc atom, the region of a core has 600 to 900 sites, and in a
microbenchmark on three 15-residue cores of the model only 0 to 10 % of
the sedd pairs are found at far points (5 of 52, 0 of 10, 3 of 9); the
near-only pass of the second view costs about 95 % of a full pass. The
construction is exact and costs nothing where it does not help; it pays
on models with sparse altlocs (6uf7: one residue), where the second view
collapses to the neighbourhood of that residue.

Where the time really goes (same microbenchmark): a 15-residue core of
the model has a sub-selection of 1300 to 2200 atoms and a box of about
2 million sedd points; the box of every cluster overlaps those of its
neighbours, so the 190 clusters evaluate the same grid points of the super
sphere many times over (the whole-model `run` in `get_clusters`, which
evaluates each point once, takes 1.3 s per view). The per-cluster
15 Bohr sub-selection exists only to make those repeated searches cheaper
and truncates the environment of the points near the box edge. One search
of the whole super sphere per pass and view, with the per-cluster pairs
taken from it by set algebra (pairs touching the core, then pairs touching
the new core), would replace about 110 s by an estimated 15 s and remove
the truncation; it changes the fragments of every clustered run, altlocs
or not, so it is a separate decision (section 6).

## 6. Out of scope and follow-up

- One whole-super-sphere search per pass and view instead of one per
  cluster (section 5): an estimated 110 s to 15 s on 1PWM, no edge
  truncation, fragments of every clustered run change once. Needs the
  user's decision; the near/far machinery of this project then applies
  once per view instead of once per cluster.

- The whole-model `run` calls in `get_clusters` (2.6 s).
- The dedup matrix of `points_and_pairs` is dense over the residue count of
  the hierarchy passed (`max(atom_in_residue)` squared bytes, the
  sub-hierarchy keeping the parent's residue indices). It is unchanged here
  and is a separate item for ribosome-scale super spheres.
- Girvan-Newman speed (follow-up (2) of the 2026-11-07 spec).
