# qr.refine runtime bottlenecks Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut the qr.refine run time (176 min on 1PWM with cctbx restraints) by fixing the weight-search stopping rules, the duplicate statistics, the per-residue OMIT maps of fit_rotatable2, and the fragmentation setup (clustering and pair interactions), with results unchanged except where the spec says so.

**Architecture:** qrefine gets a pure weight-search controller (new module `qrefine/weights.py`) that driver.refine drives; the monitor caches its statistics; the cctbx cluster route caches fragment restraints. cctbx gets a single-OMIT-map fit_rotatable2, a wfc pickle convention with an in-process cache, a prepare/run split of pair_interaction.run, and a C++ Girvan-Newman clustering with per-component recomputation that replaces the boost call bit-identically.

**Tech Stack:** Python (cctbx/qrefine), C++ (mmtbx_pair_interaction_ext, Boost.Python), libtbx tests run with `libtbx.python`; one extension rebuild with `cd /Users/pafonine/Desktop/all/phenix/build; ./bin/libtbx.scons -j8 lib/mmtbx_pair_interaction_ext.so`.

**Spec:** `qrefine/docs/superpowers/specs/2026-11-05-runtime-bottlenecks.md`

## Global Constraints

- Every qrefine test writes only into its own `qrefine_tst_XX` folder (keep the makedirs/chdir lines in run()).
- Clustering and fragment results stay bit-identical (interaction lists, clusters, fragment atoms) for Tasks 4-7; Task 2 (weight search) and Task 3 (OMIT map) change results by design.
- No occupancies in restraint targets or gradients.
- Leave changes uncommitted; the user commits.
- Timing claims only from a quiet machine; tests must not assert wall-clock times.
- PHIL: `refine.number_of_weight_search_cycles` default 50 -> 8; `refine.adjust_restraints_weight_scale_value` (default 2) is the search factor; `refine.stop_one_found_first_good_weight` is removed; `refine.r_tolerance` (0.001) is the R-free plateau tolerance.

## Review Focus

- First trial not GOOD and no GOOD trial within the cap: the search must fall back to the last trial without raising (Task 2 test `exercise_never_good`).
- `skip_weight_search=True` or `data_weight` given: the controller is not used and refinement behaves as before (Task 2, tst_44/tst_59 unchanged).
- fit_rotatable2 on a model with no rotatable H: no map is built and the sites are unchanged (Task 3 test `exercise_no_rotatable`).
- load_wfc when the mmtbx folder is not writable: the object is still returned (parsed) and no exception leaks (Task 4 test `exercise_unwritable`).
- Girvan-Newman input with duplicate or reversed pairs, self-pairs or isolated vertices: identical to the boost path, which ignores duplicates by set semantics (Task 6 test `exercise_degenerate`).

---

### Task 1: Statistics cache and elapsed time in the monitor

**Files:**
- Modify: `qrefine/results.py` (class manager: `__init__`, `update`, `get_stats`, `g_info`, `show`)
- Test: `qrefine/tests/unit/tst_60.py` (new), register in `qrefine/run_tests.py` after tst_59

**Interfaces:**
- Produces: `manager.get_stats()` returns the cached `mmtbx.model.statistics.geometry` of the current sites (computed at most once per `update`); `manager.show(prefix, suffix)` appends `time: %d s` (seconds since the monitor was created) to every line.

- [ ] **Step 1: Write the failing test**

```python
# qrefine/tests/unit/tst_60.py
def exercise_stats_cached():
  """show() twice without update() computes the statistics once."""
  model = get_model("m00_good.pdb")  # from data_files, processed with restraints
  calls = []
  orig = model.geometry_statistics
  model.geometry_statistics = lambda **kw: (calls.append(1), orig(**kw))[1]
  log = StringIO()
  monitor = results.manager(model=model,
    geometry_rmsd_manager=calculator.cctbx_geometry(model=model), log=log)
  monitor.show(prefix="a"); monitor.show(prefix="b")
  assert len(calls) == 1
  monitor.update(model=model); monitor.show(prefix="c")
  assert len(calls) == 2
  assert "time:" in log.getvalue().splitlines()[0]
```

- [ ] **Step 2: Run it: `libtbx.python qrefine/tests/unit/tst_60.py`; expected: AssertionError (2 calls)**
- [ ] **Step 3: Implement** in results.py: `self._stats = None` and `self._t0 = time.time()` in `__init__`; `update()` sets `self._stats = None`; `get_stats()` computes and stores when None; `show()` appends `"time: %d s" % (time.time()-self._t0)`.
- [ ] **Step 4: Run tst_60 (PASS) and tst_02, tst_59 (unchanged)**

### Task 2: Weight-search controller

**Files:**
- Create: `qrefine/weights.py`
- Modify: `qrefine/driver.py:72-152` (the weight loop in `refine`), `qrefine/qr.py` PHIL (`number_of_weight_search_cycles` default 8, remove `stop_one_found_first_good_weight`, help of `adjust_restraints_weight_scale_value` and `r_tolerance`), `qrefine/api.py:_auto_cust` (`stop_one_found_first_good_weight` line removed)
- Test: `qrefine/tests/unit/tst_60.py` (add), `qrefine/tests/unit/tst_56.py` (remove any stop_one_found reference)

**Interfaces:**
- Produces: `class controller(object)` in weights.py:
  - `__init__(self, start_scale=1.0, factor=2.0, max_trials=8, r_tolerance=0.001)`
  - `scale(self) -> float`: the scale of the next trial (start_scale first)
  - `record(self, scale, good, r_free, sites_cart)`: result of the trial at `scale` (good: bond rmsd within max_bond_rmsd)
  - `done(self) -> bool`
  - `choice(self) -> group_args(scale, r_free, sites_cart, reason)`; `reason` one of "best", "fallback" (no GOOD trial: the last trial)
  - `trials`: list of group_args(scale, good, r_free)
- Rules (spec 5.A): direction from trial 0 (GOOD -> down, else up) by `factor`; at the first flip one more trial at the geometric mean of the bracket, then done; done when the last two GOOD trials each failed to improve the best GOOD r_free before them by more than r_tolerance; done at max_trials. Choice: among GOOD trials with r_free <= min_good_r_free + r_tolerance, the largest scale.

- [ ] **Step 1: Write the failing tests** (tst_60):

```python
def exercise_1pwm_sequence():
  c = weights.controller(factor=1.5)   # the 1PWM log, factor 1.5
  for r in [0.1369, 0.1360, 0.1356, 0.1350, 0.1357]:
    if c.done(): break
    c.record(scale=c.scale(), good=True, r_free=r, sites_cart=None)
  assert len(c.trials) == 3, [t.r_free for t in c.trials]
  assert c.done()
  assert abs(c.choice().scale - 1/1.5) < 1e-9       # 0.1360 within 0.001 of 0.1356, larger scale
def exercise_flip_then_bisect():
  c = weights.controller()
  c.record(1.0, False, 0.20, None); assert abs(c.scale()-2.0) < 1e-9
  c.record(2.0, True, 0.19, None);  assert abs(c.scale()-math.sqrt(2.0)) < 1e-9
  c.record(c.scale(), True, 0.185, None); assert c.done()
  assert abs(c.choice().scale - math.sqrt(2.0)) < 1e-9
def exercise_never_good():
  c = weights.controller(max_trials=4)
  while not c.done(): c.record(c.scale(), False, 0.3, None)
  assert len(c.trials) == 4 and c.choice().reason == "fallback"
  assert abs(c.choice().scale - 8.0) < 1e-9
def exercise_cap():
  c = weights.controller(max_trials=3, r_tolerance=0.0)
  for r in [0.3, 0.2, 0.1]: c.record(c.scale(), True, r, None)
  assert c.done() and abs(c.choice().scale - 0.25) < 1e-9
```

- [ ] **Step 2: Run tst_60: ImportError on qrefine.weights**
- [ ] **Step 3: Implement weights.controller** (pure Python, no cctbx imports beyond libtbx.group_args).
- [ ] **Step 4: Run tst_60: PASS**
- [ ] **Step 5: Wire driver.refine**: replace the loop body: `c = weights.controller(start_scale=1.0, factor=params.refine.adjust_restraints_weight_scale_value, max_trials=params.refine.number_of_weight_search_cycles, r_tolerance=params.refine.r_tolerance)`; per trial `calculator.setw(restraints_weight_scale=c.scale())`, minimize, `c.record(scale, GOOD, monitor.r_free, monitor.sites_cart)`; keep the `monitor.show` line per trial (one line, suffix `rws: %6.3f n_fev: %d`; drop the "(NEXT)" line); after the loop `ch = c.choice()`, print `Weight search: %d trials, chosen rws %.3f (%s)`, set the scale and sites from `ch`. Remove the `up/down/r_frees/b_rmsds/a_rmsds/restraints_weight_scale` bookkeeping and the stall rule.
- [ ] **Step 6: PHIL and api changes** listed above; grep `stop_one_found` across qrefine (docs too) and remove.
- [ ] **Step 7: Run tst_02, tst_06, tst_44, tst_55, tst_56, tst_59: PASS**

### Task 3: fit_rotatable2 with one OMIT map (cctbx)

**Files:**
- Modify: `cctbx_project/mmtbx/hydrogens/__init__.py:392-500` (`fit_rotatable2`)
- Test: `cctbx_project/mmtbx/hydrogens/tst_fit_rotatable2.py` (new); register in `cctbx_project/mmtbx/run_tests.py` next to `hydrogens/tst_xh_mapper.py`

**Interfaces:**
- Consumes: `shortcut(residue, first, log)` -> list of `[axis_i_seqs, h_i_seqs]`; `map_manager.update_omit_map(omit_selection)`, `map_manager.score(sites_cart)`.
- Produces: same signature `fit_rotatable2(model, fmodel, reset_h=True, mapper=None)`; one `update_omit_map` call with the union of all rotatable H of the model; if the union is empty, no map and no change.

- [ ] **Step 1: Write the failing test**

```python
def exercise_recovers_rotation():
  # model: a SER-THR-LYS tripeptide PDB string with riding H (take the
  # fixture style of tst_xh_mapper.py); f_obs = |Fcalc| of the true model at
  # d_min 0.9 with r_free_flags; rotate HG of SER by 120 deg about CB-OG
  # with rotate_point_around_axis; run fit_rotatable2(reset_h=False)
  # assert: HG back within 0.15 A of the true site
def exercise_one_map():
  # same setup; monkeypatch mmtbx.hydrogens.map_manager.update_omit_map to count
  # calls; assert exactly 1 call, and the omit selection size equals the
  # number of rotatable H of the model (SER 1 + THR 4 + LYS 3 + N-terminal 3)
def exercise_no_rotatable():
  # GLY-GLY model: update_omit_map not called, sites unchanged
```

- [ ] **Step 2: Run: `exercise_one_map` FAILS (3 or more calls)**
- [ ] **Step 3: Implement**: first pass over residues collecting `(cl, first)` lists per residue via `shortcut`; `omit_selection` = union; if empty: restore mask params and return; `mm.update_omit_map(omit_selection)` once; second pass fits every residue's rotations in that map exactly as the current loop does (angle scan, best score). Keep the final `fmodel.update_xray_structure` + `update_all_scales()` as now.
- [ ] **Step 4: Run tst_fit_rotatable2.py and tst_xh_mapper.py: PASS**
- [ ] **Step 5: Record on 1PWM** (scratch folder): r_work/r_free after "Optimize H" old vs new from the log lines (old: 0.1347/0.1404), and the time of the step; note in the spec section 6.

### Task 4: wfc pickle convention and in-process cache (cctbx)

**Files:**
- Modify: `cctbx_project/mmtbx/pair_interaction/pair_interaction.py:24-71` (`load_wfc`)
- Rename (untracked files): `mmtbx/pair_interaction/cl_lda.wfc -> cl_wfc_obj.pkl`, `f__lda.wfc -> f__wfc_obj.pkl`, `p__lda.wfc -> p__wfc_obj.pkl` (the user commits them)
- Test: `cctbx_project/mmtbx/pair_interaction/tst_05.py` (rewrite)

**Interfaces:**
- Produces: `load_wfc(element, folder=None) -> ext.wfc`; key = element.lower() padded with "_" to two characters ("c_", "cl"); looks up `_WFC_CACHE[key]`, then `<folder>/<key>_wfc_obj.pkl` (folder default: `libtbx.env.find_in_repositories("mmtbx/pair_interaction")`), else parses `qrefine/plugin/yoink/dat/<key>_lda.wfc` and writes `<folder>/<key>_wfc_obj.pkl` (write failures ignored). Nothing is ever written into the cwd.

- [ ] **Step 1: Rewrite tst_05**

```python
def run():
  if libtbx.env.find_in_repositories("qrefine/plugin/yoink/dat") is None: return
  tmp = tempfile.mkdtemp()
  pair_interaction._WFC_CACHE.clear()
  o1 = pair_interaction.load_wfc("xe", folder=tmp)       # parsed, written
  assert os.path.isfile(os.path.join(tmp, "xe_wfc_obj.pkl"))
  assert not os.path.isfile("xe_lda.wfc")               # nothing in cwd
  pair_interaction._WFC_CACHE.clear()
  o2 = pair_interaction.load_wfc("xe", folder=tmp)       # loaded
  assert o1.ngrid == o2.ngrid and o1.zz == o2.zz
  assert pair_interaction.load_wfc("xe", folder=tmp) is o2   # cached
  ro = os.path.join(tmp, "ro"); os.mkdir(ro); os.chmod(ro, 0o500)   # exercise_unwritable
  pair_interaction._WFC_CACHE.clear()
  assert pair_interaction.load_wfc("kr", folder=ro).zz == 36
  assert pair_interaction.load_wfc("c_").zz == 6 and pair_interaction.load_wfc("cl").zz == 17
```

- [ ] **Step 2: Run: FAILS (no `folder` argument)**
- [ ] **Step 3: Implement load_wfc as specified; rename the three files**
- [ ] **Step 4: Run pair_interaction tst_00..tst_06 and qrefine tst_48: PASS; check `git status` shows no new *_lda.wfc anywhere**

### Task 5: pair_interaction.run: prepare once, run per cluster (cctbx + qrefine)

**Files:**
- Modify: `cctbx_project/mmtbx/pair_interaction/pair_interaction.py:76-148` (`run`), `qrefine/fragment.py:183-219` (`get_fragments`)
- Test: `qrefine/tests/unit/tst_61.py` (new), register after tst_60

**Interfaces:**
- Produces in pair_interaction.py:
  - `prepare(ph) -> group_args(ph, atoms, xyz (Bohr, flex.vec3_double), atom_in_residue (flex.int, 1-based), residue_atoms (list of flex.size_t per residue, 0-based i_seq), n_residues, element_wfc_dict)`; `ph` is a deep copy with Bohr coordinates and i_seqs reset (the caller's hierarchy is not changed).
  - `run_prepared(prepared, core) -> (core_atoms, fragment_atoms, fragment_residues)` identical to `run(ph.deep_copy(), core)`; `run(ph, core=None)` keeps its behaviour and becomes `prepare` + `run_prepared` (core given) or the whole-model path.
  - The 15-Bohr selection: one call `xyz.min_distance_between_any_pair`-free vectorized form: for the core atoms' xyz `c`, `selection = flex.bool(n, False); for each core atom: selection |= (xyz - a).norms() < 15` stays as is in Python unless profiling shows it matters (it is 0.1 s per cluster), otherwise a C++ `within_distance(xyz, core_xyz, 15)`.
- Consumes in fragment.get_fragments: `prepared = pair_interaction.prepare(self.pdb_hierarchy_super)` once; per cluster `pair_interaction.run_prepared(prepared, cluster)`; `two_buffers` the same way.

- [ ] **Step 1: Write the failing test** (tst_61): for 1yjp.pdb and 2olx.pdb from data_files, `il = run(ph.deep_copy())`, clusters from `qrefine.clustering.betweenness_centrality_clustering(..., maxnum 4 / 15)`; for every cluster assert `run(ph.deep_copy(), cluster) == run_prepared(prepare(ph), cluster)` (tuple of sorted lists), and that `ph` coordinates are unchanged after `prepare`.
- [ ] **Step 2: Run: AttributeError (prepare)**
- [ ] **Step 3: Implement prepare/run_prepared; rewrite run on top of them**
- [ ] **Step 4: Run tst_61, pair_interaction tst_00..06, qrefine tst_01, tst_10, tst_11, tst_48: PASS**
- [ ] **Step 5: Switch fragment.get_fragments to prepare/run_prepared; run qrefine tst_03, tst_04, tst_12, tst_29, tst_33: PASS**

### Task 6: Girvan-Newman clustering in C++ with per-component recomputation

**Files:**
- Create: `cctbx_project/mmtbx/pair_interaction/girvan_newman.h`
- Modify: `cctbx_project/mmtbx/pair_interaction/interaction_ext.cpp` (export `girvan_newman`), `qrefine/clustering.py:106-160` (`betweenness_centrality_clustering.get_clusters`, `build_graph`)
- Test: `cctbx_project/mmtbx/pair_interaction/tst_07.py` (new, register in mmtbx/run_tests.py), `qrefine/tests/unit/tst_48.py` (existing, must pass unchanged)

**Interfaces:**
- Produces: `ext.girvan_newman(n_vertices: int, edges: flex.int (2*m, pairs i,j 0-based in insertion order), threshold: double) -> flex.int (remaining edge pairs, in the input order)`. Semantics of boost `betweenness_centrality_clustering(g, bc_clustering_threshold(threshold, g, false))` on `adjacency_list<setS, vecS, undirectedS>`: edges deduplicated (an (i,j)/(j,i) repeat or a self-pair is dropped, first occurrence kept); repeat {Brandes edge betweenness, halved for undirected; e = first edge in input order with the maximum; if max < threshold stop; remove e} until no edge. Per-component recomputation: after a removal, recompute only the component(s) of the removed edge's endpoints; components whose maximum is already < threshold are final.
- Consumes in clustering.py: `get_clusters` builds the deduplicated 0-based edge list from `interaction_list`, then for threshold in 9..4: `edges = ext.girvan_newman(size, edges, threshold)`, components from the remaining edges (`connected_pieces` with residues 1..size), stop when the largest is <= maxnum; the rest (split_cluster) unchanged. `build_graph` and the boost graph imports go; keep a module-level `_boost_remaining_edges(size, edges, threshold)` only inside tst_07 (test reference), not in clustering.py.

- [ ] **Step 1: Write the failing test** (tst_07): reference function using `boost_adaptbx.graph` as clustering.py does today (copy its graph construction), returning the remaining edges as sorted pairs; cases: (a) 30 random graphs (seeded) with 8-60 vertices, 1-3 edges per vertex, thresholds 2, 4, 9, including duplicates/reversed/self pairs (`exercise_degenerate`); (b) 1yjp interaction list (pair_interaction.run on data_files/1yjp.pdb if qrefine is present) at thresholds 9..4 chained; assert `set(ext_result) == set(reference)` and the ext result keeps input order.
- [ ] **Step 2: Build fails / AttributeError girvan_newman**
- [ ] **Step 3: Implement girvan_newman.h** (adjacency as vectors, Brandes with BFS queue and dependency accumulation in vertex order, doubles; component labels by BFS; a per-component record {max_value, max_edge_index, final}); export in interaction_ext.cpp with `af::shared<int>` in/out; rebuild the extension.
- [ ] **Step 4: Run tst_07: PASS**
- [ ] **Step 5: Switch clustering.py; run qrefine tst_48, tst_03, tst_04, tst_12, tst_29, tst_33, tst_59: PASS**
- [ ] **Step 6: Time it on conformer A of 1PWM** with `scratchpad bcc/measure.py` adapted to call the new path; record old (426 s loaded) vs new in the spec section 6.

### Task 7: cctbx cluster route: cached fragment restraints, no per-evaluation super-sphere file

**Files:**
- Modify: `qrefine/restraints.py:187-230` (`from_cctbx.__init__`, `target_and_gradients`), `qrefine/cluster_restraints.py:30-44` (`target_and_gradients` file write)
- Test: `qrefine/tests/unit/tst_61.py` (add `exercise_grm_cache`), existing tst_10, tst_11

**Interfaces:**
- Produces: `from_cctbx._fragment_grm(index)` returns `self.fragment_extracts.super_sphere_geometry_restraints_manager.select(super_selection)` cached in `self._grm_cache` keyed by `(id(self.fragment_extracts.fragment_super_selections), index)`; the cache is dropped when the key's first element changes (re-clustering builds new selection lists). The super sphere file `tmp_ase.pdb` is written only when `restraints_manager` is a `from_qm`.

- [ ] **Step 1: Write the failing test**: run two evaluations of a `from_cluster` built as tst_11 does (small model, clustering=True, cctbx) and count `mmtbx.restraints.manager.select` calls via monkeypatch: second evaluation makes 0 calls; gradients identical between the two evaluations at the same sites; no `tmp_ase.pdb` in the working folder.
- [ ] **Step 2: Run: FAILS (select called again; file present)**
- [ ] **Step 3: Implement**
- [ ] **Step 4: Run tst_61, tst_10, tst_11, tst_29, tst_59: PASS**

### Task 8: Benchmark and spec update

- [ ] **Step 1: Rerun** `qr.refine 1PWM.mtz 1PWMH_charge3.pdb restraints=cctbx` in a fresh `rt3/run_<timestamp>` folder; collect the phase table from the new `time:` stamps.
- [ ] **Step 2: Full test runs**: `libtbx.python qrefine/run_tests.py` and the mmtbx tests of pair_interaction and hydrogens; list any failure against the baseline in memory `mmtbx-preexisting-failures-2026-11-02`.
- [ ] **Step 3: Update spec section 6** with the numbers (flagging machine load) and list what changed in results (weight chosen, H fit R-factors, cluster timing).
