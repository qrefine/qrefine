# Clustering with altlocs: per-fragment conformer variants. Implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** qr.refine's clustering route handles models with any number of altloc labels by clustering the altloc model once and evaluating each fragment once per distinct conformer variant, all variants in one parallel pool, with gradients verified against the existing routes with cctbx and MOPAC.

**Architecture:** `fragment.fragments` gets optional conformer selections and turns its per-fragment job list into a per-variant job list (variant = distinct atom set a global conformer takes from a fragment, with a weight and a representative label). `cluster_restraints.from_cluster` combines variant gradients as weighted means. `refine.create_restraints_manager` sends clustered altloc models to this path instead of `conformers.from_conformers`; plain and expansion routes are untouched.

**Tech Stack:** cctbx (iotbx.pdb hierarchy, cctbx.crystal.super_cell, mmtbx.pair_interaction), qrefine (conformers.py, fragment.py, cluster_restraints.py, restraints.py, refine.py), MOPAC through ASE for QM tests, libtbx test conventions.

**Spec:** `modules/qrefine/docs/superpowers/specs/2026-11-07-altloc-clustering-variants-design.md`

## Global Constraints

- Build: run everything through `<root>/build/bin/` dispatchers (`libtbx.python`, `qr.refine`, `mmtbx.reduce2`); root is `/Users/pafonine/Desktop/all/phenix`. Never source setpaths.
- Every test's `run()` starts with `os.makedirs(prefix, exist_ok=True); os.chdir(prefix)` and writes only into `qrefine_tst_XX/`; new tests are registered in `qrefine/run_tests.py`; each runs in about a minute (MOPAC tests up to two).
- Occupancies are never used in restraint targets or gradients.
- No new PHIL parameters. `altlocs.handling=conformers` is the default and selects the variant path for clustering.
- A single-conformer model must give bit-identical gradients and the same `ase/<index>/` file layout as before this change (tst_35 relies on it).
- Speed and memory are requirements: no per-conformer copies of the super sphere beyond the capper parents; one `parallel_map` per evaluation.
- The engine never sees an altloc or a dangling bond (asserted after relabeling; capper refusals kept).
- Leave the tree uncommitted at the end of every task; the user commits (commit steps below are "report the diff and stop", not git commit).
- Use fresh scratch directories (`run_$(date +%s)`) for manual runs; never `rm -rf`.

## Review Focus

1. Altlocs on waters or ligands only (HETATM residue groups): variants must form on a non-protein fragment; capped and charged like any other. Test added to Task 2.
2. An altloc residue that enters a fragment only through a symmetry copy (buffer copy, master residue outside the fragment): the tiled selection must put the copy's conformer atoms into the right variant. Test added to Task 2.
3. Microheterogeneity (A is ASN, B is ASP in one residue group): variant hierarchies take the right resname per label. Test added to Task 2.
4. Neutron models with H/D exchange plus heavy-atom altlocs: hd_mapper runs first, conformers come from the de-deuterated model; the clustered route must complete. Test added to Task 7.
5. `working_folder=None` (api.refine from phenix.refine, no files) with several variants per fragment: folder naming must not touch the file system. Test added to Task 9.

---

### Task 1: conformer helpers in conformers.py

**Files:**
- Modify: `modules/qrefine/conformers.py` (`_single_conformer_hierarchy`, new `super_selections`)
- Test: `modules/qrefine/tests/unit/tst_31.py` (add `exercise_super_selections`, `exercise_single_conformer_hierarchy`)

**Interfaces:**
- Consumes: `conformer_selections(pdb_hierarchy)` (existing), `cctbx.crystal.super_cell.manager` (`super_cell_hierarchy`, `keep_selection`, `super_sphere_hierarchy`).
- Produces:
  - `single_conformer_hierarchy(pdb_hierarchy, iselection, label, labels) -> iotbx.pdb.hierarchy.root`: public name of today's `_single_conformer_hierarchy`, same contract (atom order of iselection preserved, every atom group blank, merged per residue group, resname rule of the docstring), `iselection` a `flex.size_t` in ascending order; `split` calls it.
  - `super_selections(conformer_selections, n_atoms, scm) -> dict`: `{label: group_args(cell=flex.bool over scm.super_cell_hierarchy, sphere=flex.bool over scm.super_sphere_hierarchy)}` for every `group_args(label, iselection, n_fallback)` of `conformer_selections`; `cell` is the master selection tiled `scm.super_cell_hierarchy.atoms_size() // n_atoms` times (assert the division is exact), `sphere = cell.select(scm.keep_selection)`.

- [ ] **Step 1: Write the failing tests in tst_31.py**

```python
def exercise_super_selections():
  # h_altconf_complete.pdb in the small box of tst_46.get_model(fn, 2.0)
  # (several images in the sphere); for every label the tiled sphere
  # selection equals the per-residue rule applied directly to the sphere
  # hierarchy (same label set, _selections on scm.super_sphere_hierarchy)
  from cctbx.crystal import super_cell
  model = tst_46.get_model("h_altconf_complete.pdb", contact=2.0)
  ph = model.get_hierarchy()
  scm = super_cell.manager(pdb_hierarchy=ph,
    crystal_symmetry=model.crystal_symmetry(), select_within_radius=10)
  sels = conformers.conformer_selections(ph)
  ss = conformers.super_selections(sels, ph.atoms_size(), scm)
  direct = dict([(s.label, s.iselection) for s in
    conformers.conformer_selections(scm.super_sphere_hierarchy)])
  assert sorted(ss.keys()) == ["A", "B"]
  for label in ss:
    assert ss[label].sphere.size() == scm.super_sphere_hierarchy.atoms_size()
    assert list(ss[label].sphere.iselection()) == list(direct[label])
    assert ss[label].cell.size() == scm.super_cell_hierarchy.atoms_size()
    # master block first
    assert list(ss[label].cell[:ph.atoms_size()].iselection()) == \
      list([s for s in sels if s.label == label][0].iselection)

def exercise_single_conformer_hierarchy():
  # public helper: same result as split() for every conformer
  pi = iotbx.pdb.input(os.path.join(qr_unit_tests_data, "h_altconf_complete.pdb"))
  model = mmtbx.model.manager(model_input=pi, log=null_out())
  ph = model.get_hierarchy()
  labels = set(conformers.conformer_labels(ph))
  for c, s in zip(conformers.split(model), conformers.conformer_selections(ph)):
    h = conformers.single_conformer_hierarchy(ph, s.iselection, s.label, labels)
    assert h.is_similar_hierarchy(c.model.get_hierarchy())
    assert list(h.altloc_indices().keys()) == [""]
```

- [ ] **Step 2: Run tst_31 to verify the new exercises fail**

Run: `libtbx.python modules/qrefine/tests/unit/tst_31.py`
Expected: AttributeError on `conformers.super_selections`.

- [ ] **Step 3: Implement `single_conformer_hierarchy` (rename, keep `_single_conformer_hierarchy` out) and `super_selections` in conformers.py; update `split` and the module docstring**

`super_selections`: build `cell` by `flex.bool` concatenation of the master bool (from `iselection`) `n_copies` times; the tiling holds because `super_cell.manager` appends chain-wise detached copies of the whole model in model order (spec 3.2).

- [ ] **Step 4: Run tst_31; expect OK**

- [ ] **Step 5: Report the diff and stop (user commits)**

---

### Task 2: variants in fragment.fragments

**Files:**
- Modify: `modules/qrefine/fragment.py` (`fragments.__init__`, `set_up_cluster_qm`, `get_fragment_hierarchies_and_charges`, `make_capper`, `get_fragment_extracts`, `get_qm_file_name_and_pdb_hierarchy`, `write_mm_charge_file`, `write_cluster_and_fragments_pdbs`)
- Create: `modules/qrefine/tests/unit/tst_64.py`
- Modify: `modules/qrefine/run_tests.py` (register tst_64)

**Interfaces:**
- Consumes: Task 1 `super_selections`, `single_conformer_hierarchy`; `conformers.conformer_selections`.
- Produces (attributes of `fragments`, lists of equal length, one entry per variant = job):
  - `fragments(..., conformer_selections=None)`: the list from `conformers.conformer_selections(pdb_hierarchy)` or None (single-conformer behaviour).
  - `variant_fragment_index` (list of int), `variant_weights` (list of int), `variant_labels` (list of str, "" for a fragment with one variant),
  - `fragment_selections`, `fragment_super_selections`, `cluster_selections`, `buffer_selections`, `fragment_cappers`, `fragment_charges`, `fragment_capped_initial`: per variant (fragment-level lists repeated per variant),
  - `variant_hierarchies` (list of hierarchies: the relabeled single-conformer variant at setup coordinates, atom order = `fragment_super_selections[i]` order),
  - `conformer_spheres` (dict label -> relabeled single-conformer super-sphere hierarchy; `{"": pdb_hierarchy_super}` without conformers),
  - `n_fragments` (int), `fragment_cluster_selections` and `fragment_full_super_selections` (per fragment: today's `cluster_selections` entry and the fragment's full super-sphere selection, the union of its variants).
  - Sorry messages of a variant's capper are re-raised with the prefix `"fragment %d, conformer %s: "` (label `"-"` for a fragment without variants).
  - `get_fragment_extracts()` carries all of the above plus the existing fields.
  - `job_folder_name(fragment_extracts, index) -> str`: `"<fragment_index>"` when the fragment has one variant, else `"<fragment_index>_<label>"`; module-level function used for `ase/` sub-folders and engine labels.

- [ ] **Step 1: Write tst_64.py (cctbx as calculator, no QM)**

Model builders (in the test file): `three_label_model()` from `altloc_scenarios.base_model()` without waters: residue 2 gets side-chain conformers A,B (`side_chain_conformers(rg, ["A","B"], [0.5,0.5], [0,120])`), residue 5 gets B,C (same helper with labels `["B","C"]`), the rest blank; crystal symmetry replaced by the P1 box of `tst_46.get_model` style with contact 30 (no images). Second model `copy_in_buffer_model()`: h_altconf_complete in the contact-2 box of tst_46 (altloc residue reaches other fragments through images). Third: `hetero_model()`: `base_model()` with waters where water 101 (`HOH`) gets `label_atoms(rg, ["O","H1","H2"], "A", 0.5)` and a shifted B copy (`backbone_conformers` is for peptides; copy the atom group by hand, shift 0.5 A along x, label B). Fourth: `micro_model()`: three_label_model where residue 2's B atom group gets `resname="ASP"` with the side chain beyond CB renamed to ASP names (CG, OD1, OD2; drop ND2's H) or, simpler, reuse tst_31's `exercise_microheterogeneity` builder if it is importable.

```python
def exercise_three_labels():
  model = three_label_model()
  ph = model.get_hierarchy()
  f = fragments(working_folder=None, make_working_folder=False,
    clustering_method="bcc", maxnum_residues_in_cluster=1,
    pdb_hierarchy=ph, crystal_symmetry=model.crystal_symmetry(),
    qm_engine_name="mopac", cap_fragments=True,
    conformer_selections=conformers.conformer_selections(ph))
  n_conf = 3
  for i in range(f.n_fragments):
    jobs = [j for j, k in enumerate(f.variant_fragment_index) if k == i]
    assert sum(f.variant_weights[j] for j in jobs) == n_conf
    labels = set(conformers.conformer_labels(
      f.pdb_hierarchy_super.select(f.fragment_full_super_selections[i])))
    assert len(jobs) <= len(labels) + 1 and len(jobs) <= n_conf
    if(len(labels) == 0):
      assert len(jobs) == 1 and f.variant_labels[jobs[0]] == ""
    for j in jobs:
      h = f.variant_hierarchies[j]
      assert list(h.altloc_indices().keys()) == [""]
      assert h.atoms_size() == f.fragment_super_selections[j].count(True)
      assert f.fragment_cappers[j] is not None
      assert f.fragment_charges[j] is not None
  # a fragment containing residue 2 ({A,B}) and residue 5 ({B,C}): A, B, C
  # give three distinct variants (A, B, C) and no fourth
  # (job lists per fragment checked against hand-derived expectations)
  ...
  # job count
  assert len(f.variant_weights) == sum(...)  # hand-derived

def exercise_variant_atoms_match_split():
  # every variant's master atoms are exactly the fragment atoms of the
  # conformer model conformers.split gives for the representative label
  ...

def exercise_capper_parent():
  # for every variant with a label, the number of caps equals the number of
  # caps of the same residues cut from the single-conformer model of that
  # label (conformers.split(model) -> fragments on that model with the same
  # clusters): no cap on a B atom next to an A fragment
  ...

def exercise_single_conformer_unchanged():
  # a model without altlocs: one variant per fragment, weights 1, labels "",
  # lists identical to fragments(conformer_selections=None)
  ...

def exercise_copy_in_buffer():   # Review Focus 2
  # copy_in_buffer_model: a fragment whose master residues have no altloc
  # but whose buffer holds an image of the altloc residue has two variants
  ...

def exercise_hetero():           # Review Focus 1
  # water with A/B: the fragment of that water has two variants, both capped
  # (zero caps) and charged 0
  ...

def exercise_micro():            # Review Focus 3
  # variant A has resname ASN, variant B ASP in residue 2
  ...

def exercise_job_folder_name():
  assert fragment.job_folder_name(fe, j) == "3" / "3_B" per the rule

def exercise_sorry_prefix():
  # three_label_model with residue 5's atom groups renamed to resname XXX
  # (no dictionary): fragments(...) raises Sorry whose message starts with
  # "fragment <n>, conformer B: Capping: no dictionary for XXX"
  ...
```

- [ ] **Step 2: Run tst_64; expect TypeError (unexpected keyword conformer_selections)**

- [ ] **Step 3: Implement in fragment.py**

- `__init__`: store `conformer_selections`; `self.labels = [s.label for s in conformer_selections]` or `[]`.
- `set_up_cluster_qm`: after `get_clusters` and `get_fragments`, compute `self._super_selections = conformers.super_selections(...)` (empty dict without conformers), build `self.conformer_cells` (dict label -> `single_conformer_hierarchy(SCM.super_cell_hierarchy, cell.iselection(), label, labels)`) and `self.conformer_spheres` (same on `pdb_hierarchy_super` with `sphere`), then `get_fragment_hierarchies_and_charges`.
- `_variants(fragment_super_selection)`: for every label `v = fragment_super_selection & sphere[label]`; group equal `v` (key: `tuple(v.iselection())`) into `[(v, weight, representative_label)]` with the smallest label as representative; without conformers `[(fragment_super_selection, 1, "")]`; a fragment whose variants are all equal (no altloc atoms) collapses to one entry of weight `n_conformers` and label "".
- `get_fragment_hierarchies_and_charges`: loop fragments, then variants; per variant: `variant_hierarchy = single_conformer_hierarchy(pdb_hierarchy_super, v.iselection(), label, labels)` (plain `select` with blank label when `label == ""`), capper from `make_capper(v, label)`, charge, `fragment_selections` entry = master atoms of `v` (positions `< pdb_hierarchy.atoms_size()`), `fragment_super_selections` entry = `v`. Fragment-level `cluster_selections`/`buffer_selections` repeated per variant. Assert the variant hierarchy has no non-blank altloc.
- `make_capper(fragment_super_selection, label)`: parent = `self.conformer_cells[label]` (the super cell itself when `label == ""`); selection over the parent: `keep.iselection().select(v.iselection())` mapped through the rank array of `cell[label].iselection()`; assert `c.fragment_hierarchy.is_similar_hierarchy(variant_hierarchy)`; a `Sorry` from the capper is re-raised with the `"fragment %d, conformer %s: "` prefix.
- `get_fragment_extracts`: add the new fields.
- `get_qm_file_name_and_pdb_hierarchy`: fragment hierarchy = `variant_hierarchies[index]` with `set_xyz(pdb_hierarchy_super.select(fragment_super_selections[index]).atoms().extract_xyz())`; sub-folder from `job_folder_name`.
- `write_mm_charge_file`: environment from `conformer_spheres[variant_labels[index]]` minus the variant (selection mapped onto that sphere through `sphere[label]`), folder from `job_folder_name`.
- `write_cluster_and_fragments_pdbs`: `<i>_cluster.pdb`, `<i>_frag.pdb` per fragment; `<job_folder_name>_capped0.pdb` per variant; `fragment_info.txt` gains a line per fragment listing variant labels and weights.

- [ ] **Step 4: Run tst_64, tst_10, tst_35 (layout), tst_46; expect OK**

- [ ] **Step 5: Report the diff and stop**

---

### Task 3: weighted combination and job plumbing

**Files:**
- Modify: `modules/qrefine/cluster_restraints.py` (`from_cluster.target_and_gradients`)
- Modify: `modules/qrefine/restraints.py` (`from_qm.target_and_gradients` clustering branch, `from_cctbx._fragment_grm` unchanged except naming)
- Create: `modules/qrefine/tests/unit/tst_65.py`
- Modify: `modules/qrefine/run_tests.py`

**Interfaces:**
- Consumes: Task 2 extracts (`variant_weights`, `variant_fragment_index`, `fragment_cluster_selections`, `n_fragments`), `job_folder_name`.
- Produces: `from_cluster.target_and_gradients(sites_cart) -> (target, gradients)` with the spec 3.4 combination; `from_cluster.n_jobs` (int) and `from_cluster.multiplicity` (`flex.double`, the per-atom weight sum `den` of the last call; 1.0 everywhere for a single-conformer model) for tests.

- [ ] **Step 1: Write tst_65.py (cctbx)**

Helper `clustered(model, handling, maxnum)` builds params as `tst_46.get_cctbx_params(expansion=False, clustering=True, handling, maxnum)` and returns `refine.create_restraints_manager(params, model)`; until Task 4 lands, the test builds the manager directly: `fragments(..., conformer_selections=...)` + `cluster_restraints.from_cluster(from_cctbx(super-sphere grm), fm, params.parallel)`. Write it against `refine.create_restraints_manager` and expect it to fail until Task 4.

```python
def exercise_one_cluster_limit():
  # maxnum=999: new route == from_conformers over expansion, 1e-6
  for fn in ["h_altconf_complete.pdb", three_label (from tst_64)]:
    model = get_model(fn, contact=2.0)      # images in the sphere
    rm_c = clustered(model, "conformers", maxnum=999)
    params_e = tst_46.get_cctbx_params(expansion=True, clustering=False,
      handling="conformers", maxnum=999)
    rm_e = refine.create_restraints_manager(params_e, model)
    sites = model.get_sites_cart()
    g_c = rm_c.target_and_gradients(sites)[1]; g_e = rm_e.target_and_gradients(sites)[1]
    assert flex.max((g_c - g_e).norms()) < 1.e-6
    assert rm_c.n_jobs == conformers.n_conformers(model.get_hierarchy())

def exercise_relation_per_atom():
  # same clustering: native clustered (a) vs conformers clustered (c)
  # altloc atoms equal; blank atoms: m*g_c == g_a + (m-1)*g_blank_only
  # where g_blank_only is the native clustered gradient of the blank-only
  # model (ph.select(blank)) with the same clusters (pass maxnum=3)
  ...

def exercise_clustering_vs_expansion():
  # altloc data files (altlocs.pdb, altlocs2.pdb, h_altconf*.pdb) in their
  # own boxes, maxnum=3, conformers handling: clustering vs expansion
  # within 1.e-4 (tst_10's one-buffer miss), reported per file
  ...

def exercise_single_conformer_bit_identical():
  # a model without altlocs (gly2_1.pdb): from_cluster over fragments built
  # with conformer_selections=None and with
  # conformer_selections=conformers.conformer_selections(ph) (an empty
  # list) give gradients with max difference exactly 0, n_jobs ==
  # n_fragments, multiplicity all 1
  ...
```

- [ ] **Step 2: Run tst_65; expect failures (n_jobs missing; combination wrong)**

- [ ] **Step 3: Implement**

`from_cluster.target_and_gradients`: job list over all variants; after `parallel_map`, `num = flex.vec3_double(n)`, `den = flex.double(n)`; for each job `j`: `w = variant_weights[j]`, master selection `fragment_selections[j]`, `num.set_selected(sel, num.select(sel) + w*g)`, `den.set_selected(sel, den.select(sel) + w)`; then zero `num` on the fragment's buffer (`buffer_selections[j]`, fragment-level) as today; `gradients = num * (1/den)` where `den > 0` (assert all cluster atoms have `den > 0`); `target = sum(w_j*e_j)/n_conformers` with `n_conformers = max(1, sum of weights of the first fragment's variants)`. `n_jobs = len(variant_weights)`; `self.multiplicity = den` (1.0 for buffer-only entries is irrelevant; set `den` to 1 where it is 0 before dividing).
`from_qm.target_and_gradients`: no logic change beyond the helpers of Task 2 (folder and label from `job_folder_name`); engine label `os.path.splitext(qm_pdb_file)[0]` stays.

- [ ] **Step 4: Run tst_65 (the parts not needing Task 4), tst_04, tst_10, tst_29, tst_35; expect OK**

- [ ] **Step 5: Report the diff and stop**

---

### Task 4: dispatch, driver, old proxy removal, tst_43 and tst_46 on the new route

**Files:**
- Modify: `modules/qrefine/refine.py` (`create_fragment_manager`, `create_restraints_manager`)
- Modify: `modules/qrefine/conformers.py` (remove `fragment_manager_proxy`, `from_conformers.fragment_manager`)
- Modify: `modules/qrefine/tests/unit/tst_31.py` (remove `exercise_proxy` and the dummies)
- Modify: `modules/qrefine/tests/unit/tst_46.py` (`frozen_environment` for variants, isinstance assertions)
- Modify: `modules/qrefine/tests/unit/tst_43.py` (assert the new route in the log)

**Interfaces:**
- Consumes: Tasks 1-3.
- Produces: `create_fragment_manager(model, params, working_folder, conformer_selections=None)`; `create_restraints_manager` returns `from_cluster` for clustered altloc models with `handling=conformers` and logs "Clustering with N conformers: per-fragment variants" (log line asserted by tst_43 and tst_68).

- [ ] **Step 1: Update tst_46 `frozen_environment` for the cluster route with variants**

For `from_cluster`: `self.jobs_of_atom[j]` = list of (job index, weight) of the variants whose master selection contains j; `energy(j, k, delta)` returns `sum(w * e_job)` over them (None if any cap moves); `exercise` uses `multiplicity = flex.double(den)` taken from the manager (`rm.multiplicity` attribute added to `from_cluster`, equal to `den` of the last call) when `rm` is a `from_cluster` of an altloc model. Change the `assert isinstance(rm, conformers.from_conformers) == (handling == "conformers")` to: expansion -> `from_conformers` iff conformers; clustering -> `from_cluster` always. `check_supercell` for the clustered conformer route: per variant is not available; skip for `from_cluster` with altlocs (the one-cluster limit of tst_65 covers it).

- [ ] **Step 2: Run tst_46; expect failure on the isinstance assertion (dispatch not changed yet)**

- [ ] **Step 3: Implement the dispatch (spec 3.1), remove the proxy, add the log line**

- [ ] **Step 4: Run tst_31, tst_43, tst_46 (cctbx and MOPAC), tst_65 fully; expect OK**

- [ ] **Step 5: Report the diff and stop**

---

### Task 5: duplicate-conformer identity through the clustered route (cctbx)

**Files:**
- Create: `modules/qrefine/tests/unit/tst_66.py`
- Modify: `modules/qrefine/run_tests.py`

**Interfaces:** consumes `tst_31.duplicate`, `tst_46.get_model`, `from_cluster.n_jobs`.

- [ ] **Step 1: Write tst_66.py**

```python
def exercise(labels):
  # h_altconf_complete without altlocs (remove_alt_confs) in the contact-30
  # box; single-conformer clustered gradients g_single (maxnum=3);
  # duplicate residue 2 into identical conformers (labels; occupancies
  # 1/len and then altered to 0.9/0.1: results unchanged): clustered
  # conformers route gives g == g_single.select(original) to 1e-9;
  # n_jobs == n_fragments + (len(labels)-1) * n_fragments_touching_residue_2
  ...
for labels in [["A","B"], ["A","B","C"]]: exercise(labels)
```

- [ ] **Step 2: Run; expect failure only if Tasks 2-4 are wrong (this test must pass on first run; if it fails, fix the code, not the test)**

- [ ] **Step 3: Register and run tst_66; expect OK**

- [ ] **Step 4: Report and stop**

---

### Task 6: MOPAC super-toy tests

**Files:**
- Create: `modules/qrefine/tests/unit/tst_67.py`
- Modify: `modules/qrefine/run_tests.py`

**Interfaces:** consumes `tst_46.get_mopac_params` (check its signature in tst_46 before use), `altloc_scenarios.base_model/side_chain_conformers/write`, `run_qrefine_in`.

- [ ] **Step 1: Write tst_67.py**

Model: `super_toy()` = `altloc_scenarios.base_model()` without waters, residue 3 side chain A/B (chi1 0 and 120, occupancies 0.5), P1 box with 30 A gap (no images), written to `super_toy.pdb`. Params: MOPAC PM7, `EPS=None` as tst_46 (`qm_engine.set(EPS=None)` on every `from_qm` of the route), `maxnum_residues_in_cluster=1`, `parallel.nproc=1`.

```python
def exercise_one_cluster_limit():
  # maxnum=999: clustered conformers route vs from_conformers over
  # expansion: same capped fragments -> gradients within 1e-4 kcal/mol/A
  # (engine noise), n_jobs == 2
def exercise_duplicate():
  # single-conformer super toy clustered (maxnum=1) vs the same with residue
  # 3 duplicated into identical A/B: gradients equal within 1e-4;
  # n_jobs == n_fragments + n_fragments_touching_residue_3
def exercise_opt():
  # qr.refine mode=opt restraints=qm quantum.engine_name=mopac
  # max_iterations_refine=3 on super_toy.pdb: clustering=False (b) and
  # clustering=True maxnum_residues_in_cluster=1 (c); final models within
  # rmsd 0.02 A; ase/<i>_A and ase/<i>_B folders exist for the fragments
  # containing residue 3 and ase/<i> for the others
```
Report the wall time of each exercise in the output.

- [ ] **Step 2: Run tst_67; must pass; if slower than two minutes, reduce `max_iterations_refine` to 2 and say so**

- [ ] **Step 3: Register; report and stop**

---

### Task 7: qr.refine end to end with cctbx: (a), (b), (c), re-clustering, H/D

**Files:**
- Create: `modules/qrefine/tests/unit/tst_68.py`
- Modify: `modules/qrefine/run_tests.py`

**Interfaces:** consumes `tst_06.write_inputs/r_work_of/check_refined/refine_fast`, `run_fmodel`, `run_qrefine_in`; the H/D-exchange-plus-heavy-altlocs fixture of tst_25 (its exercise whose docstring reads "H/D exchange (hd_mapper) together with heavy-atom altlocs", around line 251; reuse its model builder or PDB string).

- [ ] **Step 1: Write tst_68.py**

```python
def rmsd(model_1, model_2): ...   # same atom order, flex rms
def exercise_routes():
  # tst_06 data (good/poor, high_res=2): (a) native clustering=False,
  # (b) conformers clustering=False, (c) conformers clustering=True
  # maxnum_residues_in_cluster=3; refine_fast; check_refined on all three;
  # rmsd(b,c) < 0.05 A; rmsd(a,b) < 0.2 A and rmsd(a,c) < 0.2 A (print the
  # values; tighten to 2x measured if far below); (c) log contains
  # "per-fragment variants"; folders ase/<i>_A and ase/<i>_B exist and
  # ase/A does not
def exercise_reclustering():
  # (c) with re_calculate_rmsd_tolerance=0 number_of_refine_cycles=2:
  # log contains "need to redo clustering" at least twice; run completes
def exercise_hd_exchange():        # Review Focus 4
  # neutron-style model with H/D exchange plus heavy altlocs (tst_25
  # fixture), scattering_table=neutron, restraints=cctbx clustering=True:
  # completes, log "Using hd_mapper" and "per-fragment variants"
```

- [ ] **Step 2: Run tst_68; must pass (fix code if not); register; report and stop**

---

### Task 8: 6uf7 suite test

**Files:**
- Use: `modules/qrefine/tests/unit/data_files/6uf7_qmready.pdb` and `6uf7.mtz` (already in the suite; `mmtbx.reduce2 qmready=true` output of 6uf7, used by tst_63)
- Create: `modules/qrefine/tests/unit/tst_69.py`
- Modify: `modules/qrefine/run_tests.py`

- [ ] **Step 1: Invocation**: as tst_63 `exercise_6uf7_series` does (`qr.refine 6uf7.mtz 6uf7_qmready.pdb restraints=cctbx ...`, positional files, no label arguments), plus `refine_fast`.

- [ ] **Step 2: Write tst_69.py**

```python
def exercise_cctbx():
  # (a), (b), (c) as tst_68 with refine_fast on 6uf7_qmready.pdb + 6uf7.mtz;
  # R-work drops or stays within 0.002 of start in all three (ultra-high
  # resolution model is already good: assert r_end <= r_start + 0.002);
  # rmsd(b,c) < 0.05, rmsd(a,c) < 0.2; (c) log has "per-fragment variants"
  # and ase/<i>_A folders (ILE 4 A/B)
def exercise_mopac():
  # mode=opt restraints=qm mopac max_iterations_refine=2: (b) and (c)
  # maxnum_residues_in_cluster=2; rmsd(b,c) < 0.05; print wall times
```

- [ ] **Step 3: Run; register; report wall times and stop**

---

### Task 9: phenix.refine hook equivalence with clustering and altlocs

**Files:**
- Modify: `modules/phenix_regression/refinement/tst_177.py`

**Interfaces:** consumes the existing `compare(table, data_type, prefix)` flow of tst_177; `main.use_qrefine=True qrefine.restraints=cctbx qrefine.clustering=True`.

- [ ] **Step 1: Add a case**: X-ray data computed from `tst_06.pdb_str_good` (altlocs A/B), poor model `tst_06.pdb_str_poor`; phenix.refine with `qrefine.clustering=True` vs standalone `qr.refine ... clustering=True maxnum_residues_in_cluster=3` with the same settings tst_177 already aligns; equality to the tolerance tst_177 uses. The hook calls `api.refine(working_folder=None)`: assert the phenix.refine run wrote no `ase/` folder (Review Focus 5).

- [ ] **Step 2: Run tst_177; must pass; report and stop**

---

### Task 10: profiling and speed

**Files:**
- Possibly modify: `modules/qrefine/charges.py`, `modules/qrefine/restraints.py` (`from_qm.target_and_gradients`), `modules/qrefine/fragment.py`
- Create: `modules/qrefine/docs/superpowers/specs/2026-11-07-altloc-clustering-variants-design.md` section 7 "Results" (append)

- [ ] **Step 1: Profile 1PWM setup** (`rt3/1PWMH_charge3.pdb`, cctbx, clustering=True, `cProfile` around `create_restraints_manager`): record wall time, n_jobs (expect 304), top ten functions. Compare with the per-conformer path at the parent commit (`git stash` is not allowed to disturb the tree; run the old path from a second checkout or from the ledger's recorded numbers: 2 x 123 s setup).

- [ ] **Step 2: Profile one evaluation** (cctbx) and one MOPAC evaluation on 1PWM (`parallel.nproc=8`): time in capping, charges, `write_hierarchy`, engine. Refactor what exceeds 10 % of the evaluation outside the engine (candidates: `charges_class` construction per capper, PDB writes per call when the engine reads XYZ, `deep_copy` per variant per call) keeping results bit-identical (tst_64-66 unchanged).

- [ ] **Step 3: Acceptance runs**: 1PWM cctbx one cycle (`refine_fast`) before/after: setup time, evaluation time, final R-work vs the per-conformer run recorded on 2026-11-05 (0.1259 after one cycle from 0.1307); 6uf7 MOPAC refine mode (b) and (c). Append the numbers to the spec's Results section.

- [ ] **Step 4: Run the full qrefine suite** (`libtbx.run_tests_parallel module=qrefine nproc=4 run_in_unique_dirs=True` from an empty scratch directory); all OK. Report and stop.

---

### Task 11: documentation

**Files:**
- Modify: `modules/qrefine/CHANGES.md` (Unreleased / Changed: clustering with altlocs uses per-fragment conformer variants; one clustering; one parallel pool; folder layout `ase/<index>_<label>`; `from_conformers` keeps plain and expansion)
- Modify: `modules/phenix_html/rst_files/reference/AQuaRef.txt` lines 38-44: replace the "Static disorder ... is not supported" paragraph with: alternative conformations are supported on every route once the model is completed with `mmtbx.reduce2 qmready=true`; exchangeable H/D sites are handled for neutron data; unsupported entities (metal ions) still need `exclude`.
- Modify: `modules/qrefine/doc/qr/tests.html` rows for tst_64-69 if that file is maintained by hand (check for a generator first; if generated, run it).

- [ ] **Step 1: Make the edits; run `libtbx.find_clutter` on the changed qrefine files; report and stop**
