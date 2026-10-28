# Alternative conformations in qr.refine: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** qr.refine refines models with alternative conformations through
every restraints route (plain, expansion, clustering; cctbx or QM) by
splitting the model into single-conformer models and combining their
gradients.

**Architecture:** One new module, `qrefine/conformers.py`, splits a model
into complete single-conformer models and wraps one ordinary restraints
manager per conformer; `refine.create_restraints_manager` calls back into
itself per conformer with its own scratch folder. `qr.py` gains an `altlocs`
PHIL scope (`handling`, `reconcile`) and reconciles labels before restraints
are built. The old per-altloc code in `restraints.py` and `fragment.py` is
removed.

**Tech Stack:** cctbx/mmtbx (iotbx.pdb hierarchy, mmtbx.model.manager,
scitbx flex), libtbx (group_args, Sorry, easy_run), MOPAC via the ASE plugin
for QM tests. Python as used in qrefine (two-space indentation,
`from __future__` headers, `libtbx.test_utils.approx_equal`).

**Spec:** `docs/superpowers/specs/2026-10-28-altlocs-in-qr-refine-design.md`

## Global Constraints

- Changes stay uncommitted; the user commits (memory: user-commits-themselves).
- Restraints and gradients never use occupancies, including when combining
  conformers (spec 1).
- Any number of conformer labels is supported (spec 3.2): nothing may assume
  two; tests cover two, three and four labels.
- Every test runs in about a minute, at most two, and follows
  `tests/unit/tst_00.py` (PDB strings in the file where the model is new,
  `run(prefix=...)` that makes and enters `prefix`, `run_fmodel` /
  `run_qrefine` from `qrefine.tests.unit`, plain `assert`, `run()` under
  `__main__`). Existing data files in `tests/unit/data_files` may be read by
  path as `tst_04`/`tst_29` do.
- All tests listed in `run_tests.py` (`qr.test`) pass after every task; run
  at least the tests named in the task, and the whole list in Task 9.
- Tests are run as `cd <scratch dir>; qrefine.python <abs path>/tests/unit/tst_NN.py`
  from a clean directory; each test writes into its own `qrefine_tst_NN/`.
- MOPAC is the only QM engine used in tests (`quantum.engine_name=mopac`).
- Non-altloc runs keep the `ase/` layout (`ase/tmp_ase.pdb`, `ase/<i>/`).
- `energy` returned by `from_conformers` is 0, as `from_expansion` does.

## Review Focus

1. A model with three conformer labels (A/B/C) everywhere: multiplicity 3 for
   blank atoms, identity of the duplicate test still holds. Pinned in Task 2
   (tst_31, `exercise_duplicate` with three labels).
2. A model whose altloc residues have no blank atoms at all (labels only, as
   `altlocs.pdb`) through expansion and clustering: every atom has
   multiplicity 1 and each conformer model is a complete residue. Pinned in
   Task 8 (tst_42).
3. `refine.exclude` on a model with altlocs: the excluded atoms are removed
   before the split; excluding a whole residue keeps conformers complete and
   the run succeeds. Pinned in Task 4 (tst_39).
4. H/D exchange labels without heavy-atom altlocs: `n_conformers` is 1, only
   `hd_mapper` is used and the gradient is the one the previous code gave.
   Pinned in Task 6 (tst_25).
5. `altlocs.handling=native` with clustering: fragments carry both
   conformers and the gradient equals native plain, as it does without
   altlocs. Pinned in Task 5 (tst_05).

---

### Task 1: `conformers.py`: labels, completeness, split

**Files:**
- Create: `conformers.py`
- Create: `tests/unit/tst_31.py` (replaces the prototype; keep its
  `cctbx_calculator`, `exercise_vs_cctbx`, `exercise_no_altlocs`,
  `exercise_duplicate`, `exercise_incomplete_conformers` as the base)

**Interfaces:**
- Produces:
  - `conformer_labels(pdb_hierarchy) -> list[str]`: sorted altloc labels
    that have at least one atom whose `element.strip().upper()` is not
    `"H"`/`"D"`. Hydrogen-only labels are not conformers.
  - `n_conformers(pdb_hierarchy) -> int`: `len(conformer_labels(...))`.
  - `multiplicity(pdb_hierarchy) -> flex.double` of `atoms_size()`: 1.0 for
    an atom whose atom_group altloc is a conformer label, else
    `float(n_conformers)`.
  - `split(model) -> list[group_args]`, one per conformer label in
    `conformer_labels` order, each with `label` (str), `iselection`
    (`flex.size_t`, ascending, atoms with blank altloc, with this label, or
    with a hydrogen-only label), `model` (`mmtbx.model.manager`,
    single-conformer, atoms in `iselection` order). Raises
    `libtbx.utils.Sorry("Residue %s has altlocs %s but model has %s: conformers would be incomplete." % (rg.id_str(), sorted(labels), sorted(all_labels)))`
    when a residue group's conformer labels are a non-empty proper subset of
    the model's.

- [ ] **Step 1: Write the failing tests in `tests/unit/tst_31.py`**

Keep the prototype's helpers and add `exercise_split()`:

```python
def exercise_split():
  pi = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data, "h_altconf_complete.pdb"))
  model = mmtbx.model.manager(model_input=pi, log=null_out())
  ph = model.get_hierarchy()
  assert conformers.conformer_labels(ph) == ["A", "B"]
  assert conformers.n_conformers(ph) == 2
  m = conformers.multiplicity(ph)
  ci = ph.get_conformer_indices().conformer_indices
  assert m.select(ci == 0).all_eq(2.0) and m.select(ci != 0).all_eq(1.0)
  parts = conformers.split(model)
  assert [p.label for p in parts] == ["A", "B"]
  for p, c in zip(parts, [1, 2]):
    assert p.iselection.size() == ((ci == c) | (ci == 0)).count(True)
    assert list(p.iselection) == list(((ci == c) | (ci == 0)).iselection())
    h = p.model.get_hierarchy()
    assert list(h.altloc_indices().keys()) == [""]
    assert max(rg.atom_groups_size() for rg in h.residue_groups()) == 1
    assert h.atoms().extract_xyz().max_distance(
      ph.atoms().extract_xyz().select(p.iselection)) < 1.e-6
    assert list(h.atoms().extract_name()) == list(ph.atoms().select(p.iselection).extract_name())
  # hydrogen-only labels are not conformers
  rg = list(ph.residue_groups())[0]
  ag = rg.only_atom_group()
  h_atom = [a for a in ag.atoms() if a.name.strip() == "H"][0]
  ag_h = iotbx.pdb.hierarchy.atom_group(altloc="C", resname=ag.resname)
  ag_h.append_atom(h_atom.detached_copy())
  ag.remove_atom(h_atom)
  rg.append_atom_group(ag_h)
  ph.atoms().reset_i_seq()
  assert conformers.conformer_labels(ph) == ["A", "B"]
  assert conformers.multiplicity(ph)[ag_h.atoms()[0].i_seq] == 2.0
```

Keep `exercise_incomplete_conformers` but call `conformers.split(model)`
instead of `gradients_from_conformers` and keep the `"conformers would be
incomplete"` assertion.

- [ ] **Step 2: Run the test to verify it fails**

Run: `qrefine.python tests/unit/tst_31.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'qrefine.conformers'`
(import `from qrefine import conformers` at the top of the test).

- [ ] **Step 3: Implement `conformer_labels`, `n_conformers`, `multiplicity`, `split` in `conformers.py`**

`split`: tag `a.tmp = a.i_seq` on a deep copy of the hierarchy, select
`iselection` with `copy_atoms=True`, then for every residue group with more
than one atom group set every `altloc = ""` and merge into the first with
`residue_group.merge_atom_groups(primary, secondary)` +
`remove_atom_group(secondary)` (the pattern of
`hierarchy.remove_alt_confs`, which keeps atom order); `reset_i_seq()`,
`atoms_reset_serial()`; assert `extract_tmp_as_size_t()` equals
`iselection`. Build the model as
`mmtbx.model.manager(model_input=None, pdb_hierarchy=h, crystal_symmetry=model.crystal_symmetry(), restraint_objects=model.get_restraint_objects(), log=null_out())`
without processing (`get_current_pdb_interpretation_params` returns the
defaults when unprocessed, which is what `restraints.restraints` reads).

- [ ] **Step 4: Run the test to verify it passes**

Run: `qrefine.python tests/unit/tst_31.py`
Expected: `qrefine_tst_31:  OK`

---

### Task 2: `from_conformers` gradient combination

**Files:**
- Modify: `conformers.py`
- Modify: `tests/unit/tst_31.py`

**Interfaces:**
- Consumes: Task 1.
- Produces:
  - `class from_conformers(model, build, working_folder="ase", log=None)`.
    `build(conformer_model, working_folder) -> manager` is called once per
    conformer with `os.path.join(working_folder, label)`; the manager must
    have `target_and_gradients(sites_cart) -> (energy, flex.vec3_double)`.
    Attributes: `conformers` (the `split` list, each item gains
    `restraints_manager`), `multiplicity` (`flex.double`), `size` (int),
    `fragment_manager` (a `fragment_manager_proxy` when at least one
    conformer manager has a `fragment_manager` attribute, else `None`).
    Methods: `target_and_gradients(sites_cart) -> (0, flex.vec3_double)`
    with `g_i = (1/m_i) Σ g_{c,i}`; `energies_sites(sites_cart, compute_gradients=True) -> group_args(target=0, gradients=...)`.
  - `class fragment_manager_proxy(fragment_managers)`:
    `set_up_cluster_qm()` calls it on each; property `interaction_list`
    returns the concatenation.

- [ ] **Step 1: Write the failing tests**

Replace the prototype's `gradients_from_conformers` with
`from_conformers` in `exercise_vs_cctbx`, `exercise_no_altlocs`,
`exercise_duplicate` using

```python
def build_cctbx(conformer_model, working_folder):
  return restraints.from_cctbx(restraints_manager=restraints.get_cctbx_gradients(
    ph=conformer_model.get_hierarchy(), cs=conformer_model.crystal_symmetry()
    ).model.get_restraints_manager())

def gradients_from_conformers(model):
  fc = conformers.from_conformers(model=model, build=build_cctbx)
  return fc.target_and_gradients(sites_cart=model.get_sites_cart())[1]
```

(`restraints.from_cctbx.target_and_gradients` returns `(target, gradients)`.)
`exercise_no_altlocs` becomes: `n_conformers == 1`, so the caller
(`create_restraints_manager`, Task 3) never builds `from_conformers`;
assert `conformers.n_conformers(ph) == 1` and drop the gradient check.
Keep the assertions of `exercise_vs_cctbx` (altloc atoms equal cctbx, blank
atoms `g*m == g_cctbx + g_blank_only*(m-1)`, blank atoms not equal cctbx).
Add the helper `duplicate(model, labels, residue_indices) -> (dup_model, original)`
to `tst_31.py` (later tests import it): tag `a.tmp = a.i_seq`, for each
residue group index in `residue_indices` detach-copy its only atom group
once per label, give the copies the labels and occupancies
`1.0/len(labels)`, restore `tmp` on the copies (`detached_copy` resets it),
`atoms_reset_serial()`/`reset_i_seq()`, build a fresh `mmtbx.model.manager`
from the hierarchy, and return it with `original = flex.size_t([a.tmp for a in atoms])`.
Extend `exercise_duplicate(labels)` to run for `labels=["A","B"]`,
`labels=["A","B","C"]` and `labels=["A","B","C","D"]` (Review Focus 1)
with `residue_indices=[2]` and
assert `approx_equal(g, g_single.select(original))` and
`flex.max(conformers.multiplicity(ph)) == len(labels)`.
Add `exercise_proxy()`: two dummy managers with `fragment_manager` objects
having `set_up_cluster_qm` (counts calls) and `interaction_list=[(1,2)]`;
`from_conformers(...).fragment_manager.set_up_cluster_qm()` calls both and
`interaction_list == [(1,2),(1,2)]`; with managers lacking the attribute
`fragment_manager is None`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `qrefine.python tests/unit/tst_31.py`
Expected: FAIL with `AttributeError: module 'qrefine.conformers' has no attribute 'from_conformers'`

- [ ] **Step 3: Implement `from_conformers` and `fragment_manager_proxy` in `conformers.py`**

`target_and_gradients`: `g = flex.vec3_double(size, (0,0,0))`; per
conformer `g = g.set_selected(isel, g.select(isel) + g_c)`; return
`0, g * (1/self.multiplicity)`. Log one line per conformer at construction:
label, atom count, folder.

- [ ] **Step 4: Run the test to verify it passes**

Run: `qrefine.python tests/unit/tst_31.py`
Expected: `qrefine_tst_31:  OK`, under 30 s.

---

### Task 3: `create_restraints_manager` dispatch and scratch folders

**Files:**
- Modify: `qr.py` (`master_phil_str`: add the `altlocs` scope from spec 3.5
  after the `restraints` choice)
- Modify: `refine.py:63-88` (`create_fragment_manager`), `refine.py:150-212`
  (`create_restraints_manager`), `refine.py:389-390` (call site)
- Modify: `restraints.py:33-47` (`restraints.__init__`), `restraints.py:375-402`
  (`from_qm.__init__` folder), `fragment.py:294-296` (`os.mkdir` →
  `os.makedirs(..., exist_ok=True)`)
- Modify: `gtest.py:90-91`, `tests/unit/tst_04.py:35`, `tst_26.py:33`,
  `tst_27.py:68-71`, `tst_29.py:36`, `tst_32.py:50`, `tst_33.py:50` (drop
  `altlocs_present=False`)
- Create: `tests/unit/tst_11.py` (replaces the old average-method test)

**Interfaces:**
- Consumes: Task 2 `from_conformers`, `n_conformers`.
- Produces:
  - `refine.create_restraints_manager(params, model, working_folder="ase", log=None)`:
    when `conformers.n_conformers(model.get_hierarchy()) > 1` and
    `params.altlocs.handling == "conformers"` returns
    `conformers.from_conformers(model=model, build=lambda m, f: create_restraints_manager(params, m, f, log), working_folder=working_folder, log=log)`;
    otherwise the existing routes, with `working_folder` passed to
    `restraints.restraints(params, model, working_folder)` and
    `create_fragment_manager(model, params, working_folder)`.
  - `refine.create_fragment_manager(model, params, working_folder="ase")`:
    `fragments(working_folder=working_folder + "/", ...)` (keep the trailing
    slash the current code produces).
  - `restraints.restraints(params, model, working_folder="ase")`: passes
    `file_name=os.path.join(working_folder, "tmp_ase.pdb")` to `from_qm`;
    `from_qm` creates its folder with `os.makedirs(..., exist_ok=True)`.
  - PHIL `params.altlocs.handling` (`*conformers native`),
    `params.altlocs.reconcile` (bool, True).

- [ ] **Step 1: Write the failing test `tests/unit/tst_11.py`**

Model helper as in `tst_29.get_model`; `get_restraints_manager(file_name, expansion, clustering, handling)`
sets `params.restraints="cctbx"`, `params.altlocs.handling=handling`,
`params.cluster.maxnum_residues_in_cluster=999`, `params.parallel.nproc=1`
and returns `refine.create_restraints_manager(params=params, model=model)`.

```python
def exercise_routes_agree():
  """clustering vs expansion with altlocs, both handlings (old tst_25)."""
  for handling in ["conformers", "native"]:
    for fn in ["gly2_1.pdb", "gly2_2.pdb", "altlocs2.pdb", "altlocs.pdb", "h_altconf_complete.pdb"]:
      rm1 = get_restraints_manager(fn, expansion=False, clustering=True,  handling=handling)
      rm2 = get_restraints_manager(fn, expansion=True,  clustering=False, handling=handling)
      g1 = rm1.target_and_gradients(sites_cart=sites_cart)[1]
      g2 = rm2.target_and_gradients(sites_cart=sites_cart)[1]
      assert flex.max(flex.abs(g1.as_double()-g2.as_double())) < 1.e-6

def exercise_duplicate_routes():
  """A single-conformer model duplicated into A/B gives the single-conformer
  gradients through plain, expansion and clustering."""
  # single: h_altconf_complete.pdb with remove_alt_confs(always_keep_one_conformer=True)
  # dup: residue groups 2..4 duplicated (helper from tst_31) into
  # labels ["A","B"] and ["A","B","C"], occupancies 1/len(labels)
  for labels in [["A","B"], ["A","B","C"]]:
    for expansion, clustering in [(False, False), (True, False), (False, True)]:
      g_single = ...target_and_gradients(sites_cart=single.get_sites_cart())[1]
      g_dup    = ...target_and_gradients(sites_cart=dup.get_sites_cart())[1]
      assert approx_equal(g_dup, g_single.select(original))
      assert len(rm_dup.conformers) == len(labels)

def exercise_folders():
  """conformer scratch folders ase/A and ase/B are created for clustering."""
  rm = get_restraints_manager("h_altconf_complete.pdb", expansion=False, clustering=True, handling="conformers")
  assert os.path.isdir(os.path.join("ase", "A")) and os.path.isdir(os.path.join("ase", "B"))
  assert isinstance(rm, conformers.from_conformers) and rm.fragment_manager is not None
```

Import the helper: `from qrefine.tests.unit.tst_31 import duplicate`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `qrefine.python tests/unit/tst_11.py`
Expected: FAIL with `AttributeError: 'libtbx_phil_extract' object has no attribute 'altlocs'` or a `TypeError` on `altlocs_present`.

- [ ] **Step 3: Implement the PHIL scope, the dispatch and the folder plumbing**

In `create_restraints_manager` delete the `from_altlocs2` branch and the
`altlocs_present` parameter; keep everything else. Thread
`working_folder`. Update `gtest.py` and the six tests to the new signature.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `qrefine.python tests/unit/tst_11.py`, then `tst_04`, `tst_26`,
`tst_27`, `tst_29`, `tst_31`, `tst_32`, `tst_33`, `tst_35`.
Expected: all `OK`; `tst_35` still finds `ase/[0-9]*` (no altlocs: layout
unchanged).

---

### Task 4: `refine.run` composition, parameter checks, qr.refine end to end

**Files:**
- Modify: `refine.py:323-395` (`run`: hd_mapper gate, conformer checks,
  forced minimizer), `refine.py:293-321` (`check_qm_sees_crystal_neighbours`
  untouched)
- Modify: `qr.py:385-393` (`auto_cust`: delete the AIMNet2 altloc `Sorry`)
- Create: `tests/unit/tst_06.py` (replaces the old one), `tests/unit/tst_39.py`

**Interfaces:**
- Consumes: Task 3 `create_restraints_manager(params, model, working_folder, log)`.
- Produces, in `refine.run` right after `set_qm_defaults`:
  - `n_conf = conformers.n_conformers(model.get_hierarchy())`; if `n_conf > 1`
    and `params.altlocs.handling == "native"` and `params.restraints == "qm"`:
    `raise Sorry("altlocs.handling=native needs restraints=cctbx: a QM engine cannot handle alternative conformations. Use altlocs.handling=conformers.")`;
    if `n_conf > 1` and `handling == "conformers"`:
    `params.refine.minimizer = "lbfgs"; params.refine.gradient_only = True`
    and print `"  alternative conformations require minimizer=lbfgs and gradient_only=True"`.
  - hd_mapper gate: `use_hdm = len(model.get_hierarchy().exchangeable_hd_selections()) > 0`
    replaces `model.altlocs_present_only_hd()`; the `altlocs_present`
    variable disappears; `create_restraints_manager(params=params, model=model_for_restraints, log=log)`.

- [ ] **Step 1: Write the failing tests**

`tst_06.py` holds `pdb_str_good` = `data_files/h_altconf_2_complete.pdb`
and `pdb_str_poor` = `data_files/h_altconf_2_complete_modified.pdb`
(inline, as `tst_00`). Helper `r_work_of(pdb_file, r_fm)` builds an
`mmtbx.f_model.manager` as `tst_00` does and returns `r_work()` after
`update_all_scales()`.

```python
def run(prefix=...):
  """qr.refine with alternative conformations, cctbx restraints, both
  altlocs.handling routes; native with QM is refused."""
  r_fm = run_fmodel(prefix=prefix, args=[pdb_good, "high_res=2"])
  r_start = r_work_of(pdb_poor, r_fm)
  for handling in ["conformers", "native"]:
    r_qr = run_qrefine(prefix=prefix+"_"+handling, args=[pdb_poor, r_fm.mtz,
      "restraints=cctbx", "clustering=False", "altlocs.handling=%s" % handling,
      "output_folder_name=%s" % handling])
    h = r_qr.model.get_hierarchy()
    assert list(h.altloc_indices().keys()) == ["", "A", "B"]
    assert approx_equal(h.atoms().extract_occ(), poor_occupancies)
    r_end = r_work_of(r_qr.pdb, r_fm)
    assert r_end < 0.5 * r_start, (r_start, r_end)
    if handling == "conformers":
      assert r_qr.std_out_text.find("require minimizer=lbfgs") > -1
  r_qr = run_qrefine(prefix=prefix+"_sorry", sorry_expected=True, args=[pdb_poor, r_fm.mtz,
    "restraints=qm", "quantum.engine_name=mopac", "altlocs.handling=native", "dry_run=True"])
  assert r_qr.std_out_text.find("altlocs.handling=native needs restraints=cctbx") > -1
```

`run_qrefine` reads `pdb/*_refined.pdb` from the current directory, so run
each `qr.refine` in its own subdirectory (`os.makedirs(name); os.chdir(name)`
... `os.chdir("..")`) and drop the `output_folder_name` arguments; the same
convention applies to `tst_39`, `tst_41` and `tst_42`.

`tst_39.py` imports `pdb_str_good`, `pdb_str_poor`, `r_work_of` from
`tst_06` and runs the same assertions for
`["expansion=True", "clustering=False"]`,
`["clustering=True", "maxnum_residues_in_cluster=999", "parallel.nproc=1"]`
(asserts `os.path.isdir("ase/A")` after it), and
`["expansion=True", "clustering=False", "refine.exclude=resseq 99"]`
(Review Focus 3: residue 99 has no altlocs; the run succeeds and the R
factor drops).

- [ ] **Step 2: Run the tests to verify they fail**

Run: `qrefine.python tests/unit/tst_06.py`
Expected: FAIL: `assert r_qr.std_out_text.find("require minimizer=lbfgs") > -1`
(or earlier: `Sorry` text not found).

- [ ] **Step 3: Implement the `refine.run` changes and remove the AIMNet2 altloc `Sorry` in `qr.auto_cust`**

- [ ] **Step 4: Run the tests to verify they pass**

Run: `qrefine.python tests/unit/tst_06.py`; `qrefine.python tests/unit/tst_39.py`
Expected: both `OK`, each under 90 s (three `qr.refine` runs each).

---

### Task 5: conformers vs native (cctbx) comparison

**Files:**
- Create: `tests/unit/tst_05.py` (replaces the subtract test)

**Interfaces:**
- Consumes: Task 3 `create_restraints_manager`; `restraints.get_cctbx_gradients(ph, cs).gradients`;
  `calculator.sites_opt` + `scitbx.minimizers.lbfgs` as the old
  `tst_06.cctbx_opt(model, restraints_manager, max_shift=0.2, max_iterations=25)`
  did (copy that helper in).

- [ ] **Step 1: Write the test**

```python
def exercise_gradients():
  for fn in ["gly2_1.pdb", "gly2_2.pdb", "h_altconf_complete.pdb",
             "h_altconf_2_complete.pdb", "altlocs2.pdb", "altlocs.pdb"]:
    for expansion, clustering in [(False, False), (True, False), (False, True)]:
      g_c = gradients(fn, "conformers", expansion, clustering)
      g_n = gradients(fn, "native",     expansion, clustering)
      ci = ph.get_conformer_indices().conformer_indices
      altloc = ci != 0; blank = ci == 0
      assert approx_equal(g_c.select(altloc), g_n.select(altloc))
      if blank.count(True) == 0: continue
      m = conformers.multiplicity(ph).select(blank)
      g_blank_only = restraints.get_cctbx_gradients(ph=ph.select(blank), cs=cs).gradients
      assert approx_equal(g_c.select(blank) * m, g_n.select(blank) + g_blank_only * (m - 1))
      if not expansion and not clustering:
        assert approx_equal(g_n, restraints.get_cctbx_gradients(ph=ph, cs=cs).gradients)  # Review Focus 5 for clustering is the (False, True) iteration above

def exercise_minimization():
  """h_altconf_2_complete_modified.pdb minimized once per handling (plain
  route): the two results are close."""
  results = {}
  for handling in ["conformers", "native"]:
    model = get_model("h_altconf_2_complete_modified.pdb")
    rm = create_restraints_manager(params(handling), model)
    results[handling] = cctbx_opt(model=model, restraints_manager=rm)
  d = flex.sqrt((results["conformers"] - results["native"]).dot())
  print("mean/max |conformers-native|:", flex.mean(d), flex.max(d))
  assert flex.mean(d) < 0.1 and flex.max(d) < 0.3, (flex.mean(d), flex.max(d))
```

If the measured distances exceed the thresholds, report the numbers in the
task result instead of loosening them silently.

- [ ] **Step 2: Run the test**

Run: `qrefine.python tests/unit/tst_05.py`
Expected: `OK`, under 60 s. (The gradient part should pass as soon as
Task 3 is in; the minimization part is the measurement.)

---

### Task 6: H/D exchange with heavy-atom altlocs

**Files:**
- Create: `tests/unit/tst_25.py` (replaces the old clustering-with-subtract test)
- Modify: `refine.py` only if `hd_mapper` needs a fix for mixed models
  (`hd_mapper.__init__`, `refine.py:90-146`)

**Interfaces:**
- Consumes: `refine.hd_mapper(model)` with `get_single_model()`,
  `shrink(array)`, `expand(g_short)`; Task 3 `create_restraints_manager`.

- [ ] **Step 1: Write the test**

Build two models from `h_altconf_complete.pdb` in the test:

- `model_h`: unchanged (H only).
- `model_hd`: for the first three residue groups that have one blank atom
  group, move atom `H` out of the blank group into a new atom group
  `altloc="C"` (occ 0.5) and add a copy named `" D  "`, element `"D"`,
  same coordinates, in a new atom group `altloc="D"` (occ 0.5); tag
  `a.tmp` with the `model_h` i_seq (the D copy gets its H's i_seq).

```python
def exercise_mixed():
  assert conformers.conformer_labels(model_hd.get_hierarchy()) == ["A", "B"]
  hdm = refine.hd_mapper(model=model_hd)
  single = hdm.get_single_model()
  assert single.get_hierarchy().atoms().extract_element().count(" D") == 0
  rm = refine.create_restraints_manager(params=params, model=single)   # cctbx, conformers
  g = hdm.expand(rm.target_and_gradients(sites_cart=hdm.shrink(model_hd.get_sites_cart()))[1])
  rm_ref = refine.create_restraints_manager(params=params, model=model_h)
  g_ref = rm_ref.target_and_gradients(sites_cart=model_h.get_sites_cart())[1]
  original = flex.size_t([a.tmp for a in model_hd.get_hierarchy().atoms()])
  assert approx_equal(g, g_ref.select(original))

def exercise_hd_only():
  """Review Focus 4: H/D labels without heavy altlocs: n_conformers == 1."""
  # same construction on remove_alt_confs(always_keep_one_conformer=True) model
  assert conformers.n_conformers(ph) == 1
  assert refine.create_restraints_manager(...) is not a conformers.from_conformers instance
  # and g == g_ref as above
```

Also run `qr.refine` once on `model_hd` written to a file with
`restraints=cctbx mode=opt number_of_micro_cycles=1 max_iterations_refine=2`
and assert the output keeps the `C`/`D` atom groups (the log says
`Using hd_mapper`).

- [ ] **Step 2: Run the test**

Run: `qrefine.python tests/unit/tst_25.py`
Expected: `OK`. If `hd_mapper.expand` asserts (a gradient of exactly
`(0,0,0)`) or `exchangeable_hd_selections` pairs nothing, fix
`hd_mapper` for the mixed case and record what was wrong in the task
result.

---

### Task 7: MOPAC through every route

**Files:**
- Create: `tests/unit/tst_40.py`

**Interfaces:**
- Consumes: Task 3 `create_restraints_manager`; `duplicate` helper from
  `tst_31`; `data_files/m2_complete_box_large.pdb` (7 residues, QM-ready,
  box so large that the super sphere is the model itself).

- [ ] **Step 1: Write the test**

```python
def run(prefix=...):
  """Duplicate-conformer identity with MOPAC: plain, expansion, clustering.
  Conformer scratch folders ase/A and ase/B are used."""
  single = get_model("m2_complete_box_large.pdb")   # model_from_hierarchy-style processing
  dup, original = duplicate(single, labels=["A", "B"], residue_indices=[2, 3, 4])
  dup3, original3 = duplicate(single, labels=["A", "B", "C"], residue_indices=[2, 3, 4])
  # three labels: plain route only (3 MOPAC runs); asserted like the A/B case
  for expansion, clustering in [(False, False), (True, False), (False, True)]:
    params = qr.get_default_params(); params.restraints = "qm"
    params.quantum.engine_name = "mopac"; params.quantum.method = "PM7"; params.quantum.basis = ""
    params.expansion = expansion; params.cluster.clustering = clustering
    params.cluster.maxnum_residues_in_cluster = 2; params.parallel.nproc = 1
    g_single = refine.create_restraints_manager(params, single).target_and_gradients(single.get_sites_cart())[1]
    g_dup    = refine.create_restraints_manager(params, dup).target_and_gradients(dup.get_sites_cart())[1]
    assert approx_equal(g_dup, g_single.select(original), 1.e-4)
  assert os.path.isfile(os.path.join("ase", "A", "tmp_ase.pdb"))
  assert os.path.isdir(os.path.join("ase", "B", "0"))   # clustering fragments of conformer B
  assert os.path.isfile(os.path.join("ase", "C", "tmp_ase.pdb"))  # from the three-label plain run
```

Use `refine.set_qm_defaults(params, log)` instead of setting method/basis
by hand if that is what `qr.py` does for mopac.

- [ ] **Step 2: Run the test**

Run: `qrefine.python tests/unit/tst_40.py`
Expected: `OK` in under 90 s (21 MOPAC runs on 107 atoms). If slower,
drop the expansion route from the loop (plain already covers the
whole-model QM input) and say so in the task result.

---

### Task 8: reconcile labels in qr.py; labels-only model

**Files:**
- Modify: `qr.py:300-330` (`Program.run`: before `self.model.process(...)`)
- Create: `tests/unit/tst_41.py`, `tests/unit/tst_42.py`

**Interfaces:**
- Consumes: `iotbx.pdb.hierarchy.root.reconcile_altlocs(crystal_symmetry, log=...)`
  returning `group_args(n_sites, n_components, swaps, unresolved)`;
  `residue_group.relabel_altlocs(mapping)`; `params.altlocs.reconcile`.
- Produces: in `qr.Program.run`, after the model is fetched and before
  `process`: if `params.altlocs.reconcile` and
  `conformers.n_conformers(hierarchy) > 1`: `self.header("Reconciling altlocs")`,
  `r = hierarchy.reconcile_altlocs(crystal_symmetry=self.model.crystal_symmetry(), log=self.logger)`,
  then print `"  %d residues relabeled" % len(r.swaps)` and, if
  `r.unresolved`, `"  WARNING: altlocs of %d groups of residues cannot be made consistent; refinement continues with the labels as they are." % len(r.unresolved)`.

- [ ] **Step 1: Write `tst_41.py` (reconcile)**

In the test: load `h_altconf_complete.pdb`, assert
`reconcile_altlocs(cs, check_only=True).swaps == []`; then find the first
altloc residue group `rg` such that after `rg.relabel_altlocs({"A":"B","B":"A"})`
`reconcile_altlocs(cs, check_only=True)` reports exactly one swap on
`rg.id_str()` (undo and try the next one otherwise; assert one is found).
Write that model to `pdb_swapped`. Then:

```python
  common = [pdb_swapped, "restraints=cctbx", "mode=opt", "clustering=False",
            "number_of_micro_cycles=1", "max_iterations_refine=1"]
  r_qr = run_qrefine(prefix=prefix+"_on", args=common)        # in subdirectory "on"
  assert r_qr.std_out_text.find("Reconciling altlocs") > -1
  assert r_qr.std_out_text.find("1 residues relabeled") > -1
  assert labels_of(r_qr.pdb, rg_id) == labels_of(original_pdb, rg_id)   # restored
  r_qr = run_qrefine(prefix=prefix+"_off", args=common+["altlocs.reconcile=False"])  # in "off"
  assert r_qr.std_out_text.find("Reconciling altlocs") == -1
  assert labels_of(r_qr.pdb, rg_id) == swapped_labels            # untouched
```

`labels_of(pdb_file, rg_id)` returns the per-atom altloc list of that
residue group in atom order.

- [ ] **Step 2: Write `tst_42.py` (labels only, old tst_47 revived)**

`pdb_str` = `data_files/altlocs.pdb` inline. `run_fmodel(high_res=2)` on
it, then `run_qrefine([pdb, mtz, "restraints=cctbx", "clustering=False"])`:
asserts `r_work_of(r_qr.pdb) < 0.05` and output altlocs `["A", "B"]`. Then
in-process (Review Focus 2): `create_restraints_manager` with
`expansion=True` and with `clustering=True` on the same model returns a
`conformers.from_conformers` whose `multiplicity.all_eq(1.0)` and whose
`target_and_gradients` runs; both gradients agree with the plain route
to 1e-6 (one residue in a box: no neighbours inside the sphere; if the
cell makes neighbours, assert only that the call succeeds and say so).

- [ ] **Step 3: Run both to verify they fail**

Run: `qrefine.python tests/unit/tst_41.py`
Expected: FAIL at `"Reconciling altlocs"` not found.

- [ ] **Step 4: Implement the `qr.py` change**

- [ ] **Step 5: Run both to verify they pass**

Run: `qrefine.python tests/unit/tst_41.py`; `qrefine.python tests/unit/tst_42.py`
Expected: both `OK`.

---

### Task 9: remove the old altloc machinery

**Files:**
- Modify: `restraints.py` (delete `from_altlocs2` 208-272, `from_cctbx_altlocs`
  274-315, `from_altlocs` 317-335, the `rm_only` branch of
  `get_cctbx_gradients` 192-206 if nothing else uses it, the
  `gradients_scale` lines in `from_cctbx.target_and_gradients` 377-380 and
  `from_qm.target_and_gradients` 522-579)
- Modify: `fragment.py` (delete in `fragments.__init__` the
  `altloc_method`, `bond_with_altloc_flag`, `altloc_atoms` members 75-110,
  139-140; in `get_fragments` the `phs` per-altloc construction 188-202
  (keep `phs = [self.pdb_hierarchy_super]`), the overlap block 255-338 and
  `fragment_scales`; `bond_with_altloc` 363-375, `atoms_overlap` 377-380;
  `collect_cluster_and_fragment` no longer builds `scale_list`;
  `get_fragment_extracts` drops `fragment_scales`; `write_mm_charge_file`
  keeps only the `len(altlocs)==1` branch, 533-555)
- Modify: `cluster_restraints.py:15-27` (delete `check_no_altlocs` and
  its commented call), `refine.py:63-88` (drop the two arguments),
  `qr.py:71-82` (drop `altloc_method`, `bond_with_altloc`)
- Delete: `tests/unit/broken_to_fix/tst_43.py`, `tests/unit/broken_to_fix/tst_47.py`
- Modify: `run_tests.py` (add `tst_39`–`tst_42` to `tst_list`)

**Interfaces:**
- Consumes: everything above; no new interfaces.

- [ ] **Step 1: Grep for every use before deleting**

Run: `grep -rn "altloc_method\|bond_with_altloc\|fragment_scales\|gradients_scale\|from_altlocs\|from_cctbx_altlocs\|check_no_altlocs\|altloc_atoms\|atoms_overlap\|rm_only" --include=*.py . | grep -v broken_to_fix`
Expected after the deletions: no matches outside `docs/`.

- [ ] **Step 2: Delete the code and parameters listed above**

- [ ] **Step 3: Run the whole suite**

Run: `qrefine.python run_tests.py` (or `qr.test`), from a clean scratch
directory, with `MOPAC_COMMAND` set.
Expected: every test `OK`; report the per-test times and flag any test over
120 s.

---

### Task 10: profile the altloc path

**Files:**
- Modify: only files on the altloc path if a hot spot is found
  (`conformers.py`, `fragment.py:get_fragments` `copy.deepcopy(ph)` per
  cluster, `restraints.from_qm.target_and_gradients` PDB write per call)

- [ ] **Step 1: Profile one gradient evaluation per route**

Run with `qrefine.python -m cProfile -s cumtime` a script that calls
`create_restraints_manager` + `target_and_gradients` on
`h_altconf_complete.pdb` (cctbx, conformers) for plain, expansion and
clustering, and on the duplicated `m2_complete_box_large.pdb` with mopac
clustering.
Expected: a table of the five most expensive functions per route in the
task result.

- [ ] **Step 2: Fix hot spots inside the altloc path only**

Candidates from the spec (3.10): the per-cluster `copy.deepcopy(ph)` in
`fragments.get_fragments` (replace by one copy outside the loop if
`pair_interaction.run` does not mutate it beyond coordinates, verified by
`tst_04`/`tst_29`/`tst_36` still passing); the PDB file written on every
`from_qm` call when the engine takes XYZ. Anything outside the altloc path
is listed, not changed.

- [ ] **Step 3: Run the whole suite again**

Run: `qrefine.python run_tests.py`
Expected: every test `OK`; times reported.
