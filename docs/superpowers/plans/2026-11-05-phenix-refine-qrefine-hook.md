# phenix.refine qr.refine Hook Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `phenix.refine main.use_qrefine=True` refines coordinates in the last macro-cycle through `qrefine.api.refine`, refuses unsupported setups at startup, and leaves every other phenix.refine step and every `use_qrefine=False` run unchanged.

**Architecture:** A self-contained PHIL scope in phenix.refine is translated into a finalized qrefine extract by `xyz_qrefine.make_params`. Startup checks live in `Program.validate` (import, engine, joint mode) and in `fsr.preprocess.qrefine_checks` (twin, ias, den, qi, amber, selections, solvent auto-off, `api.check_model`). The task table gets a `qrefine_cycle` flag that forces the riding-H and cctbx-geometry coordinate steps off in the last macro-cycle, and `macro_cycle.xyzrec` branches to `xyz_qrefine.run`, which calls the API with `working_folder=None`, applies the returned sites, re-scales with the same call qrefine used and asserts the R factors (D13).

**Tech Stack:** phenix (phenix/refinement, phenix/programs/phenix_refine.py, __init__.params), qrefine.api, phenix_regression/refinement, phenix_html.

**Spec:** `qrefine/docs/superpowers/specs/2026-11-05-phenix-refine-qrefine-hook-design.md` (umbrella `2026-11-04-qrefine-in-phenix-refine-design.md`)

## Global Constraints

- Repositories: `modules/phenix` (code, PHIL), `modules/phenix_regression` (tests), `modules/phenix_html` (docs). qrefine is imported lazily inside functions only (phenix.refine must start without it when `use_qrefine` is False). The user commits; leave changes uncommitted. Never `rm -rf`; scratch in fresh `rt/run_$(date +%s)` folders under `/Users/pafonine/Desktop/all/phenix`.
- Python: `/Users/pafonine/Desktop/all/phenix/build/bin/libtbx.python`; dispatchers by full path from `/Users/pafonine/Desktop/all/phenix/build/bin/`. After editing `__init__.params` no refresh is needed (read at import); after adding a new module no refresh either.
- Runs with `use_qrefine=False` are untouched: the only shared code changed is the task table (one flag, nine `and not self.qrefine_cycle`), `macro_cycle.xyzrec` (one early branch), `preprocess.run` (one call guarded by `use_qrefine`), `Program.validate` (one guarded block) and the PHIL.
- PHIL exactly as spec section 3 (names, choices, defaults, help). `aquaref` -> qrefine engine `aimnet2`. `clustering`/`max_bond_rmsd` None means qrefine's default; `max_bond_rmsd` is applied after `finalize_params`.
- The API call: `api.refine(model, qp, fmodel=fmodels.fmodel_xray(), working_folder=None, intermediate_folder=None, bulk_solvent_params=params.bulk_solvent_and_scale, log=log)`; afterwards sites applied, `update_xray_structure(update_f_calc=True, update_f_mask=True, force_update_f_mask=True)`, `model.set_xray_structure`, `update_all_scales(params=params.bulk_solvent_and_scale, remove_outliers=False, refine_hd_scattering=False)`, `create_target_functors()`, `approx_equal(r_work, result.r_work, 1e-4)` and r_free likewise (D13).
- Gated off in the qrefine cycle: setrh, regHxyz, ncsrot, realsrg, realsrl, realsrl2, sacart, tardy, den. Everything else keeps its rule.
- Refusals (Sorry) and the solvent auto-off exactly as spec section 5; every Sorry names the parameter to change.
- Tests write only into their own folder named after the prefix (tst_173 style: recreate the folder, chdir into it); MOPAC-dependent tests skip with a printed note when `qrefine.api.check_engine` raises.
- cctbx style (2-space indent, 80 columns where feasible, `print(..., file=log)`, `Sorry` for user errors).

## Review Focus

1. A user who sets `use_qrefine=True` without qrefine installed must get the one-line Sorry before any data processing, not an ImportError traceback. Test `test_no_qrefine_import` in Task 1 (simulated by a `sys.modules["qrefine"] = None` injection in a unit-level call of the validate block, or by an `importlib` stub).
2. `number_of_macro_cycles=1`: the only cycle is the qrefine cycle; no riding manager is ever created; the output H lengths are X-ray. Covered by tst_177 (Task 4) which runs one cycle.
3. A `refine.sites.individual` selection: only the selected atoms move; the exclude string built from several selection strings parses. Test `test_selection_moves_only_selected` in tst_178 (Task 4).
4. Neutron data with riding H (`hydrogens.refine=riding`, exchangeable H/D): the qrefine cycle must not fail on `refine.exclude` (none given) and must not re-idealize H afterwards. Covered by tst_177 neutron (Task 4) with `hydrogens.refine=riding` forced.
5. `main.ordered_solvent=True` with `use_qrefine`: auto-off with the boxed message and no water picking in any cycle (the log has no "Update (add/remove/refine) ordered solvent" run). Test `test_ordered_solvent_auto_off` in tst_178 (Task 2).

---

### Task 1: PHIL, `xyz_qrefine.make_params`, `Program.validate` check

**Files:**
- Modify: `phenix/phenix/refinement/__init__.params` (`main` scope ~line 291: add `use_qrefine` after `nqh_flips`; new `qrefine` scope after `bulk_solvent_and_scale` ~line 745)
- Create: `phenix/phenix/refinement/xyz_qrefine.py`
- Modify: `phenix/phenix/programs/phenix_refine.py` (`Program.validate`, ~line 360, first statements)
- Create: `phenix_regression/refinement/tst_178.py` (started here with the validate-level tests; extended in Tasks 2 and 4)

**Interfaces:**
- Produces: `xyz_qrefine.make_params(params, model, fmodel, log=None) -> qrefine extract` (finalized; spec section 4); `xyz_qrefine.engine_name(params) -> str` (`aquaref` -> `aimnet2`); `xyz_qrefine.raw_params_for_engine_check(params) -> qrefine extract` (unfinalized, `restraints`, `quantum.engine_name`, `quantum.server_url` set); `xyz_qrefine.check_qrefine_available()` raising `Sorry("main.use_qrefine=True needs the qrefine module, which is not installed in this Phenix")`.
- Consumes: `qrefine.qr.get_default_params`, `qrefine.api.finalize_params`, `qrefine.api.check_engine`, `qrefine.conformers.n_conformers`.

- [ ] **Step 1: Write the failing tests** (tst_178, first part; run through `run_phenix_refine` with `sorry_expected=True` where a Sorry is expected)

```python
def test_phil_present():
  # phenix.refine --show-defaults (or master_params()) contains
  # "use_qrefine = False" and a "qrefine" scope with restraints, qm_engine,
  # clustering, max_bond_rmsd
def test_make_params():
  # unit level: params = phenix.refinement.master_params().extract();
  # params.refinement.main.use_qrefine = True; qm_engine=aquaref -> qp.quantum.engine_name == "aimnet2",
  # qp.finalized, qp.refine.max_bond_rmsd == the aquaref preset for d_min 1.5 (0.02);
  # with qrefine.max_bond_rmsd=0.05 -> 0.05 after finalize; clustering=None keeps the preset (False);
  # sites.individual=["chain A"] on a model with chains A and S -> qp.refine.exclude == "not ((chain A))"
def test_joint_refused():
  # run_phenix_refine with an x-ray model/data and a neutron model/data (joint) and use_qrefine=True:
  # Sorry "joint"
def test_no_qrefine_import():
  # xyz_qrefine.check_qrefine_available() with sys.modules["qrefine"] = None (restored in finally) -> Sorry naming qrefine
def test_engine_unavailable():
  # run_phenix_refine with use_qrefine=True qrefine.qm_engine=server qrefine.restraints=qm and
  # an unreachable server url (qrefine's default http://127.0.0.1:8000 unless something listens; use
  # a fake engine route documented in the test): Sorry "not available"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/phenix_regression/refinement/tst_178.py`
Expected: FAIL at `test_phil_present` (unknown parameter `use_qrefine`).

- [ ] **Step 3: Implement**

PHIL per spec section 3. `xyz_qrefine.py` with the three functions (imports of qrefine inside the functions). `Program.validate`: at its top, when `self.params.refinement.main.use_qrefine` (single mode; in joint mode the flag sits under `self.params.xray.refinement.main`, test both): `check_qrefine_available()`, `if self.joint_xn: raise Sorry(...)`, `api.check_engine(raw_params_for_engine_check(params))`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; then `libtbx.python .../phenix_regression/refinement/tst_173.py` from a scratch folder (unchanged behaviour with the flag off). Expected: OK, OK.

- [ ] **Step 5: Checkpoint**

`git -C modules/phenix status --short` shows `__init__.params`, `phenix_refine.py`, `xyz_qrefine.py`; phenix_regression shows tst_178.py. Do not commit.

---

### Task 2: `preprocess.qrefine_checks`

**Files:**
- Modify: `phenix/phenix/refinement/fsr/preprocess.py` (new function; call in `run` after `set_refinement_flags`, before `set_TLSplusADP_params`)
- Modify: `phenix_regression/refinement/tst_178.py`

**Interfaces:**
- Produces: `preprocess.qrefine_checks(params, model, fmodels, log)` (spec section 5 item 2), called only when `params.main.use_qrefine`.
- Consumes: Task 1 `make_params`; `qrefine.api.check_model`; `qrefine.refine.exchangeable_hd_pairs`; `mmtbx.geometry_restraints.quantum_interface.is_quantum_interface_active`.

- [ ] **Step 1: Write the failing tests** (tst_178, second part; each a `run_phenix_refine(..., sorry_expected=True)` on the X-ray fixture with `use_qrefine=True qrefine.restraints=cctbx` unless stated)

```python
def test_twin_refused():          # twin_law="-h,-k,l" (any valid law for the P1 cell: use a P2 fixture
                                  # or pass the data manager's twin_law parameter): Sorry "twin"
def test_ias_refused():           # main.ias=True: Sorry "ias"
def test_den_refused():           # strategy=individual_sites+den: Sorry "den"
def test_selection_with_clustering_refused():   # refine.sites.individual="chain A" qrefine.clustering=True: Sorry "clustering"
def test_uncapped_refused():      # stripped tst_173 fragment (remove=capping_h) with qrefine.restraints=qm qm_engine=mopac
                                  # (skip without MOPAC): Sorry "capping=true"
def test_ordered_solvent_auto_off():   # main.ordered_solvent=True: run succeeds; log has the boxed "Automatic adjustment"
                                  # message and no "ordered solvent" step executed in any cycle
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: the Task 1 Step 2 command. Expected: FAIL at `test_twin_refused` (the run proceeds instead of a Sorry).

- [ ] **Step 3: Implement**

`qrefine_checks` in the order of spec section 5; the boxed message in the style of `preprocess.hydrogens` ("*"*79 lines). The call in `run`: `if(params.main.use_qrefine): qrefine_checks(params=params, model=model, fmodels=fmodels, log=log)`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Task 1 Step 2 command; then tst_173. Expected: OK, OK.

- [ ] **Step 5: Checkpoint**

Do not commit.

---

### Task 3: task table flag, `xyz_qrefine.run`, the `xyzrec` branch; tst_175

**Files:**
- Modify: `phenix/phenix/refinement/task_table.py` (`update` sets `self.qrefine_cycle`; `run_setrh` 91, `run_ncsrot` 108, `run_realsrg` 150, `run_realsrl` 157, `run_realsrl2` 178, `run_den` 203, `run_tardy` 206, `run_sacart` 214, `run_regHxyz` 281 get `and not self.qrefine_cycle`)
- Modify: `phenix/phenix/refinement/macro_cycle.py` (`xyzrec` 429: early branch)
- Modify: `phenix/phenix/refinement/xyz_qrefine.py` (`run`)
- Create: `phenix_regression/refinement/tst_175.py`

**Interfaces:**
- Produces: `task_table.run.qrefine_cycle` (bool, set in `update`, False when `use_qrefine` is off); `xyz_qrefine.run(model, fmodels, params, prefix, log)` (spec section 7); the task-table printout line for the step reads "XYZ individual (reciprocal space) [qr.refine]" in that cycle (find where the table prints step names and append the tag when the flag is set).
- Consumes: Task 1 `make_params`; `qrefine.api.refine`, `qrefine.api.check_engine`.

- [ ] **Step 1: Write the failing test** (tst_175; fixture `in.pdb` from `mmtbx.hydrogens.tst_xh_mapper.pdb_str`, data `run_fmodel(args=[in.pdb, "high_res=1.5"], ...)`)

```python
def run():
  r = run_phenix_refine(args=[pdb, mtz, "main.use_qrefine=True", "qrefine.restraints=cctbx",
        "main.number_of_macro_cycles=3", "main.nqh_flips=False", "strategy=individual_sites+individual_adp"], prefix=prefix)
  log = open(r.log).read()
  cycles = split_by_macro_cycle(log)          # on the "macro-cycle" headers print_statistics writes
  assert "qr.refine" in cycles[3] and "qr.refine" not in cycles[1] + cycles[2]
  assert "Ensure riding H" in cycles[2] and "Ensure riding H" not in cycles[3]
  assert "Local real-space refinement" in cycles[1] and "Local real-space refinement" not in cycles[3]
  # output X-H lengths equal X-ray targets (model processed with use_neutron_distances=False) within 1e-3
  # r_work final < r_work start (from_pdb_or_log)
  # the qrefine summary lines present ("data weight", "n_fev"); no "AssertionError" / "Traceback"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/phenix_regression/refinement/tst_175.py`
Expected: FAIL at the "qr.refine" assertion (standard xyz ran in cycle 3).

- [ ] **Step 3: Implement**

Task table flag and gating; `xyzrec` branch; `xyz_qrefine.run` per spec section 7 steps 1-9 (summary header format copied from `xyz_refinement_scorer._show`; the caps note when `model.selection("capping_h").count(True) > 0`).

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; tst_178; tst_173; tst_174. Expected: all OK.

- [ ] **Step 5: Checkpoint**

Do not commit.

---

### Task 4: tst_176 (MOPAC), tst_177 (equivalence), tst_178 (rest)

**Files:**
- Create: `phenix_regression/refinement/tst_176.py`, `tst_177.py`
- Modify: `phenix_regression/refinement/tst_178.py`

**Interfaces:** none new.

- [ ] **Step 1: Write the tests**

tst_176: tst_175's run with `qrefine.restraints=qm qrefine.qm_engine=mopac qrefine.clustering=False main.number_of_macro_cycles=1`; skip with a note when `api.check_engine` raises for mopac; same assertions as tst_175 for the one cycle plus "scratch folder" absent from the log? No: with `working_folder=None` and MOPAC the API logs a temporary folder and removes it: assert the named folder does not exist afterwards.

tst_177: for table in [("n_gaussian", "xray"), ("electron", "electron"), ("neutron", "neutron")]: data with `run_fmodel(args=[pdb, "high_res=1.5", "scattering_table=%s" % table])` (check `mmtbx.command_line.fmodel` accepts that argument; else generate the X-ray data and pass the table only to the refinements); phenix.refine args: `main.use_qrefine=True qrefine.restraints=cctbx main.number_of_macro_cycles=1 strategy=individual_sites main.scattering_table=<table> main.nqh_flips=False main.fix_cablam=off main.target=ml` plus for neutron `hydrogens.optimize_scattering_contribution=False hydrogens.real_space_optimize_x_h_orientation=False hydrogens.refine=riding`; qr.refine args (through `qrefine.tests.unit.run_qrefine_in`): `restraints=cctbx expansion=False clustering=False scattering_table=<table> number_of_macro_cycles=1` with the same `max_iterations_refine`, `number_of_refine_cycles`, `number_of_weight_search_cycles` as qrefine's defaults (phenix.refine's hook uses the defaults). Compare `model.get_sites_cart()` of the two output models (atom order identical: both come from the same file; assert names equal) within 1e-4 A and the final r_work/r_free (phenix.refine from `from_pdb_or_log`; qr.refine from its log's "Best Rwork, Rfree (after refinement)") within 1e-5. Print the measured maxima. If the gate fails, do NOT loosen it silently: report with the numbers (the controller rules).

tst_178 additions: `test_selection_moves_only_selected` (`refine.sites.individual="resseq 4"` with clustering off: atoms outside resseq 4 unchanged within 1e-6 A, inside moved), `test_caps_kept` (capped tst_173 model with `qrefine.restraints=cctbx`: output keeps HCG at occupancy 0, log has the strip note), `test_one_macro_cycle` (`number_of_macro_cycles=1`: "qr.refine" in cycle 1, no "Ensure riding H" anywhere).

- [ ] **Step 2: Run the new tests**

Run each from a fresh scratch folder. Expected: tst_176 OK or skipped with the note; tst_177 OK with printed maxima; tst_178 OK.

- [ ] **Step 3: Checkpoint**

Do not commit.

---

### Task 5: Documentation and the existing-suite sample

**Files:**
- Modify: `phenix_html/rst_files/reference/refinement.pre` (new subsection after "Refinement of coordinates" content, before "Refinement of Atomic Displacement Parameters")
- Modify: `phenix_html/rst_files/reference/AQuaRef.txt` (paragraph "Inside phenix.refine")

- [ ] **Step 1: Write the documentation** per spec section 9: purpose, the last-macro-cycle rule and what runs in it, preparation (`mmtbx.reduce2 qmready=true [capping=true]`), parameters with defaults, refusals, solvent auto-off, caps note, engine availability (aquaref needs the AIMNet2 install; MOPAC from conda-forge), example command lines (`phenix.refine model.pdb data.mtz main.use_qrefine=True` and with `qrefine.qm_engine=mopac`).

- [ ] **Step 2: Run the existing-suite sample**

Run from a fresh scratch folder: tst_173, tst_174, and the first ten refinement tests that exist among `tst_001`..`tst_020` (skip any that needs data not in this checkout, and say so). Expected: all OK with the flag off.

- [ ] **Step 3: Report**

Files changed in the three repositories; test results including tst_177's measured maxima; anything skipped.
