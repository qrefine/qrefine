# qrefine Callable API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** qr.refine becomes `qrefine.api.refine(model, params, fmodel, map_data, ...)` returning a result object, with `finalize_params`, `check_engine` and `check_model` callable before any heavy work; the command line keeps its behaviour and output layout by calling the same functions.

**Architecture:** `refine.run` keeps the orchestration but takes a working folder and an optional intermediate folder, writes nothing else, and returns a `group_args`; `api.refine` wraps it with deep copies of the caller's model and fmodel. Parameter rules move out of `qr.Program` and `refine.run` into `api.finalize_params`; the checks into `api.check_engine` / `api.check_model`. qr.py splits into `prepare_inputs` (data manager to objects), the API calls, and output writing through `model_as_pdb_or_mmcif`.

**Tech Stack:** qrefine (api.py new; qr.py, refine.py, driver.py, calculator.py, results.py, restraints.py, fragment.py, cluster_restraints.py, super_cell.py, conformers.py), cctbx (mmtbx.model, iotbx.pdb hierarchy writers, libtbx), qrefine tests/unit.

**Spec:** `qrefine/docs/superpowers/specs/2026-11-05-qrefine-api-design.md` (umbrella: `2026-11-04-qrefine-in-phenix-refine-design.md`)

## Global Constraints

- Repository: `modules/qrefine` only (a chem_data monomer file is read, never written). The user commits; leave changes uncommitted. Never run `rm -rf`; scratch goes in a fresh `rt/run_$(date +%s)` folder under `/Users/pafonine/Desktop/all/phenix`.
- Python: `/Users/pafonine/Desktop/all/phenix/build/bin/libtbx.python`; dispatchers by full path from `/Users/pafonine/Desktop/all/phenix/build/bin/`.
- Command-line numbers must not move: tst_55, tst_00, tst_01, tst_53 and the rest of the suite give the same results after every task; the `pdb/` and `ase/` layout of `qr.refine` stays.
- `api.refine` never modifies the caller's model or fmodel; all its writes stay under `working_folder` (engine scratch, fragments, expansion, conformer sub-folders, debug pickles, dump_gradients) and, when given, `intermediate_folder` (per-cycle models). `log=None` means `null_out()`.
- Every `update_all_scales` call inside the API: `params=bulk_solvent_params, remove_outliers=False, refine_hd_scattering=False` and no back-trace.
- `params.finalized` is set only by `finalize_params`; `refine.run` and `api.refine` assert it. No function other than `finalize_params` mutates params.
- `rst_file` and `max_atoms` are removed from the PHIL and from every use (D16, D17); `dump_gradients` ends with `Sorry` after dumping (D18); aimnet2 species read from a fresh `AIMNet2Calculator` (D19); the result carries the API's own refined model (D20).
- `working_folder=None` (D21): qrefine writes no file at all; in-process engines (aimnet2, aimnet2-old, torchani, fairchem, pyscf) and cctbx restraints touch the disk nowhere; external-program engines (mopac, xtb, orca, terachem, turbomole, gaussian) run in a `tempfile.mkdtemp()` folder removed in a `finally` when `refine` returns. `api.ENGINES = {name: group_args(needs_files=bool, lookup=callable, hint=str, elements=set or None)}` is the one table for engine facts.
- Sorry messages name the offenders (up to ten) and the preparation command (`mmtbx.reduce2 qmready=true [capping=true]`) or the install hint.
- Every hierarchy or model write goes through `write_pdb_or_mmcif_file` / `model_as_pdb_or_mmcif`; the returned file name is the one used afterwards.
- cctbx style: 2-space indent, 80 columns where feasible, `print(..., file=log)`, `Sorry` for user errors, docstrings on new public functions. Every qrefine test writes only into its own `qrefine_tst_XX` folder and is registered in `run_tests.py`.

## Review Focus

1. Two consecutive `api.refine` calls in one process with the same `working_folder` (sub-project 3 may call per run on the same prefix): identical results, no stale engine scratch interfering. Test `test_called_twice` in Task 4 (tst_56).
2. A caller whose fmodel shares the model's `xray_structure` object (phenix.refine does): after `api.refine` the caller's model sites are unchanged. Test `test_inputs_untouched` in Task 4 builds the fmodel from `model.get_xray_structure()` itself.
3. A caller's model with ligand restraint objects: the API copy must carry them, or processing fails on the ligand. Covered by tst_57 (Task 5) whose ligand has its own CIF, run through the API path (the command line uses the API).
4. `intermediate_folder` that does not exist yet: created, cycle files inside, nothing elsewhere. Test `test_intermediate_folder` in Task 4.
5. A water with one H (partial completion) must be refused by `check_model` with QM restraints, not silently accepted as "has H". Test `test_water_one_h` in Task 3 (tst_58).

---

### Task 1: Plumbing — log, working folder, intermediates, debug outputs

**Files:**
- Modify: `qrefine/driver.py` (module-level `log = sys.stdout` line 13; `refine` line 56, `opt` line 207; prints at lines 68-234; `monitor.write_pdb_file` calls lines 195, 225; debug pickle line 239)
- Modify: `qrefine/calculator.py` (`sites.target_and_gradients` dump_gradients block lines ~447-457)
- Modify: `qrefine/results.py` (`write_pdb_file` 93; keep `finalize` for now)
- Modify: `qrefine/refine.py` (`create_fragment_manager` 65, `create_restraints_manager` 227, `run` 406: thread `working_folder`, `intermediate_folder`; drop the `rst_file` path code lines ~456-465 but keep the PHIL parameter until Task 2)
- Modify: `qrefine/restraints.py` (`from_qm` default `file_name` 185/232; keep behaviour)
- Modify: `qrefine/fragment.py` (`expansion_file` line 149-153: into `self.working_folder`)
- Modify: `qrefine/cluster_restraints.py:30`, `qrefine/super_cell.py:71` (callers pass a path under the working folder)
- Modify: `qrefine/qr.py` (`run`: pass `working_folder="ase"`, `intermediate_folder=self.params.output_folder_name`)
- Modify: `qrefine/tests/unit/tst_36.py` (debug pickle path now under `ase/`)

**Interfaces:**
- Produces: `driver.refine(fmodel, params, monitor, calculator, geometry_rmsd_manager, log, intermediate_folder=None)`, `driver.opt(params, monitor, calculator, log, working_folder, intermediate_folder=None)`; `refine.run(model, fmodel, map_data, params, prefix, log, working_folder="ase", intermediate_folder=None)` (still writes the final file through `monitor.finalize` until Task 4); `calculator.sites(..., dump_gradients=None, working_folder=None)`: dumps `<working_folder>/<dump_gradients>_{dg,rg,g}` then raises `Sorry("dump_gradients: gradients dumped to %s_{dg,rg,g}, stopping.")`; `driver.opt` debug pickle at `<working_folder>/debug_optimization_clustering<bool>.pkl`; `fragment.fragments` writes `expansion.pdb` into its working folder.
- Consumes: nothing new.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/tst_36.py`: load the two debug pickles from `os.path.join(r_qr.folder, "ase", ...)` instead of `r_qr.folder`. Add a new test file `tests/unit/tst_59.py` (plumbing; numbers 56-58 are reserved for Tasks 2-5; register every new test in numeric order):

```python
def test_no_cwd_writes():
  # qr.refine on the tst_xh_mapper fixture, restraints=cctbx, expansion=False,
  # clustering=False, 1 macro cycle, from folder "run": afterwards the folder
  # contains exactly: the input files, "pdb/" and "ase/" (sorted listing).
def test_dump_gradients_sorry():
  # qr.refine with dump_gradients=g.pkl restraints=cctbx: rc != 0, "Sorry"
  # and "gradients dumped to" in the output; ase/g.pkl_dg, _rg, _g exist.
def test_expansion_file_in_working_folder():
  # qr.refine with clustering=True restraints=cctbx expansion=False (cctbx
  # fragments, no QM) on the fixture: "ase/expansion.pdb" exists and no
  # "expansion.pdb" in the run folder.
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/qrefine/tests/unit/tst_59.py`
Expected: FAIL at `test_expansion_file_in_working_folder` (expansion.pdb in the run folder) and `test_dump_gradients_sorry` (rc 0, no Sorry).

- [ ] **Step 3: Implement**

Remove the module-level log from driver.py and add the `log` argument to `refine` and `opt`; replace `params.output_folder_name` in the two `write_pdb_file` calls by `intermediate_folder` guarded by `if intermediate_folder is not None`; thread `working_folder` through `refine.run` to `create_restraints_manager`, `create_fragment_manager`, `conformers.from_conformers`, `calculator.sites`, `driver.opt`; move the debug pickle and the dump files; in `fragment.fragments._expand` write `expansion.pdb` under `self.working_folder`; delete the restart-path block of `refine.run` (the PHIL parameter stays until Task 2). `results.manager.write_pdb_file` creates the folder if missing.

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; then `tst_36.py`, `tst_55.py`, `tst_53.py`, `tst_00.py` each from a fresh scratch folder; then `grep -n "open(\|write_pdb_file\|write_pdb_or_mmcif_file\|easy_pickle.dump\|os.getcwd\|makedirs\|mkdir" *.py` from `qrefine/` and confirm every hit uses a path derived from a working or intermediate folder (list the hits in the report). Expected: all OK; grep list justified.

- [ ] **Step 5: Checkpoint**

`git status --short` shows the files above plus tst_59.py and run_tests.py. Do not commit.

---

### Task 2: `api.finalize_params`, `rst_file` removal, `Auto` target

**Files:**
- Create: `qrefine/api.py`
- Modify: `qrefine/qr.py` (`master_phil_str`: `refinement_target_name = *Auto ml mlhl ls_wunit_k1` line 130, delete `rst_file` lines 272-275 and `max_atoms` lines ~19-21; `Program.validate` 301-328 keeps the data-manager counts only; `Program.auto_cust` 330-390 deleted; `run` calls `api.finalize_params` and resolves `Auto` before `set_target_name`)
- Modify: `qrefine/refine.py` (`set_qm_defaults` 311-361 moves to api.py; `run` 406+ loses its param mutations and the `max_atoms` warning; asserts `params.finalized`)
- Create: `qrefine/tests/unit/tst_56.py` (test (e) now; extended in Task 4)
- Modify: `qrefine/run_tests.py`

**Interfaces:**
- Produces: `api.finalize_params(params, d_min=None, scattering_table=None, has_map=False, has_data=False, log=None) -> params` (same object, mutated, `params.finalized = True` injected via `params.__inject__("finalized", True)`); rules 1-6 of spec section 4 in that order; `Auto` target resolved by the caller: `api.resolve_target_name(params, fmodel) -> str` ("mlhl" when `fmodel.hl_coeffs() is not None` and the parameter is Auto, "ml" for Auto otherwise, the explicit value else).
- Consumes: Task 1 signatures.

- [ ] **Step 1: Write the failing tests** (tst_56)

```python
def test_finalize_auto_cust():
  for d_min, b, a in [(1.0, 0.025, 2.5), (2.0, 0.02, 2.0), (3.5, 0.01, 1.7)]:
    p = qr.get_default_params(); p.quantum.engine_name = "aimnet2"
    api.finalize_params(p, d_min=d_min, has_data=True, log=null_out())
    assert p.finalized
    assert (p.refine.max_bond_rmsd, p.refine.max_angle_rmsd) == (b, a)
    assert p.expansion and not p.cluster.clustering
    assert p.refine.minimizer == "lbfgs" and p.refine.gradient_only
    assert p.refine.number_of_weight_search_cycles == 10
  p = qr.get_default_params(); p.quantum.engine_name = "aimnet2"
  api.finalize_params(p, has_map=True, log=null_out())
  assert p.refine.minimizer == "lbfgsb" and not p.expansion and p.refine.max_bond_rmsd == 0.009
def test_finalize_refuses():
  p = qr.get_default_params(); p.expansion = True; p.cluster.clustering = True
  try: api.finalize_params(p, d_min=2.0, has_data=True)
  except Sorry: pass
  else: raise AssertionError
  p = qr.get_default_params(); p.refine.minimizer = "lbfgsb"; p.refine.gradient_only = True
  try: api.finalize_params(p, d_min=2.0, has_data=True)
  except Sorry: pass
  else: raise AssertionError
def test_finalize_modes():
  p = qr.get_default_params(); api.finalize_params(p)
  assert p.refine.mode == "opt"
  p = qr.get_default_params(); p.cluster.clustering = True; api.finalize_params(p, d_min=2, has_data=True)
  assert p.refine.gradient_only and p.refine.minimizer == "lbfgs"
  p = qr.get_default_params(); p.quantum.engine_name = "mopac"; api.finalize_params(p)
  assert p.quantum.method == "PM7" and p.quantum.basis == ""
def test_rst_file_and_max_atoms_gone():
  assert "rst_file" not in qr.master_phil_str and "max_atoms" not in qr.master_phil_str
  assert "max_atoms" not in open(refine.__file__).read()
def test_resolve_target_name():
  # fmodel-like stub with hl_coeffs() -> None / not None (use a small class);
  # Auto -> "ml"/"mlhl"; explicit "ls_wunit_k1" stays.
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/qrefine/tests/unit/tst_56.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'qrefine.api'`.

- [ ] **Step 3: Implement**

`api.py` module docstring: the four public functions and that qr.py is their first caller. Move the code bodies as named in Files; `finalize_params` prints the auto-cust table qr.py printed today (same text) when `log` is given. `refine.run` top becomes: `assert getattr(params, "finalized", False), "call api.finalize_params first"`. qr.py's `run` computes `d_min`, `has_data`, `has_map` and calls `finalize_params` right after the model is processed (where `auto_cust` was), then `resolve_target_name` and `set_target_name`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; then `tst_55.py`, `tst_53.py`, `tst_46.py`, `tst_36.py` from scratch folders. Expected: OK, OK, OK, OK, OK; tst_55's logged parameters still show the aimnet2-free defaults.

- [ ] **Step 5: Checkpoint**

Register tst_56. Do not commit.

---

### Task 3: `api.check_engine` and `api.check_model`; tst_58

**Files:**
- Modify: `qrefine/api.py`
- Modify: `qrefine/qr.py` (`run`: call both after `finalize_params`; the `conformers.check_completed` call lines ~468-470 moves into `check_model`)
- Modify: `qrefine/refine.py` (`check_qm_sees_crystal_neighbours` 363-386 called from `check_model` instead of `run`)
- Create: `qrefine/tests/unit/tst_58.py`
- Modify: `qrefine/run_tests.py`

**Interfaces:**
- Produces: `api.check_engine(params)` and `api.check_model(model, params, has_map=False)`, both returning None or raising `Sorry`; `api.ENGINES` (Global Constraints) with `needs_files` True for mopac, xtb, orca, terachem, turbomole, gaussian (False for the rest, server included) and these lookups: mopac `Mopac().get_command()`, xtb `GFNxTB().get_command()`, orca `Orca().get_command()` (`ORCA_COMMAND`/`Orca_COMMAND`), terachem `TeraChem().get_command()` (`TeraChem_COMMAND`/`TERACHEM_COMMAND`), turbomole `Turbomole().get_command()` (`TM_COMMAND`), gaussian `shutil.which("g16") or shutil.which("g09")` (the plugin keeps no lookup), aimnet2/aimnet2-old `import aimnet2calc`, torchani `import torchani`, pyscf `import pyscf`, fairchem `import fairchem`, server GET with 3 s timeout (`urllib.request.urlopen`); `elements` `{"H","C","N","O"}` for torchani, None for the others (aimnet2 species read at check time); `hint` the install hint per engine.
- Consumes: Task 2 `finalize_params` (check_* assert `params.finalized`).

- [ ] **Step 1: Write the failing tests** (tst_58; fixtures: the tst_173 ARG/SER fragment through `mmtbx.reduce2 qmready=true capping=true` as tst_53 does, giving `capped.pdb`; variants made by editing the hierarchy in the test and writing with `write_pdb_or_mmcif_file`)

```python
def test_check_engine():
  p = finalized(restraints="cctbx"); api.check_engine(p)            # no raise
  p = finalized(engine="mopac"); api.check_engine(p)                # build has MOPAC
  if shutil.which("xtb") is None and "XTBHOME" not in os.environ:
    p = finalized(engine="xtb"); assert_sorry(api.check_engine, p, "not available")
  p = finalized(engine="server"); p.quantum.server_url = "http://127.0.0.1:1"
  assert_sorry(api.check_engine, p, "server")
def test_check_model_capped_ok():      # capped.pdb, restraints=qm, mopac: no raise
def test_check_model_uncapped():       # cap removed (remove=capping_h): Sorry mentions "ARG" and "3" and "capping=true"
def test_check_model_missing_h():      # all H removed: Sorry mentions "reduce2 qmready=true"
def test_water_no_h():                 # water O without H added: Sorry mentions the water
def test_water_one_h():                # water with one H: Sorry
def test_check_model_elements():       # if aimnet2 importable: a ZN atom added -> Sorry naming "Zn"; else skip with a print
def test_altloc_completion():          # the uncompleted altloc fixture `tst_50.pdb_str` (import it; qr.refine refuses it today) with handling=conformers: Sorry mentions "Altloc completion"
def test_check_model_cctbx_lenient():  # uncapped model with restraints=cctbx: no raise
def test_crystal_neighbours():         # P1 fixture in a tiny cell with restraints=qm, expansion=False, clustering=False, has_map=False: Sorry "Symmetry-related"; has_map=True: no raise
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/qrefine/tests/unit/tst_58.py`
Expected: FAIL with `AttributeError: module 'qrefine.api' has no attribute 'check_engine'`.

- [ ] **Step 3: Implement**

`check_model` per spec section 6, in that order; the residue listing format `"<chain> <resname> <resseq><icode>"`; caps per residue via `model.selection("capping_h")` intersected with the residue's atoms; water: atom groups whose resname is a common water name with fewer than two H/D atoms. Elements: `set(e.strip().capitalize() for e in hierarchy.atoms().extract_element())`, D -> H.

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; then `tst_53.py`, `tst_55.py`, `tst_52.py` (altloc completion path through qr.py) from scratch folders. Expected: all OK.

- [ ] **Step 5: Checkpoint**

Register tst_58. Do not commit.

---

### Task 4: `api.refine`, `refine.run` returning a result, qr.py split; tst_56 extended

**Files:**
- Modify: `qrefine/refine.py` (`run`: returns `group_args`; real-space branch returns too; no `monitor.finalize`)
- Modify: `qrefine/api.py` (`refine`)
- Modify: `qrefine/results.py` (delete `finalize`, `write_final_pdb_files`)
- Modify: `qrefine/qr.py` (`prepare_inputs` module function; `Program.run` = prepare, finalize, checks, `api.refine`, write)
- Modify: `qrefine/tests/unit/tst_56.py`, `qrefine/tests/unit/__init__.py` (`run_qrefine_in` also finds `pdb/*_refined.cif`)

**Interfaces:**
- Produces: `api.refine(model, params, fmodel=None, map_data=None, working_folder="ase", intermediate_folder=None, bulk_solvent_params=None, log=None) -> group_args(sites_cart, sites_cart_nuclear, r_work, r_free, data_weight, restraints_weight_scale, b_rmsd, a_rmsd, n_fev, n_cycles, cc_mask, mode, model)`; `qr.prepare_inputs(data_manager, params, log) -> group_args(model, fmodel, map_data, d_min, prefix)`; `qr.write_output(result, params, prefix, log) -> file name`.
- Consumes: Tasks 1-3.

- [ ] **Step 1: Write the failing tests** (tst_56, fixture `in.pdb` from `tst_xh_mapper.pdb_str`, data from `run_fmodel`; a helper `objects_like_cli(args)` builds a `DataManager`, loads the files, extracts params with `qr.master_phil_str` from the same argument strings as the command line, and calls `qr.prepare_inputs`)

```python
def test_api_equals_cli_cctbx():
  cli = run_qrefine_in(args + ["restraints=cctbx", "expansion=False", "clustering=False"], prefix, folder="cli")
  o = objects_like_cli(same args); api.finalize_params(...); api.check_engine; api.check_model
  r = api.refine(o.model, o.params, fmodel=o.fmodel, working_folder="api_ase", log=log)
  assert approx_equal(r.r_work, cli_r_work_from_log, 1e-6)
  assert r.sites_cart.rms_difference(cli.model.get_sites_cart()) < 1e-6   # same convention: X-ray lengths
def test_api_equals_cli_mopac():    # same with restraints=qm quantum.engine_name=mopac clustering=False, tolerance 1e-6
def test_inputs_untouched():        # fmodel built on o.model.get_xray_structure() (shared object): after api.refine,
                                    # o.model sites, o.model hierarchy as_pdb_string and o.fmodel.r_work() unchanged
def test_no_writes_outside():       # listing of the run folder before/after: only working_folder appears; with
                                    # intermediate_folder=None no cycle files anywhere
def test_intermediate_folder():     # intermediate_folder="cycles" (absent before): created, contains *_refine_cycle.pdb
def test_called_twice():            # two api.refine calls, same working_folder: identical r_work and sites
def test_no_files_at_all():         # working_folder=None, restraints=cctbx, expansion=False, clustering=False,
                                    # intermediate_folder=None: recursive listing of the run folder and the set of
                                    # entries of tempfile.gettempdir() identical before and after
def test_no_files_at_all_clustering():  # same with clustering=True (cctbx fragments, save_clusters left at its
                                    # default): still nothing written
def test_mopac_tempdir_removed():   # restraints=qm mopac clustering=False working_folder=None: succeeds; the log
                                    # names the temporary folder ("scratch folder: <path> (removed on return)");
                                    # that path does not exist afterwards
def test_opt_mode():                # fmodel=None: mode == "opt", r_work is None, sites_cart == sites_cart_nuclear
def test_raw_params_refused():      # api.refine with qr.get_default_params() unfinalized -> AssertionError
def test_auto_target_mlhl():        # fmodel with abcd (HL from the model's phases) -> resolve_target_name "mlhl"; then
                                    # api.refine runs with it (target set on the fmodel by the test as qr.py does)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/qrefine/tests/unit/tst_56.py`
Expected: FAIL with `AttributeError: module 'qrefine.api' has no attribute 'refine'`.

- [ ] **Step 3: Implement**

`api.refine`: `log = log or null_out()`; assert finalized; resolve the folder: `working_folder` as given, or None, or (`restraints=qm` with an engine whose `needs_files` is True) `tempfile.mkdtemp(prefix="qrefine_")`, logged and removed with `shutil.rmtree` in a `finally` (the one permitted recursive delete: it only removes the folder this call created); thread `working_folder=None` down so every write in restraints.py, fragment.py, cluster_restraints.py, super_cell.py, conformers.py, driver.py, calculator.py is skipped (`if working_folder is not None:`); build the model copy as spec section 3 says (`mmtbx.model.manager(model_input=None, pdb_hierarchy=model.get_hierarchy().deep_copy(), crystal_symmetry=model.crystal_symmetry(), restraint_objects=model.get_restraint_objects(), log=null_out())`, process as qr.py, `setup_scattering_dictionaries(scattering_table=model.get_scattering_table(), d_min=...)`), `fmodel = fmodel.deep_copy()` with `bulk_solvent_params` stored for the driver's calls (add the argument to `driver.refine` and pass it to both `update_all_scales` calls), `refine.run(...)`, then the result: `sites_cart` from `mapper.shrink` (or nuclear), `model` the API copy. qr.py: `prepare_inputs` from today's `run` body; `write_output` through `model_as_pdb_or_mmcif(target_format="pdb")` with the extension chosen by `can_be_output_as_pdb()`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; then `tst_55.py`, `tst_53.py`, `tst_00.py`, `tst_36.py`, `tst_46.py` from scratch folders. Expected: all OK, same numbers as before Task 1 (compare the R factors in tst_55's log with the Task 1 run).

- [ ] **Step 5: Checkpoint**

Do not commit.

---

### Task 5: mmCIF identifiers audit; tst_57

**Files:**
- Modify: `qrefine/restraints.py` (`from_qm.target_and_gradients` write at line ~354 and the `file_name` label derivation), `qrefine/fragment.py` (writes at lines ~354-449), `qrefine/cluster_restraints.py:30`, `qrefine/super_cell.py:71`, `qrefine/charges.py` (`write_pdb_hierarchy_*` if reached from fragments), `qrefine/results.py` (`write_pdb_file`)
- Create: `qrefine/tests/unit/tst_57.py`
- Modify: `qrefine/run_tests.py`

**Interfaces:**
- Produces: every write listed by the Task 1 grep goes through `write_pdb_or_mmcif_file(target_filename=<path>.pdb, target_format="pdb")` and keeps the returned name; `from_qm` derives the engine label from the returned name's stem.

- [ ] **Step 1: Write the failing test** (tst_57)

Fixture: `tst_xh_mapper.pdb_str` with chain `A` renamed `AA` and the water chain `S` renamed `SW`, plus one ligand residue `ETHA` (an ethanol: copy `chem_data/mon_lib/e/EOH.cif`, replace every `EOH` by `ETHA` in a temporary file in the test folder, place its atoms 6 A from the SER OG) written as mmCIF with `write_pdb_or_mmcif_file`; data from `run_fmodel` on that file (mmCIF input is fine for `mmtbx.command_line.fmodel`).

```python
def test_cctbx():   # qr.refine model.cif data.mtz etha.cif restraints=cctbx expansion=False clustering=False:
                    # output pdb/*_refined.cif exists, chain ids {"AA","SW"} and resname "ETHA" present, atom count equal
def test_mopac_no_clustering():   # restraints=qm quantum.engine_name=mopac clustering=False expansion=True: same checks;
                                  # ase/ contains a .cif scratch model and the MOPAC .mop input has the right atom count
def test_mopac_clustering():      # clustering=True: same checks; fragment files under ase/ are .cif
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.python /Users/pafonine/Desktop/all/phenix/modules/qrefine/tests/unit/tst_57.py`
Expected: FAIL (empty scratch PDB files or a crash in the engine path, or a `.pdb` output with truncated identifiers).

- [ ] **Step 3: Implement**

Replace each `write_pdb_file(file_name=X)` with `name = h.write_pdb_or_mmcif_file(target_filename=X, target_format="pdb")` and use `name` afterwards; in `from_qm`, `label = os.path.splitext(name)[0]` from the returned name; `results.manager.write_pdb_file` likewise (file name without extension in, returned name logged).

- [ ] **Step 4: Run the tests to verify they pass**

Run: the Step 2 command; then `tst_53.py`, `tst_36.py`, `tst_55.py` from scratch folders. Expected: all OK; the `pdb/` names of PDB-format runs unchanged (`*_refined.pdb`).

- [ ] **Step 5: Checkpoint**

Register tst_57. Do not commit.

---

### Task 6: Whole-suite run

**Files:** none expected (report only; a failing test is a bug to fix in the task that caused it).

- [ ] **Step 1: Run the qrefine suite**

Run: `cd /Users/pafonine/Desktop/all/phenix/rt/run_$(date +%s) && libtbx.run_tests_parallel module=qrefine nproc=4 run_in_unique_dirs=True > suite.log 2>&1; tail -20 suite.log`
Expected: every test OK, including tst_56-59. A number that moved is a bug (nothing in this sub-project may change command-line results): stop and report it.

- [ ] **Step 2: Report**

Files changed in qrefine; the suite summary; any pre-existing oddity met (not fixed).
