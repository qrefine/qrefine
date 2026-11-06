# phenix.refine hook for qr.refine: design

Date: 2026-11-05. Status: approved 2026-11-05 (D22-D25 confirmed). Sub-project 3 of
`2026-11-04-qrefine-in-phenix-refine-design.md` (umbrella D1, D7-D13).
Repositories: phenix (phenix/refinement), phenix_regression, phenix_html.
qrefine is consumed through `qrefine.api` (sub-project 2, committed
a7177ba); no qrefine change is expected.

## 1. Goal

`phenix.refine ... main.use_qrefine=True` refines coordinates in the last
macro-cycle with qr.refine instead of the cctbx coordinate minimization and
does everything else as today. A user gets QM-restrained coordinates plus
phenix.refine's ADP, occupancy and anomalous refinement, its output files
and its statistics. Unsupported setups stop at startup with a Sorry; nothing
changes for runs with `use_qrefine=False`.

## 2. Findings

- `phenix.programs.phenix_refine.Program.validate` (phenix_refine.py ~360)
  gathers the models and builds the fmodel(s) from the data manager;
  `self.joint_xn` marks joint X-ray/neutron mode (two models, two data
  sets). `Program.run` then calls `driver.run0` (single) or `run0_joint`.
- `driver.run0` -> `refinement_objects.__init__` -> `fsr.preprocess.run(
  params_global, master_params, model, fmodels, call_back_handler, log)`,
  "order of calls is important": `check_params` (IAS, ions), `tls1`,
  `finalize_model` (pdb_interpretation with CDL, neutron distances for
  neutron/electron tables, `model.process(make_restraints=True)`,
  scattering dictionaries at the data d_min), `tidy_us`, ...,
  `hydrogens` (resolves `hydrogens.refine=Auto`: X-ray -> riding unless
  d_min < 0.9 or ias; neutron -> riding if d_min > 1.5), 
  `set_refinement_flags` (`model.refinement_flags.sites_individual` from
  `refine.sites.individual`), ..., `setup_ncs`, monitors. It returns the
  monitors. This is where the model is complete and processed and where
  startup checks that need it belong.
- `driver.run.run()` loops `macro_cycle.run(...)` for
  `main.number_of_macro_cycles` cycles; afterwards `clean_ions`,
  `clean_water`, outputs. No H idealization after the loop (grep).
- `macro_cycle.run.__init__` executes the step list in a fixed order
  through `call(func, prefix)`, each gated by `task_table.status(func)`;
  `task_table.update(params, macro_cycle, model, fmodel,
  neutron_refinement)` recomputes every `run_<step>` flag per macro-cycle.
  Relevant flags: `run_setrh` (riding only, mode "last"), `run_regHxyz`
  (riding), `run_ncsrot`, `run_realsrg`, `run_realsrl`, `run_realsrl2`,
  `run_sacart`, `run_tardy`, `run_den`, `run_xyzrec`, `run_fitrh`
  (neutron or individual H only), `run_regHadp` (riding ADP of H).
- Riding H: no manager exists in macro-cycle 1 (H refine freely there);
  `regHxyz` ("Ensure riding H", `model.idealize_h_riding()`) creates the
  manager at the end of cycle 1, so cycles 2.. ride; `setrh` re-creates it
  in the last cycle. With `use_qrefine` both are off in the last cycle, so
  the riding manager of cycles 2..n-1 is never applied to QM positions.
- `macro_cycle.xyzrec` (macro_cycle.py 429) is the only caller of
  `xyz_reciprocal_space.run`; `call` asserts model/fmodel consistency
  before and after every step (`assert_or_set_expected_state`:
  `assert_xray_structures_equal`, grad flags reset, `model.set_xray_structure`).
- phenix.refine's bulk-solvent step passes `params.bulk_solvent_and_scale`
  (mmtbx `bulk_solvent_and_scaling.master_params`), `remove_outliers` from
  the data manager, `update_f_part1`, `apply_back_trace` (only with ADP
  strategies), `optimize_mask`, `refine_hd_scattering` (neutron with riding
  H and `hydrogens.optimize_scattering_contribution`).
- Twinning: `fmodels.fmodel_xray().twin` is truthy for a twin fmodel. The
  `qi` scope exists only when injected (`hasattr(params, "qi")`);
  `mmtbx.geometry_restraints.quantum_interface.is_quantum_interface_active`.
  amber: `hasattr(params, "amber") and params.amber.use_amber`.
- `main.ordered_solvent` defaults to False (the umbrella assumed True);
  `main.place_ions` defaults to None; `main.ias` False; `den` is a strategy
  choice; `main.fix_cablam` defaults to off; `main.nqh_flips` True.
- `hydrogens.high_resolution_limit_to_include_scattering_from_h` is
  defined but unused; H always scatter.
- qrefine API facts (sub-project 2 final review notes): `api.check_engine`
  works on unfinalized params; `finalize_params` once per extract
  (second call Sorry); `api.refine` uses the fmodel's target as set, takes
  `bulk_solvent_params`, re-scales three times with it; with
  `working_folder=None` writes nothing (aquaref/cctbx) or a tempdir
  removed on return (external engines); the model copy keeps hierarchy,
  symmetry, restraint objects, SS annotation, scattering table, exact
  sites, and loses pdb_interpretation customizations (CDL, edits, links
  by LINK records, reference restraints), the riding manager and ias;
  `refine.exclude` is a selection string; the model's H must be in the
  convention of its scattering table; `result.sites_cart` is the shrunk
  (caller-convention) view; `check_model` requires finalized params and a
  processed model.
- phenix_regression: tests live in `phenix_regression/refinement/tst_N.py`
  (next free number 175) and are discovered by directory; helpers
  `run_phenix_refine(args, prefix, sorry_expected=...)` (attributes `.pdb`,
  `.log`, `.geo`, `.r_work_start`, `.r_free_start`), `run_fmodel`,
  `from_pdb_or_log`. tst_173 (capped model through phenix.refine) is the
  closest template. qrefine's `tests/unit` has `run_qrefine_in` and the
  tst_xh_mapper fixture.
- Docs: `phenix_html/rst_files/reference/refinement.pre`, section
  "Refinement of coordinates" (line 313), and `AQuaRef.txt`.

## 3. User-facing parameters

In `phenix/refinement/__init__.params`:

```
main {
  use_qrefine = False
    .type = bool
    .help = Refine coordinates in the last macro-cycle with qr.refine (QM \
            or cctbx restraints through the qrefine module) instead of the \
            standard coordinate minimization; everything else is unchanged. \
            Requires a model prepared with mmtbx.reduce2 qmready=true \
            (capping=true for truncated residues). See the qrefine scope.
    .short_caption = Coordinates with qr.refine (QM restraints)
    .expert_level = 1
}
qrefine
  .help = Settings of the qr.refine coordinate step (main.use_qrefine=True)
  .short_caption = qr.refine coordinate refinement
  .expert_level = 1
{
  restraints = cctbx *qm
    .type = choice(multi=False)
    .help = qm: restraints from the QM engine; cctbx: standard restraints \
            through qrefine (fast; for testing)
  qm_engine = *aquaref mopac xtb orca pyscf torchani terachem turbomole \
              gaussian server
    .type = choice(multi=False)
    .help = QM engine; aquaref is the AIMNet2 model of AQuaRef
  clustering = None
    .type = bool
    .help = Fragment the model (qrefine cluster.clustering); None: the \
            engine's default (off for aquaref)
  max_bond_rmsd = None
    .type = float
    .help = Bond rmsd accepted during qr.refine's weight search; None: \
            qrefine's default (resolution-dependent for aquaref)
}
```
The scope lives at the level of `bulk_solvent_and_scale`. Nothing from
qrefine's PHIL is included: phenix.refine starts without qrefine installed
as long as `use_qrefine` is False. GUI exposure is not designed here
(umbrella deferral); the expert level keeps it out of the main dialog.

## 4. Translation to qrefine parameters

`phenix/refinement/xyz_qrefine.py`, `make_params(params, model, fmodel)`
-> finalized qrefine extract:
- `qp = qrefine.qr.get_default_params()`; `qp.restraints =
  params.qrefine.restraints`; `qp.auto_cust = (qp.restraints == "qm")`
  (the AQuaRef preset is a QM preset: with cctbx restraints qrefine's
  plain defaults apply, lbfgsb and no expansion, and `clustering` defaults
  to False unless the user sets `qrefine.clustering`); `qp.quantum.engine_name = "aimnet2" if
  params.qrefine.qm_engine == "aquaref" else params.qrefine.qm_engine`;
  `qp.cluster.clustering = params.qrefine.clustering` when not None;
  then, AFTER `finalize_params` (whose aquaref preset rewrites
  `restraints`, `cluster.clustering`, `expansion` and `max_bond_rmsd`),
  the user's explicit choices are re-applied: `qp.restraints` always,
  `qp.cluster.clustering` when not None (with `qp.expansion = not
  clustering`), `qp.refine.max_bond_rmsd` when not None; the preset
  supplies only what the user did not choose. `qp.write_nuclear_h_positions
  = False`; `qp.refine.exclude` built
  from `params.refine.sites.individual` when it does not select every atom:
  `"not (%s)" % " or ".join("(%s)" % s for s in strings)`.
- `api.finalize_params(qp, n_conformers=conformers.n_conformers(
  model.get_hierarchy()), d_min=fmodel.f_obs().d_min(),
  scattering_table=params.main.scattering_table, has_data=True,
  has_map=False, log=log)`. The `scattering_table` argument is informative;
  qrefine follows the model's table (sub-project 2).
- `qp.parallel.nproc` and everything else stay at qrefine's defaults (D11).
`make_params` is pure parameter logic and is called twice: at startup for
`check_model`, and in the last macro-cycle for `refine`.

## 5. Startup checks

Two places, both before any macro-cycle:

1. `Program.validate` (phenix_refine.py), first thing when
   `params.refinement.main.use_qrefine` (or the joint `xray`/`neutron`
   copies) is True: `import qrefine` else `Sorry("main.use_qrefine=True
   needs the qrefine module, which is not installed in this Phenix")`;
   joint X-ray/neutron mode -> `Sorry("main.use_qrefine is not available for
   joint X-ray/neutron refinement")`; `api.check_engine(qp)` with a raw
   extract carrying `restraints`/`engine_name`/`server_url` (works
   unfinalized). This runs before the data manager builds fmodels.
2. `fsr.preprocess.qrefine_checks(params, model, fmodels, log)`, called
   from `preprocess.run` right after `set_refinement_flags` (model
   processed, H mode resolved, selections known) and before `setup_ncs`:
   - Sorry: `fmodels.fmodel_xray().twin` ("twinned data"); `params.main.ias`;
     `"den" in params.refine.strategy`; `hasattr(params, "qi") and
     is_quantum_interface_active(params)`; `hasattr(params, "amber") and
     params.amber.use_amber`; any `refine.sites.individual` selection that
     does not cover every atom (`xyz_qrefine._exclude_selection(params,
     model)` not None: "refine.sites.individual selections are not
     supported yet with main.use_qrefine=True ..."), because qrefine's
     `refine.exclude` path diverges (it restrains only the selected atoms
     while the weight search targets the whole model's bond rmsd); the
     refusal is lifted when that is fixed. Each message names the
     parameter to change.
   - Auto-off with a boxed log message ("Automatic adjustment:
     use_qrefine=True disables ordered_solvent and place_ions: waters added
     during refinement would reach the QM engine without hydrogens"):
     `params.main.ordered_solvent = False`, `params.main.place_ions =
     None`. (D8; a Sorry would be wrong only for users who turned
     ordered_solvent on explicitly, and the message tells them.)
   - `qp = make_params(...)`; `api.check_model(model, qp, has_map=False)`
     on phenix.refine's processed model (completeness, caps, waters,
     elements, altloc completion, crystal neighbours when neither expansion
     nor clustering).
   - With `restraints=qm` and `hydrogens.refine == "individual"` nothing
     special; with riding nothing special (D5: the last cycle is free of
     riding steps).

## 6. The macro-cycle

`task_table.update` sets `self.qrefine_cycle = params.main.use_qrefine
and macro_cycle == params.main.number_of_macro_cycles`. In that cycle
these flags are forced False: `run_setrh`, `run_regHxyz`, `run_ncsrot`,
`run_realsrg`, `run_realsrl`, `run_realsrl2`, `run_sacart`, `run_tardy`,
`run_den`. `run_xyzrec` keeps its rule (the step runs; its body branches).
Every other step keeps its rule, so bulk solvent, target, restraint
updates, NQH flips, cablam, rigid body, NCS operators, weights (still
needed for the ADP weight), the coordinate step, ADP, riding H ADP
regularization, H/D recentering, occupancies, fp/fdp run as today. The
step header reads "XYZ individual (reciprocal space) [qr.refine]" in that
cycle (the task table has no printout of its own).

`macro_cycle.xyzrec`: `if self.tt.qrefine_cycle:
xyz_qrefine.run(model, fmodels, params, prefix, log); return` before the
existing code.

## 7. `xyz_qrefine.run(model, fmodels, params, prefix, log)`

1. `print_statistics.make_header(prefix + " [qr.refine]", out=log)`.
2. `fmodel = fmodels.fmodel_xray()`; `assert_xray_structures_equal(fmodel
   .xray_structure, model.get_xray_structure())` (as every step does).
3. `qp = make_params(params, model, fmodel)`; `api.check_engine(qp)`.
4. `result = api.refine(model, qp, fmodel=fmodel, working_folder=None,
   intermediate_folder=None, bulk_solvent_params=
   params.bulk_solvent_and_scale, log=log)`. No file is written by
   qrefine (D21; aquaref and cctbx write nothing, an external engine uses
   a tempdir removed on return). The model goes in with phenix.refine's H
   convention and comes back in it (`result.sites_cart`).
5. Apply: `xrs = fmodel.xray_structure; xrs.set_sites_cart(
   result.sites_cart); fmodels.update_xray_structure(xray_structure=xrs,
   update_f_calc=True, update_f_mask=True, force_update_f_mask=True);
   model.set_xray_structure(xrs)`.
6. Re-scale as qrefine did: `fmodel.update_all_scales(params=
   params.bulk_solvent_and_scale, remove_outliers=False,
   refine_hd_scattering=False)`; `fmodels.create_target_functors()`.
7. D13 assert: `approx_equal(fmodel.r_work(), result.r_work, 1e-4)` and the
   same for r_free. If the regression data of section 10 show a persistent
   gap, the documented fallback is for the hook to adopt the API's final
   scales (sub-project 2 can expose them); the assert stays.
8. Summary, in the style of `xyz_refinement_scorer._show`: one table line
   with r_work, r_free, delta, bond and angle rmsd (from phenix.refine's own
   `model.geometry_statistics()` after the step, X-ray geometry, so it is
   comparable with the other cycles), and a second line with qrefine's
   data weight, restraints weight scale, `n_fev`, `n_cycles`.
9. The riding manager, if one exists from earlier cycles, is left as is;
   nothing idealizes H afterwards in this cycle (section 6), and no
   standard cycle follows (D1).

Caps stay in the model (D12); when `model.selection("capping_h")` is non
empty the log says, once, how to strip them before deposition
(`phenix.pdbtools model.pdb remove=capping_h`).

## 8. Neutron and electron data

`scattering_table=neutron`: the model's H/D are nuclear; the API's X-H
mapper is the identity; exchangeable H/D are handled by qrefine's
hd_mapper; `hydrogens.refine` is whatever phenix.refine chose, and the
riding steps are off in the qrefine cycle like for X-ray. `fitrh` (fit
rotatable H) keeps its rule (neutron or individual H) and runs before the
coordinate step as today. `electron`: nuclear convention, identity mapper;
phenix.refine's own charge handling is kept (the API does not neutralize;
`setup_scattering_dictionaries` with the electron table neutralizes the
API's copy as it does phenix.refine's model).

## 9. Documentation

`refinement.pre`, new subsection under "Refinement of coordinates":
"Coordinate refinement with QM restraints (qr.refine, AQuaRef)": what it
does, the last-macro-cycle rule, model preparation with reduce2, the
parameters, the refusals, the caps note, engine availability, and an
example command line. `AQuaRef.txt` gains a paragraph pointing to it.

## 10. Tests (phenix_regression/refinement)

Fixtures: the tst_xh_mapper fragment (protein, water, altloc SER; nuclear
H from reduce2) written as `in.pdb`; data from `run_fmodel` at 1.5 A with
`n_gaussian`, and separately with `scattering_table=electron` and
`neutron`; the capped tst_173 fragment for the caps case. MOPAC presence
checked with `qrefine.api.check_engine` (skip with a printed note).

- tst_175: `use_qrefine=True qrefine.restraints=cctbx` with three
  macro-cycles, `main.nqh_flips=False`: the log shows "qr.refine" in the
  coordinate step of cycle 3 only; "Ensure riding H", "Local real-space
  refinement" lines present in cycles 1-2 and absent in cycle 3; the
  output model's X-H lengths equal the X-ray targets within 1e-3 A (model
  processed with X-ray distances); R-work decreased from the start; the
  D13 assert passed (no AssertionError; the log shows the qrefine summary).
- tst_176: the same with `qrefine.restraints=qm qrefine.qm_engine=mopac
  qrefine.clustering=False` and one macro-cycle; skipped without MOPAC.
- tst_177: one macro-cycle, `strategy=individual_sites`, `use_qrefine=True`
  with cctbx restraints, against X-ray, electron and neutron data (three
  runs), each compared with standalone `qr.refine` on the same files with
  the same settings (`restraints=cctbx expansion=False clustering=False`,
  same cycle counts): refined coordinates equal within 1e-4 A, R factors
  within 1e-5 (umbrella). Settings pinned on the phenix.refine side:
  `main.nqh_flips=False`, `main.fix_cablam=off`, `main.bulk_solvent_and_scale`
  default, `main.target=ml`, and for neutron
  `hydrogens.optimize_scattering_contribution=False` and
  `hydrogens.real_space_optimize_x_h_orientation=False` (the two steps
  qr.refine does not have). If a residual difference remains, the test
  states which setting causes it and the gate is loosened only with that
  reason written in the test.
- tst_178: the refusals and adjustments: twin_law set; `main.ias=True`;
  `den` strategy; `qrefine.qm_engine=xtb` when xtb is absent (else a fake
  `server`); a `refine.sites.individual` selection of part of the model
  (one residue with clustering off, a chain with clustering on); an
  uncapped truncated model with `restraints=qm`; each a Sorry with the
  expected words. Plus: `main.ordered_solvent=True` and `main.place_ions`
  auto-off with the boxed message; the capped model refining with the
  caps kept at occupancy 0 and the strip note in the log;
  `number_of_macro_cycles=1`.
- Existing suite: tst_173, tst_174 and a sample of refinement tests
  (tst_000-tst_010 or the ones the implementer can run in reasonable time)
  unchanged with `use_qrefine` off; the task table and the PHIL are the
  only touched shared code.

## 11. Implementation outline

1. PHIL (section 3) and `xyz_qrefine.make_params`; Program.validate check.
2. `preprocess.qrefine_checks`; tst_178 refusals.
3. `task_table.qrefine_cycle` and the gated flags; `xyz_qrefine.run` and
   the `xyzrec` branch; tst_175.
4. tst_176 (MOPAC), tst_177 (equivalence, three tables), tst_178 rest.
5. Documentation; sample of the existing suite.

## 12. Decisions to confirm

- D22. `ordered_solvent`/`place_ions` auto-off with a boxed message rather
  than a Sorry (ordered_solvent defaults to False, so only users who asked
  for it see the message).
- D23. `make_params` called twice (startup and last cycle) rather than
  stored; stateless and cheap.
- D24. The summary line reports phenix.refine's X-ray geometry statistics
  (comparable with the other cycles), with qrefine's own numbers on a
  second line.
- D25. `qrefine.clustering` and `max_bond_rmsd` are `None`-able PHIL values
  meaning "qrefine's default" (the GUI-friendly form of Auto).
