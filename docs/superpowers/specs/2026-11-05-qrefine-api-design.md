# qrefine callable API: design

Date: 2026-11-05. Status: approved 2026-11-05 (D16-D20 confirmed, D21 proposed); implemented 2026-11-05. Sub-project 2 of
`2026-11-04-qrefine-in-phenix-refine-design.md` (umbrella D3, D4, D7, D10,
D13, D14, D15). Repository: qrefine (plus one qrefine test fixture from
chem_data). No cctbx change is expected.

## 1. Goal

qr.refine becomes a function that takes objects in and returns objects out,
so phenix.refine (sub-project 3) can call it for the coordinate step of its
last macro-cycle. The command line `qr.refine` keeps its behaviour and
output layout (`pdb/` for models, `ase/` for engine scratch) by calling the
same function. Checks that today fail late (engine missing, model
incomplete) become functions a caller can run before any heavy work.

## 2. Findings

- `qr.Program.run` (qr.py ~387-489) does everything in one method: fmodel
  from the data manager plus `update_all_scales(refine_hd_scattering=False)`
  and `set_target_name`; map normalisation; model reconciliation of
  altlocs; `conformers.check_completed`; `neutralize_scatterers` for the
  electron table; `process` with neutron distances, no CDL, no atom
  sorting, no GRM normalisation; scattering dictionaries at the data
  resolution; `auto_cust`; `refine.run`. `validate` holds the
  expansion/clustering exclusivity and minimizer rules; `auto_cust` the
  aimnet2 presets by resolution; `refine.set_qm_defaults` the engine
  method/basis defaults; `refine.run` mutates params again (conformers and
  clustering force `minimizer=lbfgs`, `gradient_only=True`).
- `refine.run(model, fmodel, map_data, params, rst_file, prefix, log)`
  returns nothing. It creates `params.output_folder_name` and the restart
  pickle path (`rst_file`, removed at start, never written since the
  restart feature was dropped), builds everything, and ends with
  `monitor.finalize` writing `<prefix>_refined.pdb` into the output folder.
  Real-space mode writes `real_space_refined.pdb` into the current
  directory and returns. `driver.refine` writes `<n>_refine_cycle.pdb` per
  cycle, `driver.opt` `<n>_opt_cycle.pdb`, both through the monitor.
- `driver.py` binds `log = sys.stdout` at import and uses it in 18 prints;
  `refine.py`, `results.py`, `calculator.py` take a log. `driver.opt` with
  `calculator.debug` pickles to the current directory.
- `calculator.sites.target_and_gradients` with `dump_gradients` pickles
  three arrays and calls `sys.exit(0)`; `tests/unit/tst_36.py` relies on
  that exit.
- Engine scratch: `restraints.restraints(params, model, working_folder="ase")`
  and `from_qm(file_name="./ase/tmp_ase.pdb")`; `from_qm` writes
  `tmp_ase.pdb` with `pdb_hierarchy.write_pdb_file` on every evaluation
  when not clustering, then builds ASE atoms from the hierarchy in memory;
  the engines get the atoms and a label derived from the file name, so the
  PDB file is a record, not an input. `fragment.fragments` writes
  `expansion.pdb` into the current directory and per-fragment files into
  its working folder (`save_clusters`); `cluster_restraints` writes the
  super hierarchy; `super_cell.write_super_cell_selected_in_sphere` writes
  a PDB. `conformers.from_conformers` builds one working folder per
  conformer under the given one.
- `iotbx.pdb.hierarchy.as_pdb_string` returns an empty string, without
  error, when the hierarchy does not fit the PDB format (chain IDs longer
  than two characters, residue names longer than three, ...);
  `write_pdb_file` is deprecated in favour of `write_pdb_or_mmcif_file`,
  which returns the file name it wrote. `mmtbx.model.manager.model_as_pdb`
  has the same limitation; `model_as_pdb_or_mmcif` and
  `can_be_output_as_pdb` exist. `results.manager` uses `model_as_pdb`.
- Engine availability today is discovered by the plugins when first used:
  `Mopac.get_command()` (mopac_manager, `which("mopac")`, `MOPAC_COMMAND`),
  `GFNxTB.get_command()` (`XTBHOME`, `which("xtb")`), `aimnet2_qr` imports
  `aimnet2calc` at module import, `torchani_qr`/`pyscf_qr`/`fairchem`
  import their packages; orca, terachem, turbomole, gaussian, server have
  no check before the first run.
- `aimnet2calc.AIMNet2ASE` exposes `implemented_species` (atomic numbers)
  when the loaded model provides it and refuses atoms outside it at
  calculation time.
- Model completeness: `mmtbx.model.manager.missing_atoms_by_residue(
  include_hydrogens)` (residue id -> missing atom names, heavy or heavy+H,
  water H never required, HIS tautomers and disulfide HG handled);
  `model.selection("capping_h")` and `capping.capping_h_selection`
  identify caps; `conformers.check_completed` raises for altloc
  completion.
- `fmodel.update_all_scales(params=...)` passes the
  `mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params` extract on to
  `f_model_all_scales` (bulk_solvent and anisotropic_scaling switches
  among others); `None` means defaults. phenix.refine passes its
  `bulk_solvent_and_scale` scope.
- phenix.refine's target rule (`phenix.refinement.misc.set_data_target_type`):
  `ls` when the test set is too small, `mlhl` when experimental phases are
  present and allowed, else `ml`. qr.py offers `*ml ls_wunit_k1` and
  ignores phases.
- `tests/unit/__init__.py`: `run_fmodel`, `run_qrefine_in` (expects
  `pdb/*_refined.pdb`). No qrefine test uses an mmCIF model.

## 3. The API: `qrefine/api.py`

```
refine(model, params, fmodel=None, map_data=None, working_folder="ase",
       intermediate_folder=None, bulk_solvent_params=None, log=None)
```
- `model`: an `mmtbx.model.manager` with its restraint objects (ligand
  CIFs) attached, scattering table set, H positions in the convention of
  its scattering table (sub-project 1). It is never modified: the API
  builds its own model from a deep copy of the hierarchy with the same
  crystal symmetry and restraint objects, processed as qr.py does today
  (neutron distances, CDL off, sort_atoms False, grm_normalization False)
  and with scattering dictionaries at the data resolution.
- `fmodel`: an `mmtbx.f_model.manager` whose scatterers are in the model's
  atom order (checked, `Sorry` otherwise), with its target set, scaled by
  the caller. Deep-copied; the caller's fmodel is never modified. `None`
  with `map_data=None` means geometry optimisation (`refine.mode=opt`).
- `map_data`: a normalised real-space map (as qr.py prepares it) for
  real-space mode; unchanged behaviour, included so the command line can
  use the same function.
- `params`: qrefine's PHIL extract after `finalize_params` (section 4).
  The API asserts `params.finalized` (a flag `finalize_params` sets) so a
  raw `get_default_params()` is refused.
- `working_folder`: engine scratch, fragments, expansion file, conformer
  sub-folders. Created if needed. Nothing is written anywhere else.
  `None` (what phenix.refine passes) means qrefine writes no file at all:
  no scratch model, no expansion or fragment files, no charge files, no
  debug output; the engines that run in this process (aimnet2,
  aimnet2-old, torchani, fairchem, pyscf) and cctbx restraints then touch
  the disk nowhere. An engine that is an external program (mopac, xtb,
  orca, terachem, turbomole, gaussian) cannot run without its input and
  output files: with `working_folder=None` it gets a private temporary
  folder (`tempfile.mkdtemp`) that is removed when `refine` returns, so
  nothing persists (D21). `api.ENGINES[name].needs_files` records which
  engines are which.
- `intermediate_folder`: per-cycle model files (`<n>_refine_cycle`,
  `<n>_opt_cycle`) when not None.
- `bulk_solvent_params`: passed to every `update_all_scales` call the API
  makes on its fmodel copy (inside the refine cycles, after re-scaling at
  the end); `None` means cctbx defaults. The other arguments of those
  calls stay fixed: `remove_outliers=False`, `refine_hd_scattering=False`,
  no back-trace.
- `log`: every print goes there; `None` means `null_out()`.

Returns a `libtbx.group_args`:
`sites_cart` (caller's convention), `sites_cart_nuclear`, `r_work`,
`r_free` (None without data), `data_weight`, `restraints_weight_scale`,
`b_rmsd`, `a_rmsd` (nuclear geometry), `n_fev`, `n_cycles`, `cc_mask`
(real-space mode only, else None), `mode` ("refine", "opt" or
"real_space"), `model` (the API's own refined model, nuclear sites, for
callers that want the hierarchy; the caller must not assume its atom
order differs from the input: it does not).

Also in `api.py`:
```
finalize_params(params, d_min=None, scattering_table=None, has_map=False,
                has_data=False, log=None)
check_engine(params)
check_model(model, params, has_map=False)
```
described in sections 4 to 6.

`refine.run` becomes the implementation behind `api.refine`: same
signature as `api.refine` minus the copying, returns the result; qr.py and
the API never call it directly except through `api.refine`.

## 4. `finalize_params`

One place for every parameter rule, called once before `check_*` and
`refine`:
1. From `qr.Program.validate`: `expansion` and `cluster.clustering` are
   exclusive (`Sorry`); `minimizer=lbfgsb` refuses `gradient_only`;
   expansion or clustering force `minimizer=lbfgs`, `gradient_only=True`.
2. From `qr.Program.auto_cust` (unchanged logic, now a function of
   `d_min`, `has_data`, `has_map`): the aimnet2 presets.
3. From `refine.set_qm_defaults`: engine method and basis defaults.
4. From `refine.run`: mode forced to `opt` without data and map, `refine`
   with a map; `altlocs.handling=native` with QM restraints is a `Sorry`;
   alternative conformations with `conformers` handling and clustering
   force `minimizer=lbfgs`, `gradient_only=True`.
5. `refinement_target_name` gains `Auto` as default: `ml`, or `mlhl` when
   the fmodel carries experimental phases (D15b); explicit choices win.
   qr.py resolves `Auto` once it has the fmodel and sets the target.
6. Sets `params.finalized = True` (injected attribute).
After this, no function mutates params. `refine.run` asserts the flag.

The `rst_file` parameter is removed (the restart feature it served is
broken and may be rebuilt later, if ever), and so is `max_atoms` (a
leftover with no purpose); `output_file_name_prefix`,
`output_folder_name`, `write_nuclear_h_positions` stay as command-line
concerns handled by qr.py.

## 5. `check_engine(params)`

Nothing for `restraints=cctbx`. Otherwise, by `quantum.engine_name`:
- `mopac`: `Mopac().get_command()` not None.
- `xtb`: `GFNxTB().get_command()` not None.
- `orca`, `terachem`, `turbomole`, `gaussian`: the executable the plugin
  runs (`orca`, `terachem`, `define`/`ridft`, `g16`/`g09`) found with
  `shutil.which` or the environment variable the plugin honours; the
  plan lists each plugin's lookup after reading it.
- `aimnet2`, `aimnet2-old`, `torchani`, `pyscf`, `fairchem`: import of the
  package (`aimnet2calc`, `torchani`, `pyscf`, `fairchem`), with the
  warning filter aimnet2 uses.
- `server`: an HTTP GET of `quantum.server_url` with a 3 s timeout; a
  connection error is a `Sorry`.
A failure raises `Sorry("QM engine '<name>' is not available: <what was
looked for>. <install hint>")`. The hints: MOPAC and xTB from conda-forge,
AIMNet2 via `config/update_phenix.sh aimnet2`, the others from their
vendors. `refine.set_qm_defaults`'s aimnet2 import check moves here.

## 6. `check_model(model, params)`

Runs on the caller's model (not the API copy) so it can precede the copy.
Raises `Sorry` listing up to ten offenders, with the preparation command:
1. (removed together with `max_atoms`)
2. With `altlocs.handling=conformers` and more than one conformer:
   `conformers.check_completed` (moved from qr.py).
3. With `restraints=qm`:
   a. Missing heavy atoms (`missing_atoms_by_residue(False)`) in a residue
      that carries no capping H (`model.selection("capping_h")` restricted
      to the residue is empty): "truncated residues without caps: run
      mmtbx.reduce2 qmready=true capping=true".
   b. Missing H (`missing_atoms_by_residue(True)`), and waters with fewer
      than two H or D in any conformer of the residue group (the blank
      group's atoms plus each altloc group's; a water with O in the blank
      group and D1/D2 in altloc groups passes): "run mmtbx.reduce2
      qmready=true".
   c. Elements outside the engine's set, when the engine declares one:
      aimnet2 reads `implemented_species` from a freshly constructed
      `AIMNet2Calculator(params.quantum.method)` (one extra model load,
      seconds); torchani: H, C, N, O; the others declare none. D counts as
      H.
4. Crystal data (fmodel, no map) with QM restraints and neither expansion
   nor clustering: `check_qm_sees_crystal_neighbours` (moved from
   refine.run; needs `has_map`, so the signature is
   `check_model(model, params, has_map=False)`).

## 7. qr.py after the change

`Program.run`:
1. `prepare_inputs(data_manager, params, log) -> group_args(model,
   fmodel, map_data, d_min, prefix)`: a module-level function holding the
   data-manager work of today's `run` (fmodel with explicit
   `update_all_scales(remove_outliers=True, update_f_part1=True,
   refine_hd_scattering=False)`; map normalisation; altloc
   reconciliation; electron-table neutralisation; model processing and
   scattering dictionaries). Tests reuse it to build the same objects the
   command line uses.
2. `api.finalize_params(params, d_min=..., ...)`; the target is set on the
   fmodel from the resolved `refinement_target_name`.
3. `api.check_engine(params)`, `api.check_model(model, params, has_map)`.
4. `result = api.refine(model, params, fmodel, map_data,
   working_folder="ase", intermediate_folder=params.output_folder_name,
   log=self.logger)`.
5. Output: `pdb/<prefix>_refined.<ext>` written from `result.model` with
   the sites of `result.sites_cart` (or `sites_cart_nuclear` when
   `write_nuclear_h_positions`), through `model_as_pdb_or_mmcif` so the
   extension is `.cif` when the identifiers do not fit the PDB format;
   real-space mode writes `real_space_refined.<ext>` as today plus the CC
   line. The final log lines ("See ... in ... folder") move here.

`validate` keeps only the data-manager counts. `auto_cust` and the
altloc/completion code leave `Program`.

## 8. Cleanup that makes the API honest

- `driver.py`: no module-level log; `refine(...)` and `opt(...)` take
  `log`; `run_minimize` unchanged. The debug pickle of `opt` goes into the
  working folder.
- `calculator.sites`: `dump_gradients` dumps into the working folder and
  raises `Sorry("dump_gradients: gradients dumped to ..., stopping.")`
  instead of `sys.exit`; tst_36 expects the Sorry.
- `results.manager`: `finalize` and `write_final_pdb_files` are removed;
  `write_pdb_file(folder, name)` stays for intermediates and is called only
  when `intermediate_folder` is set (the driver receives the folder, None
  skips).
- `restraints.py`, `fragment.py`, `cluster_restraints.py`,
  `super_cell.py`, `conformers.py`: every file path derives from the
  working folder passed down; `expansion.pdb` moves into it. No
  `os.getcwd()`-relative write remains (grep in the plan). Every such
  write is skipped when the working folder is None: the scratch model
  `from_qm` writes before each evaluation is a record, not an engine
  input (the engines receive ASE atoms built in memory); fragment,
  cluster, expansion and super-cell files are records (`save_clusters`);
  the point-charge file of charge embedding is an input only for engines
  that read it, which are external-program engines and so always have a
  folder.
- `refine.run` neither creates output folders nor builds a restart path.

## 9. mmCIF identifiers (D14)

Every hierarchy or model write in qrefine goes through
`write_pdb_or_mmcif_file` / `model_as_pdb_or_mmcif`, and the name those
return is the name used afterwards (engine labels are derived from it). The
engines receive ASE atoms built from the hierarchy, so a `.cif` scratch
file changes nothing for them. `charges_class`, `fragment_capping`,
`super_cell`, `conformers` work on hierarchies and need no change beyond
the writes. A test (section 10, tst_57) runs a model with two-character
chain IDs and a four-character residue name (a ligand whose restraints CIF
is a renamed copy of a chem_data monomer) through cctbx restraints and
MOPAC, clustering off and on, and checks the output is mmCIF with the
identifiers intact.

## 10. Tests (qrefine/tests/unit, registered)

- `tst_56` (API): with the tst_xh_mapper fixture and `run_fmodel` data:
  (a) `api.refine` with cctbx restraints equals the command line
  (`run_qrefine_in` on the same files, same params) in `r_work`, `r_free`
  and sites (1e-6 and 1e-6 A) — the objects come from `qr.prepare_inputs`
  on a data manager built from the same files; (b) the same with MOPAC,
  clustering off; (c) the caller's model and fmodel are unchanged after
  the call (sites, hierarchy string, r_work); (d) nothing is written
  outside `working_folder`, and nothing at all with
  `intermediate_folder=None` (directory listing before and after); with
  `working_folder=None` and cctbx restraints, clustering off and on, the
  recursive listing of the run folder and the entries of
  `tempfile.gettempdir()` are identical before and after; with MOPAC and
  `working_folder=None` the run succeeds, the log names the temporary
  folder, and that folder no longer exists afterwards; (e)
  `finalize_params` reproduces qr.py's old `auto_cust` on three
  resolutions (bonds/angles thresholds, expansion, minimizer) and refuses
  expansion with clustering; a raw params extract is refused by
  `api.refine`; (f) `Auto` target resolves to `mlhl` when the fmodel has
  HL coefficients (fmodel built with `abcd` from `run_fmodel` data plus
  phases from the model) and `ml` otherwise; (g) mode=opt (no fmodel)
  returns nuclear sites and `r_work=None`.
- `tst_58` (checks): `check_engine` raises for `quantum.engine_name=xtb`
  when `which("xtb")` is None (skip if present) and a fake server URL;
  passes for cctbx and mopac; `check_model`: the capped tst_173 fragment
  passes with QM restraints; the same with the cap removed raises naming
  ARG 3; a model with its H stripped raises pointing to reduce2; a water
  without H raises; an element outside aimnet2's set raises only when
  aimnet2 is importable (skip otherwise); an altloc-incomplete model
  raises with `conformers` handling.
- `tst_57` (mmCIF): section 9.
- `tst_59` (plumbing): a command-line run leaves nothing in the current
  directory but `pdb/` and `ase/`; `dump_gradients` ends with the Sorry and
  its files under `ase/`; clustering writes `expansion.pdb` under `ase/`.
- `tst_36`: the debug pickles of `driver.opt` are read from `ase/`.
- Existing suite unchanged otherwise: the command line produces the same
  numbers (tst_55 and the X-ray tests re-baselined in sub-project 1 stay
  as they are), the same `pdb/` and `ase/` layout.

## 11. Implementation outline

1. Plumbing: driver log, working folder everywhere, intermediates optional,
   dump_gradients Sorry, debug pickle path; tst_36 adapted.
2. `api.finalize_params` from validate/auto_cust/set_qm_defaults/run;
   `rst_file` removed; `Auto` target; tests (e).
3. `api.check_engine`, `api.check_model`; tst_58.
4. `refine.run` returns a result; `api.refine` with copies; qr.py split
   into `prepare_inputs`/run/write; results.manager without finalize;
   tst_56.
5. mmCIF audit and tst_57.
6. Full suite; re-baseline nothing (numbers must not move).

## 12. Decisions to confirm

- D16. `rst_file` removed from the PHIL (unknown-parameter error for old
  scripts that pass it). Confirmed 2026-11-05.
- D17. `max_atoms` removed: the parameter and every check. Confirmed
  2026-11-05.
- D18. `dump_gradients` ends with a `Sorry` (non-zero exit) instead of
  `sys.exit(0)`; a `debug`-level command-line feature phenix.refine never
  reaches. Confirmed 2026-11-05.
- D21. `working_folder=None`: no file at all from qrefine; in-process
  engines write nothing; external-program engines run in a temporary
  folder removed on return. Proposed 2026-11-05 (the alternative is to
  refuse external-program engines without a folder).
- D19. `check_model` loads the aimnet2 model once more to read its species.
- D20. The API returns its own refined model object in addition to the
  site arrays.
