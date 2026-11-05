# qr.refine inside phenix.refine: umbrella design

Date: 2026-11-04. Status: decisions approved in conversation 2026-11-04.
This note records the whole picture; each sub-project below gets its own
spec and plan. Repositories: cctbx_project (mmtbx), qrefine, phenix,
phenix_regression.

## 1. Goal

phenix.refine gets an optional mode, `main.use_qrefine=True`, in which the
"XYZ individual (reciprocal space)" step of the last macro-cycle is done by
qr.refine (any QM engine, AQuaRef included) instead of the cctbx coordinate
minimization. Everything else phenix.refine does (bulk solvent, ADP,
occupancy, anomalous, output) stays as it is. Users get the geometry of
QM-restrained coordinates together with the rest of phenix.refine.

qr.refine refines coordinates only, today only from the command line. It
becomes callable from Python with objects in and objects out.

## 2. Decisions

D1. **Switch and placement.** `main.use_qrefine=False` in phenix.refine. When
True, the coordinate step of the last macro-cycle (macro_cycle ==
main.number_of_macro_cycles) calls qrefine; all earlier macro-cycles run
exactly as today, riding H included. The step runs where xyzrec runs now:
after weights, before ADP and occupancy refinement. Running qrefine in more
macro-cycles is a later option, not built now.

D2. **Scope.** Reciprocal-space phenix.refine only. phenix.real_space_refine
later, if ever.

D3. **Model preparation is external.** The input is prepared by
`mmtbx.reduce2 qmready=true` (with `capping=true` for truncated residues).
phenix.refine checks the model and refuses with a Sorry that names the
preparation step; it never rewrites the model (no altloc reconciliation, no
H placement inside phenix.refine).

D4. **Hand-off convention.** qrefine's API returns refined sites and a
summary. It never mutates the caller's model or fmodel; it deep-copies both.
phenix.refine applies the sites to its own objects. Standalone qr.refine
uses the same API and writes its files afterwards.

D5. **Hydrogens.** No riding model ever runs inside qrefine. Inside qrefine
the working coordinates are nuclear at all times; the only consumer of
shortened X-H positions is the scattering side (fmodel: Fcalc, R factors,
maps, data target and gradients), and only for X-ray scattering tables.
Neutron and electron tables: nuclear everywhere, as today. On output
(hand-off to phenix.refine or file) H positions follow the scattering-table
convention: X-ray lengths for X-ray tables, nuclear otherwise. A parameter
forces nuclear output for comparison with old runs. In phenix.refine the
model therefore goes in with X-ray H positions, is expanded once to nuclear,
refined, and comes back shrunk; phenix.refine never learns about nuclear
positions. H ADPs stay riding, H occupancies untouched.
Sub-project 1 (`2026-11-04-xh-mapper-design.md`) builds this.

D6. **Mapping direction.** Nuclear positions are the refined parameters and
the data gradient is chain-ruled back. The alternative (X-ray parameters,
elongated H for QM) leaves the parameter X-H length unrestrained and it
drifts. See the sub-project 1 spec for the mapping and its Jacobians.

D7. **Weights.** phenix.refine computes no wxc for the qrefine step; qrefine
runs its own protocol (gradient-norm data weight, restraint-scale search
accepted on bond rmsd, refine cycles). One qrefine macro-cycle means the
full qr.refine protocol once. The weights step of phenix.refine still runs
because it also produces the ADP weight.

D8. **Ordered solvent and ions.** Never updated when use_qrefine is on:
`main.ordered_solvent` and `main.place_ions` are set False at startup with a
prominent log line (not a Sorry: ordered_solvent defaults to True). Reason:
waters added as bare O would reach the QM engine without H. Adding H to new
waters before the qrefine cycle is a possible later extension.

D9. **Steps gated off in the qrefine macro-cycle** (task table flag
`qrefine_cycle`): setup riding H model (setrh), ensure riding H (regHxyz),
NCS rotamer correction (ncsrot), global and local real-space refinement
(realsrg, realsrl), post-minimization rotamer correction (realsrl2),
simulated annealing torsion and Cartesian (tardy, sacart), DEN (den).
Everything else runs in today's order: bss, settarget, restraint updates,
QI (refused, see D10), NQH flips, cablam, update SS, rotamer and RDL
updates, rigid body, NCS operators, weights, fit rotatable H, the QM
coordinate step, then ADP, riding H ADP regularization, H/D recentering,
occupancies, IAS steps (refused), fp/fdp. NCS, secondary structure,
Ramachandran, reference model and rotamer restraints need no global
disabling: they act in the standard cycles and are simply not consulted by
qrefine. The riding manager built in standard cycles keeps stale torsion
parameters after the qrefine step; no standard cycle follows, and the
implementation verifies that the output path performs no final H
idealization.

D10. **Refused at startup** (Sorry, before any data processing): qrefine
not importable; the requested engine unavailable (`check_engine`);
twin_law set; joint X-ray plus neutron data; the `qi` quantum-interface
scope active; amber; `main.ias`; `den` in strategy; a
`refine.sites.individual` selection combined with cluster.clustering or
with exchangeable H/D sites (qrefine's exclude selection supports
neither). After the model is loaded and processed: `check_model` (altloc
completion needed, uncapped missing atoms, elements the engine cannot
handle, more than max_atoms atoms). `individual_sites_real_space` in the
strategy is not refused: its steps are gated off in the qrefine cycle and
run before.

D11. **PHIL in phenix.refine.** Self-contained, no include of qrefine's
PHIL, so phenix.refine starts without qrefine installed:

    main.use_qrefine = False
    qrefine {
      restraints = cctbx *qm        # cctbx: fast path for tests and smoke runs
      qm_engine = *aquaref mopac xtb orca pyscf torchani terachem turbomole gaussian server
      clustering = Auto             # Auto: engine auto-cust default
      max_bond_rmsd = Auto          # Auto: by resolution, as qr.refine auto-cust
    }

`aquaref` maps to qrefine engine `aimnet2` plus its auto-cust settings.
Expansion is not a user choice: QM restraints imply expansion unless
clustering is on; cctbx restraints handle symmetry natively and never need
it. `parallel.nproc` and every other qrefine parameter keep qrefine's
defaults (exposing nproc or an expert pass-through string is a later
option).

D12. **Capping H** (as implemented 2026-11-04, see the capping-H spec):
occupancy 0, so no Fcalc contribution, zero data gradient on position and
ADP, excluded from the bulk-solvent mask (ignore_zero_occupancy_atoms),
from occupancy refinement and from probe/clashscore (minimum_occupancy
0.02). Restraints: bond and angle proxies with origin id `capping` plus
nonbonded as H; they act in the cctbx coordinate target, are excluded
from bond/angle statistics and listed under "Capping H" in the .geo file.
Not riding: the riding manager cannot parameterize them, so they are free
atoms held by their own restraints. H ADP regularization covers them
through the X-H connectivity table (verify in sub-project 3). In the
qrefine step a cap is an ordinary H to the QM engine. The final
phenix.refine output keeps the caps; the log says how to strip them
(`phenix.pdbtools remove=capping_h`). `check_model` treats a capped
truncation as complete and an uncapped one as an error pointing to
reduce2 capping=true.

D13. **R-factor consistency at hand-off.** After applying the returned
sites, phenix.refine updates f_calc and f_mask (forced), re-scales with the
same call qrefine's driver uses after its last cycle
(`update_all_scales(params=<bulk_solvent_and_scale>, remove_outliers=False,
refine_hd_scattering=False)`, no back-trace), and asserts r_work and r_free
equal to the returned values within 1e-4 (the tolerance qrefine's driver
itself uses for the same comparison). Because the API re-scales with the
caller's own bulk-solvent parameters (D15), the values are expected to
agree to minimizer noise; if the regression data show otherwise, the
fallback is for the API to also return the scale parameters it ended with
and for the hook to apply them; the assert stays.

D14. **mmCIF identifiers.** phenix.refine accepts mmCIF models with
two-character chain IDs and residue names longer than three characters.
qrefine round-trips hierarchies through PDB text in places (engine input
files, cluster files, charge calculation). Sub-project 2 audits every such
place and adds a test; qrefine is fixed where it breaks.

D15. **fmodel setup: one convention.** Both programs build the fmodel with
the same `data_manager.get_fmodel`; the differences were only in what each
did afterwards, with no need behind them. Two measures remove them:
(a) called from phenix.refine, the API uses the caller's fmodel as is (deep
copy): target, scales, mask, accuracy and outlier state; the only thing the
API does to an fmodel on its own is re-scale between its cycles, for which
it takes the bulk-solvent parameters (`bulk_solvent_and_scale` scope) from
the caller and holds the rest fixed (remove_outliers=False,
refine_hd_scattering=False, no back-trace); (b) standalone qr.py sets up
its fmodel as phenix.refine does, at the cctbx level with no phenix
dependency: default scattering table n_gaussian (was wk1995), scattering
dictionaries at the data resolution (was a fixed 1.0 A), the target chosen
by the same rule (ml; mlhl when experimental phases are present; ls on
request), explicit `update_all_scales` arguments equal to phenix.refine's
defaults. (b)'s table and resolution change is done in sub-project 1 so
the X-ray tests are re-baselined once; the rest of (b) and (a) are
sub-project 2. tNCS correction stays phenix.refine-only (off by default).

## 3. Sub-projects, in order

1. **X-H mapper** (`mmtbx.hydrogens.xh_mapper`, name open) and its use in
   qrefine's reciprocal-space calculator, driver and output. Spec:
   `2026-11-04-xh-mapper-design.md`, which also switches qr.py's default
   scattering table to n_gaussian with dictionaries at the data resolution
   (D15b). Done and tested first; it changes qr.refine's X-ray results
   (start R factors, scattering table, refined H positions), re-baselined
   once.
2. **qrefine callable API** (`qrefine.api`): `refine(model, fmodel, params,
   working_folder, log, bulk_solvent_params=None, write_intermediate=False)`
   (D15a; None means cctbx defaults) returning a group_args
   (sites_cart in the caller's convention, sites_cart_nuclear, r_work,
   r_free, data_weight, restraints_weight_scale, bond and angle rmsd,
   n_fev, n_cycles); `finalize_params(params, d_min, scattering_table,
   has_map, log)` absorbing qr.py's validate, auto_cust and
   set_qm_defaults; `check_engine(engine_name)`; `check_model(model,
   params, engine_name)`. Cleanup it requires: driver.py takes a log
   argument instead of its module-level sys.stdout; results.manager writes
   intermediate PDBs only when asked; the restart pickle goes; dump_gradients
   raises instead of sys.exit; the working folder is an argument wherever
   "ase" and "pdb" are hard-coded. qr.py becomes: data manager in,
   finalize_params, check_engine, check_model, api.refine, write outputs.
   fit_rotatable2 and the conditional pre-optimization stay inside the API
   as part of the protocol. qr.py's fmodel setup follows D15b (target rule,
   explicit scaling arguments). The mmCIF audit (D14) is part of this
   sub-project.
3. **phenix.refine hook**: PHIL (D11), startup checks (D10), task-table
   flag and gating (D9), `phenix/refinement/xyz_qrefine.py` called from
   macro_cycle.xyzrec (translates the scope into qrefine params, builds the
   exclude selection from refinement_flags.sites_individual, calls
   api.refine with a working folder named after the output prefix, applies
   the sites, re-scales and asserts per D13, syncs model and fmodels,
   recreates target functors, prints a summary line in the table style),
   regression tests, documentation.

## 4. Tests

Sub-project 1 (details in its spec): mmtbx tst_xh_mapper (finite
differences of the chain rule; shrink/expand lengths; identity for
neutron and electron; caps, water, D, altlocs, unbonded H); qrefine
finite-difference test of the calculator with the mapper; end-to-end X-ray
run checking output H lengths and R-factor consistency; the neutron and
electron tests of the existing suite unchanged.

Sub-project 2: command line and API give identical results on the same
inputs (cctbx restraints and MOPAC); refine() leaves its inputs untouched
and writes nothing with write_intermediate off; check_engine and
check_model Sorry cases, caps accepted, uncapped truncation refused;
finalize_params reproduces qr.py's auto-cust; mmCIF model with
two-character chain IDs and four-character residue names through qr.refine
with cctbx restraints and MOPAC.

Sub-project 3 (phenix_regression):
- use_qrefine with cctbx restraints over three macro-cycles: only the last
  cycle used qrefine (log), gated steps absent in that cycle and present
  before, output H lengths are X-ray, R factors sane;
- the same with MOPAC, skipped when MOPAC is absent;
- use_qrefine with X-ray, electron and neutron data;
- one macro-cycle, strategy=individual_sites, use_qrefine, against X-ray,
  electron and neutron data, compared with standalone qr.refine on the same
  inputs: refined coordinates equal within 1e-4 A and R factors within
  1e-5, numerical noise only. With D15 the two fmodel setups coincide; the
  test disables the phenix.refine steps that move atoms before the
  coordinate step (NQH flips, cablam) and passes identical qrefine
  settings to both. If the gate ever has to be loosened, the test states
  which difference remains and why;
- the in-code R-factor assert of D13 exercised by every run above;
- a capped model; a sites.individual selection; number_of_macro_cycles=1;
- every Sorry of D10; ordered_solvent and place_ions auto-off logged;
- the full existing phenix.refine suite unchanged with use_qrefine off.

## 5. Deferred

H on waters added by ordered solvent; qrefine in more than the last
macro-cycle and weight reuse across cycles; real-space refinement; GUI
exposure of the scope; exposing parallel.nproc or an expert pass-through
string; the real-space mode of standalone qr.refine with an X-ray table
(stays nuclear).
