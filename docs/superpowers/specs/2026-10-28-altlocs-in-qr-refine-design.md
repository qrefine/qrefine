# Alternative conformations in qr.refine: design

Date: 2026-10-28. Status: approved 2026-10-28 (all decisions settled).

## 1. Goal

qr.refine refines models with alternative conformations (altlocs) using a
restraints source that only understands single-conformer models (a QM engine:
AQuaRef, xTB, and MOPAC for development). It must work with every restraints
route qr.refine has: plain (asymmetric unit), `expansion=True` (super sphere)
and `cluster.clustering=True` (fragments), for both `restraints=qm` and
`restraints=cctbx`, in `mode=opt` and `mode=refine` (reciprocal and real space).

The approach was agreed on 2026-10-22 and prototyped in
`tests/unit/tst_31.py` (`gradients_from_conformers`):

1. Relabel altlocs of interacting residues so that the same label means the
   same physical state (`iotbx.pdb.hierarchy.root.reconcile_altlocs`, done and
   committed in cctbx).
2. Split the model into complete single-conformer models, one per label:
   conformer X = atoms with blank altloc + atoms with altloc X.
3. Compute restraint gradients for every conformer model with the ordinary
   machinery (plain / expansion / clustering, cctbx or QM).
4. Combine per atom: g_i = (1/m_i) Σ_{c ∋ i} g_{c,i}, where m_i is the number
   of conformer models that contain atom i (altloc atoms: 1; blank atoms: the
   number of conformers). Only gradients are combined. There is no energy.
   Occupancies are never used.

## 2. What exists today

- `refine.run` computes `altlocs_present` and, when heavy-atom altlocs are
  present, `create_restraints_manager` returns `restraints.from_altlocs2`
  ("developmental, NOT IN PRODUCTION"). It implements `altloc_method=subtract`
  (needs a gradient of the blank-only model, impossible with QM) and
  `altloc_method=average` (marked "does not work"), each conformer through
  `from_expansion` only. Clustering is not reachable with altlocs.
- `fragment.py` carries a second, independent altloc scheme ("Min's code"):
  per-altloc hierarchies, overlap clusters, `fragment_scales` of -1 or
  1/frequency, `bond_with_altloc` (an O(N·M) Python distance loop). The capper
  refuses fragments with two or more altlocs, so this path cannot cap.
- H/D exchange sites (neutron) are handled by `refine.hd_mapper`: the model is
  de-deuterated into a single protonated model for restraints and gradients are
  expanded back. It is only used when the model has *no* other altlocs.
- `calculator.sites` and `calculator.sites_opt` hand the full-model
  `sites_cart` to `restraints_manager.target_and_gradients` and expect
  full-model gradients back. `driver` touches the restraints manager only for
  re-clustering (`fragment_manager.set_up_cluster_qm()`, `interaction_list`).
- QM scratch files live in `./ase/` (`ase/tmp_ase.pdb`, `ase/<i>/` per
  fragment); file names double as MOPAC job labels.
- Tests of the old scheme: `tst_05` (subtract vs cctbx), `tst_06`
  (from_altlocs2 minimization equals from_cctbx), `tst_11` (average, no
  assertions), `tst_25` (clustering vs expansion with subtract,
  bond_with_altloc=False). `tests/unit/broken_to_fix/tst_43` (fragments with
  altloc_method) and `tst_47` (qr.refine on altlocs.pdb) are altloc related;
  the rest of `broken_to_fix` is not (TorchANI, PySCF, pair_interaction,
  charges, qr.finalise, qr.gtest).

## 3. Design

### 3.1 One layer above the existing routes

A new module `qrefine/conformers.py` owns everything altloc specific. The
rest of qrefine keeps seeing single-conformer models.

```
refine.run
  model_for_restraints  (exclude_selection applied, then hd_mapper applied)
  create_restraints_manager(params, model, working_folder="ase")
    if conformers.n_conformers(model) > 1 and params.altlocs.handling == "conformers":
      return conformers.from_conformers(
        model   = model,
        build   = lambda conformer_model, folder:
                    create_restraints_manager(params, conformer_model, folder),
        log     = log)
    ... existing plain / expansion / clustering code, unchanged ...
```

`from_conformers` builds one ordinary restraints manager per conformer by
calling back into `create_restraints_manager` with a single-conformer model and
its own scratch folder (`ase/A`, `ase/B`, ...). Every route (plain, expansion,
clustering; cctbx or QM) is therefore supported without touching it, and new
routes get altloc support for free.

`from_conformers.target_and_gradients(sites_cart)`:

```
g = 0 (full size); m = 0 (full size)
for conformer c:
  e_c, g_c = manager_c.target_and_gradients(sites_cart.select(isel_c))
  g[isel_c] += g_c ; m[isel_c] += 1
return 0, g / m          # energy undefined, as from_expansion
```

It also provides `energies_sites` (real-space refinement calls it), and the
re-clustering hooks `driver` needs (see 3.6).

### 3.2 Conformer models

`conformers.split(model)` returns, for each conformer label X, a
`group_args(label, model, iselection)` where

- *conformer labels* are the altloc labels that have at least one non-hydrogen
  atom (the rule `reconcile_altlocs` uses: hydrogen-only labels such as
  neutron H/D exchange sites are not conformers);
- `iselection` selects blank atoms, atoms with label X, and atoms of
  hydrogen-only labels, in model order;
- `model` is `model.select(iselection)` with every remaining atom group relabeled
  to blank and merged into one atom group per residue group. The conformer model
  is a genuine single-conformer model: no downstream code (fragments, capper,
  charges, QM engines) ever sees an altloc label. Atom order is preserved
  (verified with `remove_alt_confs`-style merging on `h_altconf_complete.pdb`).

Completeness: every residue group that has conformer labels must have *all*
of the model's conformer labels; otherwise conformer models would lack atoms
(a residue with A/B/C in a model where the others have A/B leaves conformer C
without those residues). This raises `Sorry` naming the residue, as agreed on
2026-10-22 (decision D3 below). Planned extension, not in this version: fill
the labels such a residue lacks from its first label (the residue then
contributes the same atoms to the conformers it does not distinguish). The
`split` function is the single place to add it.

Multiplicity m_i therefore is 1 for atoms with a conformer label and the
number of conformers for everything else. The number of conformer labels is
not limited: A/B is the common case, but A/B/C/D... is handled the same way
(one conformer model and one restraints manager per label), and the tests
exercise two, three and four labels.

### 3.3 Reconciling labels first

`qr.py` calls `hierarchy.reconcile_altlocs(crystal_symmetry, log)` on the
input model before `model.process(make_restraints=True)` (restraints must be
built after relabeling) and prints the relabeled residues. Unresolved groups
(no relabeling removes the clash) are printed as a warning and refinement
continues (decision D4). Controlled by `altlocs.reconcile` (3.5). It applies
to both handling routes. The written output model carries the reconciled
labels.

### 3.4 H/D exchange together with heavy-atom altlocs

Composition order in `refine.run`: `exclude_selection` → `hd_mapper` →
conformer split. `hd_mapper` already pairs H/D atom groups regardless of other
altlocs; after it, the single protonated model has only conformer labels, and
3.2 applies. The `altlocs_present_only_hd` gate becomes "exchangeable H/D
present" (use `hd_mapper`) and `n_conformers > 1` (use `from_conformers`),
independently. A small neutron-style test file (H/D on a few sites of
`h_altconf_complete.pdb`) covers the mixed case.

### 3.5 Two ways to handle altlocs with cctbx restraints

A new PHIL scope selects the route:

```
altlocs {
  handling = *conformers native
    .type = choice(multi=False)
    .help = conformers: split the model into single-conformer models, compute \
            gradients for each and combine them (works with any restraints \
            source). native: hand the model with altlocs to the restraints \
            source as is. Only restraints=cctbx can do this.
  reconcile = True
    .type = bool
    .help = Relabel altlocs of interacting residues so that the same label \
            means the same conformer (iotbx.pdb.hierarchy.reconcile_altlocs) \
            before restraints are built.
}
```

`conformers` (default) is the path of 3.1 and is the same for cctbx and QM.
With cctbx it means the model is "pretended" single-conformer-only, which
makes cctbx an exact, fast stand-in for QM in tests.

`native` keeps cctbx's own altloc-aware restraints (one restraints manager on
the whole model) and is the reference the conformer path is compared against.
It is available for every route: plain and expansion already work on models
with altlocs; clustering works once the old per-altloc fragment code is
removed (fragments then contain all conformers of their residues, and the
super-sphere restraints manager selects them like any other atoms; the
capper is not involved with cctbx). `handling=native` with `restraints=qm`
raises `Sorry`.

Expected relation between the two (quantified by `tst_31`): altloc atoms
get identical gradients; for blank atoms m·g_conformers = g_native +
(m−1)·g_blank_only, so they differ only where blank atoms interact with
altloc atoms, and the refined models are very similar. A test (T2b) checks
both: gradient agreement per atom class and the closeness of two refinements
of the same model, one per route.

The energy is undefined on the conformer path, so `handling=conformers` with
altlocs forces `minimizer=lbfgs` and `gradient_only=True`, as expansion and
clustering already do (`qr.validate`). `native` keeps the user's minimizer.

### 3.6 Driver and calculator

- `driver.refine` / `driver.opt` reach into `restraints_manager.fragment_manager`
  for re-clustering and a log line. `from_conformers` exposes a small
  `fragment_manager` proxy with `set_up_cluster_qm()` (re-clusters every
  conformer) and `interaction_list` (concatenated). No change to the driver
  logic; the `try/except` on `restraints_manager.clustering` stays as is.
- `calculator.sites`, `sites_opt`, `sites_real_space`, `compute_weight`: no
  change. They already pass the full `sites_cart` and consume full gradients.
- `results.manager`, `cctbx_geometry` (bond/angle rmsd), `fit_rotatable2`,
  pre-optimization: operate on the full model with altlocs as today.

### 3.7 Scratch folders and QM job labels

`restraints.restraints` / `from_qm` / `fragments` get the scratch folder from
one `working_folder` argument (default `"ase"`, so non-altloc runs and
`tst_35`'s `ase/[0-9]*` layout are unchanged). Conformer X uses
`ase/<X>/`, so MOPAC labels and fragment files never collide between
conformers.

### 3.8 Removal of the old altloc machinery

Removed with their tests (decision D2):

- `restraints.from_altlocs2`, `from_cctbx_altlocs`, `from_altlocs`,
  `get_cctbx_gradients(rm_only=...)` leftovers that only they use;
- `cluster.altloc_method`, `cluster.bond_with_altloc` PHIL parameters and
  their plumbing (`refine.create_fragment_manager`, `fragments.__init__`);
- in `fragment.py`: the per-altloc `phs` loop, `atoms_overlap`,
  `bond_with_altloc`, `altloc_atoms`, `fragment_scales` other than 1.0 (the
  `gradients_scale` multiplication in `from_qm` / `from_cctbx` becomes a
  no-op and is removed), `check_no_altlocs` in `cluster_restraints.py`;
- `write_mm_charge_file`'s altloc branches (the super sphere is
  single-conformer now).

The `altlocs_present` argument of `create_restraints_manager` goes away; the
function decides from the model and `altlocs.handling`. After the removal,
`fragments` treats a model with altlocs like any other (all conformers of a
residue belong to the same cluster and fragment), which is what
`handling=native` with clustering needs.

### 3.9 Error handling

`Sorry` with a precise message for: a residue whose conformer labels are a
proper subset of the model's (3.2); `exclude_selection` with H/D (existing);
a QM run that fails for one conformer propagates as today. Unresolved
`reconcile_altlocs` groups: warning only.

### 3.10 Performance

Setup and each gradient evaluation cost `n_conformers` times a single-conformer
run (two for the usual A/B case). Clustering builds `n_conformers` fragment
managers (each runs `pair_interaction` on its super sphere). Conformers are
processed one after the other in this version (decision D5).

While implementing, profile one altloc run per route and fix what is in the
altloc path. Already identified: `fragments.bond_with_altloc` (O(N·M) Python
loop, removed); `copy.deepcopy(ph)` per cluster in `fragments.get_fragments`
before each `pair_interaction.run`; `from_qm.target_and_gradients` writing a
PDB on every call even when the engine reads XYZ. Anything outside the altloc
path is reported, not changed.

## 4. Tests

All new tests follow `tests/unit/tst_00.py`: models as PDB strings inside the
test file (`pdb_str_...`), `run(prefix=...)` that creates and enters its own
folder, `run_fmodel` / `run_qrefine` from `qrefine.tests.unit` for
end-to-end runs, plain `assert`s, `run()` under `__main__`, registered in
`run_tests.py`. They run in about a minute, at most two, so they use small
models, in-process API calls (`refine.create_restraints_manager`) where the
check is about gradients, and at most two or three `qr.refine` invocations per
test (each costs ~10 s of start-up). MOPAC is the QM engine.

Unit tests of `conformers.py` (cctbx as the calculator, from `tst_31`):

- T1 split: labels, selections, multiplicities, single-conformer models in
  model order; H-only labels are not conformers; completeness `Sorry`.
- T2 combination vs cctbx on the six altloc data files: altloc atoms equal
  cctbx; blank atoms satisfy m·g = g_cctbx + (m−1)·g_blank_only; no altlocs:
  one call on the whole model (`tst_31`, moved to the new module's test).
- T2b conformers vs native (cctbx): on the altloc models, altloc atoms get
  identical gradients and blank atoms satisfy the relation of 3.5; then the
  same model is minimized once per `altlocs.handling` (in-process, as
  `tst_06` did) and the two results agree to a small rmsd. Both with and
  without `expansion`, and with `clustering`.
- T3 duplicate-conformer identity: a single-conformer model whose residues are
  duplicated into identical A and B conformers gives exactly the gradients of
  the original model. This is the key end-to-end check and is run through
  every route: plain, expansion, clustering (cctbx) and plain, expansion,
  clustering (MOPAC on a small peptide, one gradient evaluation each). It
  also shows occupancy independence (occupancies of the copies are changed).
- T4 clustering vs expansion agreement with real altlocs (cctbx; the altloc
  data files), the analogue of `tst_04`/`tst_29`/`tst_25`, for both
  `altlocs.handling` values.
- T5 H/D exchange with heavy altlocs: hd_mapper then conformer split; gradient
  equals the gradient of the de-deuterated model with conformers (new small
  neutron-style file).

Integration tests through `qr.refine`:

- T6 `restraints=cctbx` refinement of `h_altconf_2_complete_modified.pdb`
  (reciprocal space, `tst_06` replacement): runs, output keeps altlocs and
  occupancies, bond rmsd and R drop; the same via `expansion=True` and via
  `clustering=True`; `altlocs.handling=native` once as well, and
  `handling=native restraints=qm` gives the `Sorry`.
- T7 `restraints=qm quantum.engine_name=mopac mode=opt` on a two-residue altloc
  model (`gly2_1.pdb`), `number_of_micro_cycles=1`, with `clustering=True` and
  with `expansion=True`; scratch folders `ase/A`, `ase/B` exist and the
  gradient of the duplicated-conformer model equals the single-conformer run.
- T8 `altlocs.reconcile`: a two-residue model with swapped labels is relabeled
  by `qr.refine` (log line) and the output model carries the new labels;
  `reconcile=False` leaves them.
- T9 `tst_47` (qr.refine on `altlocs.pdb`, labels only, no blank atoms) revived.

Existing tests: `tst_04`, `tst_26`, `tst_29`, `tst_32`–`tst_38` keep passing
unchanged apart from the removed `altlocs_present` argument. `tst_05`,
`tst_06`, `tst_11`, `tst_25` and `broken_to_fix/tst_43`, `tst_47` are replaced
by T2, T6, T4, T9. The unrelated files in `broken_to_fix` are left alone.
`tst_31` becomes the unit test of `conformers.py`.

## 5. Implementation outline

1. `conformers.py` with `split`, `n_conformers`, `from_conformers` and its
   test (T1–T3 cctbx part, T2b gradients), driven by `tst_31`'s prototype.
   No qr.refine changes yet.
2. `working_folder` plumbing through `restraints.restraints`, `from_qm`,
   `create_fragment_manager`; `create_restraints_manager` recursion; remove
   `altlocs_present`; T3 clustering/expansion (cctbx), T4.
3. `refine.run` composition (exclude → hd_mapper → conformers),
   `altlocs.handling` dispatch and the `native` + QM `Sorry`, `qr.validate`
   forcing lbfgs/gradient_only, driver proxy hooks; T2b minimizations, T5,
   T6.
4. MOPAC: T7.
5. `reconcile_altlocs` in `qr.py`; T8, T9.
6. Remove the old machinery (3.8) and its tests; update `run_tests.py`;
   profile and fix hot spots in the altloc path; run the whole `qr.test`.

Each step ends with the full test suite green. Changes stay uncommitted.

## 6. Decisions to confirm

- D1 (settled 2026-10-28) both routes are available through
  `altlocs.handling`; the conformer path is the default and a test compares
  the two (3.5, T2b).
- D2 (settled) the old altloc machinery and its tests are removed (3.8).
- D3 (settled) incomplete conformers raise `Sorry`; filling the missing
  labels from the first label is the planned extension noted in 3.2.
- D4 (settled: warn and continue) `reconcile_altlocs` tries every relabeling of a group of interacting
  residues and keeps the one with the fewest clashes between paired labels.
  When even the best one leaves two paired conformers clashing badly (the
  labels cannot be made consistent: for example residue 10 A clashes with
  residue 20 A, and B with B, whichever way they are paired), the residues
  are reported as unresolved and left unchanged. The choice is whether
  qr.refine then stops (`Sorry`) or prints the report and continues; with
  QM the clash is then simply part of one conformer's geometry. Recommended:
  warn and continue.
- D5 (settled: sequential) In `from_conformers.target_and_gradients` the conformers are processed
  one after the other: the QM gradient of conformer A, then of conformer B.
  They are independent, so they could run at the same time (two QM jobs at
  once, halving the wall time when the machine has the resources). This
  version runs them one after the other; running them in parallel would
  later change only that loop (e.g. with `libtbx.easy_mp`), nothing else.
