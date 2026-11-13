# Clustering with alternative conformations: per-fragment conformer variants

Date: 2026-11-07. Status: approved in conversation 2026-11-07, implemented
(Tasks 1-14); sections 3.2, 3.3, 3.4, 3.7, 4.1, 5, 6 and 7 corrected
2026-11-12 from the task reports (where a report and the first draft
disagree, the report is right).

## 1. Goal

qr.refine must refine a model with any number of altloc labels on any
residues with a QM engine through `cluster.clustering=True`, at a cost that
scales with how much of the model has altlocs, not with the number of
labels. Today the clustering route reaches altloc models only through
`conformers.from_conformers`: every conformer is a whole single-conformer
model, clustered and evaluated on its own, one after the other. Setup,
memory and QM jobs are multiplied by `n_conformers` even when a dozen
residues have altlocs.

Target: a clustered run on a model with few altloc residues costs about a
clustered run on the same model without altlocs; a model with altlocs
everywhere (1PWM: 168 of 912 residue groups, 149 of 155 fragments touched)
costs at most `n_conformers` times, the physical minimum since a QM energy
is not additive. Setup is paid once, all QM jobs of all conformers run in
one parallel pool. Everything must be runtime and memory efficient enough
not to become the obstacle when ribosome-scale models are attacked later
(a separate project).

Measured 2026-11-07 on 1PWM (8452 atoms, labels A/B): one clustering of the
altloc model takes 123 s and gives 155 fragments; per-fragment variants
give 304 QM jobs per evaluation against 310 for the per-conformer scheme,
so the saving there is the second clustering and the second sequential
pass, not the job count (the implementation found 373 jobs, section 7). On
a typical 1.5-2.5 A model with a handful of altloc residues the job count
approaches that of a single-conformer run.

## 2. What exists today

- `conformers.py`: conformer labels, per-residue-group rule (blank atoms and
  H-only-label atoms to every conformer; own label; fallback to the first
  own label for a residue group lacking the label), `conformer_selections`,
  `multiplicity`, `check_links`, `check_completed`, `split` (relabel to
  blank and merge atom groups), `from_conformers` (one restraints manager
  per conformer, gradient g_i = (1/m_i) sum over conformers), and
  `fragment_manager_proxy` for re-clustering of per-conformer fragment
  managers.
- `refine.create_restraints_manager`: with more than one conformer and
  `altlocs.handling=conformers` it returns `from_conformers`, which calls
  back into `create_restraints_manager` per conformer model; with
  `clustering=True` each call builds a `fragments` object and a
  `from_cluster`.
- `fragment.fragments`: super sphere (`super_cell.expand`), interaction list
  (`mmtbx.pair_interaction`), clusters (`clustering.py`), fragments (cluster
  plus interaction buffer), per-fragment capper (`fragment_capping.capper`,
  cut bonds against the super cell) and charge, selection lists consumed by
  `cluster_restraints.from_cluster`, `restraints.from_qm` and
  `restraints.from_cctbx`.
- `cluster_restraints.from_cluster.target_and_gradients`: one `parallel_map`
  over the fragments; gradient of a cluster atom from its fragment, buffer
  atoms zero, fragment energies summed.
- `handling=native` with clustering (cctbx only): `fragments` on the altloc
  model, no conformer awareness; pdb_interpretation's altloc-aware
  restraints of the super sphere selected per fragment.
- Tests of the current clustering with altlocs: tst_43 (qr.refine, cctbx,
  conformers with clustering and native plain), tst_46 (frozen-environment
  finite differences, expansion and clustering routes, cctbx both handlings
  and MOPAC conformers), tst_31 (relation between standard and conformer
  gradients).

## 3. Design

### 3.1 Architecture

`fragment.fragments` gains optional conformer awareness. Given the model's
conformer selections (`conformers.conformer_selections`, the existing rule),
it builds variants per fragment; without them it behaves exactly as today.
`conformers.py` remains the single owner of the altloc rules and gains two
helpers: tiling a master conformer selection onto the super cell and super
sphere, and relabeling an arbitrary hierarchy selection to a
single-conformer hierarchy (the body of today's
`_single_conformer_hierarchy`, generalized). `from_conformers` is unchanged
and keeps serving the plain and expansion routes (whole-model routes have
one fragment, so their variants are the conformers, which is what
`from_conformers` computes).

Dispatch in `refine.create_restraints_manager`:

```
if n_conformers > 1 and handling == "conformers":
  if clustering:
    conformers.check_links on the model (same Sorry as split)
    fragment_manager = fragments(..., conformer_selections = conformers.conformer_selections(ph))
    return from_cluster(restraints_source.restraints_manager, fragment_manager, ...)
  else:
    return from_conformers(...)        # plain and expansion, as today
... existing code ...
```

`handling=native` builds `fragments` without conformer selections (today's
behaviour, cctbx only). The phenix.refine hook (`xyz_qrefine.make_params`)
passes `clustering` through and needs no change.

### 3.2 Setup, once per (re-)clustering

One super sphere of the altloc model and one clustering; the interaction
list, the bonded residue pairs and the buffer search run per single-conformer
view and are united (3.4). Then, for every fragment F (a
selection on the super sphere) and every global conformer c:

- `S_c`: conformer c's master selection (flex.bool over the model) tiled
  onto the super cell. The super cell is the master followed by chain-wise
  copies of the whole model in model order, so atom `k*N + i` of the super
  cell is master atom `i`; the tiling is index arithmetic. The super-sphere
  selection is `S_c` cut by `SCM.keep_selection`.
- variant atom set `V_c(F) = sel_F & S_c` (on the super sphere).

Equal atom sets collapse into one variant carrying a weight (the number of
global conformers producing it) and a representative label (the smallest
label among them). A fragment without altloc atoms has one variant of
weight `n_conformers`. The weights of the variants containing an atom sum
to the atom's conformer multiplicity by construction. The number of
variants of a fragment is at most the number of distinct labels among its
residue groups plus one, never more than `n_conformers`.

Residue groups that exist in some conformers only (Task 12, every route):
a residue group with no blank atoms (atoms whose only label is on hydrogens
count as blank), whose own conformer labels are a proper subset of the
model's and which is bonded to no other residue group, belongs only to the
conformers of its own labels; `S_c` of every other conformer omits it (no
fallback), so its atoms' multiplicity is the number of its own labels (a
water labelled B only, 1PWM HOH 1684/1704/1322: multiplicity 1, absent from
A; 1PWM has 10 such residue groups absent from A and 19 from B, CL 1321
among them, and no capping failure). Bonded is by distance
(`cctbx.crystal.distance_based_connectivity`: 1.3 x the tabulated bond
length, else 1.3 x the larger van der Waals radius; pairs tabulated at 0,
metal-O among them, are never bonded; the criterion of
`fragment.bonded_residue_pairs`)
between non-hydrogen atoms that can meet in one of the residue group's own
conformers (the partner's atoms by the fallback rule, or by its own labels
if it is absent itself; resolved by iteration), plus a peptide link to a
chain neighbour under any labels (`iotbx.pdb.altloc_completion.bonded`).
Bonds through hydrogens are ignored and bonds across symmetry are not seen
(the test is on the model alone). The rule is about existence, the bond
test only protects covalently linked groups: a B-only water on a blank
metal at 2.1 A is not a `distance_based_connectivity` bond and is absent
from A, which is physically right. Computed once per model in
`conformers._selections`, which reports `n_absent`. Residue groups with
blank atoms or a bonded partner keep the fallback and `check_links`.
Candidates bonded to each other by atoms of the same label are one entity
(connected components: a B-only two-residue ligand, glycan chain or
peptide), absent as a whole unless a member is peptide-linked or bonded
to a residue group outside the entity that the conformer contains.
The single-conformer views of the interaction and buffer searches lack the
absent residue groups (their residue indices are mapped back to the
model's), and a conformer giving a fragment no cluster atom gives it no
variant, so such a fragment's weights sum to less than `n_conformers`.

The fragment manager's job list becomes the list of variants:
`fragment_selections` (master atoms of the variant, on the model),
`fragment_super_selections` (on the full super sphere), `fragment_cappers`,
`fragment_charges`, `cluster_selections` and `buffer_selections` (those of
the fragment, repeated per variant), plus `variant_weights`,
`variant_labels` and `variant_fragment_index`. A single-conformer model
yields one variant per fragment with weight 1 and the same lists as before
the change, so nothing changes for it, including file layout.

### 3.3 Single-conformer view, capping, charges

A variant is evaluated as a genuine single-conformer fragment: its atoms
selected from the super sphere and relabeled to blank with the atom groups
of each residue group merged (resname rule as in `split`: the atom group
with the representative label, else the parent label, else the first).
Capper, charge computation and engines never see an altloc; this is
asserted after relabeling (no atom group with a non-blank altloc). The
relabeled hierarchy is built once at setup; per evaluation only coordinates
are pushed.

Capping: `fragment_capping.capper` finds cut bonds by distance against the
super cell. With the altloc super cell as parent it would see a B atom next
to an A fragment as a cut bond and cap it wrongly, so the capper of a
variant gets as parent the conformer's super cell: `SCM.super_cell_hierarchy`
selected by the tiled `S_c`, relabeled. One parent per conformer label,
shared by all variants of that label; the parents exist only while the
cappers are built (`get_fragment_hierarchies_and_charges`, once per
(re-)clustering) and are released after them (Task 14). A capper drops its
processed fragment model once the caps are placed: `cap()` needs the
dangling-atom records only.
Cut-bond detection, dangling atoms, cap placement and the capper's Sorry
conditions (unsatisfied valence that is no cut, unknown residue, cap
clash) are unchanged: no dangling bond ever reaches an engine. For a
single-conformer model the parent is the super cell itself, as today.

Charges: the charge of the capped variant comes from the capper, so every
variant has its own charge (conformers of a residue may differ in
protonation). Charge embedding (engines other than MOPAC) takes the point
charges from the conformer's super sphere (the super sphere selected by
`S_c`, relabeled) minus the variant, within `charge_cutoff`, so the
environment is single-conformer too.
Known limitation: a fragment whose variants merge into one job (label "")
takes its point charges from the conformer sphere of the smallest label,
so the field of a merged job is conformer A's even where the environment
within `charge_cutoff` differs between conformers.

### 3.4 Evaluation

`from_cluster.target_and_gradients` runs one `parallel_map` over all
variants of all fragments. Combination per atom i:

```
num[i] += w_v * g_v[i]   for every variant v containing i (master atoms)
den[i] += w_v
g[i] = num[i] / den[i]   on cluster atoms; 0 on buffer atoms
```

`den[i]` equals the conformer multiplicity of i, so this is the current
per-atom formula of `from_conformers` applied per fragment; for a
single-conformer model (weights 1, one variant per fragment) it is the
sum computed today. Energy: `sum_v w_v e_v / n_conformers`, the mean over
conformers of the sum of fragment energies, documented as undefined for
refinement (gradient_only is already forced with clustering). With absent
residue groups (3.2) the weights of a fragment sum to less than
`n_conformers`, so `den[i]` (not `n_conformers`) is the divisor of the
gradient, and `n_conformers` is the number of labels. The target equals the
mean of the conformer energies only without symmetry images in the sphere
(tst_65, one-cluster limit): with images the altloc super sphere includes
edge residues that only one conformer brings within the radius, so a
variant's buffer is a superset of the split model's sphere and buffer-buffer
terms enter the fragment energy (17 % in the contact-2 box, gradients agree
to 5.7e-14). The target is not a model energy.

QM (`restraints.from_qm.target_and_gradients`): receives a job index, looks
up the variant's super selection, relabeled hierarchy, capper, charge and
folder, pushes the current coordinates, caps, runs the engine, drops caps
and buffer atoms, returns the gradient over the variant's master atoms.
Folder and engine label: a fragment with one variant of full weight (all
conformers; always so for a single-conformer model) keeps `ase/<index>/` and
label `<index>` (single-conformer runs and tst_35's layout unchanged);
every other job, including a single variant of partial weight (a fragment
whose cluster exists in some conformers only), uses `ase/<index>_<label>/`
and label `<index>_<label>`, so MOPAC jobs never collide.
`fragment.job_folder_name` is the one place that decides.

Parallelism (Task 10): QM jobs run in one `parallel_map` over all variants;
cctbx fragment restraints are evaluated in process, whatever
`parallel.nproc` (cctbx is the test path; a process per job re-pickled the
super-sphere restraints at every evaluation). A cctbx exception propagates
as it is: no Sorry wrapper, no `keep_error_copy` folder.

cctbx (`restraints.from_cctbx`): the super-sphere restraints manager is
built once on the altloc super sphere by pdb_interpretation, which is
altloc-aware (no bonds or nonbonded terms between A and B); selecting it to
a variant gives exactly the single-conformer restraints of that conformer.
Cached per variant as the per-fragment cache does today; invalidated by new
selection lists after a re-clustering.

Interaction search, bonded residue pairs and buffer search run on
single-conformer views (Task 4, fix round 1), not on the altloc hierarchy:
`mmtbx.pair_interaction` assigns a grid point near an altloc residue to the
A and B copies of that residue as its two nearest atoms and discards the
same-residue pair, so interactions between residues that both carry altlocs
were lost (differences of 49-131 kcal/mol/A in the duplicate-conformer
identity, even in the one-cluster limit). There is one view per conformer
label (`conformers.single_conformer_hierarchy` of the model, and of the super
sphere with the tiled selections; absent residue groups are missing from a
view and its residue indices are mapped back to the model's); the
interaction lists, the bonded pairs and the per-cluster buffers are united
over the views (atoms mapped back to full-sphere positions through the
view's iselection), and a fragment is whole residues of any view. Cost:
`n_conformers` times the interaction-search part of setup (1PWM: 191
clusters and 373 jobs from the union at the time, 155 and 304 from the
altloc hierarchy; section 7); the views live during setup only
(`conformer_spheres` are kept for charge embedding). A shortcut running one
view where no altloc atom is near a cluster never applies on 1PWM. Not
done for `altlocs.handling=native` with clustering (section 6).

### 3.5 Parameters, driver, order of operations

No new PHIL. `altlocs.handling=conformers` (default) means: plain and
expansion routes split into whole conformer models (as today), the
clustering route uses per-fragment variants. `handling=native` unchanged.
`api.finalize_params` already forces lbfgs with gradient_only for
clustering and for conformers.

`driver` re-clusters through `fragment_manager.set_up_cluster_qm()`; on
this route the fragment manager is a plain `fragments` object on the altloc
model, so one call rebuilds clusters, buffers, variants, conformer parents
and cappers. `fragment_manager_proxy` became dead code and was removed.
`interaction_list` for the driver's log line is the single list.

`refine.run` order unchanged: exclude selection, then hd_mapper (the
de-deuterated single model), then `create_restraints_manager`. Conformer
selections are computed on `model_for_restraints`, so H/D exchange sites
are gone before conformers are counted, as now. `api.check_model` keeps
its completion check.

### 3.6 Errors and files

Sorry messages name fragment and variant: "fragment 12, conformer B: ...".
A failed QM job keeps today's behaviour (working folder copied to `error/`,
Sorry). `conformers.check_links` runs once on the model before
fragmentation.

`ase/expansion.pdb` is the altloc super sphere (one file). `save_clusters`
writes `<index>_cluster.pdb`, `<index>_frag.pdb` and one capped fragment
per variant, `<index>_<label>_capped0.pdb`; `fragment_info.txt` lists
variants and weights per fragment.

### 3.7 Performance

Setup: one expansion, one interaction list, one clustering, `n_conformers`
selections of the super cell for the capper parents, one capper per
variant (today: one per fragment per conformer, about the same count).
Evaluation: one `parallel_map` over all variants; the relabeled hierarchies
are prepared at setup. Speed is a requirement throughout: profile the
1PWM setup and one evaluation (cctbx; one MOPAC evaluation for the QM
path) and refactor what shows up, in particular the charge computation and
`from_qm.target_and_gradients`, both known to carry overhead, as long as
results are unchanged (section 7). Memory: no per-conformer copies of the
super sphere beyond the relabeled conformer super cells used as capper
parents; these exist only while the cappers are built (Task 14).

## 4. Verification

Meaning of "match", agreed 2026-11-07:

- (a) cctbx standard whole model vs (b) cctbx conformers whole model: not
  equal by design (altloc atoms identical; blank atoms obey
  m*g_conformers = g_standard + (m-1)*g_blank_only); unchanged here.
- (b) vs (c) conformers clustered: identical to 1e-6 in the one-cluster
  limit; with real clusters within the frozen-buffer tolerance clustering
  vs expansion has today for single-conformer models (tst_04, tst_29).
- (a) clustered standard vs (c): same clustering, so the relation holds
  exactly per fragment and is tested per atom.
- MOPAC: (b) vs (c) in the one-cluster limit to engine noise; variant
  gradients by frozen-environment finite differences; duplicate-conformer
  identity.

Test model for MOPAC: a super-toy peptide (tri- or tetra-ALA, A/B on one or
two residues built in the test by copying and shifting a residue as tst_41
and altloc_scenarios do), `maxnum_residues_in_cluster=1`, in a P1 box large
enough that no image is within the nonbonded cutoff. MOPAC jobs of 20-40
atoms take about a second; each test runs in about a minute.

### 4.1 Unit tests (tests/unit, registered in run_tests.py)

1. Variants: hand-made model with residues carrying {A,B}, {B,C} and none,
   global labels A, B, C: expected variant atom sets, weights and labels per
   fragment; tiled super-sphere selections equal the per-residue rule
   applied directly to the sphere hierarchy; a single-conformer model gives
   one variant of weight 1 per fragment and the same selection lists as
   before the change.
2. One-cluster limit: `maxnum_residues_in_cluster=999`, new route vs
   `from_conformers` over expansion: 1e-6 (cctbx) on the super-toy peptide
   and h_altconf_complete; engine noise (MOPAC) on the super-toy peptide.
3. Clustered standard vs clustered conformers (cctbx, same clustering):
   altloc atoms identical, blank atoms satisfy the relation per atom.
4. Duplicate-conformer identity through the clustered route: a
   single-conformer peptide copied into identical A and B (and A, B, C with
   altered occupancies) gives exactly the single-conformer clustered
   gradients, cctbx and MOPAC; job counts show the dedupe (one variant per
   fragment where only one residue is duplicated and its fragment does not
   reach, two where it does).
5. Frozen-environment finite differences per variant (tst_46 pattern) on
   the super-toy peptide, cctbx and MOPAC.
6. Clustering vs expansion with real altlocs on the altloc data files,
   cctbx, within tst_29's tolerance.
7. Re-clustering: `re_calculate_rmsd_tolerance=0` on an altloc model forces
   a rebuild every cycle; the run completes and the variant lists are
   rebuilt.
8. qr.refine end to end on the tst_06 toy data, cctbx: (a) standard, (b)
   conformers whole model, (c) conformers clustered; R-work drops in all
   three, (b) vs (c) within a small rmsd, (a) within tst_43's rmsd; folder
   layout `ase/<index>_<label>` checked. MOPAC: (b) vs (c) optimization of
   the super-toy peptide for a few iterations.
9. phenix_regression tst_177 extended: an altloc model with
   `qrefine.clustering=True`; phenix.refine equal to standalone qr.refine
   (1e-4; the case uses a float32-exact cell, 18.375 20.875 31.375: with the
   tst_06 cell, which MTZ storage rounds, the two programs adopt the data's
   cell through different conversions, start 1-4 ULP apart and the
   gradient-only minimizer amplifies that to 0.01 A; Task 9b).
10. 6uf7 as a standard suite test (real data, 95 atoms, A/B altlocs,
   macrocycle closed through symmetry, refinements under a minute; an ultra
   high resolution niche case covered for free): the reduce2 qmready model
   and the data go into tests/unit/data_files; qr.refine (a), (b), (c) with
   cctbx and (b), (c) with MOPAC (opt mode, few iterations), compared as in
   test 8. As built (tst_69): cctbx (a), (b), (c) only. The MOPAC exercise
   was first gated by an environment variable and removed on 2026-11-12
   (user ruling: a gate for runtime is not worth keeping; one evaluation is
   10 jobs of 0.3-5 s and an optimization needs many, so it never fits the
   suite); MOPAC with altlocs and clustering is tst_67, and 6uf7 with MOPAC
   is a manual acceptance run. Route (b) with expansion stops in the capper
   (N PRO 1 closed through a symmetry mate, an old defect of the expansion
   route).

As built: tst_64 (variants, weights, labels, capper parents, job folder
names, charge embedding, absent residue groups, interactions of altloc
residues), tst_65 (weighted combination, one-cluster limit, relation per
atom, clustering vs expansion), tst_66 (duplicate-conformer identity), tst_67
(MOPAC super toy), tst_68 (qr.refine on the toy data, re-clustering, H/D),
tst_69 (6uf7); tst_31 and tst_49-52 cover the conformer rules. Times are
between 11 s (tst_66) and about 80 s.

Existing tests: tst_43 and tst_46 run on the new route (assertions stand);
tst_31 unchanged.

### 4.2 Acceptance runs (manual, outside the suite)

- 6uf7 in refine mode with MOPAC (the suite test keeps MOPAC in opt mode
  for runtime).
- 1PWM with cctbx clustering, one cycle: 304 jobs (373 as built, section 7),
  setup once, evaluation
  time before and after, final R-work against the per-conformer run. One
  MOPAC gradient evaluation on 1PWM fragments to confirm the variant jobs
  run and to profile charges, capping and `from_qm` overhead.

## 5. Also in this project

The AQuaRef documentation sentence "Static disorder (i.e., alternative
conformations) is not supported" (phenix_html
rst_files/reference/AQuaRef.txt) is replaced by the current state
(alternative conformations supported on every route once the model is
completed by `mmtbx.reduce2 qmready=true`; exchangeable H/D sites handled
for neutron data; unsupported entities such as metal ions still need
`exclude`). Done in Task 11, with CHANGES.md and the tests table
(doc/qr/tests.html).

## 6. Out of scope

Ribosome-scale setup (super-sphere memory, pair_interaction and
Girvan-Newman on 100k+ atoms) is a separate project; this design must not
make it harder (section 3.7). MOPAC's cost and SCF failures above about 700
atoms per fragment are an engine property (MOPAC is the test engine only)
handled by `maxnum_residues_in_cluster`.

Known limitations of this design:

- `altlocs.handling=native` with clustering (cctbx only) still runs the
  interaction search on the altloc hierarchy and can miss interactions
  between altloc residues (3.4).
- The one-cluster target is the mean of the conformer energies only without
  symmetry images; it is not a model energy (3.4).
- Bonds across symmetry are not seen by the test that a residue group exists
  in some conformers only (3.2).
- Fragment weight sums can be below `n_conformers` (3.2).

## 7. Results (2026-11-07, Tasks 10-14 to 2026-11-12)

Measured on the 16-core development Mac, 1PWM = rt3/1PWMH_charge3.pdb
(8452 atoms, labels A/B), default `maxnum_residues_in_cluster=15`.

**Jobs.** Task 10: 191 fragments, 373 jobs (182 fragments with two
variants, 9 with one); after Task 12 (section 3.2): 190 fragments, 373 jobs.
The 304 expected in section 4.2 assumed the clustering of the native altloc
model (155 fragments); the union of the per-conformer interaction lists
(section 3.4) has more edges, and Girvan-Newman gives 191 clusters.

**cctbx setup** (`refine.create_restraints_manager`, cProfile, nproc=1):
197.6 s (190 s of it in `fragments`). Top items: `get_fragments` 119 s
(`pair_interaction.get_interactions`, C++, 108 s in 766 calls: 191
clusters x 2 conformer views x 2 passes), `get_clusters` 42 s
(Girvan-Newman 36.5 s; the two per-view `pair_interaction.run` 2.6 s and
bonded pairs 3.0 s), expansion 24.6 s (`super_cell.get_siiu` 23 s),
restraints of the model 7.8 s, variant hierarchies 3.8 s,
`_sphere_residue_atoms` 0.8 s. Nothing in qrefine's own code reaches 10 %:
the per-view searches are the cost of the conformer views. An exact
shortcut (run one view when no altloc atom is within the 15 Bohr reach of
`run_prepared` around the cluster) was implemented and measured: it never
applies on 1PWM (every cluster is near an altloc residue) and was removed.
Unchanged by this task: 191-198 s.

**cctbx evaluation.** In process: first call 20 s (super-sphere restraints
and 373 fragment selections, cached), then 0.17 s. With `parallel.nproc=4`
the jobs went through `parallel_map`, one process per job and no cache:
21.6 s per evaluation. Fixed: cctbx fragment restraints are evaluated in
process whatever nproc (`cluster_restraints.from_cluster`); results
bit-identical to the `parallel_map` path at the same coordinates and equal
to the nproc=1 path.

**QM setup** (MOPAC, capper and charge per variant): 2.35 s per capper
measured (350 ok and 23 failed in 875.8 s; 2.71 s per capper in the 12-capper
profile under cProfile). `capping.cut_bonds_by_distance` on the whole
conformer super cell (143684 atoms) took 1.70 s of it; now run on the
super-cell atoms in the fragment's box widened by the bond cutoff
(`fragment_capping.near_selection`): 0.44 s, cut bonds identical on 75 of
75 compared variants. The default pdb_interpretation scope is parsed once
(`hierarchy_utils.get_pdb_interpretation_scope`, 0.14 s per model). Now
about 1.4 s per capper (373 cappers 414-428 s, measured under load); the
rest is cctbx: `electrons.run` 0.46 s (charge), `cut_bonds_by_distance`
(`crystal.pairs_within` Python loop) 0.44 s, model processing 0.25 s.
At the time of Task 10 the full 1PWM QM setup stopped with a Sorry: 23 of
373 variants failed capping (Open, below); since Task 12 the setup completes
with no capping failure.

**MOPAC jobs** (`from_qm.__call__`, PM7 1SCF gradients with COSMO,
THREADS=1): 1PWM 28-atom variants (30 capped): 0.20 s, engine 0.19 s,
qrefine overhead 8-10 ms; 151-atom variant (161 capped): 0.86 s, overhead
24 ms (3 %). 6uf7 (`maxnum_residues_in_cluster=2`, 10 jobs, setup 18 s):
111-atom variant 0.31 s (overhead 13 ms), 400-atom variant 5.05 s
(overhead 33 ms, 0.7 %). Overhead (capping, PDB writes, ASE conversion)
is far below 10 %: `from_qm.target_and_gradients` unchanged. The 6uf7
refinement was skipped on Task 8's figure of ~5 min per job, which was
wrong: the largest job takes 5 s. An evaluation is 10 jobs and an
optimization needs many evaluations, so tst_69 has no MOPAC exercise
(removed 2026-11-12) and the 6uf7 MOPAC acceptance is manual.

**1PWM acceptance**, `qr.refine ... restraints=cctbx clustering=True
number_of_weight_search_cycles=1 number_of_refine_cycles=1
parallel.nproc=4`: before the evaluation fix 2344 s wall (start of
refinement at 269 s; weight search 54 evaluations 990 s, refinement 54
evaluations 1018 s), R-work 0.1342 -> 0.1269, R-free 0.1390 -> 0.1350.
After: 536 s wall (start at 269 s; 109 s and 91 s for the two cycles of 54
evaluations), R-work 0.1342 -> 0.1269, R-free 0.1390 -> 0.1349. R-free moved
by 1e-4 because the in-process cache keeps each fragment's updated pair
proxies over the refinement (equal to the nproc=1 path), where the old
per-process path reselected them from the super-sphere proxies at every
call. The per-conformer run of 2026-11-05: 0.1307 -> 0.1259. Peak RSS of
this run: 2.58 GB before, 3.05 GB after (the 373 cached fragment restraints
in the main process; reduced by Task 14, see Memory).

**Resolved by Task 12.** 23 1PWM variants failed capping: waters with only
label B (HOH A1684, A1704, A1322) joined conformer A by the fallback rule
and sat 0.5-1.6 A from conformer-A atoms (GLU A271 OE2, ASN A129, GLN A93,
CL A1321); `cut_bonds_by_distance` took the contact for a cut bond and the
cap clashed. The rule of section 3.2 (such residue groups are absent from
the other conformers) removes them: 0 capping failures, 373 jobs.

**Open.** cctbx setup is dominated by
`pair_interaction` (C++, two views), Girvan-Newman and the super cell; the
capper by `electrons.run` and `crystal.pairs_within`: cctbx work.

### Memory (Task 14, 2026-11-12)

1PWM (8452 atoms, `select_within_radius=10`): 16 symmetry operators, super
cell 143684 atoms (17 x 8452), super sphere 23030 atoms, 190 fragments, 373
jobs, 105222 atoms in the variant hierarchies (12.4 x the model). Measured
in one process per route (scratchpad `t14scripts/mem.py`), stage by stage:
heap in use (`mstats` bytes_used: live allocations, what a Linux RSS
roughly follows), the macOS physical footprint (current and lifetime
maximum) and RSS. RSS is not a good measure on macOS: freed pages stay in
RSS until the system needs them (after the super cell is built RSS stays at
1.8-2.0 GB while 250 MB are live), which is why the peak RSS of a run (2.0-2.3
GB here, 2.6-3.1 GB for the Task 10 qr.refine runs) overstates what is held.

| stage (heap MB added) | cctbx before | cctbx after | QM before | QM after |
|---|---|---|---|---|
| model processed (refine's, not this project's) | 71 | 71 | 71 | 71 |
| `restraints.restraints` (model copy; cctbx: + whole-model restraints) | 127 | 64 | 63 | 0 |
| super cell manager (kept: cell 41, its selection cache 18, sphere 7) | 59 | 59 | 65 | 65 |
| `set_up_conformers` (QM before: + the two conformer cells, 81) | 12 | 12 | 81 | 12 |
| fragments: prepare views (released after the search) | 14 | 14 | 14 | 14 |
| variants, selections (cctbx); + cappers (QM) | 49 | 49 | 674 | 101 |
| super-sphere restraints manager (cctbx, first evaluation) | 138 | 44 | - | - |
| 373 cached fragment restraints (cctbx, first evaluation) | 441 | 189 | - | - |
| heap in use after setup | 362 | 300 | 999 | 294 |
| heap in use after the first evaluation | 942 | 534 | - | - |
| footprint after setup / after the first evaluation | 417 / 581 | 385 / 552 | 1394 | 391 |
| peak footprint (always during the super cell build) | 1687 | 1666 | 1655 | 1603 |

The peak of every clustered run is the super cell build, not anything held
afterwards: `cctbx.crystal.super_cell.get_siiu` builds the full pair table
of the model within `select_within_radius + 1` = 11 A (5692852 (i, j,
operator) entries for 8452 atoms, 729 MB live, about 1.5 GB footprint while
the asu table and the extracted table coexist) only to find 16 distinct
operators. It is freed at once (the manager keeps 59 MB).

Reductions (qrefine only; cctbx gradients bit-identical, the 373 QM jobs
identical in charge, caps, capped coordinates and atom ids):

- `fragment_capping.capper` no longer keeps its processed fragment model
  (`_model`, about 1.5 MB per job with its restraints manager) once the
  caps are placed: `cap()` needs the dangling-atom records only. QM
  -551 MB on 1PWM.
- the relabeled conformer super cells (capper parents, 2 x ~117000 atoms)
  are built just before the cappers and released after them, instead of
  living from `set_up_conformers` to the next re-clustering. QM -68 MB.
- `fragment_capped_initial` (a capped hierarchy per job, read by nothing) is
  gone, and with it one `cap()` per job at setup.
- `restraints.restraints` copies the model only for the expansion route,
  the only user of the copy. -63 MB on every route.
- the super-sphere restraints manager is built without `plain_pairs_radius`:
  the plain pair table serves ADP restraints only and was the largest part
  of every selected fragment manager. cctbx -94 MB (sphere) and -252 MB
  (fragment cache); the 373 fragment selections of the first evaluation
  5.8 s -> 1.7 s (Task 10's 20 s first call also built the super-sphere
  restraints).

Full `qr.refine` (cctbx, clustering, one weight-search and one refinement
cycle, nproc=4, `/usr/bin/time -l`): max RSS 2786 -> 2227 MB, peak
footprint 1834 -> 1828 MB (the super cell build), R-work/R-free 0.1347/
0.1403 -> 0.1269/0.1350 (the start state of that code, not Task 10's
0.1342/0.1390) and the refined coordinates identical in both.

Not changed, to keep cctbx bit-identical: the whole-model cctbx restraints
manager that `restraints.restraints` builds on the clustering route is
unused there (64 MB), but building it runs pdb_interpretation on the
caller's hierarchy, which rewrites the coordinates of 3643 of its atoms by
up to 1.5e-14 A (fractional round trip); skipping it changes the 1PWM
gradients by 5e-10. Dropping it is a one-line change once that is accepted.

**What scales with model size (ribosome project).** Per model atom at
`select_within_radius` R (1PWM: R = 10, about 0.1 atoms per cubic A with H):
(1) the operator search: about 670 pair-table entries per atom at R + 1 = 11
A, growing with (R + 1)^3, about 180 kB of transient footprint per atom: some
50 GB for a 300000-atom ribosome at R = 10, 20 GB at R = 7. This is the first
wall and is in cctbx. (2) The super cell: (number of operators + 1) x the
whole model, about 290 B per atom plus 125 B for its atom selection cache,
kept for the whole run because `super_cell.manager.update()` moves the sphere
through the cell; with 10-16 operators a ribosome cell holds 3-5 million
atoms, 1-2 GB, plus the same again transiently for each relabeled conformer
cell (QM). (3) The super sphere: the model plus every symmetry atom within R,
2.7 x the model for 1PWM's small cell, nearer 1.3-1.5 x for a large one; the
variant hierarchies and each capper's fragment hierarchy hold the fragment
atoms (12.4 x the model on 1PWM: ~30 MB each here, ~1 GB each for a ribosome).
(4) Jobs: the selections are flex.bool arrays over the whole sphere or model,
per job (`fragment_super_selections`, `fragment_selections`) and per fragment
(cluster, buffer, full super selection): 19 MB here, but quadratic in model
size: a ribosome with ~13000 jobs and ~6700 fragments over a 400000-atom
sphere needs ~15 GB. Index selections (8 bytes per selected atom) would make
this linear; this is the largest qrefine-side item for the ribosome. Each job
also holds its capper records. (5) cctbx only: the super-sphere restraints
(426939 nonbonded pairs for 23030 atoms, now ~1.9 kB per sphere atom) and the
cached fragment restraints (~1.8 kB per fragment atom, 12.4 x the model: ~7 GB
for a ribosome); selecting a fragment costs 0.7 ms in process, so a ribosome
cctbx run should select per evaluation instead of caching (results then equal
those of the former nproc>1 path, not bit-identical to the cache). QM runs
hold neither.

**Proposed for cctbx** (not made): (a) `get_siiu` for the manager needs the
distinct operators only: stream the pairs (`neighbors_fast_pair_generator`
over the asu mappings, memory O(atoms)) and keep the first-seen order of
the operators so the super cell keeps its atom order; a cheaper but not
exact alternative searches with the non-hydrogen atoms and a radius widened
by the X-H length (half the atoms, a superset of the operators). (b) let
`super_cell.manager` keep, per sphere atom, its operator and master index
and update the sphere directly, so the super cell and its selection cache
can be dropped after construction (`update()` today rebuilds the sphere by
selection from the whole cell at every evaluation); qrefine would then cut
capper parents from the cell atoms within the bond cutoff of the sphere
(Task 10's `near_selection` margin), which the manager could return as a
"sphere plus margin" hierarchy instead of the cell.

### Three routes side by side (2026-11-12, 1PWM, cctbx, one CPU)

Identical settings: `restraints=cctbx number_of_weight_search_cycles=1
number_of_refine_cycles=1 parallel.nproc=1`, run one after another on an
idle machine, wall time and peak RSS from `/usr/bin/time -l`.

| route | total wall | to first monitor line | one macro cycle | peak RSS | R-work/R-free end |
|---|---|---|---|---|---|
| native (clustering=False, handling=native) | 286 s | 72 s | 76 s | 1.31 GB | 0.1275/0.1351 |
| conformers (clustering=False, handling=conformers) | 274 s | 74 s | 75 s | 1.30 GB | 0.1269/0.1349 |
| clustering (clustering=True, handling=conformers) | 474 s | 246 s | 98 s | 2.72 GB | 0.1269/0.1350 |

All three start at R-work/R-free 0.1347/0.1403 and take 52-54 function
evaluations. The per-conformer whole-model route costs the same as native
(cctbx geometry is cheap next to the X-ray target) and gives the clustered
answer to 1e-4 in R-free; with cctbx restraints clustering buys nothing on
a model of this size, its purpose is the QM engines. Its extra cost is
174 s of setup (pair_interaction on two views, Girvan-Newman, super cell)
plus ~20 s on the first evaluation (super-sphere restraints), then 23 s per
cycle. Peak RSS 1.3 GB on the whole-model routes is the X-ray machinery; of
the clustered route's 2.7 GB, 1.5 GB is the transient super-cell operator
search (a standalone cctbx-only script, `super_cell_peak.py`, kept outside
the repository in the session scratchpad and handed to the user on
2026-11-12, reproduces it: pair table of 5.7 million entries, 743 MB live,
1.50 GB peak, to find 16 operators; the finished manager holds 83 MB).

### Known limitations and follow-ups

- The peak memory of every clustered or expansion run is cctbx's
  `super_cell` operator search (~1.6-1.8 GB on 1PWM, proportional to
  atoms * (R + 1)^3); proposal (a) above: stream the pairs.
- The whole-model cctbx restraints that `restraints.restraints` builds on the
  clustering route and nothing reads (64 MB on 1PWM) are kept for bit
  identity (building them moves the caller's coordinates by up to 1.5e-14 A);
  dropping them is a decision for the user.
- Per-job `flex.bool` selections over the sphere or model are quadratic in
  model size; index arrays (`flex.size_t`) would make them linear.
- The cctbx fragment-restraints cache (373 managers, 189 MB on 1PWM, ~7 GB
  for a ribosome) could be replaced by selecting per evaluation (0.7 ms per
  fragment), at the price of results that equal the former nproc>1 path
  instead of the cache.
- Adopting the unit cell of the data the way phenix.refine does (Task 9b: the
  two programs convert it differently, which the gradient-only minimizer
  amplifies to ~0.01 A on a cell that is not exact in float32) is a
  follow-up, its own project if large.
- `altlocs.handling=native` with clustering keeps the lost-interaction
  defect of `pair_interaction` on altloc hierarchies (section 6).
- Setup speed (2026-11-12, 1PWM cctbx one CPU: 174 s from model to the first
  monitor line, against 72 s without clustering; the items above are memory
  items and would recover ~20 s of it). The two stages that hold the time
  are not on the memory list: (1) the `mmtbx.pair_interaction` grid search,
  108 s of 198 s in the Task 10 profile, run in full on every
  single-conformer view because each view's grid origin is the minimum
  coordinate of its own sub-hierarchy, so a second view restricted to the
  altloc neighbourhoods is not exact (an exact skip was tried in Task 10 and
  never applied on 1PWM). Follow-up, in cctbx: let the C++ grid search take
  a fixed grid origin and spacing from the caller; the first view then
  provides the grid, the other views search only the residues within the
  interaction cutoff of an altloc residue group, and the union is identical
  by construction (seconds instead of ~57 s per extra view). Running the
  views on separate CPUs is the stop-gap. (2) Girvan-Newman clustering in
  `clustering.py`, 37 s of pure-Python betweenness on the residue
  interaction graph. Follow-up: a C++ betweenness or a cheaper community
  algorithm (seconds); it changes nothing physical but may move cluster
  boundaries, so it must pass the frozen-buffer checks of tst_65/tst_68.
  Setup is paid once per clustering and re-clustering is triggered only by
  a coordinate shift beyond the tolerance, so over a multi-cycle refinement
  the fixed cost amortizes; the per-cycle difference is 98 s against 75 s
  (section 7, three routes side by side).
