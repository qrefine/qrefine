# qr.refine runtime bottlenecks (restraints=cctbx, 1PWM) and plan

Status: analysis done 2026-11-05; decisions taken 2026-11-05 (section 4);
implementation plan to be written from section 5.

## 1. The run

    qr.refine 1PWM.mtz 1PWMH_charge3.pdb restraints=cctbx

8452 atoms (4609 H), two conformers (6959 atoms each), 193,485 reflections,
d_min 0.92 A. Wall clock 176 min (log in rt3/). Defaults that matter:
cluster.clustering=True, number_of_weight_search_cycles=50,
max_iterations_weight=50, number_of_refine_cycles=5.

Context (user): cctbx restraints in qr.refine are for testing only; the
production path is QM with clustering (only AQuaRef runs without it). So the
fragmentation setup and the number of gradient evaluations (weight-search
trials) are what matter in production; the per-evaluation cctbx cost matters
for the tests.

Timeline reconstructed from file mtimes in rt3/ (no timestamps in the log):

| phase | wall | evidence |
|---|---|---|
| start to conformer A expansion written | 1 min | ase/A/expansion.pdb 10:02 |
| conformer A fragmentation | 9.5 min | ase/B/expansion.pdb 10:11 |
| conformer B fragmentation + super-sphere restraints + fit_rotatable2 + first stats | 16 min | ase/*/tmp_ase.pdb first written 10:27 |
| weight search, 19 cycles, ~1007 gradient evaluations | 135 min | pdb/0_refine_cycle.pdb 12:42 |
| 5 refinement cycles, 156 evaluations | 15 min | pdb/4_refine_cycle.pdb 12:57 |

## 2. Profiles

Two profiled runs (cProfile, scratchpad prof_cluster/ and prof_nocluster/),
1 weight cycle x 5 iterations + 1 refine cycle x 5 iterations, i.e. 16
gradient evaluations. The machine was loaded (load average ~7), so absolute
times are inflated; the proportions and call counts are what the plan rests
on. Final before/after timing must be done on a quiet machine.

### 2a. Default (clustering=True): 33 min for 16 evaluations

| item | time | notes |
|---|---|---|
| create_restraints_manager (both conformers) | 1155 s | all of it fragmentation |
| - clustering.get_clusters (betweenness centrality, boost graph) | 568 s | 284 s per conformer |
| - fragment.get_fragments: pair_interaction.run per cluster, 345 calls | 533 s | 1.5 s per cluster on the 26k-atom super sphere |
|   - of which Python per-atom loops in run() (6.4M set.add/strip/lower) | 374 s | atom_in_residue/element bookkeeping rebuilt for the same super sphere 172 times |
|   - of which load_wfc, 2776 calls | 60 s | the 8 wfc pickles unpickled on every call |
|   - of which C++ points_and_pairs | 94 s | |
| - expansion (cctbx super_cell, get_siiu) | 33 s | |
| - super-sphere GRM build (pdb_interpretation on 26k atoms), 2 | 18 s | |
| fit_rotatable2 | 443 s | see 2b |
| per gradient evaluation | 12.4 s | restraints 10.1 s (5.5 s per conformer: 172 fragments x 27 ms, each a fresh grm.select + pair proxies + nonbonded asu mapping, update_plain_pair_sym_table 19 ms of it; plus tmp_ase.pdb write 0.18 s) + X-ray term 1.8 s |
| monitor.show (full MolProbity), 6 calls | 83 s | 14 s each |
| update_all_scales, 6 calls | 95 s | 16 s each |

The super-sphere file ase/<X>/tmp_ase.pdb (2 MB) is rewritten on every
evaluation; the super-sphere GRM is built lazily on the first evaluation.
Re-clustering (rmsd > 0.5 A) never happened in the original run.

Clustering graph of conformer A (scratchpad bcc/measure.py): 912 residues
(593 waters), 4264 interaction edges (water-water 1106, water-protein 1597,
protein-protein 1561). One boost betweenness_centrality_clustering pass at
threshold 9 took 426 s (loaded machine) and left 138 components (largest 17);
the threshold-8 continuation took 0.3 s. The boost algorithm (Girvan-Newman)
recomputes edge betweenness on the whole graph after every single edge
removal, O(V E) each; here nearly all 4264 edges are removed (unnormalized
betweenness > 9 is exceeded by almost every edge of a connected 912-vertex
graph), hence ~4000 x 4M operations.

### 2b. clustering=False: 9 min for 16 evaluations

| item | time | notes |
|---|---|---|
| prepare_inputs | 43 s | model processed 5 times in total over the run (data manager, prepare_inputs, _model_copy, xh_mapper, ...) 28 s |
| fit_rotatable2 (mmtbx.hydrogens) | 369 s | 249 residues; per residue: fmodel.deep_copy 0.35 s + full FFT f_calc 0.46 s + fft_map 0.25 s + map coeffs/scaling 0.13 s = ~1.3 s |
| per gradient evaluation (calculator.update) | 1.3 s | f_calc 0.46 s + X-ray gradients 0.69 s + xh_mapper.gradients 0.076 s (Python loop over 4609 H) + restraints of both conformers 0.01 s |
| monitor.show, 6 calls | 54 s | 9 s each: clashscore 2.7, rotamers 1.8, rama_z 1.5, cablam 0.8, rama, omega, cbeta |
| update_all_scales, 6 calls | 89 s | outlier rejection 70 s of it |
| create_restraints_manager | 9 s | |

## 3. Root causes, in order of impact on the 176 min

1. Clustering with cctbx restraints: 19 min fragmentation setup and ~10 s of
   the ~12 s per evaluation (1160 evaluations). The setup cost is the same in
   QM runs (production); the per-evaluation cost is the cctbx test path only.
2. Weight search never stops. GOOD tests bond rmsd only (the angle test is
   commented out); bond rmsd started at 0.015 < max_bond_rmsd 0.03 and crept
   to 0.029 without ever crossing it, so no flip happened; the only other stop
   is "two consecutive identical rounded rmsds", which took 19 cycles. R-free
   stopped improving at cycle 5 (0.1347) and cycles 6-18 (13 x 53
   evaluations, ~95 min) changed nothing. stop_one_found_first_good_weight is a
   documented parameter the code ignores (commented out). In a QM run every
   trial is 53 QM evaluations of all fragments: this is the main lever there.
3. Full MolProbity statistics twice per weight cycle: the "(NEXT)" line
   recomputes the identical statistics. 46 show() calls in the run, 9-14 s
   each: 7-10 min.
4. fit_rotatable2 recomputes the entire structure factor set and copies the
   fmodel once per residue (249 x): 6-7 min. cctbx code
   (mmtbx/hydrogens/__init__.py), called by qrefine when d_min < 1.2 and by
   phenix.fit_h.
5. Inherent per-evaluation X-ray cost, ~1.1-1.3 s on the loaded machine (FFT
   f_calc + FFT gradients at 193k reflections). Not avoidable.

Side findings:
- wfc files. mmtbx/pair_interaction/<el>_wfc_obj.pkl (H, C, N, O, S, D, Zn)
  are pickled mmtbx_pair_interaction_ext.wfc objects added to cctbx on
  2019-04-21 (commit b8fd0abbb6 "Move wv to c++; add pkl for wfc"). The text
  tables for all elements are qrefine/plugin/yoink/dat/<el>_lda.wfc. load_wfc
  looks for <el>*.pkl in the mmtbx folder; if absent it parses the text table
  into a wfc object (0.02 s) and pickles it into the current directory under
  the text file's name (<el>_lda.wfc); tst_05 asserts that write. Nothing ever
  reads that file: the loader never looks in the cwd and in the mmtbx folder
  only at .pkl names. So the on-the-fly cache was never closed, and the three
  files moved into mmtbx/pair_interaction/ (cl_lda.wfc, f__lda.wfc,
  p__lda.wfc; verified to be wfc pickles, e.g. Cl: ngrid 5776, Z 17) are
  ignored and regenerated into the cwd on every run.
- xh_mapper shrink/gradients are Python loops over all H: ~0.13 s per
  evaluation, ~2.5 min over the run.
- The model is processed (pdb_interpretation) 5 times before refinement
  starts (~30 s).

## 4. Decisions (user, 2026-11-05)

- cctbx restraints are a test path; clustering stays the default and its
  speed is what matters (QM production path). No clustering=Auto.
- Weight search: the two goals are lowest R-free under max_bond_rmsd. Do not
  restore the angle test (it is an obstacle in other cases). Sample the weight
  more efficiently and detect the R-free plateau: 0.1369 vs 0.1365 is the
  same.
- Statistics: keep full MolProbity on the lines if affordable (clashes and
  rama Z are useful); the duplicate computation goes.
- fit_rotatable2: one OMIT map, once, for the whole model.
- Waters stay in the clustering (no water-free clustering).

## 5. Plan

Expected result for the 1PWM cctbx run: setup from 26 min to about 2 min,
weight search from 19 trials to 3-5, per cctbx evaluation on the cluster
route from 12 s to about 2 s (X-ray term plus cached fragment restraints).
That is roughly 15-20 min on the loaded machine for the same refinement
cycles, less on a quiet one; in QM runs the gain is the trial count and the
setup.

A. Weight search (qrefine driver.py, medium). Replace the loop:
   1. Trial 0 at rws=1 (gradient-norm data weight as now); GOOD = bond rmsd
      <= max_bond_rmsd (unchanged, no angle test).
   2. Direction from trial 0: GOOD -> scale rws down, else up, by factor 2
      (adjust_restraints_weight_scale_value, default 2 already; the code
      hard-codes 1.5 and ignores it).
   3. Stop at the first flip (GOOD <-> not GOOD): the bracket is found; one
      bisection trial at the geometric mean of the bracket, then stop.
   4. Plateau stop: while GOOD, if r_free is within refine.r_tolerance
      (0.001) of the best GOOD r_free so far for 2 consecutive trials, stop.
      On the 1PWM log this stops after trial 2 (0.1369, 0.1360, 0.1356);
      with factor 2 it needs 3-4 trials.
   5. Choice: among GOOD trials, the largest rws whose r_free is within
      r_tolerance of the minimum (equal R-free -> better geometry), instead
      of the raw minimum.
   6. number_of_weight_search_cycles default 50 -> 12 (hard cap; the log
      warns when the cap, not a rule, stopped the search); drop the dead
      stop_one_found_first_good_weight parameter; keep
      max_iterations_weight=50 (unchanged; make it a knob to lower for QM).
   8. (added 2026-11-06 after the 6uf7 test, where R-free kept falling at
      every halving from a start weight 100x too large and the cap of 8
      stopped the search) the step accelerates: factor for the first move,
      then factor^2, factor^3, ... at every trial that keeps the direction
      and makes progress (going up: every not-GOOD trial; going down: every
      GOOD trial improving R-free by more than r_tolerance); after the flip
      the bracket is bisected until its ratio is at most sqrt(factor).
      6uf7: 7 trials (1, 1/2, 1/8, 1/64, 1/1024, 1/32768, 1/2^20), chosen
      1/1024 at the plateau, R-free 0.0745 (the cap-limited run: 8 trials,
      0.0753 and still falling).
   7. Log: timestamp (elapsed seconds) on every cycle line.
   Tests: the existing weight-search tests plus a synthetic run on 1PWM
   (cctbx) checking trial count and chosen rws.

B. Statistics (qrefine results.py, small). Cache the geometry statistics in
   the monitor, invalidated by update(); show() reuses it, so the "(NEXT)"
   line costs nothing. Full MolProbity stays on every line. Saves half of
   item 3 (the rest goes away with fewer trials).

C. fit_rotatable2 (cctbx mmtbx/hydrogens/__init__.py, small). Collect the
   rotatable H selections of all residues first (shortcut per residue, as
   now), build one OMIT map with all of them omitted (one f_calc, one FFT),
   then fit every residue in that map. Per-residue fmodel copies and f_calc
   disappear; 249 x 1.3 s -> a few seconds. Also used by phenix.fit_h.
   Verify on 1PWM: fitted H positions against the current code (expect small
   differences where neighbouring rotatable H were in the per-residue maps).

D. Fragmentation setup (qrefine + cctbx mmtbx/pair_interaction +
   boost_adaptbx, medium; results identical, snapshot harness of the
   2026-10-31 review):
   D1. Girvan-Newman with per-component recomputation, in C++ (new function
       in boost_adaptbx/graph/clustering_algorithm_ext.cpp or in
       mmtbx_pair_interaction_ext): after removing the max-betweenness edge,
       recompute betweenness only in the component that lost the edge;
       components whose max betweenness is already <= threshold are final.
       Same removal sequence as boost (global max across components; keep
       boost's first-max tie-break by iterating edges in the same order).
       Expected: 284-426 s -> seconds (the only expensive steps are the
       first few hundred removals on the 912-vertex component). The
       threshold loop 9..4 in clustering.py stays.
   D2. pair_interaction.run(ph, core): compute atom_in_residue, element
       types, Bohr coordinates and the per-residue atom lists once per
       hierarchy (the super sphere is the same for all 172 clusters) and
       pass them in; per cluster only the 15-Bohr selection and the C++
       call remain. ~370 s of 530 s. The selection loop over core atoms can
       be one C++/flex call too.
   D3. load_wfc: module-level cache of wfc objects (60 s), and the file
       convention: read <el>_wfc_obj.pkl from mmtbx/pair_interaction/; when
       absent, parse the text table and write <el>_wfc_obj.pkl there (cwd
       only if the folder is not writable). Rename the three moved files
       accordingly; fix tst_05 to the new convention. Add Cl, F, P pickles
       to the repository.
   D4. Skip writing ase/<X>/tmp_ase.pdb on every evaluation when the engine
       does not read it (cctbx, in-process engines); fragments are written
       separately anyway.

E. Cluster route with cctbx restraints (tests, qrefine restraints.py /
   cluster_restraints.py, small): select the fragment GRMs once per
   (re)clustering and keep them (pair proxies cached, atoms do not change
   between evaluations) -> 27 ms -> ~2 ms per fragment. Makes the clustering
   tests about 5 x faster per evaluation.

F. Small items (each 1-3 min on this run, later): xh_mapper shrink/gradients
   in flex/C++; process the model once; expansion get_siiu (31 s).

Not changed: the X-ray term per evaluation (FFT), the floor at ~1 s here.

## 6. Verification (results, 2026-11-06)

Benchmark: the same command, all changes in, rt3/run_1762449360 (machine
loaded: the qrefine test suite and another agent's jobs ran alongside).

| phase | before (2026-11-05) | after |
|---|---|---|
| start to first statistics line ("Start:") | ~26 min | 6.4 min (381 s; conformer A fragmentation 2.3 min, B 2.3 min, fit_rotatable2 seconds) |
| weight search | 19 trials, 135 min | 3 trials, 5.7 min (343 s; ~105 s per trial of 53 evaluations) |
| 5 refinement cycles | 15 min | 8.9 min (533 s; 264 evaluations) |
| total wall clock | 176 min | 22 min 17 s |

Per evaluation on the cctbx cluster route: ~1.3-1.5 s (was ~12 s): the
X-ray term, with the fragment restraints cached. Results: R-work/R-free
0.1263/0.1349 (before 0.1253/0.1346), bond rmsd 0.016 (before 0.021), chosen
restraints scale 0.5 (before 0.132: the plateau rule keeps the larger scale
whose R-free is within 0.001 of the best). H fit: r_work/r_free after
"Optimize H" 0.1347/0.1403 (before 0.1347/0.1404).

fit_rotatable2 old vs new on the full 1PWM model (scratchpad hfit/compare.py,
same fmodel, reset_h=False): 430 s -> 63 s; r_work/r_free after the fit
0.1334/0.1392 (old) vs 0.1334/0.1390 (new); 350 of 8452 atoms differ by more
than 0.01 A, 69 by more than 0.3 A (max 2.1 A): a few rotatable groups pick
another rotation in the whole-model OMIT map, at equal or slightly better
R-free.

Clustering: ext.girvan_newman is identical to boost (tst_07; order-included
identity of clustering.get_clusters on 1yjp, 2olx, 2lvr, m00_good); 1PWM
conformer A get_clusters 43 s (was 284-426 s), 148 clusters as before.

Tests: qrefine tst_60 (monitor cache, weight controller), tst_61
(prepare/run_prepared, fragment restraints cache); cctbx
pair_interaction/tst_05 (wfc convention), tst_07 (Girvan-Newman),
hydrogens/tst_fit_rotatable2; tst_11 updated (no tmp_ase.pdb on the cctbx
route). Suites: qrefine run_tests.py, mmtbx pair_interaction and hydrogens
tests (see the plan's ledger for the runs).

Deferred (review minors): weights.controller asserts on bad PHIL values
instead of Sorry; the no-GOOD-trial fallback refines at the last scale (the
old loop stepped up once more); tst_05 leaves its temp folder; tst_fit_rotatable2
takes ~90 s with unseeded free flags. number_of_weight_search_cycles (now 8)
also caps the real-space weight search (documented in its help).
