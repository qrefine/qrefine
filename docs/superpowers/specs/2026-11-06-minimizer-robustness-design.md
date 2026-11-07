# Minimizer robustness in qr.refine: gradient-only line search, trial convergence, step cap

Status: proposal 2026-11-06, for review before an implementation plan is
written. Nothing is implemented. Evidence from one case (6uf7); the design
targets the general case and section 7 says how that is checked.

## 1. Problem

    qr.refine 6uf7.mtz 6uf7H_charge0.pdb restraints=cctbx [r_tolerance=0.0001]

95 atoms (5-residue peptide, A/B altlocs), cell 12 x 11 x 26 A, R-work
0.067 at start. The weight search picks restraints scale 1.0 with the
default r_tolerance and 0.002 with r_tolerance=0.0001. The trial series
(R-free per trial) is not a function of the weight: 0.0645, 0.0642,
0.0638, 0.0643, 0.0633, 0.0642, 0.0642 for scales 1, 1/2, 1/8, 1/64,
1/512, 1/8192, 1/131072.

Instrumented runs (scratchpad qr_probe.py, qr_probe2.py: every target
evaluation of a trial logged with f, |g|, R-work and the largest atom
shift) show that the weight controller (qrefine/weights.py) does what its
spec says; the trials it compares are unconverged minimizations whose end
point is set by where the minimizer happened to be cut off.

Per-evaluation trace of trial 0 (scale 1, default stpmax=0.2):

| evaluation | f | &#124;g&#124; | R-work | largest atom shift in the step |
|---|---|---|---|---|
| 1 (start) | 834.4 | 1108 | 0.0670 | 0 |
| 2 | 5124.1 | 24808 | 0.2530 | 0.200 |
| 3 | 5124.1 | 24808 | 0.2530 | 0 (same point again) |
| 4 to 30 | falls to 829.7 | falls to 185 | falls to 0.0665 | 0.04 to 0.001 |
| 31 to 52 | rises to 831.3 | ~80 | 0.0665 to 0.0666 | 0.003 |

Same trial with stpmax=0.02: evaluation 2 raises f from 834.4 to 858.8
(R-work 0.0700), every later iteration moves the largest atom by 0.0003 A
and evaluates each point twice; after 92 evaluations f is 835.3, above
the start.

Variants of the full run, R-free per trial:

| variant | series | chosen |
|---|---|---|
| as shipped, r_tolerance 0.0001 | 0.0645 0.0642 0.0638 0.0643 0.0633 0.0642 0.0642 | 0.002 |
| early-stop callback disabled | 0.0645 0.0642 0.0637 0.0632 0.0633 0.0633 | 0.016 |
| callback disabled, 200 iterations | 0.0643 0.0639 0.0629 0.0626 0.0620 0.0616 0.0624 0.0615 | 7e-9 |
| stpmax=0.02 | 0.0650 0.0651 0.0649 | 1.0 |
| clustering=False altlocs.handling=native minimizer=lbfgsb (function-value line search) | 0.0644 0.0641 0.0639 0.0642 0.0642 | 0.125 |
| clustering=False altlocs.handling=native minimizer=lbfgs gradient_only=False | crash: "Line search failed: The step is at the upper bound stpmax()" | - |

With a function-value line search every trial converges (14 to 51
evaluations), the five refinement cycles after the search reproduce the
same R-work to 5 decimals, and the series is smooth.

## 2. Mechanisms

Three independent defects, each general (none depends on 6uf7):

A. The gradient-only line search (scitbx/lbfgs.h, the "Insertion" blocks
   in mcsrch and mcstep_g) has no descent safeguard. Per iteration it
   evaluates the trial point, computes a secant step from the two
   gradients (Snyman), evaluates that, and accepts it whatever happened
   (nfev == 2). Consequences:
   - the first trial step of the minimization is the L-BFGS convention
     stp = 1/|g| (1 A of total displacement over all coordinates), then
     clamped so the largest atom shift is in [0.0001, stpmax] A. For
     8000 atoms that is 0.01 A per atom; for 95 atoms it hits the 0.2 A
     cap. The step size depends on the model size, not on the target;
   - a step that overshoots is kept, so the L-BFGS memory receives
     (s, y) pairs with s.y <= 0 and the quasi-Newton matrix stops being
     positive definite: the trace shows f rising again after evaluation 30;
   - when both the trial step and the secant step hit the same clamp the
     same point is evaluated twice (half of all evaluations in the
     stpmax=0.02 run, evaluation 3 in every default run);
   - refine.stpmax (PHIL help: "maximum step length") is the per-atom
     per-step cap only in this mode.
   qr.refine forces this mode whenever cluster.clustering=True (the
   default), expansion=True, or the model has alternative conformations
   with altlocs.handling=conformers (api.py finalize_params): the cluster
   and conformer routes return target 0 (restraints.py, conformers.py),
   there is no function value. So this is the production QM path for
   every user, and it is what the 6uf7 cctbx run used.

B. The early stop in calculator.sites.callback_after_step stops the
   minimizer when R-work rounded to 5 decimals is identical over the last
   three evaluations (not iterations: r_works has one entry per
   evaluation, so a duplicated evaluation from A counts twice). It fires
   on false plateaus during the walk of A (trials 3, 5, 6 in the second
   log: 28 to 31 evaluations, R-free 0.0642; trial 4 ran to the cap and
   reached 0.0633), and a 1e-5 absolute change in R-work means different
   things at R 0.06 and at R 0.30.

C. driver.run_minimize passes refine.stpmax (0.2) as the line-search
   multiplier bound (core_parameters.stpmax) of the function-value L-BFGS.
   After the first iteration the Moré-Thuente search starts at stp = 1 >
   0.2 and throws "The step is at the upper bound stpmax()". The
   function-value minimizer is unusable in qr.refine (restraints=cctbx,
   clustering=False, no conformers), which is also the only mode where a
   descent property can be asserted in tests.

The weight controller itself (choice within r_tolerance, plateau stop) is
not changed by this proposal; its rules can only be judged once the
trials it compares are converged.

## 3. Scope

In scope: A, B, C above, the log line of a trial, the tests and the
validation of section 7.

Out of scope, listed for later decisions:
- the controller's choice rule (largest scale within r_tolerance of the
  best R-free) and plateau rule, discussed 2026-11-06;
- the weight search running to zero restraints weight when the data alone
  hold the geometry (6uf7 with 200 iterations: bond rmsd stays 0.012 at
  scale 7e-9), i.e. whether max_bond_rmsd is a sufficient GOOD criterion
  at sub-A resolution;
- sites_opt (geometry optimization) and sites_real_space calculators and
  their own callbacks;
- lbfgsb (bounded) mode: its box of +-stpmax on the total shift is a
  different semantics and is left as is.

## 4. Design A: gradient-only line search with a curvature condition

Location: scitbx/lbfgs.h (the gradient_only branches of mcsrch and
mcstep_g are replaced by one routine), scitbx/lbfgs/__init__.py
(core_parameters), scitbx/lbfgs/tst_func_free_line_search.py (kept) plus
a new test. The only production callers with gradient_only=True are
qrefine and the phenix xyz_qrefine hook (mmtbx has it commented out,
maptbx and diffBragg pass False).

Approaches considered:

1. Minimal patch: keep the Snyman secant; fix the first step; reject a
   secant point when s.y <= 0 and halve; skip duplicate evaluations.
   Smallest diff, but descent still rests on a two-point secant with no
   bracketing, and the acceptance rule stays "whatever the second
   evaluation gives".
2. Gradient-only Wolfe curvature line search (recommended): accept a step
   only when the directional derivative satisfies the strong curvature
   condition; bracket and interpolate on overshoot, expand on undershoot.
   This is the standard way to run BFGS without function values; it
   guarantees s.y > 0 at every accepted step, hence a positive definite
   L-BFGS matrix and descent directions. A few more evaluations per
   iteration in the worst case, usually the same two.
3. Function-value line search using the data term plus the restraint
   energy where it exists. Not available on the cluster and conformer
   routes (target is 0 there by construction), so it cannot be the
   production path.

Design (approach 2):

- Notation: p the L-BFGS direction, d(stp) = g(x + stp p).p the
  directional derivative, d0 = d(0) < 0, dr_max(p) the largest per-atom
  norm of p, so the largest atom shift of a step is stp * dr_max(p).
- Descent check: if d0 >= 0 the L-BFGS memory is reset and p = -g (today
  this case throws only in function mode). With the s.y guard below it
  should not occur; the reset is the safety net.
- Step cap: every trial stp is clamped to stp <= max_shift / dr_max(p),
  max_shift a new core parameter (A per atom per step; qrefine passes
  refine.stpmax, default 0.2; None = no cap). This is the semantics
  refine.stpmax already has in gradient-only mode, made explicit and
  applied in both modes (section 6).
- First trial step of a minimization: stp such that the largest atom
  shift is first_step_max_shift (new core parameter, default 0.05 A;
  qrefine does not expose it). Later iterations start at stp = 1 as
  L-BFGS expects. The 1/|g| convention goes for gradient-only mode (it
  stays for function mode, where the Moré-Thuente search corrects it).
  The value matters little once the search is correct: an overshoot is
  bracketed on the next evaluation; an undershoot doubles.
- Acceptance: |d(stp)| <= gtol * |d0| (gtol = 0.9, existing parameter).
  This implies d(stp) > d0 and therefore s.y = stp (d(stp) - d0) > 0.
- Overshoot, d(stp) > 0: bracket [lo, stp] where lo is the last point
  with negative derivative (0 at first); next stp is the secant zero of
  d on the bracket, safeguarded into the middle 80 % of the bracket.
- Undershoot, d(stp) < gtol * d0 and no bracket yet: stp <- 2 stp,
  clamped by the cap; if the clamp leaves stp unchanged the point is
  accepted (the cap is the user's limit, not a failure). With a bracket,
  the secant/bisection continues with lo <- stp.
- Never evaluate the same point twice: a proposed stp within 1e-10
  relative of an evaluated one accepts the evaluated point.
- Exhaustion: after maxfev (20) evaluations the evaluated point with the
  smallest |d| is accepted; if its s.y <= 0 the pair is not stored.
- L-BFGS update guard (gradient-only mode only): a pair with
  s.y <= 1e-10 |s||y| is not stored. Function mode keeps the Wolfe
  guarantees it has.
- line_search=False (existing flag) keeps its meaning: the first trial
  point is accepted; the s.y guard still applies.
- Behaviour change for current users: every QM refinement with clustering
  takes a different path through coordinate space (first step 0.05 A
  instead of up to 0.2 A on the largest-gradient atom; no uphill
  acceptances). Expected: fewer evaluations to the same or lower target,
  never a trial ending above its start. This is checked in section 7, not
  assumed.

## 5. Design B: trial convergence and early stop (qrefine)

Location: qrefine/calculator.py (class sites), qrefine/driver.py.

- Remove the R-work plateau stop from sites.callback_after_step.
- Replace it with a step-based stop evaluated once per accepted L-BFGS
  iteration (the callback is called after each iteration; the calculator
  records the largest atom shift of the accepted step): stop when that
  shift is below 0.001 A for two consecutive iterations. With the line
  search of A, tiny accepted steps mean the line minimum is at the
  current point in two consecutive directions, which is convergence; with
  the current search they only meant the clamp floor. The value is a
  constant in the calculator (decision 8.3: constant or PHIL).
- The stop applies in both the weight-search trials and the refinement
  cycles; max_iterations_weight and max_iterations_refine stay the caps.
- The trial line gains the stop reason, e.g. "n_fev: 34 stop: shift" or
  "stop: cap", so a cap-limited trial is visible in the log (today an
  unconverged trial is indistinguishable from a converged one).
- The weight controller is told nothing new; it keeps comparing R-free of
  trials, which are now converged to a common criterion.
- r_tolerance: the PHIL help still says "R factor delta between
  refinement cycles to determine convergence"; after this change it is
  used by the weight controller only, and the help is corrected.

## 6. Design C: stpmax semantics in function-value mode (qrefine)

Location: qrefine/driver.py run_minimize, qrefine/qr.py PHIL help.

- The function-value L-BFGS gets the default line-search bounds
  (stpmin, stpmax of core_parameters untouched) and the per-atom cap
  through max_shift (section 4), so refine.stpmax means the same thing in
  both L-BFGS modes: largest atom shift per minimizer step.
- The "step at the upper bound" exception then only occurs at the real
  default bound (1e20) and is left as is.
- PHIL help of refine.stpmax rewritten: "Largest shift of any atom in one
  minimizer step (A) for minimizer=lbfgs; for lbfgsb the bound on the
  total shift of each coordinate during one minimization."
- No change to lbfgsb.

## 7. Generality and validation

One case drove the analysis; the fixes are argued general because each
removes a model-size or R-factor-scale dependence (section 2). The
validation is designed to catch where that argument fails.

Unit tests (written first, each failing before its fix):
- scitbx: a new test next to tst_func_free_line_search.py running the
  gradient-only minimizer on (a) a quadratic with condition number 1e4,
  (b) the Rosenbrock function, (c) a sum-of-Gaussians "crystallographic"
  target (periodic in one coordinate, so an overshoot is a real risk),
  asserting: no two consecutive evaluations at the same point, s.y > 0 at
  every stored pair, f at the end below f at the start, and convergence
  to the known minimum. The existing toy test keeps passing.
- qrefine: tst_60 (controller) unchanged; a new test on the 6uf7 files
  (copied into the test data, as other tests do) with restraints=cctbx,
  asserting for every weight-search trial that no two consecutive
  evaluations are at the same coordinates, that the trial line reports
  its stop reason, and that the trial series is unimodal within 0.0003
  in R-free (a property of this case: refined start model, data that
  hold the geometry). A second part runs the function-value path
  (clustering=False altlocs.handling=native minimizer=lbfgs
  gradient_only=False), checks it no longer raises, and asserts that the
  total target of every trial is non-increasing over the accepted
  iterates and ends below its start.
- qrefine suite: tests/unit via run_tests.py, with particular attention
  to the minimizer tests tst_05, tst_24, tst_35, tst_36, tst_56 and the
  weight-search tests tst_02, tst_06, tst_44, tst_51, tst_52, tst_55,
  tst_59; phenix_regression refinement tst_175 and tst_178 (the
  qrefine-in-phenix hook, gradient-only).

Cases (before/after on the same machine, same inputs; before = current
code, after = all three fixes):

| case | path exercised | what is compared |
|---|---|---|
| 6uf7 cctbx (this report) | gradient-only, conformers | trial series, chosen scale, final R-work/R-free, bond rmsd, evaluations |
| 1PWM cctbx (2026-11-05 benchmark, 8452 atoms, d_min 0.92 A) | gradient-only, clustering, conformers | the same, plus wall clock against the 22 min 17 s of 2026-11-06 |
| helix example of the test suite with MOPAC PM7 (MOPAC_COMMAND is set on this machine), restraints=qm, clustering on | gradient-only with a real QM engine | per-trial evaluations, final target, R-factors; the QM path must not get slower |
| a mid-resolution protein with cctbx restraints (1US0 or 1R6J from rt3/, or a 2 A phenix_regression model) | gradient-only, clustering, larger shifts per step | trial series shape, no trial ending above its start |
| 6uf7 function-value path | lbfgs with Moré-Thuente | runs; trial series agrees with the gradient-only series within 0.0005 |

Acceptance criteria. What is general and what is case-specific is kept
apart: the minimizer's objective is the total target (data term times
weight plus restraints), not R-work; R-work may legitimately rise in a
trial at a large restraints scale on a distorted start, and on the
cluster and conformer routes no function value exists at all, so descent
cannot be asserted from a run there.
- function-value mode: the total target is non-increasing over the
  accepted iterates of every trial and ends below its start (exact,
  guaranteed by the Wolfe conditions);
- gradient-only mode: every stored L-BFGS pair has s.y > 0, no two
  consecutive evaluations at the same coordinates, every trial line
  reports its stop reason; descent of the method itself is verified on
  the scitbx toy targets, where f is known;
- R-factors are case-specific expectations only: the trial series of
  each case is smooth (unimodal within 0.0003), and before/after results
  are compared as below, never "R-work not above start";
- 1PWM: R-work/R-free and bond rmsd within 0.002 / 0.003 A of the
  2026-11-06 result or better, wall clock not more than 10 % longer;
- QM case: evaluations per trial not more than before, final target not
  higher;
- all listed tests pass; libtbx.find_clutter clean on changed files.

If a case fails a criterion the fix is revisited, not the criterion.

## 8. Decisions needed

1. Approach for A: 2 (curvature line search, recommended) or 1 (minimal
   patch).
2. First-step largest atom shift: 0.05 A proposed. Alternatives: 0.02 A
   (safer for sub-A data, one or two more evaluations at mid resolution)
   or scaled from the data d_min (not recommended: the QM term has no
   resolution).
3. Step-based stop tolerance 0.001 A: constant in the calculator
   (proposed) or a PHIL parameter.
4. Keep the old gradient-only search reachable during validation (an
   internal argument of scitbx.lbfgs.run, removed afterwards) so before
   and after can be run from one build, or run "before" from a separate
   checkout. Proposed: separate checkout, no transitional code.
5. The validation set of section 7: which mid-resolution case; whether a
   second QM engine is needed.

## 9. Implementation order

Three separate changes, each with its tests, in this order; the 6uf7
series is rerun after each:

1. C (qrefine driver/PHIL) together with the max_shift core parameter in
   scitbx, since C depends on it. Small; makes the function-value path
   usable as a reference.
2. A (scitbx line search). The core change.
3. B (qrefine early stop and log line).
Then the validation matrix of section 7, then the follow-up decisions of
section 3 (controller rules, zero-weight gate) in a separate design.

cctbx_project rules (AGENTS.md): one problem per change, regression test
first, run scitbx and mmtbx suites after the lbfgs.h change, no
co-authorship trailers in cctbx commits, the user commits and pushes.

## 10. Results (2026-11-06, implementation and validation)

Implemented as planned (docs/superpowers/plans/2026-11-06-minimizer-robustness.md)
with two deviations, both ledgered: the step-based stop of section 5 was
dropped (it fired after 6-7 evaluations on unconverged 6uf7 trials; a
per-evaluation trace showed accepted steps of 0.001-0.03 A for 50
iterations, so a per-step shift says nothing about convergence on an
ill-conditioned target); there is no early stop now, the iteration budget
is the common criterion and every line says stop: cap|conv. And
driver.refine recomputes the bulk-solvent mask at the end of every
refinement cycle, because the stale mask made the phenix hook's R-factor
consistency check (1e-4) fail once the model actually moved. The fresh
review added: the real-space calculator and phenix.deviation now pass the
per-atom cap as max_shift (they lost it under the new stpmax semantics);
the gradient-only search keeps the last point instead of a zero step when
no trial reduced |d|; a zero gradient after a step converges instead of
throwing.

Tests: scitbx tst_max_shift, tst_gradient_only_line_search (both
registered), scitbx suite (one pre-existing unrelated failure), mmtbx suite
(the 6 known failures plus a registered-but-missing
ligands/tst_electrons_his.py, all pre-existing), qrefine tst_02 05 06 24 35
36 44 51 52 55 56 59 60 61 62 63, phenix_regression tst_175 and tst_178:
all pass.

Before = unmodified build (rt3/minimizer_baseline), after = final code
(rt3/minimizer_after), same inputs, the four cctbx/MOPAC cases run
simultaneously in both passes, 1ubq alone.

| case | R-free per trial, before | after | chosen | final Rw/Rf before -> after | bond | wall |
|---|---|---|---|---|---|---|
| 6uf7, default | 0.0645 0.0642 0.0640 | 0.0643 0.0638 0.0637 | 1.0 / 1.0 | 0.0667/0.0645 -> 0.0661/0.0653 | 0.009 | 53 s -> 57 s |
| 6uf7, r_tolerance 1e-4 | 0.0645 0.0642 0.0638 0.0643 0.0633 0.0642 0.0642 | 0.0643 0.0638 0.0635 0.0627 0.0633 0.0634 | 0.002 / 0.016 | 0.0634/0.0633 -> 0.0596/0.0594 | 0.012 -> 0.011 | 68 s -> 71 s |
| 1PWM | 0.1367 0.1357 0.1349 | 0.1370 0.1354 0.1347 0.1347* 0.1344 | 0.5 / 0.5 | 0.1263/0.1349 -> 0.1258/0.1349 | 0.016 | 20 min 13 s -> 25 min 04 s |
| helix, MOPAC PM7, clustering | 0.0264 0.0242 0.0210 0.0197 0.0207 0.0206 | 0.0264 0.0239 0.0209 0.0194 0.0197 0.0200 | 0.016 / 0.016 | 0.0146/0.0147 -> 0.0144/0.0146 | 0.012 -> 0.011 | 5 min 19 s -> 4 min 45 s |
| 1ubq, 1.8 A, generated flags | 0.2139 0.2181 0.2160 | 0.2194 0.2229 0.2204 | 1.0 / 1.0 | 0.1793/0.2220 -> 0.1795/0.2285 | 0.022 | 71 s -> 76 s |

(* bond rmsd 0.031, not GOOD; the bracket was bisected.) 1ubq starts at
R-work/R-free 0.1920/0.1929: the deposited model was refined against the
reflections now flagged free, so every movement raises R-free (both
builds), and its R-free column measures that de-biasing, not the
minimizer. Its series shape, chosen scale and bond rmsd agree.

Acceptance criteria of section 7:
- function-value mode: total target non-increasing over accepted iterates
  and below the start: holds (tst_63 exercise_function_mode_descends; the
  6uf7 function path runs).
- gradient-only: every stored pair has s.y > 0 (tst_gradient_only_line_search);
  no two consecutive evaluations at the same point in the per-evaluation
  traces of 6uf7 trial 0 and of the helix MOPAC trial 0 (the one repeat is
  qrefine's apply_x re-evaluation at the end of a minimization); every
  line has a stop reason.
- series unimodal within 0.0003: holds for all five cases.
- 1PWM: R-work/R-free within 0.002 and bond rmsd within 0.003 A: holds
  (0.1258/0.1349 vs 0.1263/0.1349, 0.016). Wall clock not more than 10 %
  longer: FAILS, +24 %. Cause: 561 target evaluations instead of 423 (five
  trials instead of three, since R-free now keeps improving down to scale
  0.0625, and no early stop in the refinement cycles); the cost per
  evaluation fell (2.68 s vs 2.87 s). This is the price of the dropped
  early stop and of a search that explores what the controller asks for;
  max_iterations_weight is the knob.
- helix MOPAC: evaluations per trial unchanged (52-53), R-free lower in
  every trial and at the end; the target comparison was replaced by the
  R-factors (a baseline trace needs the old build).

Observed and left open:
- With the restraints scale fixed at 1.0 on an already refined model (6uf7
  default, 1ubq), the refinement cycles raise R-free a little (6uf7
  0.0643 -> 0.0653 over five cycles) while R-work and bond rmsd do not
  change: the per-evaluation trace shows the total target creeping up by
  0.2 % after its minimum at evaluation 13, i.e. accepted steps along flat
  directions (H atoms move up to 0.4 A). This is inherent to a line search
  without function values: the strong curvature condition cannot reject a
  step that is a line minimum beyond a hump. Tightening gtol from 0.9 to
  0.1 triples the evaluations (155 vs 53) and does not remove it (f end
  831.58 vs 831.82, minimum 830.28). The old search had the same drift
  after its initial 0.2 A jump. A function-value line search where the
  target exists (cctbx restraints without clustering, altlocs native) does
  not drift; the production QM path has no target.
- The weight controller's choice rule picked scale 0.5 on 1PWM over 0.0625
  (R-free 0.1344 vs 0.1354, within r_tolerance) and 1.0 on 6uf7 over 0.25
  (0.0637 vs 0.0643): the rule of section 3, unchanged, still to be
  decided.
- The early stop: the trial budget is now the only stop. A relative
  R-work/R-free plateau measured per iteration is the candidate if the
  cost matters (1PWM +24 %).
