# Minimizer Robustness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every qr.refine minimization (weight-search trial and refinement cycle) a converged, descent-safe minimization in both L-BFGS modes, so that the weight search compares converged trials.

**Architecture:** Three separate changes. (1) scitbx L-BFGS core: a per-atom step cap `max_shift` honoured by both line searches, and a new gradient-only line search with the strong curvature condition, bracketing, a first step of fixed largest atom shift and a positive-curvature guard on stored pairs, replacing the Snyman two-point secant. (2) qrefine driver: the function-value L-BFGS gets the cap instead of a misused line-search bound. (3) qrefine calculator: the R-work plateau early stop becomes a step-based stop, and every trial line reports why it stopped. A baseline is recorded before any change and a validation matrix is run after all three.

**Tech Stack:** C++ (scitbx/lbfgs.h, boost.python ext), Python 3 (scitbx, qrefine), libtbx build (`make` in `<root>/build`), qr.refine runs with cctbx restraints and MOPAC.

**Spec:** `docs/superpowers/specs/2026-11-06-minimizer-robustness-design.md` (qrefine). Decisions taken 2026-11-06: approach 2 (curvature line search); first step largest atom shift 0.05 A; stop tolerance a constant 0.001 A; "before" numbers come from a baseline run of the unmodified build (Task 0), not a second checkout.

Paths below are relative to `<root>` = `/Users/pafonine/Desktop/all/phenix`; `modules/cctbx_project` and `modules/qrefine` are separate git repositories. `python` means `<root>/build/bin/libtbx.python`; `qr.refine` means `<root>/build/bin/qr.refine`.

## Global Constraints

- cctbx_project AGENTS.md: one problem per change, smallest contained diff, regression test written and failing first, no `from __future__` in new files, 2-space indent, 80 columns, numpydoc docstrings on new public functions, `libtbx.find_clutter` clean, no `Co-Authored-By` trailer in cctbx commits; the user commits and pushes (nothing is committed by the implementer in either repository; each task ends with the diff ready and listed).
- qrefine: every test writes only into its own `qrefine_tst_XX` folder (keep the `os.makedirs(prefix)` / `os.chdir(prefix)` lines in `run()`); never use occupancy in restraints; no `rm -rf` in commands.
- Function-value L-BFGS behaviour must be bit-identical to today for callers that do not set `max_shift` (mmtbx, maptbx, diffBragg): Task 2 proves it.
- After any `lbfgs.h` change: `cd <root>/build && make` (rebuilds the scitbx_lbfgs extension), then the scitbx lbfgs tests and the mmtbx suite.
- Values from the spec: `first_step_max_shift` 0.05 A (passed by qrefine, default 0 = off in scitbx); `max_shift` = `refine.stpmax` (default 0.2 A); gtol 0.9; maxfev 20; curvature guard `s.y <= 1e-10 |s||y|` not stored; duplicate point tolerance 1e-10 relative in stp; bracket safeguard: new stp inside the middle 80 % of the bracket; shift stop 0.001 A for two consecutive accepted steps.

## Review Focus

Inputs the spec implies but no task's test would exercise without the lines added here:

1. `n` not divisible by 3 (non-atomic parameter vectors through the gradient-only path): the cap must apply per coordinate, not crash or silently ignore the remainder. Test added to Task 1.
2. A direction with `d0 >= 0` (non-descent, after a skipped pair or an exact minimum): the gradient-only search must reset to steepest descent, not loop. Test added to Task 2.
3. A zero gradient at the start (already at the minimum): first step `1/|g|` is infinite; the minimizer must return without evaluating. Test added to Task 2.
4. `line_search=False` with the new search: the first trial point is accepted, pair guard still applied. Test added to Task 2.
5. A weight-search trial that hits the cap (unconverged): the trial line must say `stop: cap`, and the weight controller still records it. Test added to Task 4.

---

### Task 0: Baseline of the unmodified build

**Files:**
- Create (outside the repositories): `rt3/minimizer_baseline/<case>/` run folders and `rt3/minimizer_baseline/README.md` listing the commands.

**Interfaces:**
- Produces: the "before" column of the validation matrix (Task 7): per-trial lines, `Best restraints scale`, `Best Rwork, Rfree (after refinement)`, wall clock from the `time:` field.

- [ ] **Step 1: 6uf7, both tolerances**

Run, each in its own folder under `rt3/minimizer_baseline/6uf7_default` and `rt3/minimizer_baseline/6uf7_tol1e-4`:
```bash
qr.refine <root>/rt3/6uf7/6uf7.mtz <root>/rt3/6uf7/6uf7H_charge0.pdb restraints=cctbx > run.log 2>&1
qr.refine <root>/rt3/6uf7/6uf7.mtz <root>/rt3/6uf7/6uf7H_charge0.pdb restraints=cctbx r_tolerance=0.0001 > run.log 2>&1
```
Expected: the series of section 1 of the spec (0.0645 0.0642 0.0640 chosen 1.0; and 0.0645 0.0642 0.0638 0.0643 0.0633 0.0642 0.0642 chosen 0.002).

- [ ] **Step 2: 1PWM**

In `rt3/minimizer_baseline/1PWM`:
```bash
qr.refine <root>/rt3/1PWM.mtz <root>/rt3/1PWMH_charge3.pdb restraints=cctbx > run.log 2>&1
```
Expected: about 22 min on a quiet machine; record wall clock, trial series, chosen scale, final R-work/R-free, bond rmsd (2026-11-06 values: 0.1263/0.1349, 0.016 A, scale 0.5).

- [ ] **Step 3: helix with MOPAC (gradient-only, real QM engine)**

In `rt3/minimizer_baseline/helix_mopac`, build the data as `qrefine/tests/unit/run_tests.py:setup_helix_example` does (write `m00_good.mtz` from `m00_good.pdb` at d_min 1.5 A), then:
```bash
qr.refine m00_good.mtz <root>/modules/qrefine/tests/unit/data_files/m00_poor.pdb restraints=qm quantum.engine_name=mopac clustering=True parallel.nproc=1 > run.log 2>&1
```
Expected: completes; record per-trial `n_fev`, series, chosen scale, final R-factors, wall clock.

- [ ] **Step 4: mid-resolution protein (decision: proposed 1UBQ, 1.8 A, 660 atoms)**

In `rt3/minimizer_baseline/1ubq`: `phenix.fetch_pdb --mtz 1ubq`, add H with `phenix.ready_set` or `reduce2` as the other cases were prepared (`6uf7H_charge0.pdb` came from the qmready pipeline; use `qr.finalise` the same way), then:
```bash
qr.refine 1ubq.mtz 1ubqH.pdb restraints=cctbx > run.log 2>&1
```
Expected: completes; record the same quantities. If the user prefers another entry, substitute it here and in Task 7.

- [ ] **Step 5: write `rt3/minimizer_baseline/README.md`** with the four commands, the build's `git log -1 --oneline` of cctbx_project and qrefine, and the extracted numbers in one table.

---

### Task 1: `max_shift` step cap in the scitbx L-BFGS core

**Files:**
- Modify: `modules/cctbx_project/scitbx/lbfgs.h` (class `minimizer` constructor and members near lines 1040-1075, 1296-1310; `generic_run` lines 1430-1446 where `mcsrch_instance.run` is called; `mcsrch::run` signature lines 382-400 and 490-520, the throw at lines 639-642)
- Modify: `modules/cctbx_project/scitbx/lbfgs/ext.cpp` (the `init<optional<...>>` at line 294-295; new accessors)
- Modify: `modules/cctbx_project/scitbx/lbfgs/__init__.py` (`core_parameters`, `run_c_plus_plus` construction at lines 104-112 and the log lines 96-102)
- Test: `modules/cctbx_project/scitbx/lbfgs/tst_max_shift.py` (new; registered in `modules/cctbx_project/scitbx/run_tests.py` next to line 98)

**Interfaces:**
- Produces (C++): `minimizer(n, m=5, maxfev=20, gtol=0.9, xtol=1e-16, stpmin=1e-20, stpmax=1e20, max_shift=0, first_step_max_shift=0)`; accessors `FloatType max_shift() const`, `FloatType first_step_max_shift() const`; `static FloatType largest_shift(SizeType n, const FloatType* p)` on `minimizer` (largest Euclidean norm over consecutive triplets of `p` when `n % 3 == 0`, else the largest `|p_j|`).
- Produces (Python): `scitbx.lbfgs.core_parameters(m=5, maxfev=20, gtol=0.9, xtol=1e-16, stpmin=1e-20, stpmax=1e20, max_shift=0, first_step_max_shift=0)`; `ext.minimizer` constructed with the two extra positional arguments.
- Semantics: when `max_shift > 0`, every line search of this minimizer (both modes) runs with `stpmax_eff = min(stpmax, max_shift / largest_shift(n, p))`, where `p` is the current search direction; in function mode a step at `stpmax_eff < stpmax` that satisfies the sufficient decrease condition is accepted (no "step is at the upper bound" error). `first_step_max_shift` is stored here and used in Task 2 only.
- `mcsrch::run` loses the `gradient_only` and `line_search` arguments (Task 2 deletes the gradient-only branches; this task only threads `accept_at_stpmax`) — do it in this task to keep one signature change: new signature `run(gtol, stpmin, stpmax, accept_at_stpmax, n, x, f, g, s, is0, stp, ftol, xtol, maxfev, info, nfev, wa)`; the two removed arguments are passed as members until Task 2 removes the branches (`mcsrch_instance.gradient_only`, `.line_search` already exist as members at lines 259-260).

- [ ] **Step 1: Write the failing test** `scitbx/lbfgs/tst_max_shift.py`

```python
def exercise_cap_function_mode():
  # quadratic in 2 "atoms" (n=6), minimum at (10,0,0, 0,0,0) from the origin
  e = evaluator(target=quadratic, x0=[0]*6)
  m = scitbx.lbfgs.run(target_evaluator=e,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.5))
  # no accepted step moved any atom by more than 0.5 (steps recorded in
  # e.callback_after_step from successive x)
  assert max(e.largest_step_per_iteration) <= 0.5 + 1e-12
  assert approx_equal(e.x, [10,0,0, 0,0,0], 1e-4)

def exercise_cap_per_coordinate():
  # n=4: not divisible by 3, cap applies per coordinate
  e = evaluator(target=quadratic4, x0=[0]*4)
  scitbx.lbfgs.run(target_evaluator=e,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.1))
  assert max(e.largest_coordinate_step_per_iteration) <= 0.1 + 1e-12

def exercise_no_cap_identical():
  # max_shift=0 (default) reproduces today's iterates exactly
  e1 = evaluator(target=rosenbrock3, x0=[-1.2, 1.0, 0.5]); scitbx.lbfgs.run(target_evaluator=e1)
  e2 = evaluator(target=rosenbrock3, x0=[-1.2, 1.0, 0.5]); scitbx.lbfgs.run(target_evaluator=e2,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0))
  assert e1.history == e2.history   # list of (x, f) per evaluation

def exercise_accessors():
  mm = scitbx.lbfgs.ext.minimizer(6, 5, 20, 0.9, 1e-16, 1e-20, 1e20, 0.5, 0.05)
  assert mm.max_shift() == 0.5 and mm.first_step_max_shift() == 0.05
```
`evaluator` is a small class in the test with `x`, `compute_functional_and_gradients()`, `callback_after_step(minimizer)` recording `x` per iteration. Register the file in `scitbx/run_tests.py` after `tst_func_free_line_search.py`.

- [ ] **Step 2: Run it to verify it fails**

Run: `python modules/cctbx_project/scitbx/lbfgs/tst_max_shift.py`
Expected: `TypeError` (unexpected keyword `max_shift`) from `core_parameters`.

- [ ] **Step 3: Add the parameters to `core_parameters` and the `ext.minimizer` construction in `__init__.py`**, log them next to `stpmax` when `log` is given.

- [ ] **Step 4: Implement in `lbfgs.h`**: the two constructor parameters and members `max_shift_`, `first_step_max_shift_` (validate `>= 0`, else `error_improper_input_parameter`); `largest_shift`; in `generic_run`, before `mcsrch_instance.run`, compute `stpmax_eff` and `accept_at_stpmax = (stpmax_eff < stpmax_)` when `max_shift_ > 0` (else `stpmax_`, `false`); change `mcsrch::run` to the new signature; in the throw block at line 639 replace the `error_line_search_failed("The step is at the upper bound stpmax().")` by `if (accept_at_stpmax) { info = 1; break; } throw ...`. Keep `mcsrch` members `gradient_only`, `line_search` set from the minimizer's `generic_run` arguments so the gradient-only branches behave as before until Task 2.

- [ ] **Step 5: Bind in `ext.cpp`**: extend `init<optional<std::size_t, std::size_t, std::size_t, double, double, double, double, double, double> >()` and add `.def("max_shift", &w_t::max_shift)`, `.def("first_step_max_shift", &w_t::first_step_max_shift)`.

- [ ] **Step 6: Rebuild and run the test**

Run: `cd <root>/build && make 2>&1 | tail -3 && python modules/cctbx_project/scitbx/lbfgs/tst_max_shift.py`
Expected: build without errors; test prints `OK`.

- [ ] **Step 7: Run the existing scitbx lbfgs tests** (bit-identical function mode)

Run: `for t in tst_ext tst_curvatures tst_lbfgs_fem tst_func_free_line_search tst_mpi_split_evaluator; do python modules/cctbx_project/scitbx/lbfgs/$t.py || echo FAIL $t; done; python modules/cctbx_project/scitbx/examples/lbfgs_recipe.py; python modules/cctbx_project/scitbx/tests/tst_minimizer_lbfgs_general.py`
Expected: all `OK`, no `FAIL`.

- [ ] **Step 8: `libtbx.find_clutter`** on the three changed files; show the diff (`git -C modules/cctbx_project diff --stat`) for the user to commit as "scitbx.lbfgs: per-atom step cap max_shift for both line searches".

---

### Task 2: Gradient-only line search with the curvature condition

**Files:**
- Modify: `modules/cctbx_project/scitbx/lbfgs.h`: delete the "Insertion" blocks in `mcsrch::run` (lines 529, 563-583, 623-628, 663-680) and `mcstep_g` (lines 462-491 declaration, 902-940 definition) and the members `gradient_only`, `line_search`, `gy` of `mcsrch`; add class `detail::gradient_only_line_search`; change `generic_run` (first step at lines 1358-1360 and 1431, the pair storage at lines 1446-1453, the `iter_`-based bookkeeping at lines 1364-1380 and 1395-1400).
- Modify: `modules/cctbx_project/scitbx/lbfgs/ext.cpp`: add `.def("n_pairs", &w_t::n_pairs)`.
- Modify: `modules/cctbx_project/scitbx/lbfgs/__init__.py`: nothing functional; the `run_c_plus_plus` loop keeps passing `gradient_only, line_search` to `minimizer.run`.
- Test: `modules/cctbx_project/scitbx/lbfgs/tst_gradient_only_line_search.py` (new; registered in `scitbx/run_tests.py`); `tst_func_free_line_search.py` unchanged and still passing.

**Interfaces:**
- Consumes: `max_shift_`, `first_step_max_shift_`, `largest_shift` from Task 1.
- Produces (C++): `detail::gradient_only_line_search<FloatType, SizeType>` with
  `void run(FloatType gtol, FloatType stpmax, bool line_search, SizeType maxfev, SizeType n, FloatType* x, const FloatType* g, const FloatType* s, SizeType is0, FloatType& stp, int& info, SizeType& nfev, FloatType* wa)` using the reverse-communication protocol of `mcsrch` (`info == -1`: caller evaluates `g` at `x` and calls again; `info == 1`: done, `x` is the accepted point, `stp` its step, `g` as passed in the last call is its gradient). `minimizer::n_pairs() const` (number of stored curvature pairs, at most `m`).
- Produces (behaviour, gradient-only mode): first trial step of a minimization `stp = first_step_max_shift / largest_shift(n, p)` when `first_step_max_shift > 0`, else `1/|g|`; later iterations `stp = 1`; all clamped by `stpmax_eff` (Task 1). Acceptance `|d(stp)| <= gtol * |d0|`. Overshoot (`d > 0`): bracket `[lo, stp]`, next `stp` = secant zero of `d` on the bracket, pulled into `[lo + 0.1 w, hi - 0.1 w]`, `w = hi - lo`. Undershoot (`d < gtol * d0`, no bracket): `stp <- 2 stp` clamped; if the clamp leaves `stp` unchanged, accept. With a bracket, undershoot moves `lo`. A proposed `stp` within `1e-10` relative of an evaluated one accepts that evaluated point. After `maxfev` evaluations the evaluated point with the smallest `|d|` is chosen: if it is the last evaluated point accept it; otherwise set `x` to it and request one more evaluation (`info = -1`), then accept. `line_search == false`: the first trial point is accepted. If `d0 >= 0` the minimizer discards its pairs (`n_pairs_ = 0`) and uses `p = -g` for this iteration. If `|g| == 0` at the first call, `generic_run` returns `false` (converged) without a line search.
- Pair guard (gradient-only mode only): after the line search, with `s = stp p`, `y = g_new - g_old`, the pair is stored only when `s.y > 1e-10 |s||y|`; storage advances `point` and `n_pairs_`; the top-of-iteration scaling (`ys/yy`) and `bound` use the last stored pair and `n_pairs_` instead of `iter_`; in function mode pairs are always stored, so `n_pairs_ == min(iter_-1, m)` and the indices are those of today.

- [ ] **Step 1: Write the failing test** `scitbx/lbfgs/tst_gradient_only_line_search.py`

Targets (each a function of a flex.double returning `(f, g)`, with `f` used by the test only): `ill_quadratic` (n=6, Hessian diag `[1e4, 1, 1, 1e4, 1, 1]`, minimum at `[1,2,3,4,5,6]`); `rosenbrock3` (n=3, Rosenbrock in `x0,x1` plus `(x2-1)^2`); `gaussian_wells` (n=3, `f = -sum_k exp(-|x-c_k|^2/(2*0.3^2))` with `c_k` at `(0,0,0)` and `(1.2,0,0)`, start `(0.5,0.1,0)`); `abs_target` (n=3, `f = sum |x_i|`, non-smooth, to exhaust the search).

```python
def run_go(target, x0, **core):
  e = evaluator(target, x0)   # records every evaluation (x, f, g) and, in
                              # callback_after_step, the per-iteration x and g
  m = scitbx.lbfgs.run(target_evaluator=e, gradient_only=True,
    core_params=scitbx.lbfgs.core_parameters(**core),
    termination_params=scitbx.lbfgs.termination_parameters(max_iterations=500))
  return e, m

def exercise_converges():
  for target, x0, xmin in [(ill_quadratic, [0]*6, [1,2,3,4,5,6]),
                           (rosenbrock3, [-1.2,1,0], [1,1,1]),
                           (gaussian_wells, [0.5,0.1,0], [0,0,0])]:
    e, m = run_go(target, x0, max_shift=0.2, first_step_max_shift=0.05)
    assert approx_equal(e.x, xmin, 1e-3)
    assert e.f_history[-1] < e.f_history[0]

def exercise_no_duplicate_evaluation():
  for target, x0 in [(ill_quadratic, [0]*6), (gaussian_wells, [0.5,0.1,0])]:
    e, m = run_go(target, x0, max_shift=0.2, first_step_max_shift=0.05)
    for a, b in zip(e.x_history, e.x_history[1:]):
      assert (a - b).norm() > 0, "same point evaluated twice"

def exercise_curvature_pairs():
  e, m = run_go(ill_quadratic, [0]*6, max_shift=0.2, first_step_max_shift=0.05)
  for (x0,g0),(x1,g1) in zip(e.iterates, e.iterates[1:]):
    s = x1 - x0; y = g1 - g0
    assert s.dot(y) > 0

def exercise_first_step():
  e, m = run_go(ill_quadratic, [0]*6, max_shift=0.2, first_step_max_shift=0.05)
  first = e.x_history[1] - e.x_history[0]
  assert approx_equal(largest_atom_shift(first), 0.05, 1e-9)

def exercise_exhaustion_and_guard():
  e, m = run_go(abs_target, [0.3,-0.2,0.1], max_shift=0.2, first_step_max_shift=0.05)
  assert m.n_pairs() <= 5           # never more than m
  assert e.f_history[-1] <= e.f_history[0]

def exercise_non_descent_reset():
  # a target whose gradient flips sign abruptly produces d0 >= 0 once the
  # memory is wrong: f = sum(min(x^2, (x-2)^2)); start at 0.9 (between wells)
  e, m = run_go(two_wells, [0.9, 0.9, 0.9], max_shift=0.5, first_step_max_shift=0.05)
  assert e.f_history[-1] < e.f_history[0]

def exercise_zero_gradient_start():
  e, m = run_go(ill_quadratic, [1,2,3,4,5,6], max_shift=0.2, first_step_max_shift=0.05)
  assert len(e.x_history) == 1

def exercise_no_line_search_flag():
  e = evaluator(ill_quadratic, [0]*6)
  scitbx.lbfgs.run(target_evaluator=e, gradient_only=True, line_search=False,
    core_params=scitbx.lbfgs.core_parameters(max_shift=0.2, first_step_max_shift=0.05))
  # exactly one evaluation per iteration
  assert len(e.x_history) == len(e.iterates)

def exercise_toy_unchanged_in_function_mode():
  # function mode: iterates identical to Task 1's run (history compare)
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python modules/cctbx_project/scitbx/lbfgs/tst_gradient_only_line_search.py`
Expected: `exercise_no_duplicate_evaluation` fails (the Snyman search evaluates the clamped point twice) and `exercise_first_step` fails (first shift is 0.2, the cap).

- [ ] **Step 3: Implement `detail::gradient_only_line_search` in `lbfgs.h`** per the Interfaces block (state: `d0`, `lo`, `d_lo`, `hi`, `d_hi`, `has_hi`, `best_stp`, `best_abs_d`, `evaluated` vector of stp, `nfev`, `final_eval` flag, start point copy in `wa`).

- [ ] **Step 4: Wire it into `generic_run`**: when `gradient_only` call the new class instead of `mcsrch_instance`; first-step rule; `|g| == 0` return; `d0 >= 0` reset; pair guard and `n_pairs_` bookkeeping. Delete the Insertion blocks and `mcstep_g`. Bind `n_pairs` in `ext.cpp`.

- [ ] **Step 5: Rebuild and run the new test**

Run: `cd <root>/build && make 2>&1 | tail -3 && python modules/cctbx_project/scitbx/lbfgs/tst_gradient_only_line_search.py`
Expected: `OK`.

- [ ] **Step 6: Run the scitbx lbfgs tests of Task 1 Step 7 plus `tst_max_shift.py`**
Expected: all `OK`. `tst_func_free_line_search.py` (gradient_only=True on a quartic, `stpmax=500`) must still reach `[2, -5, 500]`.

- [ ] **Step 7: Run the scitbx and mmtbx suites** from an empty scratch directory

Run: `cd <scratch>/scitbx_$(date +%s) && <root>/build/bin/libtbx.run_tests_parallel module=scitbx nproc=8 run_in_unique_dirs=True` and the same with `module=mmtbx`.
Expected: no new failures against the 2026-11-02 baseline list in memory (`mmtbx-preexisting-failures-2026-11-02`); quote any failure in full.

- [ ] **Step 8: `libtbx.find_clutter`** on changed files; show the diff for the user to commit as "scitbx.lbfgs: gradient-only line search with curvature condition".

---

### Task 3: qrefine driver: cap instead of line-search bound (spec section 6)

**Files:**
- Modify: `modules/qrefine/driver.py:31-50` (`run_minimize`)
- Modify: `modules/qrefine/qr.py:170-172` (`stpmax` help), `qr.py:200-203` (`r_tolerance` help)
- Test: `modules/qrefine/tests/unit/tst_63.py` (new; part 1), registered in `modules/qrefine/run_tests.py` after `tst_62.py`

**Interfaces:**
- Consumes: `scitbx.lbfgs.core_parameters(max_shift=, first_step_max_shift=)` (Task 1), the descent property of the function-value search (Task 2 unchanged for that mode).
- Produces: `run_minimize(calculator, params, mode=None)` unchanged signature; in the `lbfgs` branch `core_params = scitbx.lbfgs.core_parameters(max_shift=params.refine.stpmax, first_step_max_shift=0.05)` (stpmin/stpmax defaults); constant `FIRST_STEP_MAX_SHIFT = 0.05` at module level of `driver.py`.
- PHIL help, verbatim: `stpmax`: "Largest shift of any atom in one minimizer step (A) for minimizer=lbfgs; for lbfgsb the bound on the total shift of each coordinate during one minimization." `r_tolerance`: "R-free changes within this value count as a plateau in the weight search (see weights.controller)."

- [ ] **Step 1: Write the failing test** `tests/unit/tst_63.py`, part 1, using the `fixture` of `tests/unit/tst_54.py` (`from qrefine.tests.unit.tst_54 import fixture, get_model, get_fmodel, get_restraints_manager`, `from qrefine.tests.unit import run_fmodel`):

```python
def exercise_function_mode_descends():
  """lbfgs with gradient_only=False: no stpmax crash; the total target is
  non-increasing over accepted iterates and ends below its start."""
  r = run_fmodel(args=["in.pdb", "high_resolution=1.5"], prefix="fm")  # as tst_54 does
  fx = fixture(r)
  c = fx.calculator()
  params = qr.get_default_params()
  params.refine.minimizer = "lbfgs"; params.refine.gradient_only = False
  params.refine.max_iterations_weight = 30
  fs = []
  orig = c.callback_after_step
  def rec(minimizer=None):
    fs.append(c.f); return orig(minimizer=minimizer)
  c.callback_after_step = rec
  driver.run_minimize(calculator=c, params=params, mode="weight")
  assert len(fs) > 2
  assert all(b <= a + 1e-9 for a, b in zip(fs, fs[1:])), fs
  assert fs[-1] < fs[0]

def exercise_function_mode_cap():
  """No accepted step moves an atom by more than refine.stpmax."""
  ... same setup, params.refine.stpmax = 0.05, record largest atom shift
  between consecutive callbacks; assert max <= 0.05 + 1e-9
```
`run()` of the test creates and enters `qrefine_tst_63` (pattern of tst_59 lines 65-71).

- [ ] **Step 2: Run it to verify it fails**

Run: `python modules/qrefine/tests/unit/tst_63.py`
Expected: `RuntimeError: lbfgs error: Line search failed: The step is at the upper bound stpmax().`

- [ ] **Step 3: Change `run_minimize`** as in Interfaces; rewrite the two PHIL help strings.

- [ ] **Step 4: Run the test**
Run: `python modules/qrefine/tests/unit/tst_63.py`
Expected: `OK`.

- [ ] **Step 5: 6uf7 quick check** (gradient-only path now uses the new search; the stop is still the R-work plateau):
Run in a fresh scratch folder: `qr.refine <root>/rt3/6uf7/6uf7.mtz <root>/rt3/6uf7/6uf7H_charge0.pdb restraints=cctbx r_tolerance=0.0001 > run.log 2>&1; grep -n "^[0-9]: Rw\|Weight search\|Best" run.log`
Expected: every trial's first evaluations no longer jump to R-work 0.25 (check with the scratchpad `qr_probe2.py` if in doubt: evaluation 2 moves the largest atom by 0.05 A). Record the series for Task 7.

- [ ] **Step 6: Show the qrefine diff** (`git -C modules/qrefine diff --stat`) for the user to commit as "driver: refine.stpmax is the per-atom step cap in both lbfgs modes".

---

### Task 4: Step-based stop and stop reason on the trial line (spec section 5)

**Files:**
- Modify: `modules/qrefine/calculator.py:237-420` (class `sites`: `_initialize` line 342-343, `update` lines 398-401, `callback_after_step` lines 413-419)
- Modify: `modules/qrefine/driver.py` (`run_minimize` return; the trial line at line 118-119 and the refinement line at lines 165-167)
- Test: `tests/unit/tst_63.py` part 2 (unit test of the callback) and the log-line assertion

**Interfaces:**
- Produces: `sites.SHIFT_TOLERANCE = 0.001` (class attribute, A); `sites.callback_after_step(self, minimizer=None) -> bool` returns True when the largest atom shift between the sites at this call and at the previous call was below `SHIFT_TOLERANCE` for two consecutive calls; it sets `self.stop_reason = "shift"` when it returns True. `sites.reset_stop()` sets `stop_reason = None` and forgets the previous sites (called by `run_minimize` before each minimization). `r_works` and its assertion are removed from `sites` (no other user: grep shows calculator.py only).
- `driver.run_minimize(calculator, params, mode=None)` returns `group_args(minimized=<minimizer wrapper>, stop_reason=<"shift"|"cap"|"conv">)`: `"shift"` when the calculator set it, else `"cap"` when `minimized.minimizer.iter() >= max_iterations` (lbfgs) or `minimized.minimizer.n_calls >= max_iterations` (lbfgsb), else `"conv"`.
- Trial line suffix becomes `"rws: %.3g n_fev: %d stop: %s"`; the refinement-cycle line the same.

- [ ] **Step 1: Write the failing tests** in `tst_63.py`:

```python
def exercise_shift_stop():
  c = fixture(r).calculator()
  c.reset_stop()
  sites = c.sites_cart_nuclear.deep_copy()
  def move(dx):
    s = sites.deep_copy(); s[0] = (s[0][0] + dx, s[0][1], s[0][2])
    c.set_sites(s)        # sets sites_plus_x without minimizing
  for dx, expected in [(0.0, False), (0.01, False), (0.0105, False),
                       (0.0108, True)]:          # last two steps 0.0005, 0.0003
    move(dx); assert c.callback_after_step() is expected, dx
  assert c.stop_reason == "shift"
  c.reset_stop()
  for dx, expected in [(0.0, False), (0.0005, False), (0.02, False)]:
    move(dx); assert c.callback_after_step() is expected
  assert c.stop_reason is None

def exercise_no_rwork_plateau_stop():
  """Three evaluations with identical R-work but large shifts do not stop."""
  ... as above with shifts 0.05, 0.05 and the fmodel's r_work patched to a
  constant; assert callback returns False twice

def exercise_trial_line_has_stop_reason():
  r_qr = run_qrefine_in(folder="run_stop", prefix="stop", args=[mtz, pdb,
    "restraints=cctbx", "max_iterations_weight=3", "number_of_weight_search_cycles=1",
    "number_of_refine_cycles=1", "max_iterations_refine=3"])
  lines = [l for l in r_qr.std_out_text.splitlines() if l.startswith("0:")]
  assert lines and all(" stop: cap" in l for l in lines), lines
```
(mtz/pdb: the helix example data built as in `run_tests.setup_helix_example`, in the test folder.)

- [ ] **Step 2: Run to verify they fail**
Run: `python modules/qrefine/tests/unit/tst_63.py`
Expected: `AttributeError: 'sites' object has no attribute 'reset_stop'`.

- [ ] **Step 3: Implement** in `calculator.sites`: `SHIFT_TOLERANCE`, `reset_stop`, the new `callback_after_step` (largest atom shift via `flex.max(flex.sqrt((sites_plus_x - prev).dot()))`), remove `r_works`; in `driver.run_minimize`: `calculator.reset_stop()` before minimizing, the return value, and the two log suffixes.

- [ ] **Step 4: Run the test**
Run: `python modules/qrefine/tests/unit/tst_63.py`
Expected: `OK`.

- [ ] **Step 5: Run the qrefine tests that read trial lines or the minimizer**
Run: `for t in 02 05 06 24 35 36 44 51 52 55 56 59 60 61 62; do python modules/qrefine/tests/unit/tst_$t.py > tst_$t.log 2>&1 || echo FAIL $t; done` from an empty scratch folder.
Expected: no `FAIL`; a test asserting on the old suffix format or on `r_works` is updated to the new line.

- [ ] **Step 6: 6uf7 quick check** as Task 3 Step 5; expected: no trial stops at 28-31 evaluations with R-free 0.0642; each line ends with `stop: shift` or `stop: cap`.

- [ ] **Step 7: Show the diff** for the user to commit as "calculator: step-based stop, stop reason on cycle lines".

---

### Task 5: 6uf7 regression test (spec section 7)

**Files:**
- Create: `modules/qrefine/tests/unit/data_files/6uf7.mtz` (copy of `rt3/6uf7/6uf7.mtz`, 155 kB; `6uf7_qmready.pdb` already there is atom-identical to `6uf7H_charge0.pdb`)
- Modify: `tests/unit/tst_63.py` part 3

**Interfaces:**
- Consumes: `run_qrefine_in(args, prefix, folder)` from `qrefine/tests/unit/__init__.py`; trial lines of the form `N: Rw: x Rf: y ... rws: s n_fev: k stop: reason`.

- [ ] **Step 1: Write the test**

```python
def parse_trials(text):
  """[(scale, r_free, n_fev, stop)] from the 'Optimal weight search' block."""

def exercise_6uf7_series():
  r_qr = run_qrefine_in(folder="run_6uf7", prefix="w", args=[
    data("6uf7.mtz"), data("6uf7_qmready.pdb"), "restraints=cctbx",
    "r_tolerance=0.0001", "number_of_refine_cycles=1"])
  trials = parse_trials(r_qr.std_out_text)
  assert len(trials) >= 3
  assert all(t.stop in ("shift", "cap") for t in trials)
  # unimodal within 0.0003: after the minimum, no later trial is lower than
  # the minimum by more than 0.0003, and before it no rise above 0.0003
  rf = [t.r_free for t in trials]; k = rf.index(min(rf))
  assert all(rf[i] >= rf[i+1] - 0.0003 for i in range(k))
  assert all(rf[i] <= rf[i+1] + 0.0003 for i in range(k, len(rf)-1))

def exercise_6uf7_function_path():
  r_qr = run_qrefine_in(folder="run_6uf7_f", prefix="f", args=[
    data("6uf7.mtz"), data("6uf7_qmready.pdb"), "restraints=cctbx",
    "clustering=False", "altlocs.handling=native", "minimizer=lbfgs",
    "gradient_only=False", "r_tolerance=0.0001", "number_of_refine_cycles=1"])
  trials = parse_trials(r_qr.std_out_text)
  assert len(trials) >= 3 and all(t.stop in ("shift","cap","conv") for t in trials)
```

- [ ] **Step 2: Run it**
Run: `python modules/qrefine/tests/unit/tst_63.py`
Expected: `OK` in under 3 minutes. If `exercise_6uf7_series` fails on unimodality, this is a finding about the fixes (spec: the fix is revisited, not the criterion): stop and report the series.

- [ ] **Step 3: Register `tst_63.py` in `modules/qrefine/run_tests.py`** (after `tst_62.py`) and run the whole qrefine suite once: `python modules/qrefine/run_tests.py` from an empty scratch folder.
Expected: all pass.

- [ ] **Step 4: Show the diff** (test, data file, run_tests.py) for the user to commit as "tst_63: minimizer stop, function path, 6uf7 weight-search series".

---

### Task 6: phenix hook tests

**Files:** none modified; verification only.

- [ ] **Step 1: Run** `python modules/phenix_regression/refinement/tst_175.py` and `tst_178.py` from an empty scratch folder.
Expected: `OK` for both. A failure that quotes the old trial-line format or `r_works` is fixed in the test (phenix_regression repo, separate diff shown to the user); any other failure is reported with its full output before going on.

---

### Task 7: Validation matrix (spec section 7)

**Files:**
- Create (outside the repositories): `rt3/minimizer_after/<case>/`, and `docs/superpowers/specs/2026-11-06-minimizer-robustness-design.md` section 10 "Results" (qrefine repo) with the before/after table.

- [ ] **Step 1: Rerun the four Task 0 cases** with the final build, same commands, into `rt3/minimizer_after/`.

- [ ] **Step 2: Extract the same quantities** and fill the table: per case, trial series, chosen scale, final R-work/R-free, bond rmsd, evaluations per trial, wall clock, before and after.

- [ ] **Step 3: Check the acceptance criteria of spec section 7** and write each verdict under the table:
  - function-value path on 6uf7: trial targets non-increasing (from tst_63 part 1 and the `qr_probe2.py` trace of one trial: `f` column);
  - gradient-only: every trial line has a stop reason; no duplicate consecutive evaluations in a `qr_probe2.py` trace of trial 0 of 6uf7 and of the helix MOPAC case;
  - each case's series unimodal within 0.0003;
  - 1PWM: R-work/R-free within 0.002 and bond rmsd within 0.003 A of baseline or better; wall clock not more than 10 % longer;
  - helix MOPAC: evaluations per trial not more than baseline; final target not higher (compare `f` from a one-trial trace, since the log has no target);
  - all tests of Tasks 1-6 pass.

- [ ] **Step 4: Report.** For any criterion that fails, name the case and the numbers and stop: the fix is revisited (spec), the criterion is not. Otherwise list the three diffs ready to commit (cctbx_project: two commits; qrefine: three commits) and the follow-ups left open by the spec (controller rules, zero-weight gate).
