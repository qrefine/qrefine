from __future__ import absolute_import, division, print_function
import os, math
from io import StringIO
from mmtbx.hydrogens.tst_xh_mapper import pdb_str, get_model
from qrefine import results, calculator, weights

def exercise_stats_cached():
  """show() twice without update() computes the statistics once; every line
  carries the elapsed time."""
  model = get_model(pdb_str, "n_gaussian", use_neutron_distances=False)
  calls = []
  orig = model.geometry_statistics
  def counting(**kw):
    calls.append(1)
    return orig(**kw)
  model.geometry_statistics = counting
  log = StringIO()
  monitor = results.manager(model=model,
    geometry_rmsd_manager=calculator.cctbx_geometry(model=model), log=log)
  monitor.show(prefix="a")
  monitor.show(prefix="b")
  assert len(calls) == 1, len(calls)
  monitor.update(model=model)
  monitor.show(prefix="c")
  assert len(calls) == 2, len(calls)
  lines = log.getvalue().splitlines()
  assert len(lines) == 3, lines
  for l in lines:
    assert "time:" in l, l

def exercise_1pwm_sequence():
  """The 1PWM log (factor 1.5, all GOOD): the R-free plateau stops the search
  after three trials; the chosen scale is the largest within r_tolerance of
  the best R-free."""
  c = weights.controller(factor=1.5)
  for r in [0.1369, 0.1360, 0.1356, 0.1350, 0.1357]:
    if(c.done()): break
    c.record(scale=c.scale(), good=True, r_free=r, sites_cart=None)
  assert len(c.trials) == 3, [t.r_free for t in c.trials]
  assert c.done()
  ch = c.choice()
  assert abs(ch.scale - 1/1.5) < 1e-9, ch.scale
  assert ch.reason == "best", ch.reason

def exercise_flip_then_bisect():
  """Not GOOD at 1: up by the factor; GOOD at 2 is the flip; one trial at the
  geometric mean, then done."""
  c = weights.controller()
  assert abs(c.scale() - 1.0) < 1e-9
  c.record(1.0, False, 0.20, None)
  assert abs(c.scale() - 2.0) < 1e-9, c.scale()
  c.record(2.0, True, 0.19, None)
  assert not c.done()
  assert abs(c.scale() - math.sqrt(2.0)) < 1e-9, c.scale()
  c.record(c.scale(), True, 0.185, None)
  assert c.done()
  assert abs(c.choice().scale - math.sqrt(2.0)) < 1e-9

def exercise_never_good():
  """No GOOD trial within the cap: the step accelerates (2, 4, 8) and the
  fallback is the last trial."""
  c = weights.controller(max_trials=4)
  while(not c.done()):
    c.record(c.scale(), False, 0.3, None)
  assert [t.scale for t in c.trials] == [1.0, 2.0, 8.0, 64.0], \
    [t.scale for t in c.trials]
  ch = c.choice()
  assert ch.reason == "fallback", ch.reason
  assert abs(ch.scale - 64.0) < 1e-9, ch.scale
  assert c.stopped_by_cap()

def exercise_cap():
  """Improving every trial: the cap stops the search; the best is chosen
  and the controller reports that it was still improving."""
  c = weights.controller(max_trials=3, r_tolerance=0.0)
  for r in [0.3, 0.2, 0.1]:
    assert not c.done()
    c.record(c.scale(), True, r, None)
  assert c.done()
  assert abs(c.choice().scale - 0.125) < 1e-9, c.choice().scale
  assert c.stopped_by_cap()

def exercise_accelerates_then_bisects():
  """
  The 6uf7 case: R-free keeps falling at every halving from a start weight
  100 times too large. The step grows while trials improve (1, 1/2, 1/8,
  1/64, 1/1024), a not-GOOD trial brackets the scale, the bracket is bisected
  until its ratio is at most sqrt(factor), the best GOOD trial is chosen.
  """
  r_free = {1.0: 0.2434, 0.5: 0.2212, 0.125: 0.2054, 1/64.: 0.0772,
    1/1024.: 0.0700}
  c = weights.controller(max_trials=12)
  for k in range(5):
    sc = c.scale()
    assert any(abs(sc - x) < 1e-12 for x in r_free), sc
    c.record(sc, True, r_free[[x for x in r_free if abs(sc - x) < 1e-12][0]],
      None)
  assert not c.done()
  bad = c.scale()
  assert abs(bad - 1/32768.) < 1e-15, bad
  c.record(bad, False, 0.0690, None)
  assert not c.done()
  n_bisect = 0
  lo, hi = 1/32768., 1/1024.
  while(not c.done()):
    sc = c.scale()
    assert lo < sc < hi, (lo, sc, hi)
    assert abs(sc - math.sqrt(lo*hi)) < 1e-15
    good = sc > 3e-4
    c.record(sc, good, 0.0695 if good else 0.0690, None)
    if(good): hi = sc
    else:     lo = sc
    n_bisect += 1
  assert hi / lo <= math.sqrt(2.0) + 1e-9, (lo, hi)
  assert n_bisect == 4, n_bisect
  assert c.choice().reason == "best"
  assert c.choice().r_free <= 0.0700
  assert not c.stopped_by_cap()

def exercise_sites_carried():
  """choice() carries the sites of the chosen trial."""
  c = weights.controller(r_tolerance=0.0)
  c.record(1.0, True, 0.2, "s1")
  c.record(0.5, True, 0.21, "s2")
  c.record(0.25, True, 0.22, "s3")
  assert c.done()
  assert c.choice().sites_cart == "s1"

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_stats_cached()
  exercise_1pwm_sequence()
  exercise_flip_then_bisect()
  exercise_never_good()
  exercise_cap()
  exercise_accelerates_then_bisects()
  exercise_sites_carried()

if(__name__ == "__main__"):
  run()
  print("OK")
