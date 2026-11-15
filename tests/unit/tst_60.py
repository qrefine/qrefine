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

def feed(c, r_free, bad=()):
  """Run the controller on the table r_free (scale -> R-free); scales in bad
  are not GOOD. Returns the scales tried."""
  def lookup(sc, table):
    keys = [x for x in table if abs(sc - x) < 1e-9 * x]
    assert len(keys) == 1, (sc, sorted(table))
    return keys[0]
  while(not c.done()):
    sc = c.scale()
    key = lookup(sc, r_free)
    c.record(sc, key not in bad, r_free[key], None)
  return [t.scale for t in c.trials]

def exercise_start_is_best():
  """Both neighbours of the start are worse: three trials, the start kept."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.200, 0.5: 0.210, 2.0: 0.210})
  assert scales == [1.0, 0.5, 2.0], scales
  assert c.choice().scale == 1.0 and c.choice().reason == "best"
  assert not c.stopped_by_cap()

def exercise_minimum_above_start():
  """Trial 0 is GOOD but R-free is lower at tighter restraints: the search
  goes up (it went down only before) and stops past the minimum."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.250, 0.5: 0.260, 2.0: 0.240, 4.0: 0.230,
    8.0: 0.235})
  assert scales == [1.0, 0.5, 2.0, 4.0, 8.0], scales
  assert c.choice().scale == 4.0, c.choice().scale

def exercise_minimum_below_start():
  """An interior minimum below the start, away from the bond rmsd limit."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.250, 0.5: 0.240, 0.25: 0.230, 0.125: 0.236,
    2.0: 0.260})
  assert scales == [1.0, 0.5, 2.0, 0.25, 0.125], scales
  assert c.choice().scale == 0.25, c.choice().scale

def exercise_plateau():
  """R-free equal within r_tolerance does not move the search: flat R-free
  at ever tighter restraints must not drive the scale up without end."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.2000, 0.5: 0.1995, 2.0: 0.1995})
  assert scales == [1.0, 0.5, 2.0], scales
  assert c.choice().scale == 1.0, c.choice().scale
  # flat above a real improvement: stops one step past it
  c = weights.controller(max_trials=12)
  scales = feed(c, {1.0: 0.2500, 0.5: 0.2600, 2.0: 0.2400, 4.0: 0.2395})
  assert scales == [1.0, 0.5, 2.0, 4.0], scales
  assert c.choice().scale == 2.0, c.choice().scale

def exercise_both_sides_of_first_good():
  """The 4xcr case: the scale below the start is not GOOD and the halved gap
  improves a little; the scale above the start is tried all the same."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.2507, 0.5: 0.2526, 2.0: 0.2450,
    4.0: 0.2470}, bad=(0.5,))
  assert scales == [1.0, 0.5, 2.0, 4.0], scales
  assert c.choice().scale == 2.0, c.choice().scale

def exercise_not_good_then_bisect():
  """Not GOOD at 1: up by the factor; GOOD at 2; the scale above it is
  tried, then the gap to the not-GOOD trial is halved once (ratio
  sqrt(factor) left), then done."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.20, 2.0: 0.19, 4.0: 0.195,
    math.sqrt(2.0): 0.185}, bad=(1.0,))
  assert len(scales) == 4, scales
  assert scales[:3] == [1.0, 2.0, 4.0], scales
  assert abs(scales[3] - math.sqrt(2.0)) < 1e-9, scales
  assert abs(c.choice().scale - math.sqrt(2.0)) < 1e-9
  assert not c.stopped_by_cap()

def exercise_not_good_then_up():
  """Not GOOD at 1 and 2, GOOD at 8 (the 2yhj case): the scale above (16),
  where R-free is lowest, is tried; the search went no further than the
  first GOOD scale before."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.3435, 2.0: 0.3425, 8.0: 0.3455, 16.0: 0.3400,
    32.0: 0.3420}, bad=(1.0, 2.0))
  assert scales == [1.0, 2.0, 8.0, 16.0, 32.0], scales
  assert c.choice().scale == 16.0, c.choice().scale

def exercise_bond_limit_binds():
  """R-free falls down to the bond rmsd limit: the limit is approached from
  the GOOD side until the ratio to the not-GOOD scale is sqrt(factor)."""
  c = weights.controller()
  scales = feed(c, {1.0: 0.300, 0.5: 0.290, 2.0: 0.310, 0.25: 0.280,
    math.sqrt(0.5*0.25): 0.285}, bad=(0.25,))
  assert len(scales) == 5, scales
  assert abs(c.choice().scale - math.sqrt(0.125)) < 1e-9, c.choice().scale
  assert not c.stopped_by_cap()

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
  assert c.choice().scale == c.trials[-1].scale
  assert abs(c.choice().r_free - 0.1) < 1e-12
  assert c.stopped_by_cap()

def exercise_accelerates():
  """
  The 6uf7 case: R-free keeps falling at every halving from a start weight
  far too large. From the third improving move in a row the step grows, so
  the minimum near 1/300 is passed within nine trials; the gap the long step
  leaves is halved next.
  """
  def r_free(scale):
    return 0.07 + 0.02 * abs(math.log(scale * 300.))
  c = weights.controller(max_trials=9)
  while(not c.done()):
    sc = c.scale()
    c.record(sc, True, r_free(sc), None)
  scales = [t.scale for t in c.trials]
  assert scales[:3] == [1.0, 0.5, 2.0], scales
  down = [s for s in scales if s < 1]
  ratios = [down[i] / down[i+1] for i in range(len(down) - 1)]
  assert max(ratios) > 4, ratios # accelerated
  assert min(scales) < 1/300., scales
  assert 1/2400. < c.choice().scale < 1/37., c.choice().scale
  assert c.stopped_by_cap()

def exercise_sites_carried():
  """choice() carries the sites of the chosen trial."""
  c = weights.controller(r_tolerance=0.0)
  c.record(1.0, True, 0.2, "s1")
  c.record(0.5, True, 0.21, "s2")
  c.record(2.0, True, 0.22, "s3")
  assert c.done()
  assert c.choice().sites_cart == "s1"

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  exercise_stats_cached()
  exercise_start_is_best()
  exercise_minimum_above_start()
  exercise_minimum_below_start()
  exercise_plateau()
  exercise_both_sides_of_first_good()
  exercise_not_good_then_bisect()
  exercise_not_good_then_up()
  exercise_bond_limit_binds()
  exercise_never_good()
  exercise_cap()
  exercise_accelerates()
  exercise_sites_carried()

if(__name__ == "__main__"):
  run()
  print("OK")
