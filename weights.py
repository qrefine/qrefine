"""
Weight-search controller: which restraints weight scale to try next, when to
stop, and which trial to keep (spec docs/superpowers/specs/
2026-11-05-runtime-bottlenecks.md, section 5.A, with the accelerated step).

A trial is one minimization at a restraints weight scale (rws); it is GOOD
when the bond rmsd is within max_bond_rmsd. The rules:

  - trial 0 at start_scale; the direction comes from it: GOOD -> scale down,
    else up;
  - the step is factor for the first move and grows by factor at every
    trial that keeps the direction and makes progress (going up: every
    not-GOOD trial; going down: every GOOD trial that improves R-free by
    more than r_tolerance), so a scale far from the start is reached in a
    few trials (1, 1/2, 1/8, 1/64, ...);
  - the first flip (GOOD <-> not GOOD) brackets the acceptable scale; the
    bracket is bisected (geometric mean) until its ratio is at most
    sqrt(factor), then done;
  - R-free plateau: done when the last two GOOD trials each failed to improve
    the best GOOD r_free before them by more than r_tolerance;
  - max_trials caps the search; stopped_by_cap() tells when the cap, not a
    rule, ended it.

The choice is, among the GOOD trials whose r_free is within r_tolerance of
the best GOOD r_free, the one with the largest scale (equal R-free: the
better geometry); without a GOOD trial, the last trial (fallback).
"""
from __future__ import absolute_import, division, print_function
import math
from libtbx import group_args

class controller(object):

  def __init__(self, start_scale=1.0, factor=2.0, max_trials=8,
               r_tolerance=0.001):
    assert factor > 1, factor
    assert max_trials >= 1, max_trials
    self.start_scale = start_scale
    self.factor = factor
    self.max_trials = max_trials
    self.r_tolerance = r_tolerance
    self.trials = []
    self._direction = None # +1: scale up, -1: scale down
    self._multiplier = factor # the step of the next move
    self._bracket = None   # (good scale, bad scale) after the flip
    self._bracket_done = False
    self._n_plateau = 0

  def scale(self):
    """The scale of the next trial."""
    if(len(self.trials) == 0): return self.start_scale
    if(self._bracket is not None):
      return math.sqrt(self._bracket[0] * self._bracket[1])
    last = self.trials[-1].scale
    if(self._direction > 0): return last * self._multiplier
    return last / self._multiplier

  def _best_good_r_free(self):
    goods = [t.r_free for t in self.trials if t.good]
    if(len(goods) == 0): return None
    return min(goods)

  def record(self, scale, good, r_free, sites_cart):
    """The result of the trial at scale."""
    improved = False
    if(good):
      best = self._best_good_r_free()
      improved = best is None or r_free < best - self.r_tolerance
      self._n_plateau = 0 if improved else self._n_plateau + 1
    self.trials.append(group_args(scale=scale, good=good, r_free=r_free,
      sites_cart=sites_cart))
    if(self._bracket is not None):
      # a bisection trial narrows the bracket on its side
      g, b = self._bracket
      self._bracket = (scale, b) if good else (g, scale)
      ratio = max(self._bracket) / min(self._bracket)
      if(ratio <= math.sqrt(self.factor) * (1 + 1.e-9)):
        self._bracket_done = True
    elif(self._direction is None):
      self._direction = -1 if good else 1
    elif((good and self._direction > 0) or (not good and self._direction < 0)):
      prev = self.trials[-2].scale
      self._bracket = (scale, prev) if good else (prev, scale)
    elif((not good) or improved):
      # the direction is kept and the trial made progress: longer step
      self._multiplier *= self.factor

  def _done_by_rule(self):
    return self._bracket_done or self._n_plateau >= 2

  def done(self):
    return self._done_by_rule() or len(self.trials) >= self.max_trials

  def stopped_by_cap(self):
    """True when max_trials ended the search before a stopping rule did."""
    return len(self.trials) >= self.max_trials and not self._done_by_rule()

  def choice(self):
    """group_args(scale, r_free, sites_cart, reason): the trial to keep."""
    goods = [t for t in self.trials if t.good]
    if(len(goods) == 0):
      t = self.trials[-1]
      return group_args(scale=t.scale, r_free=t.r_free,
        sites_cart=t.sites_cart, reason="fallback")
    r_min = min(t.r_free for t in goods)
    candidates = [t for t in goods if t.r_free <= r_min + self.r_tolerance]
    t = max(candidates, key=lambda t: t.scale)
    return group_args(scale=t.scale, r_free=t.r_free,
      sites_cart=t.sites_cart, reason="best")
