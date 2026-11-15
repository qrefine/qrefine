"""
Weight-search controller: which restraints weight scale to try next, when to
stop, and which trial to keep.

A trial is one minimization at a restraints weight scale (rws); it is GOOD
when the bond rmsd is within max_bond_rmsd. The aim is the GOOD trial with
the lowest R-free; both directions from the start scale are searched, since
the R-free minimum may be at a tighter or at a looser scale than the start
and need not be at the bond rmsd limit. The rules:

  - trial 0 at start_scale;
  - no GOOD trial yet: scale up, the step growing by factor at every trial
    (1, 2, 8, 64, ...), until a trial is GOOD;
  - the best trial is the first GOOD one, replaced by a later GOOD trial
    only when that one lowers R-free by more than r_tolerance (equal R-free
    never moves the search: at tight restraints the model hardly moves and
    R-free is flat, which must not drive the scale up without end);
  - both sides of the first GOOD trial are tried (down, then up), then the
    neighbours of the best trial: down (best/factor) and up (best*factor).
    A neighbour that becomes the best is followed in the same direction
    first; from the third such move in a row the step grows by factor per
    move, so a scale far from the start is reached in a few trials;
  - a gap next to the best trial is halved (geometric mean): towards a
    not-GOOD trial below while the ratio exceeds sqrt(factor) (the bond rmsd
    limit is approached from the GOOD side), towards any tried scale while
    the ratio exceeds factor (a long step went past the minimum);
  - done when neither side of the best trial has a scale left to try;
  - max_trials caps the search; stopped_by_cap() tells when the cap, not a
    rule, ended it.

The choice is the best trial; without a GOOD trial, the last trial
(fallback).
"""
from __future__ import absolute_import, division, print_function
import math
from libtbx import group_args

_EPS = 1.e-9

class controller(object):

  def __init__(self, start_scale=1.0, factor=2.0, max_trials=8,
               r_tolerance=0.001):
    assert factor > 1, factor
    assert max_trials >= 1, max_trials
    assert start_scale > 0, start_scale
    self.start_scale = start_scale
    self.factor = factor
    self.max_trials = max_trials
    self.r_tolerance = r_tolerance
    self.trials = []
    self._incumbent = None # the best trial (module doc)
    self._anchor = None    # the first GOOD trial
    self._run_dir = 0 # direction of the last move that gave a new best
    self._run_len = 0 # such moves in a row, in that direction

  def _best(self):
    """The GOOD trial to keep (None without one): see the module doc."""
    return self._incumbent

  def _candidates(self):
    """(down, up): the scales left to try next to the best trial, or None."""
    f = self.factor
    best = self._best()
    s = best.scale
    step = f
    if(self._run_len >= 3): step = f**(self._run_len - 1)
    below = [t for t in self.trials if t.scale < s*(1 - _EPS)]
    down = None
    if(len(below) == 0):
      down = s / (step if self._run_dir < 0 else f)
    else:
      t = max(below, key=lambda t: t.scale)
      limit = f if t.good else math.sqrt(f)
      if(s / t.scale > limit*(1 + _EPS)): down = math.sqrt(s * t.scale)
    above = [t for t in self.trials if t.scale > s*(1 + _EPS)]
    up = None
    if(len(above) == 0):
      up = s * (step if self._run_dir > 0 else f)
    else:
      t = min(above, key=lambda t: t.scale)
      if(t.scale / s > f*(1 + _EPS)): up = math.sqrt(s * t.scale)
    return down, up

  def _next(self):
    """The scale of the next trial; None: nothing left to try."""
    if(len(self.trials) == 0): return self.start_scale
    if(self._best() is None):
      # no GOOD trial yet: up, accelerating
      return max(t.scale for t in self.trials) * self.factor**len(self.trials)
    down, up = self._candidates()
    # both sides of the first GOOD trial
    a = self._anchor.scale
    if(len([t for t in self.trials if t.scale > a*(1 + _EPS)]) == 0 and
       len([t for t in self.trials if t.scale < a*(1 - _EPS)]) > 0):
      return a * self.factor
    if(self._run_dir > 0 and up is not None): return up
    if(down is not None): return down
    return up

  def scale(self):
    """The scale of the next trial."""
    result = self._next()
    assert result is not None
    return result

  def record(self, scale, good, r_free, sites_cart):
    """The result of the trial at scale."""
    previous = self._best()
    t = group_args(scale=scale, good=good, r_free=r_free,
      sites_cart=sites_cart)
    self.trials.append(t)
    if(good and (previous is None or
                 r_free < previous.r_free - self.r_tolerance)):
      self._incumbent = t
    if(good and self._anchor is None): self._anchor = t
    if(self._best() is t and previous is not None):
      direction = 1 if scale > previous.scale else -1
      if(direction == self._run_dir): self._run_len += 1
      else: self._run_dir, self._run_len = direction, 1
    else:
      self._run_dir, self._run_len = 0, 0

  def _done_by_rule(self):
    return self._next() is None

  def done(self):
    return len(self.trials) >= self.max_trials or self._done_by_rule()

  def stopped_by_cap(self):
    """True when max_trials ended the search before a stopping rule did."""
    return len(self.trials) >= self.max_trials and not self._done_by_rule()

  def choice(self):
    """group_args(scale, r_free, sites_cart, reason): the trial to keep."""
    t = self._best()
    reason = "best"
    if(t is None):
      t = self.trials[-1]
      reason = "fallback"
    return group_args(scale=t.scale, r_free=t.r_free,
      sites_cart=t.sites_cart, reason=reason)
