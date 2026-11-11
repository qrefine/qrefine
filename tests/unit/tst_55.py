from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import mmtbx.f_model
import qrefine.qr as qr
import mmtbx.hydrogens.xh_mapper
from libtbx.test_utils import approx_equal
from mmtbx.hydrogens.tst_xh_mapper import pdb_str, get_model, xh_lengths, dist
from qrefine.tests.unit import run_qrefine_in, run_fmodel, fast_args

# fast_args, with 20 refinement iterations: H must reach their nuclear X-H
# lengths (write_nuclear_h_positions, neutron table)
common_args = ["in.pdb", "data.mtz", "restraints=cctbx", "expansion=False",
  "clustering=False"] + [a for a in fast_args
  if not a.startswith("max_iterations_refine=")] + ["max_iterations_refine=20"]

def atom_id(atom):
  """Identifier of an atom that does not depend on the atom order."""
  ag = atom.parent()
  return (ag.parent().parent().id, ag.parent().resseq.strip(), ag.resname,
    ag.altloc, atom.name.strip())

def targets(use_neutron_distances):
  """
  {H id: (parent id, ideal X-H length, parent element)} of the fixture
  processed with X-ray (False) or nuclear (True) X-H lengths.
  """
  m = get_model(pdb_str, "n_gaussian",
    use_neutron_distances=use_neutron_distances)
  atoms = m.get_hierarchy().atoms()
  result = {}
  for i_h, (i_x, d) in xh_lengths(m, m.get_sites_cart()).items():
    result[atom_id(atoms[i_h])] = (atom_id(atoms[i_x]), d,
      atoms[i_x].element.strip().upper())
  assert len(result) > 10, result
  return result

def lengths(model, targets):
  """{H id: X-H length in model} for the H of targets."""
  atoms = model.get_hierarchy().atoms()
  index = dict((atom_id(a), i) for i, a in enumerate(atoms))
  sites = model.get_sites_cart()
  return dict((h, dist(sites, index[h], index[x]))
    for h, (x, d, e) in targets.items())

def run_xray():
  return run_qrefine_in(common_args + ["scattering_table=n_gaussian"],
    "in", folder="xray")

def run_nuclear():
  return run_qrefine_in(common_args + ["scattering_table=n_gaussian",
    "write_nuclear_h_positions=True"], "in", folder="nuclear")

def test_output_lengths(out, nuc):
  """
  Default output: H at the X-ray X-H lengths (exactly, up to the PDB
  coordinate rounding); write_nuclear_h_positions=True: H at refined nuclear
  positions, longer than the X-ray targets.
  """
  tx = targets(use_neutron_distances=False)
  assert "X-H mapper: active" in out.std_out_text
  assert "Starting r-factors (H at X-ray positions)" in out.std_out_text
  for h, l in lengths(out.model, tx).items():
    assert approx_equal(l, tx[h][1], 1.e-3), (h, l, tx[h])
  assert "X-H mapper: active" in nuc.std_out_text
  for h, l in lengths(nuc.model, tx).items():
    x, d_x, e = tx[h]
    if(e == "C"): assert l > d_x + 0.05, (h, l, d_x)
    else:         assert l > d_x + 0.03, (h, l, d_x)
  # the default output is exactly the mapper's shrink of the nuclear result
  with open(nuc.pdb) as fi:
    nm = get_model(fi.read(), "n_gaussian", use_neutron_distances=True)
  assert nm.size() == out.model.size()
  expected = mmtbx.hydrogens.xh_mapper.mapper(nm).shrink(nm.get_sites_cart())
  diff = max(max(abs(a-b) for a, b in zip(s1, s2)) for s1, s2 in
    zip(expected, out.model.get_sites_cart()))
  assert diff < 2.e-3, diff

def rebuilt_r_factors(out, r):
  """r_work, r_free of the written model, fmodel rebuilt as qr.refine does."""
  m = out.model
  m.setup_scattering_dictionaries(scattering_table="n_gaussian",
    d_min=r.f_obs.d_min())
  fmodel = mmtbx.f_model.manager(f_obs=r.f_obs, r_free_flags=r.flags,
    xray_structure=m.get_xray_structure())
  fmodel.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
  return fmodel.r_work(), fmodel.r_free()

def reported_r_factors(out):
  """The last r_work, r_free qr.refine reports."""
  key = "Best Rwork, Rfree (after refinement):"
  lines = [l for l in out.std_out_text.splitlines() if l.startswith(key)]
  assert len(lines) > 0
  return [float(v) for v in lines[-1].split(":")[1].split()]

def test_r_factors_consistent(out, nuc, r):
  """
  R-factors of the written model (X-ray H positions) equal the last ones
  qr.refine reports: the output is what fmodel saw. The tolerance covers
  the PDB rounding of coordinates and B (2.5e-4 here) and the printed
  precision (5e-5). The nuclear output of the same refinement
  (write_nuclear_h_positions=True) is not what fmodel saw: its R factors
  differ from the reported ones by several times more than the X-ray
  output's do.
  """
  eps = 6.e-4
  diffs_xray = [abs(rb - rp) for rb, rp in
    zip(rebuilt_r_factors(out, r), reported_r_factors(out))]
  assert max(diffs_xray) < eps, diffs_xray
  diffs_nuclear = [abs(rb - rp) for rb, rp in
    zip(rebuilt_r_factors(nuc, r), reported_r_factors(nuc))]
  assert min(diffs_nuclear) > 3 * max(diffs_xray), (diffs_xray, diffs_nuclear)

def test_default_table():
  """The default scattering table is n_gaussian, as in phenix.refine."""
  assert qr.get_default_params().scattering_table == "n_gaussian"
  out = run_qrefine_in(common_args, "in", folder="default")
  assert "scattering_table = wk1995 it1992 *n_gaussian" in out.std_out_text
  assert "X-H mapper: active" in out.std_out_text

def test_mapper_independent_of_table():
  """The X-ray X-H lengths do not depend on the X-ray table (wk1995)."""
  wk = run_qrefine_in(common_args + ["scattering_table=wk1995"], "in",
    folder="wk1995")
  assert "X-H mapper: active" in wk.std_out_text
  tx = targets(use_neutron_distances=False)
  for h, l in lengths(wk.model, tx).items():
    assert approx_equal(l, tx[h][1], 1.e-3), (h, l, tx[h])

def test_neutron_nuclear():
  """Neutron table: identity mapper, H at nuclear X-H lengths."""
  out = run_qrefine_in(common_args + ["scattering_table=neutron"], "in",
    folder="neutron")
  assert "X-H mapper: identity (neutron scattering table)" in out.std_out_text
  assert "X-H mapper: active" not in out.std_out_text
  tn = targets(use_neutron_distances=True)
  for h, l in lengths(out.model, tn).items():
    assert approx_equal(l, tn[h][1], 0.02), (h, l, tn[h])

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine with the X-H mapper on the tst_xh_mapper model (X-ray data at
  1.5 A): X-ray tables write H at X-ray X-H lengths (nuclear with
  write_nuclear_h_positions=True), the written model reproduces the final
  R-factors, the neutron table maps nothing.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  # The fixture in a 25 A box instead of 50 A: 1.5 A data (H positions
  # visible in R, test_r_factors_consistent) with 8x fewer reflections.
  cryst1 = "CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1"
  assert pdb_str.count(cryst1) == 1
  with open("in.pdb", "w") as fo:
    fo.write(pdb_str.replace(cryst1,
      "CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1"))
  r = run_fmodel(args=["in.pdb", "high_res=1.5"], prefix="data")
  out = run_xray()
  nuc = run_nuclear()
  test_output_lengths(out, nuc)
  test_r_factors_consistent(out, nuc, r)
  test_default_table()
  test_mapper_independent_of_table()
  test_neutron_nuclear()

if(__name__ == "__main__"):
  run()
  print("OK")
