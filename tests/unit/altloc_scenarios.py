"""
Made-up altloc scenarios for the qr.refine altloc tests (tst_49 onward),
built in process from the seven-residue peptide with H of
tst_00.pdb_str_good: side-chain conformers are chi1 rotations of the
original, backbone conformers are rotations of one residue's backbone about
the C(i-1)...N(i+1) axis.

Convention: the editing helpers (side_chain_conformers, label_atoms,
backbone_conformers) edit the hierarchy of the model in place and reset its
i_seq and serial numbers; they return nothing. The X-ray structure of that
model manager is then stale, so every helper here that needs one (shaken,
r_work, f_calc_rel_diff) first rebuilds a manager from (a copy of) the
hierarchy; the others read the hierarchy only.
"""
from __future__ import absolute_import, division, print_function
import math
import random
import iotbx.pdb
import iotbx.pdb.hierarchy
import mmtbx.model
import mmtbx.f_model
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.utils import null_out
from qrefine.tests.unit import tst_00

not_beyond_cb = set(["N", "CA", "C", "O", "CB", "OXT", "H", "HA", "HB2",
  "HB3", "H1", "H2", "H3"])
backbone = ["N", "CA", "C", "O", "H", "HA", "HA2", "HA3"]

def _root(residue_group):
  return residue_group.parent().parent().parent()

def _reset(residue_group):
  root = _root(residue_group)
  root.atoms().reset_i_seq()
  root.atoms_reset_serial()

def _rebuilt(model):
  """
  New manager of a deep copy of the hierarchy of model.
  """
  ph = model.get_hierarchy().deep_copy()
  ph.atoms().reset_i_seq()
  ph.atoms_reset_serial()
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=model.crystal_symmetry(), log=null_out())

def _rotated(xyz, origin, axis, degrees):
  """
  xyz rotated by degrees about the axis through origin (both matrix.col).
  """
  if(degrees == 0): return xyz
  r = (matrix.col(xyz) - origin).rotate_around_origin(
    axis=axis.normalize(), angle=math.radians(degrees))
  return (r + origin).elems

def _atom_group(residue_group, label):
  """
  Atom group of residue_group with altloc label, appended if absent (resname
  of the first atom group).
  """
  for ag in residue_group.atom_groups():
    if(ag.altloc == label): return ag
  ag = iotbx.pdb.hierarchy.atom_group(altloc=label,
    resname=residue_group.atom_groups()[0].resname)
  residue_group.append_atom_group(ag)
  return ag

def _find(residue_group, name):
  """
  Atom name of residue_group: the blank one, else that of the first label.
  """
  found = []
  for ag in residue_group.atom_groups():
    for a in ag.atoms():
      if(a.name.strip() == name): found.append((ag.altloc.strip() != "", a))
  assert len(found) > 0, name
  return sorted(found, key=lambda x: x[0])[0][1]

def _move(atoms, residue_group, label, occupancy):
  """
  Move blank atoms of residue_group into the atom group of label.
  """
  ag_to = _atom_group(residue_group, label)
  for a in atoms:
    a.parent().remove_atom(a)
    a.occ = occupancy
    ag_to.append_atom(a)

def base_model():
  """
  Unprocessed model of tst_00.pdb_str_good (P 1 21 1; GLY 1, ASN 2, ASN 3,
  GLN 4, GLN 5, ASN 6, TYR 7 with H and OXT, 7 waters).
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=tst_00.pdb_str_good)
  return mmtbx.model.manager(model_input=pdb_inp, log=null_out())

def residue(model, resseq):
  """
  The residue group with resseq (int) of model.
  """
  rgs = [rg for rg in model.get_hierarchy().residue_groups()
    if rg.resseq_as_int() == resseq]
  assert len(rgs) == 1, resseq
  return rgs[0]

def side_chain_conformers(rg, labels, occupancies, chi1_degrees):
  """
  In place: the blank atoms of rg beyond CB (CG onward and their H) become
  one atom group per label, rotated about the CA-CB axis by chi1_degrees
  (0 keeps the original), with occupancies. labels[0] reuses the existing
  atoms; the other labels get rotated copies. Other atoms stay blank.
  """
  assert len(labels) == len(occupancies) == len(chi1_degrees)
  ca = matrix.col(_find(rg, "CA").xyz)
  cb = matrix.col(_find(rg, "CB").xyz)
  atoms = [a for ag in rg.atom_groups() if ag.altloc.strip() == ""
    for a in ag.atoms() if a.name.strip() not in not_beyond_cb]
  assert len(atoms) > 0
  copies = []
  for label, occupancy, degrees in zip(labels[1:], occupancies[1:],
                                       chi1_degrees[1:]):
    atoms_c = []
    for a in atoms:
      a_c = a.detached_copy()
      a_c.xyz = _rotated(a.xyz, cb, cb - ca, degrees)
      atoms_c.append(a_c)
    copies.append((label, occupancy, atoms_c))
  for a in atoms:
    a.xyz = _rotated(a.xyz, cb, cb - ca, chi1_degrees[0])
  _move(atoms, rg, labels[0], occupancies[0])
  for label, occupancy, atoms_c in copies:
    ag = _atom_group(rg, label)
    for a_c in atoms_c:
      a_c.occ = occupancy
      ag.append_atom(a_c)
  _reset(rg)

def label_atoms(rg, names, label, occupancy):
  """
  In place: move the blank atoms of rg named in names into the atom group of
  label (appended if absent) and set their occupancy.
  """
  names = set(names)
  atoms = [a for ag in rg.atom_groups() if ag.altloc.strip() == ""
    for a in ag.atoms() if a.name.strip() in names]
  assert len(atoms) == len(names), names
  _move(atoms, rg, label, occupancy)
  _reset(rg)

def backbone_conformers(model, resseq, label, occupancy, degrees, names=None):
  """
  In place: add a copy of the atoms names (default N, CA, C, O, H, HA, HA2,
  HA3; those present) of residue resseq under label, rotated by degrees
  about the axis C(resseq-1) -> N(resseq+1), so both peptide bond lengths
  are kept, with occupancy. Blank originals are moved to the first existing
  label of the residue (or "A") with occupancy 1 - occupancy; labelled
  originals are copied from that first label and keep their occupancies.
  """
  names = backbone if names is None else names
  rg = residue(model, resseq)
  c_prev = matrix.col(_find(residue(model, resseq-1), "C").xyz)
  n_next = matrix.col(_find(residue(model, resseq+1), "N").xyz)
  labels = sorted([ag.altloc for ag in rg.atom_groups()
    if ag.altloc.strip() != ""])
  first = (labels + ["A"])[0]
  assert label != first
  blank = [a for ag in rg.atom_groups() if ag.altloc.strip() == ""
    for a in ag.atoms() if a.name.strip() in names]
  if(len(blank) > 0):
    _move(blank, rg, first, 1 - occupancy)
  originals = [a for a in _atom_group(rg, first).atoms()
    if a.name.strip() in names]
  assert len(originals) > 0
  ag = _atom_group(rg, label)
  for a in originals:
    a_c = a.detached_copy()
    a_c.xyz = _rotated(a.xyz, c_prev, n_next - c_prev, degrees)
    a_c.occ = occupancy
    ag.append_atom(a_c)
  _reset(rg)

def write(model, file_name):
  """
  Write the hierarchy of model with its crystal symmetry as PDB.
  """
  model.get_hierarchy().write_pdb_file(file_name=file_name,
    crystal_symmetry=model.crystal_symmetry())

def shaken(model, mean_distance=0.3):
  """
  New model: a copy of model with all sites shaken (fixed random seeds).
  """
  result = _rebuilt(model)
  xrs = result.get_xray_structure()
  flex.set_random_seed(0)
  random.seed(0)
  xrs.shake_sites_in_place(mean_distance=mean_distance)
  result.set_sites_cart(xrs.sites_cart())
  return result

def r_work(f_obs, flags, model):
  """
  R-work of model (with bulk solvent and scaling, update_all_scales) against
  f_obs, with the qr.refine default scattering table.
  """
  from qrefine import qr
  m = _rebuilt(model)
  m.setup_scattering_dictionaries(
    scattering_table=qr.get_default_params().scattering_table)
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = flags,
    xray_structure = m.get_xray_structure())
  fmodel.update_all_scales()
  return fmodel.r_work()

def f_calc_rel_diff(model_1, model_2, d_min=1.2):
  """
  max |F1 - F2| / max |F1| of the complex structure factors (direct
  summation) of the two models on identical index sets.
  """
  fs = []
  for m in [model_1, model_2]:
    xrs = _rebuilt(m).get_xray_structure()
    fs.append(xrs.structure_factors(d_min=d_min,
      algorithm="direct").f_calc())
  f1, f2 = fs
  assert f1.indices().all_eq(f2.indices())
  return flex.max(flex.abs(f1.data() - f2.data())) / \
    flex.max(flex.abs(f1.data()))

def labels_and_occ(model):
  """
  {(resseq, atom name, altloc): occupancy} of the non-hydrogen atoms (H
  placement may rename or reorder H).
  """
  result = {}
  for a in model.get_hierarchy().atoms():
    if(a.element.strip().upper() in ("H", "D")): continue
    ag = a.parent()
    key = (ag.parent().resseq_as_int(), a.name.strip(), ag.altloc)
    assert key not in result, key
    result[key] = a.occ
  return result

def assert_scattering_neutral(model_in, model_done, copies_expected=True):
  """
  Completion is scattering neutral: structure factors of model_done equal
  those of model_in (f_calc_rel_diff < 1e-6). With copies_expected,
  model_done must have more atoms than model_in, so the check is not
  vacuous.
  """
  d = f_calc_rel_diff(model_in, model_done)
  assert d < 1.e-6, d
  if(copies_expected):
    n_in   = model_in.get_hierarchy().atoms_size()
    n_done = model_done.get_hierarchy().atoms_size()
    assert n_done > n_in, (n_in, n_done)
