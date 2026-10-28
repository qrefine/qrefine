from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import os
import iotbx.pdb
import libtbx.load_env
from libtbx.utils import Sorry
from scitbx.array_family import flex
from qrefine import fragment_capping
from qrefine.charges import charges_class

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine, "tests", "unit", "data_files")

def fragment_of(file_name, sel_str):
  pi = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data, file_name))
  ph = pi.construct_hierarchy()
  return ph, ph.atom_selection_cache().selection(sel_str), pi.crystal_symmetry()

def cap_names(r):
  return sorted(a.name.strip()
    for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection))

def check_fragment_atoms(r, ph, sel):
  frag = ph.select(sel)
  kept = r.pdb_hierarchy.atoms().select(r.fragment_selection)
  assert kept.size() == frag.atoms().size()
  assert list(kept.extract_name()) == list(frag.atoms().extract_name())
  assert kept.extract_xyz().max_distance(frag.atoms().extract_xyz()) < 1.e-6

def test_peptide_cut():
  """
  Residues 3-5 cut from the 7-residue peptide: NH2 on N 3, one H on C 5.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  r = c.cap()
  assert (c.n_caps, c.charge) == (2, 0), (c.n_caps, c.charge)
  assert cap_names(r) == ['H2', 'HOXT'], cap_names(r)
  check_fragment_atoms(r, ph, sel)
  for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection):
    parent = {'H2': 'N', 'HOXT': 'C'}[a.name.strip()]
    p = [x for x in a.parent().atoms() if x.name.strip() == parent][0]
    assert 0.95 < a.distance(p) < 1.15, (a.name, a.distance(p))

def test_truncated_side_chain():
  """
  Same fragment with GLN 4 truncated at CD: one more cap, HNE2.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb",
    "resseq 3:5 and not (resseq 4 and (name NE2 or name HE21 or name HE22))")
  c = fragment_capping.capper(ph, sel, cs)
  assert (c.n_caps, c.charge) == (3, 0), (c.n_caps, c.charge)
  assert cap_names(c.cap()) == ['H2', 'HNE2', 'HOXT'], cap_names(c.cap())

def test_pro_cut():
  """
  PRO 3 at the N side of the cut gets one H named H (secondary amine).
  """
  ph, sel, cs = fragment_of("2lvr.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  r = c.cap()
  assert (c.n_caps, c.charge) == (2, 0), (c.n_caps, c.charge)
  assert cap_names(r) == ['H', 'HOXT'], cap_names(r)
  h = [a for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection)
       if a.name.strip() == 'H'][0]
  assert h.parent().resname == 'PRO'

def test_disulfide_cut():
  """
  A CYS cut from its disulfide partner: HG on SG.
  """
  from mmtbx.hydrogens.tst_capping import pdb_disulfide, model_with_h
  model = model_with_h(pdb_disulfide)
  ph = model.get_hierarchy()
  sel = ph.atom_selection_cache().selection("resseq 232")
  c = fragment_capping.capper(ph, sel, model.crystal_symmetry())
  assert (c.n_caps, c.charge) == (1, 0), (c.n_caps, c.charge)
  assert cap_names(c.cap()) == ['HG'], cap_names(c.cap())

def test_whole_model_no_caps():
  """
  The whole complete model: nothing to cap, the charge is qr.charges'.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "all")
  c = fragment_capping.capper(ph, sel, cs)
  assert c.n_caps == 0, c.n_caps
  expected = charges_class(pdb_hierarchy=ph.deep_copy(),
    crystal_symmetry=cs).get_total_charge()
  assert c.charge == expected, (c.charge, expected)
  r = c.cap()
  assert r.fragment_selection.all_eq(True)

def test_cap_moves_with_sites():
  """
  Re-capping a rigidly shifted fragment shifts the caps with it.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  r0 = c.cap()
  shift = (1.0, 2.0, 3.0)
  r1 = c.cap(sites_cart=c.fragment_hierarchy.atoms().extract_xyz() + shift)
  assert cap_names(r0) == cap_names(r1)
  d = r1.pdb_hierarchy.atoms().extract_xyz() - \
      r0.pdb_hierarchy.atoms().extract_xyz()
  assert (d - shift).max_distance(flex.vec3_double(d.size(), (0,0,0))) < 1.e-6

def test_cap_twice():
  """
  Caps do not accumulate; the cached fragment is not changed.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  c = fragment_capping.capper(ph, sel, cs)
  n = c.fragment_hierarchy.atoms().size()
  a, b = c.cap(), c.cap()
  assert a.pdb_hierarchy.atoms().size() == b.pdb_hierarchy.atoms().size() == n+2
  assert c.fragment_hierarchy.atoms().size() == n

def test_altloc_sorry():
  """
  Two alternate conformations are refused; a single conformer that still
  carries its altloc label (what qrefine's per-conformer hierarchies look
  like) is capped.
  """
  from mmtbx.hydrogens.tst_capping import pdb_ser_altloc, model_with_h
  model = model_with_h(pdb_ser_altloc)
  ph = model.get_hierarchy()
  cs = model.crystal_symmetry()
  try:
    fragment_capping.capper(ph, flex.bool(ph.atoms().size(), True), cs)
  except Sorry as e:
    assert 'alternate conformations' in str(e), str(e)
  else:
    raise AssertionError('no Sorry')
  one = ph.select(ph.atom_selection_cache().selection("altloc A or altloc ' '"),
                  copy_atoms=True)
  assert sorted(a for a in one.altloc_indices().keys() if a.strip()) == ['A']
  c = fragment_capping.capper(one, flex.bool(one.atoms().size(), True), cs)
  assert c.charge is not None

def test_parent_untouched():
  """
  The parent hierarchy keeps its i_seqs, and the cached fragment does not
  follow later coordinate changes of the parent.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "resseq 3:5")
  i_seqs = list(ph.atoms().extract_i_seq())
  xyz = ph.atoms().extract_xyz()
  c = fragment_capping.capper(ph, sel, cs)
  assert list(ph.atoms().extract_i_seq()) == i_seqs
  ph.atoms().set_xyz(xyz + (5.0, 0.0, 0.0))
  assert c.fragment_hierarchy.atoms().extract_xyz().max_distance(
    xyz.select(sel)) < 1.e-6

def single_conformer_1bdw():
  pi = iotbx.pdb.input(file_name=os.path.join(qr_unit_tests_data,
    "1bdw_ala_refine_001_complete_minimized.pdb_modified.pdb"))
  ph0 = pi.construct_hierarchy()
  ph = ph0.select(ph0.atom_selection_cache().selection(
    "altloc A or altloc ' '"), copy_atoms=True)
  for ag in ph.atom_groups(): ag.altloc = ''
  ph.atoms().reset_i_seq()
  ph.atoms_reset_serial()
  return ph, pi.crystal_symmetry()

def test_sphere_boundary_fragment():
  """
  A fragment with a symmetry-copy residue at the edge of the super sphere:
  its peptide partner lies outside the sphere, in the super cell. Cut bonds
  must be found against the super cell (fragments.make_capper), where the N
  gets its H2; against the sphere the N looks like an unsatisfied valence.
  """
  from qrefine.fragment import fragments
  from mmtbx.hydrogens import capping
  ph, cs = single_conformer_1bdw()
  fm = fragments(pdb_hierarchy=ph, crystal_symmetry=cs, qm_engine_name="mopac",
    maxnum_residues_in_cluster=5, select_within_radius=10.0,
    cap_fragments=False, make_working_folder=True)
  SCM = fm.expansion.SCM
  boundary = capping.cut_bonds_by_distance(SCM.super_cell_hierarchy,
                                           SCM.keep_selection)
  satoms = fm.pdb_hierarchy_super.atoms()
  n_boundary = [k for k in boundary if satoms[k].name.strip() == 'N']
  assert n_boundary, "no sphere-boundary N in this expansion"
  tested = 0
  for sel in fm.fragment_super_selections:
    hits = [k for k in n_boundary if sel[k]]
    if not hits: continue
    try:
      fragment_capping.capper(fm.pdb_hierarchy_super, sel, fm.expansion.cs_box)
    except Sorry as e:
      assert 'not QM-ready' in str(e), str(e)
    else:
      raise AssertionError('parent=sphere should miss the cut bond')
    c = fm.make_capper(sel)
    assert c.fragment_hierarchy.is_similar_hierarchy(
      fm.pdb_hierarchy_super.select(sel))
    r = c.cap()
    capped = [a for a in r.pdb_hierarchy.atoms().select(~r.fragment_selection)]
    for k in hits:
      n_atom = satoms[k]
      h2 = [a for a in capped if a.name.strip() == 'H2' and
            a.parent().parent().resseq == n_atom.parent().parent().resseq and
            a.parent().parent().parent().id == n_atom.parent().parent().parent().id]
      assert len(h2) == 1, (n_atom.id_str(), [a.id_str() for a in capped])
    tested += 1
    break
  assert tested == 1

def test_write_fragment_pdbs_without_cappers():
  """
  write_cluster_and_fragments_pdbs works for cctbx-restraints clustering,
  where the fragments carry no cappers (cap_fragments=False).
  """
  from qrefine.fragment import fragments, write_cluster_and_fragments_pdbs
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb", "all")
  fm = fragments(pdb_hierarchy=ph, crystal_symmetry=cs, qm_engine_name="mopac",
    maxnum_residues_in_cluster=2, cap_fragments=False, save_clusters=True,
    make_working_folder=True)
  write_cluster_and_fragments_pdbs(fragments=fm.get_fragment_extracts(),
                                   directory="frag_pdbs_nocap")
  assert os.path.exists(os.path.join("frag_pdbs_nocap", "0_frag.pdb"))

def test_not_qmready_sorry():
  """
  A backbone N short of bonds that is not a cut: the model was not prepared
  with reduce2 qmready=True.
  """
  ph, sel, cs = fragment_of("m2_complete_box_large.pdb",
    "not (resseq 1 and (name H1 or name H2 or name H3))")
  frag = ph.select(sel)
  try:
    fragment_capping.capper(frag, flex.bool(frag.atoms().size(), True), cs)
  except Sorry as e:
    assert 'qmready' in str(e), str(e)
  else:
    raise AssertionError('no Sorry')

def test_nucleic_cut():
  """
  Cut through a phosphodiester bond: the O3' side gets HO3', the P side is
  hypervalent and refused.
  """
  from mmtbx.hydrogens.tst_capping import pdb_dna, model_with_h
  model = model_with_h(pdb_dna)
  ph = model.get_hierarchy()
  cs = model.crystal_symmetry()
  resseqs = [rg.resseq.strip() for rg in ph.residue_groups()]
  first, last = resseqs[0], resseqs[-1]
  c = fragment_capping.capper(ph,
    ph.atom_selection_cache().selection("resseq %s" % first), cs)
  assert "HO3'" in cap_names(c.cap()), cap_names(c.cap())
  try:
    fragment_capping.capper(ph,
      ph.atom_selection_cache().selection("resseq %s" % last), cs)
  except Sorry as e:
    assert ' P ' in str(e) or 'P ' in str(e), str(e)
  else:
    raise AssertionError('no Sorry')

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Exercise fragment capping (fragment_capping.capper).
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  for f in [test_peptide_cut, test_truncated_side_chain, test_pro_cut,
            test_disulfide_cut, test_whole_model_no_caps,
            test_cap_moves_with_sites, test_cap_twice, test_altloc_sorry,
            test_not_qmready_sorry, test_nucleic_cut,
            test_parent_untouched, test_sphere_boundary_fragment,
            test_write_fragment_pdbs_without_cappers]:
    f()
    print(f.__name__, "OK")

if(__name__=='__main__'):
  run()
