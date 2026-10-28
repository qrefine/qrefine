"""
Capping of QM fragments for qr.refine with mmtbx.hydrogens.capping, the
machinery of mmtbx.reduce2 capping=True.

A fragment is parent_hierarchy.select(selection): whole residues cut out of
the super sphere (clusters and their buffers) or the super sphere cut out of
the super cell (expansion). Every bond the cut severed gets one H on the
inside atom (capping.cut_bonds_by_distance, rule 9 of capping.py: NH2 on a
backbone N, one H on a carbonyl C, HG on a cut SG, HO3' on a cut O3'), and
every atom with an unsatisfied valence by the dictionary (a truncated side
chain) gets its H as reduce2 would. Nothing else is protonated: the master
model is expected to come from reduce2 qmready=True.

The fragment is interpreted once (pdb_interpretation, for the bond table,
the dangling-atom records and the total charge of the capped fragment,
computed before the caps exist as electrons.run(extra_sigma_bonds=)); the
records are cached, and cap() re-places the caps on new coordinates on an
unprocessed copy, without interpreting again. The capped hierarchy keeps
the caps inside their residues; fragment_selection says which atoms are the
fragment's. A capped hierarchy cannot be interpreted (invented cap names).
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx import group_args
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex
from mmtbx.hydrogens import capping
from qrefine.utils.hierarchy_utils import model_from_hierarchy

AMINO_ACID_CLASSES = ('common_amino_acid', 'modified_amino_acid',
                      'd_amino_acid')

class capper(object):
  """
  Caps the QM fragment parent_hierarchy.select(selection).

  Attributes: fragment_hierarchy (the uncapped fragment at the setup
  coordinates, i_seq reset, never changed), crystal_symmetry, charge (total
  charge of the capped fragment), n_caps, lines (one report line per capped
  atom), dangling (the cached find_dangling_atoms result).

  Sorry: a fragment with two or more alternate conformations; a residue without a
  dictionary (supply its restraint CIF); a backbone N short of a bond that
  is no cut (the model is not QM-ready); an internal gap, a hypervalent P or
  S, a cap clash, an atom with no bonded neighbour (from capping.py).
  """
  def __init__(self, parent_hierarchy, selection, crystal_symmetry,
               cif_objects=None, log=None):
    self.crystal_symmetry = crystal_symmetry
    # a copy: the parent keeps its i_seqs and later coordinate changes
    self.fragment_hierarchy = parent_hierarchy.select(selection,
                                                      copy_atoms=True)
    self.fragment_hierarchy.atoms().reset_i_seq()
    # one conformer may keep its altloc label (qrefine's per-conformer
    # hierarchies do); two conformers cannot be charged or capped as one
    altlocs = [a for a in self.fragment_hierarchy.altloc_indices().keys()
               if a.strip()]
    if len(altlocs) > 1:
      raise Sorry('Capping: the fragment has alternate conformations (%s); '
        'QM fragments must be single-conformer' % ', '.join(sorted(altlocs)))
    cut = capping.cut_bonds_by_distance(parent_hierarchy, selection)
    self._model = model_from_hierarchy(
      pdb_hierarchy    = self.fragment_hierarchy.deep_copy(),
      crystal_symmetry = crystal_symmetry,
      cif_objects      = cif_objects,
      log              = log)
    self.dangling = capping.find_dangling_atoms(self._model, cut_bonds=cut,
                                                log=log)
    if self.dangling.unknown_residues:
      raise Sorry('Capping: no dictionary for %s in a QM fragment; supply '
        'its restraint CIF to qr.refine' %
        ', '.join(self.dangling.unknown_residues))
    self._check_backbone_n(cut)
    r = capping.cap_model(self._model, dangling=self.dangling, log=log)
    self.charge = r.total_charge
    self.n_caps = r.capped.n_added
    self.lines  = r.capped.lines
    self._mon_lib_srv = self._model.get_mon_lib_srv()

  def _check_backbone_n(self, cut):
    # An amino-acid backbone N with fewer than 3 bonds counting the cut: the
    # dictionary pass cannot see it (it misses no heavy partner) and the
    # charge would silently be off by one. The master is expected QM-ready.
    get_class = iotbx.pdb.common_residue_names_get_class
    table = self.dangling.bond_table
    for ag in self._model.get_hierarchy().atom_groups():
      if get_class(name=ag.resname) not in AMINO_ACID_CLASSES: continue
      for atom in ag.atoms():
        if atom.name.strip() != 'N': continue
        n_bonds = len(table.get(atom.i_seq, ())) + cut.get(atom.i_seq, 0)
        if n_bonds < 3:
          raise Sorry('Capping: %s has an unsatisfied valence that is not a '
            'fragment cut; the model is not QM-ready, run mmtbx.reduce2 '
            'qmready=True first' % atom.id_str())

  def cap(self, sites_cart=None):
    """
    The fragment at sites_cart (its current coordinates, in fragment order;
    None: the setup coordinates) with the caps placed on them.
    Returns group_args(pdb_hierarchy, fragment_selection): the capped
    hierarchy (caps inside their residues, i_seq and serial reset) and a
    flex.bool over its atoms, True for the fragment atoms, in fragment order.
    """
    h = self.fragment_hierarchy.deep_copy()
    if sites_cart is not None:
      assert sites_cart.size() == h.atoms().size(), (
        sites_cart.size(), h.atoms().size())
      h.atoms().set_xyz(sites_cart)
    for a in h.atoms(): a.tmp = 1
    m = mmtbx.model.manager(
      model_input      = None,
      pdb_hierarchy    = h,
      crystal_symmetry = self.crystal_symmetry,
      log              = null_out())
    capping.cap_model(m, dangling=self.dangling, calculate_charge=False,
                      mon_lib_srv=self._mon_lib_srv)
    h = m.get_hierarchy()
    fragment_selection = h.atoms().extract_tmp_as_size_t() == 1
    return group_args(
      pdb_hierarchy      = h,
      fragment_selection = fragment_selection)
