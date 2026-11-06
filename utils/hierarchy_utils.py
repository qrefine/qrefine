from __future__ import print_function
import iotbx.pdb
import mmtbx.model
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from libtbx.utils import null_out

_mon_lib_srv = None

def get_mon_lib_srv():
  """
  The one monomer library server of qrefine: reading the library takes a
  few tenths of a second, and every fragment has its own model.
  """
  global _mon_lib_srv
  if _mon_lib_srv is None:
    _mon_lib_srv = mmtbx.monomer_library.server.server()
  return _mon_lib_srv

def write_hierarchy(pdb_hierarchy, file_name, **kw):
  """
  Write pdb_hierarchy into file_name in PDB format, or in mmCIF with the
  extension changed to .cif when its identifiers do not fit the PDB format
  (hierarchy.pdb_or_mmcif_string_info decides both). kw: crystal_symmetry
  and the other keywords accepted by both as_pdb_string and
  as_mmcif_string. Returns the name of the file written. Same text as
  hierarchy.write_pdb_or_mmcif_file, without the DataManager it constructs
  at every call (~30 ms; scratch models are written at every evaluation).
  """
  info = pdb_hierarchy.pdb_or_mmcif_string_info(
    target_filename = file_name,
    target_format   = "pdb",
    write_file      = False,
    **kw)
  with open(info.file_name, "w") as fo:
    fo.write(info.pdb_string)
  return info.file_name

def bond_edits_phil(bonds):
  """
  geometry_restraints.edits PHIL text adding one bond per (atom_selection_1,
  atom_selection_2, distance, sigma) of bonds.
  """
  lines = ["geometry_restraints.edits {"]
  for s1, s2, d, sigma in bonds:
    lines.extend([
      "  bond {",
      "    action = add",
      "    atom_selection_1 = %s" % s1,
      "    atom_selection_2 = %s" % s2,
      "    distance_ideal = %.3f" % d,
      "    sigma = %.3f" % sigma,
      "  }"])
  lines.append("}")
  return "\n".join(lines)

def bond_edits_params(bonds):
  """
  The pdb_interpretation parameters (mmtbx.model.manager scope) with the
  bonds of bond_edits_phil in geometry_restraints.edits.
  """
  import iotbx.phil
  scope = mmtbx.model.manager.get_default_pdb_interpretation_scope()
  return scope.fetch(iotbx.phil.parse(bond_edits_phil(bonds))).extract()

_AMINO_ACID_CLASSES = ("common_amino_acid", "modified_amino_acid",
  "d_amino_acid")

def _atom_selection(atom):
  rg = atom.parent().parent()
  return 'chain "%s" and resid %s and name %s' % (
    rg.parent().id, rg.resid().strip(), atom.name.strip())

def peptide_closures(hierarchy, cutoff=1.9):
  """
  Backbone C-N bonds between amino-acid residues of different chains of
  hierarchy (Cartesian distances, no symmetry) as [(atom_selection_1,
  atom_selection_2, distance, sigma)] for model_from_hierarchy(bonds=...)
  or bond_edits_params, with the TRANS link values (1.329 A, sigma 0.014).
  pdb_interpretation links the first and the last residue of one chain
  (cyclic peptide), never two chains; the super sphere of a macrocycle closed
  through a symmetry mate needs the bond between the master and a copy (the
  bond only, as the symmetry closure of the asymmetric unit).
  """
  from cctbx import crystal
  get_class = iotbx.pdb.common_residue_names_get_class
  atoms = hierarchy.atoms()
  candidates = flex.size_t()
  for i, a in enumerate(atoms):
    if(a.element_is_hydrogen()): continue
    if(a.name.strip() not in ("C", "N")): continue
    if(get_class(a.parent().resname) not in _AMINO_ACID_CLASSES): continue
    candidates.append(i)
  result = []
  for i, j, rt_mx, d in crystal.pairs_within(atoms.extract_xyz(), cutoff,
      i_seqs=candidates):
    if(i >= j): continue
    ai, aj = atoms[i], atoms[j]
    if(aj.element_is_hydrogen() or aj.name.strip() == ai.name.strip()):
      continue
    if(aj.name.strip() not in ("C", "N")): continue
    if(get_class(aj.parent().resname) not in _AMINO_ACID_CLASSES): continue
    if(ai.parent().parent().parent().id == aj.parent().parent().parent().id):
      continue
    result.append((_atom_selection(ai), _atom_selection(aj), 1.329, 0.014))
  return result

def model_from_hierarchy(pdb_hierarchy, crystal_symmetry, cif_objects=None,
                         log=None, bonds=None):
  """
  Processed mmtbx.model.manager (restraints built) of pdb_hierarchy with the
  interpretation settings of qr.refine (neutron distances, no CDL, atoms
  not sorted, no amino-acid flips, hydrogens not corrected). bonds: extra
  bonds as (atom_selection_1, atom_selection_2, distance, sigma)
  (bond_edits_params), e.g. the peptide closures between the chains of a
  super sphere.
  """
  if(bonds):
    params = bond_edits_params(bonds)
  else:
    params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  params.pdb_interpretation.flip_symmetric_amino_acids = False
  params.pdb_interpretation.correct_hydrogens=False
  if log is None: log = null_out()
  model = mmtbx.model.manager(
    model_input       = None,
    pdb_hierarchy     = pdb_hierarchy,
    crystal_symmetry  = crystal_symmetry,
    restraint_objects = cif_objects,
    log               = log)
  model.set_mon_lib_srv(get_mon_lib_srv())
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params,)
  return model
