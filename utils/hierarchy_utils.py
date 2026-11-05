from __future__ import print_function
import mmtbx.model
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

def model_from_hierarchy(pdb_hierarchy, crystal_symmetry, cif_objects=None, log=None):
  """
  Processed mmtbx.model.manager (restraints built) of pdb_hierarchy with the
  interpretation settings of qr.refine (neutron distances, no CDL, atoms
  not sorted, no amino-acid flips, hydrogens not corrected).
  """
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
