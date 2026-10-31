from __future__ import print_function
import mmtbx.model
from libtbx.utils import null_out

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
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params,)
  return model
