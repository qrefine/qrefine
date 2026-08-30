from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.finalise
import sys, time
from qrefine import finalise
import iotbx
import mmtbx
from mmtbx import utils
from libtbx.utils import Sorry
from libtbx.utils import null_out
from pathlib import Path

log = sys.stdout

legend = """\
Finalise a model before quantum refinement
"""

master_params_str = """
model_file_name = None
  .type = path
  .short_caption = Model file
  .multiple = False
  .help = Model file name to use as input for obtaining a complete model
  .style = file_type:pdb bold input_file
action = *model_completion capping
  .type = choice
  .help = The type of hydrogen addition requested. Model completion will \
          complete side-chains and terminii. Capping with add enough hydrogens \
          to create a stable molecule for QM convergence.
skip_validation = False
  .type = bool
  .help = Don't perform the validation of the charge after finalisation.
calculate_charge = False
  .type = bool
  .help = Will calculate total charge of molecule.
append_to_end_of_model = False
  .type = bool
reduce = True
  .type = bool
  .help = Use reduce to add hydrogens or fall back to Phenix.elbow
options
{
  neutron = *all_h all_d hd_and_h hd_and_d all_hd
    .type = choice
  hydrogen_atom_occupancies = 0.
    .type = float
}
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def print_legend_and_usage(log):
  print("-"*79, file=log)
  print("                               qr.finalise", file=log)
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  print(master_params().show(), file=log)

def get_inputs(args, log, master_params):
  inputs = mmtbx.utils.process_command_line_args(
    args                             = args,
    master_params                    = master_params,
    suppress_symmetry_related_errors = True)
  params = inputs.params.extract()
  # Check model file
  if (len(inputs.pdb_file_names) == 0 and (params.model_file_name is None)):
    raise Sorry("No model file found.")
  elif (len(inputs.pdb_file_names) == 1):
    params.model_file_name = inputs.pdb_file_names[0]
  elif (len(inputs.pdb_file_names) > 1):
  #else:
    raise Sorry("Only one model file should be given")
  return params

def run(args, log):
  if len(args)==0:
    print_legend_and_usage(log)
    return
  params = get_inputs(
    args          = args,
    log           = log,
    master_params = master_params(),
  )
  assert params.model_file_name.lower().endswith((".cif", ".pdb", ".ent"))
  del sys.argv[1:]
  model_completion=True
  # this is a pour plumbing job
  if params.action=='capping': model_completion=False
  # make model
  pi_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  pi_params.pdb_interpretation.use_neutron_distances = True
  pi_params.pdb_interpretation.restraints_library.cdl = False
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(file_name=params.model_file_name),
    log         = null_out())
  model.process(make_restraints=True, grm_normalization=True,
    pdb_interpretation_params = pi_params)
  # Run!
  model = finalise.run(
    model                     = model,
    model_completion          = model_completion,
    skip_validation           = params.skip_validation,
    append_to_end_of_model    = params.append_to_end_of_model,
    neutron_option            = params.options.neutron,
    hydrogen_atom_occupancies = params.options.hydrogen_atom_occupancies,
    use_reduce                = params.reduce,
    )
  # Calculate charge
  rc = None
  if params.calculate_charge:
    from qrefine import charges
    cc = charges.charges_class(
      pdb_hierarchy         = model.get_hierarchy(),
      crystal_symmetry      = model.crystal_symmetry(),
      ligand_cif_file_names = None,
      electrons             = True,
      verbose               = True,
    )
    rc = cc.get_total_charge(
      list_charges=False,
      assert_correct_chain_terminii=True,
    )
  # Output model with H
  if params.model_file_name.lower().endswith(".cif"):
    omo = model.model_as_mmcif()
    ext = ".cif"
  else:
    omo = model.model_as_pdb()
    ext = ".pdb"
  prefix = Path(params.model_file_name).stem
  if rc is None:
    output = "%s_complete%s"%(prefix, ext)
  else:
    output = "%s_complete_charge%s%s"%(prefix, str(rc), ext)
  with open(output, "w") as fo:
    fo.write(omo)
  print("\n  Output written to: %s" % output)

if __name__ == '__main__':
  t0 = time.time()
  print("Starting Q|R finalise", file=log)
  run(args=sys.argv[1:], log=log)
  print("Time: %6.4f" % (time.time() - t0), file=log)
