from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.charges
import sys, time
from qrefine import charges, __version__
import iotbx
import mmtbx
from mmtbx import utils

log = sys.stdout

legend = """\
A command line tool to calculate the total formal charge of a model.
"""

master_params_str = """
model_file_name = None
  .type = path
  .short_caption = Model file
  .multiple = False
  .help = Model file name
  .style = file_type:pdb bold input_file
ligand_cif_file_name = None
  .type = path
  .short_caption = Optional Ligand restraints file in CIF format
  .caption = This restraints file will be used to get the formal \
    charge of the ligand.
  .help = This restraints file will be used to get the formal \
    charge of the ligand.
  .style = file_type:cif input_file
  .multiple = True
list_charges = False
  .type = bool
  .help = Print a list of all the atoms with non-zero charge.
electrons = True
  .type = bool
  .help = debug option - do not change.
assert_correct_chain_terminii = True
  .type = bool
  .help = Tests for hydrogens on the N terminii.
verbose = False
  .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def print_legend_and_usage(log):
  print("-"*79, file=log)
  print("                               qr.charges", file=log)
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
  t0 = time.time()
  if len(args)==0:
    print_legend_and_usage(log)
    return
  params = get_inputs(
    args          = args,
    log           = log,
    master_params = master_params(),
    #validated     = validated,
  )
  if(params.verbose):
    print("Starting Q|R charge", file=log)
    print('version: ',__version__, file=log)
  del sys.argv[1:]
  cc = charges.charges_class(
    pdb_filename = params.model_file_name,
    ligand_cif_file_names=params.ligand_cif_file_name,
    electrons=params.electrons,
    verbose=params.verbose,
  )
  rc = cc.get_total_charge(
    list_charges=params.list_charges,
    assert_correct_chain_terminii=params.assert_correct_chain_terminii,
  )
  if(params.verbose):
    print('Charge: %s' % rc, file=log)
    print("Time: %6.4f" % (time.time() - t0), file=log)
  if(params.list_charges):
    print('Charges:', file=log)
    print(rc, file=log)
    print("Time: %6.4f" % (time.time() - t0), file=log)


if __name__ == '__main__':
  run(args=sys.argv[1:], log=log)
