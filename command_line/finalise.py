from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.finalise
import sys, time
from qrefine import finalise
import iotbx
import mmtbx
from mmtbx import utils
from libtbx.utils import Sorry

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
keep_alt_loc = True
  .type = bool
  .help = Retain alt loc. This is not a useful option and not tested well.
skip_validation = False
  .type = bool
  .help = Don't perform the validation of the charge after finalisation.
calculate_charge = False
  .type = bool
  .help = Will calculate total charge of molecule.
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, "                               qr.finalise"
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, master_params().show()

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
    #validated     = validated,
  )
  del sys.argv[1:]
  model_completion=True
  # this is a pour plumbing job
  if params.action=='capping': model_completion=False
  finalise.run(params.model_file_name,
               model_completion=model_completion,
               keep_alt_loc=params.keep_alt_loc,
               skip_validation=params.skip_validation,
               calculate_charge=params.calculate_charge,
               )

if __name__ == '__main__':
  t0 = time.time()
  print >> log,"Starting Q|R finalise"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
