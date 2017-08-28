from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.charges
import sys, time
from qrefine import charges
import iotbx
import mmtbx
from mmtbx import utils

log = sys.stdout

legend = """\
Finalise a model before quantum refinement
"""

master_params_str = """
model_file_name = None
  .type = path
  .short_caption = Model file
  .multiple = False
  .help = Model file name
  .style = file_type:pdb bold input_file
list_charges = False
  .type = bool
verbose = False
  .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, "                               qr.charges"
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, master_params().show()

def get_inputs(args, log, master_params):
  print dir(mmtbx)
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
  rc = charges.run(params.model_file_name,
                   list_charges=params.list_charges,
                   verbose=params.verbose,
                 )
  print rc

if __name__ == '__main__':
  t0 = time.time()
  print >> log,"Starting Q|R charge"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
