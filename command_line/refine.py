from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.refine
import os
import sys
import time
import libtbx.load_env
from libtbx import  easy_run
from libtbx.command_line import easy_qsub
import mmtbx.command_line
from qrefine import qr
from mmtbx import utils
from iotbx import reflection_file_utils
from cStringIO import StringIO

phenix_source = os.path.dirname(libtbx.env.dist_path("phenix"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
example_path = os.path.join(qrefine_path,"examples")

log = sys.stdout

legend = """
Refine a model using restraints from Quantum Chemistry
"""

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=qr.master_params_str)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, get_master_phil().show()
  
def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def run(args, log):
  print >> log,"Running refinement"
  print_legend_and_usage(log)
  cmdline = mmtbx.utils.process_command_line_args(
    args          = args,
    master_params = get_master_phil())
  cmdline.params.show(out=log, prefix="   ")
  params = cmdline.params.extract()
  # Read atomic model
  model = qr.process_model_file(
    pdb_file_name    = cmdline.pdb_file_names[0], 
    cif_objects      = cmdline.cif_objects, 
    crystal_symmetry = cmdline.crystal_symmetry)
  # Read reflection data
  fmodel = None
  if(len(cmdline.reflection_files)>0):
    rfs = reflection_file_server(
      crystal_symmetry = cmdline.crystal_symmetry,
      reflection_files = cmdline.reflection_files)
    determine_data_and_flags_result = utils.determine_data_and_flags(
      reflection_file_server  = rfs,
      keep_going              = True,
      log                     = log)
    f_obs = determine_data_and_flags_result.f_obs
    number_of_reflections = f_obs.indices().size()
    r_free_flags = determine_data_and_flags_result.r_free_flags
    test_flag_value = determine_data_and_flags_result.test_flag_value
    if(r_free_flags is None):
      r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      xray_structure = model.xray_structure,
      target_name    = params.refine.refinement_target_name)
    if(params.refine.update_all_scales):
      fmodel.update_all_scales(remove_outliers=False)
      fmodel.show(show_header=False, show_approx=False)
    print >> log, "Initial r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), 
      fmodel.r_free())
  log.flush()
  qr.run(
    model    = model, 
    fmodel   = fmodel, 
    params   = params, 
    rst_file = params.rst_file, 
    prefix   = os.path.basename(cmdline.pdb_file_names[0])[:-4],
    log      = log)

if __name__ == '__main__':
  t0 = time.time()
  print >> log,"Starting Q|R"
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
