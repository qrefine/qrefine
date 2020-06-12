from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.refine
import os
import sys
import time
import libtbx.load_env
from libtbx import  easy_run
from libtbx.command_line import easy_qsub
import mmtbx.command_line
from qrefine import qr, __version__
from mmtbx import utils
from iotbx import reflection_file_utils
from io import StringIO
from libtbx import group_args
from libtbx.utils import Sorry,Usage
from scitbx.array_family import flex

# phenix_source = os.path.dirname(libtbx.env.dist_path("phenix"))
qrefine_path = libtbx.env.find_in_repositories("qrefine")
example_path = os.path.join(qrefine_path,"examples")

log = sys.stdout

legend = """
Refine a model using restraints from Quantum Chemistry
"""

def get_help():
  print(legend)
  raise Usage("""
    qr.refine is an open-source module that carries out refinement of bio-macromolecules
    utilizing chemical restraints from ab initio calculations.

    Example:
    qr.refine model.pdb model.mtz [<param_name>=<param_value>] ...

    Options:
    qr.refine --defaults (print default parameters)
    qr.refine --version  (print version informations)
    """)
  sys.exit(0)
  return  

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=qr.master_params_str)

def print_legend_and_usage(log):
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  print(get_master_phil().show(), file=log)

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def run(args, log):
  cmdline = mmtbx.utils.process_command_line_args(
    args          = args,
    master_params = get_master_phil())
  if(len(args)==0 or '--help' in args):
    get_help()
    return
  elif('--defaults' in args or '--show' in args):
    print_legend_and_usage(log)
    return  
  elif('--version' in args):
    print(__version__)
    return
  print("Running refinement", file=log)
  cmdline.params.show(out=log, prefix="   ")
  params = cmdline.params.extract()

  # Read atomic model
  # XXX This is not Oleg's model !!!
  # need to validate the input args
  model = qr.process_model_file(
    pdb_file_name    = cmdline.pdb_file_names[0],
    cif_objects      = cmdline.cif_objects,
    crystal_symmetry = cmdline.crystal_symmetry)
  map_data = None
  fmodel = None
  if(params.refine.mode=="opt" or params.refine.mode=='gtest'):
    if (len(cmdline.reflection_files)>0 or cmdline.ccp4_map is not None):
      print("WARNING: data files not used in optimization or gradient test! ", file=log)
  elif(len(cmdline.reflection_files)>0 and params.refine.mode=="refine"):
    # Read reflection data
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
      r_free_flags=f_obs.generate_r_free_flags()
      print("WARNING: no free-R flags available in inputs. ", file=log)
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      xray_structure = model.xray_structure,
      target_name    = params.refine.refinement_target_name)
    if(params.refine.update_all_scales):
      fmodel.update_all_scales(remove_outliers=False)
      fmodel.show(show_header=False, show_approx=False)
    print("Initial r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(),
      fmodel.r_free()), file=log)
  elif(cmdline.ccp4_map is not None and params.refine.mode=="refine"): 
    # Read map
    map_data = cmdline.ccp4_map.map_data()
    # Normalize map
    map_data = map_data - flex.mean(map_data)
    sd = map_data.sample_standard_deviation()
    assert sd != 0
    map_data = map_data/sd
    model = model.model#inp.model()
    model = group_args(
      model              = model,
      processed_pdb_file = model._processed_pdb_file, # This must go, use model!
      pdb_hierarchy      = model.get_hierarchy(),     # This must go, use model!
      xray_structure     = model.get_xray_structure(),# This must go, use model!
      cif_objects        = model._restraint_objects,  # This must go, use model!
      has_hd             = model.has_hd)
  else:
    raise Sorry("Refinement requested (refine.mode==refine) but no data provided.")
  log.flush()
  qr.run(
    model    = model, # XXX This is not mmtbx.model.manager !!! (see above).
    fmodel   = fmodel,
    map_data = map_data,
    params   = params,
    rst_file = params.rst_file,
    prefix   = os.path.basename(cmdline.pdb_file_names[0])[:-4],
    log      = log)

if __name__ == '__main__':
  t0 = time.time()
  print("Starting Q|R", file=log)
  print('version: ',__version__, file=log)
  run(args=sys.argv[1:], log=log)
  print("Time: %6.4f" % (time.time() - t0), file=log)
