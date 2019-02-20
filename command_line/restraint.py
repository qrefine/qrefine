from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.restraint
import sys
import time
import os.path
import argparse
import iotbx.pdb
import libtbx.load_env
import mmtbx.command_line
import qrefine.qr as qr
from qrefine import __version__
from qrefine.restraints import from_qm, from_cctbx
from libtbx.utils import Sorry,Usage

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_path = os.path.join(qrefine_path, "core")

log = sys.stdout

legend = """\
Compute energy and gradient for a system
"""

def get_help():
  print legend
  raise Usage("""
    qr.restraint is an open-source module that ccomputes chemical restraints from ab initio calculations.

    Example:
    qr.restraint model.pdb  data.mtz [<param_name>=<param_value>] ...

    Options:
    qr.restraint --defaults (print default parameters)
    qr.restraint --version  (print version information)
    """)
  sys.exit(0)
  return

def get_master_phil():
  return mmtbx.command_line.generate_master_phil_with_inputs(
    phil_string=qr.master_params_str)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, get_master_phil().show()

def run(args, log):
  args.append('mode=opt')
  args.append('number_of_micro_cycles=1')
  args.append('number_of_macro_cycles=1')
  args.append('max_iterations_refine=1')
  args.append('clustering=false')
  cmdline = mmtbx.utils.process_command_line_args(
      args=args,
      master_params=get_master_phil())
  if (len(args) == 0 or '--help' in args):
      get_help()
      return
  elif ('--defaults' in args or '--show' in args):
      print_legend_and_usage(log)
      return
  elif ('--version' in args):
      print
      __version__
      return
  print >> log, "Computing restraint"
  cmdline.params.show(out=log, prefix="   ")
  params = cmdline.params.extract()
  print ""
  print_legend_and_usage(log)
  cmdline = mmtbx.command_line.load_model_and_data(
        args=args,
        master_phil=get_master_phil(),
        create_fmodel=False,
        out=log)
  model = qr.process_model_file(
    pdb_file_name    = cmdline.pdb_file_names[0],
    cif_objects      = cmdline.cif_objects,
    crystal_symmetry = cmdline.crystal_symmetry)
  params.refine.mode == "opt"
  qr.run( model    = model,
          fmodel   = None,
          map_data = None,
          params   = params,
          rst_file = None,
          prefix   = os.path.basename(cmdline.pdb_file_names[0])[:-4],
          log=log)


if (__name__ == "__main__"):
  print "Restraint for Q|R"
  t0 = time.time()
  print >> log, "Starting Q|R"
  print >> log,'version: ',__version__
  run(args=sys.argv[1:], log=log)
  print >> log, "Time: %6.4f" % (time.time() - t0)

