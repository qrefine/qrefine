from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.cluster
import sys
import time
import os.path
import libtbx
import iotbx.pdb
import mmtbx.utils
from qrefine.fragment import fragments

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_yoink_path =os.path.join(qrefine_path, "plugin","yoink","yoink")

log = sys.stdout

legend = """\
Cluster a system into many small pieces
"""

master_params_str ="""
maxnum_residues_in_cluster = 25
    .type = int
    .help = maximum number of residues in a cluster
bcc_threshold = 9
    .type = int
    .help = threshold value for bcc    
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def get_inputs(args, log, master_params):
  inputs = mmtbx.utils.process_command_line_args(
    args                             = args,
    master_params                    = master_params,
    suppress_symmetry_related_errors = True)
  params = inputs.params.extract()
  return params

def run(pdb_file, log):
  params = get_inputs(
    args          = args,
    log           = log,
    master_params = master_params(),
  )
  print("max number of residues in each cluster:\n", params.maxnum_residues_in_cluster, file=log)
  print("bcc threshold value:\n", params.bcc_threshold, file=log)
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  fq = fragments(
    pdb_hierarchy=ph,
    crystal_symmetry=cs,
    maxnum_residues_in_cluster=params.maxnum_residues_in_cluster,
    bcc_threshold = params.bcc_threshold,
    clusters_only = True)
  print("Residue indices for each cluster:\n", fq.clusters, file=log)
  print('# clusters  : ',len(fq.clusters), file=log)


if (__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args[0], log)
  print("Time: %6.4f" % (time.time() - t0), file=log)
