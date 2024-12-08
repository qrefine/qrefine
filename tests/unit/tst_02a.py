from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import iotbx.pdb
import libtbx.load_env
from qrefine.tests.unit import run_tests

from scitbx.array_family import flex
from qrefine import qr
from qrefine import restraints
import mmtbx.model
from libtbx import easy_run

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def run(prefix):
  """
  Exercise equivalency of from_altlocs and from_cctbx via minimization and
  qr.refine
  All three output PDBs are expected to be identical up to PDB file precision.
  """
  #
  # Run minimization
  #
  for option in ["from_altlocs", "from_cctbx"]:
    file_name = os.path.join(qr_unit_tests,"data_files","h_altconf_2_complete_modified.pdb")
    pdb_inp = iotbx.pdb.input(file_name)
    model = mmtbx.model.manager(
      model_input      = None,
      pdb_hierarchy    = pdb_inp.construct_hierarchy(),
      crystal_symmetry = pdb_inp.crystal_symmetry())

    params = qr.get_default_params()
    params.restraints="cctbx"
    params.quantum.method='PM7'
    params.quantum.basis=''
    params.cluster.clustering=False

    if option == "from_altlocs":
      rm = restraints.from_altlocs2(model = model, method="subtract",
        params = params)
      sites_cart = cctbx_opt(model=model, restraints_manager=rm)
    elif option == "from_cctbx":
      restraints_source = restraints.restraints(params = params, model = model)
      rm = restraints.from_cctbx(
        restraints_manager = restraints_source.restraints_manager)
      sites_cart = cctbx_opt(model=model,
        restraints_manager=restraints_source.restraints_manager)
    #
    model.set_sites_cart(sites_cart)
    model.get_hierarchy().write_pdb_file(
      "%s_%s.pdb"%(prefix,option),
      crystal_symmetry = pdb_inp.crystal_symmetry())
  #
  # Run qrefine
  #
  cmd = " ".join([
    "qr.refine",
    file_name,
    "clustering=false",
    "minimizer=lbfgs",
    "number_of_micro_cycles=1",
    "restraints=cctbx",
    "gradient_only=true",
    "max_iterations_refine=25",
    "output_file_name_prefix=%s"%prefix,
    "output_folder_name=%s"%prefix,
    " > %s.zlog"%prefix])
  easy_run.call(cmd)
  #
  # Check all three files are the same
  #
  cwd = os.getcwd()
  s1 = iotbx.pdb.input("%s/qrefine_tst_02a_from_altlocs.pdb"%cwd).atoms().extract_xyz()
  s2 = iotbx.pdb.input("%s/qrefine_tst_02a_from_cctbx.pdb"%cwd).atoms().extract_xyz()
  s3 = iotbx.pdb.input("%s/qrefine_tst_02a/qrefine_tst_02a_refined.pdb"%cwd).atoms().extract_xyz()
  d1 = flex.mean(flex.sqrt((s1 - s2).dot()))
  d2 = flex.mean(flex.sqrt((s1 - s3).dot()))
  assert d1 < 1.e-4, d1 # This is what we expect!
  assert d2 < 1.e-4, d2 # This is what we expect!

def cctbx_opt(model, restraints_manager, max_shift=0.2, max_iterations=25):
  from scitbx import minimizers
  from qrefine import calculator
  s1 = model.get_sites_cart().deep_copy()
  C = calculator.sites_opt(model=model, max_shift=max_shift,
    restraints_manager=restraints_manager, shift_eval="mean")
  import scitbx.lbfgs
  core_params = scitbx.lbfgs.core_parameters(
    stpmin = 1.e-9,
    stpmax = max_shift)
  minimized = minimizers.lbfgs(
    calculator     = C,
    mode           = "lbfgs",
    gradient_only  = True,
    max_iterations = max_iterations)
  C.apply_x()
  print("Moved by:", flex.mean(flex.sqrt((s1 - C.model.get_sites_cart()).dot())))
  return C.model.get_sites_cart()

if(__name__ == "__main__"):
  prefix = os.path.basename(__file__).replace(".py","")
  run_tests.runner(function=run, prefix=prefix, disable=False)
