from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os, tempfile
from six.moves import cStringIO as StringIO
from libtbx.utils import Sorry, null_out
from qrefine import qr, api, refine, conformers
from qrefine.tests.unit import run_qrefine_in, run_fmodel, fast_args

def test_finalize_auto_cust():
  for d_min, b, a in [(1.0, 0.025, 2.5), (2.0, 0.02, 2.0), (3.5, 0.01, 1.7)]:
    p = qr.get_default_params(); p.quantum.engine_name = "aimnet2"
    api.finalize_params(p, n_conformers=1, d_min=d_min, has_data=True, log=null_out())
    assert p.finalized
    assert (p.refine.max_bond_rmsd, p.refine.max_angle_rmsd) == (b, a)
    assert p.expansion and not p.cluster.clustering
    assert p.refine.minimizer == "lbfgs" and p.refine.gradient_only
    assert p.refine.number_of_weight_search_cycles == 10
  p = qr.get_default_params(); p.quantum.engine_name = "aimnet2"
  api.finalize_params(p, n_conformers=1, has_map=True, log=null_out())
  assert p.refine.minimizer == "lbfgsb" and not p.expansion
  assert p.refine.max_bond_rmsd == 0.009
  # auto_cust=False skips the presets
  p = qr.get_default_params(); p.quantum.engine_name = "aimnet2"
  p.auto_cust = False
  api.finalize_params(p, n_conformers=1, d_min=2.0, has_data=True)
  assert not p.expansion
  # aimnet2 presets with data need d_min
  p = qr.get_default_params(); p.quantum.engine_name = "aimnet2"
  try: api.finalize_params(p, n_conformers=1, has_data=True)
  except Sorry as e: assert "d_min is required" in str(e), str(e)
  else: raise AssertionError

def assert_sorry(**kw):
  p = qr.get_default_params()
  for k, v in kw.items():
    o = p
    parts = k.split("__")
    for s in parts[:-1]: o = getattr(o, s)
    setattr(o, parts[-1], v)
  try: api.finalize_params(p, n_conformers=1, d_min=2.0, has_data=True)
  except Sorry: pass
  else: raise AssertionError(kw)

def test_finalize_refuses():
  assert_sorry(expansion=True, cluster__clustering=True)
  assert_sorry(refine__minimizer="lbfgsb", refine__gradient_only=True)
  p = qr.get_default_params(); p.altlocs.handling = "native"
  p.restraints = "qm"
  try: api.finalize_params(p, d_min=2.0, has_data=True, n_conformers=2)
  except Sorry: pass
  else: raise AssertionError
  # a second call on the same extract
  p = qr.get_default_params(); api.finalize_params(p, n_conformers=1)
  try: api.finalize_params(p, n_conformers=1)
  except Sorry as e: assert "already finalized" in str(e), str(e)
  else: raise AssertionError

def test_finalize_modes():
  p = qr.get_default_params(); api.finalize_params(p, n_conformers=1)
  assert p.refine.mode == "opt"
  p = qr.get_default_params(); api.finalize_params(p, n_conformers=1, has_map=True)
  assert p.refine.mode == "refine"
  p = qr.get_default_params(); p.cluster.clustering = True
  api.finalize_params(p, n_conformers=1, d_min=2, has_data=True)
  assert p.refine.gradient_only and p.refine.minimizer == "lbfgs"
  p = qr.get_default_params(); p.quantum.engine_name = "mopac"
  api.finalize_params(p, n_conformers=1)
  assert p.quantum.method == "PM7" and p.quantum.basis == ""
  p = qr.get_default_params()
  api.finalize_params(p, d_min=2, has_data=True, n_conformers=2)
  assert p.refine.minimizer == "lbfgs" and p.refine.gradient_only

def test_rst_file_and_max_atoms_gone():
  assert "rst_file" not in qr.master_phil_str
  assert "max_atoms" not in qr.master_phil_str
  src = open(refine.__file__.replace(".pyc", ".py")).read()
  assert "max_atoms" not in src and "rst_file" not in src

class _fmodel(object):
  def __init__(self, hl): self.hl = hl
  def hl_coeffs(self): return self.hl

def test_resolve_target_name():
  p = qr.get_default_params()
  assert p.refine.refinement_target_name == "Auto" or \
    str(p.refine.refinement_target_name) == "Auto"
  assert api.resolve_target_name(p, _fmodel(None)) == "ml"
  assert api.resolve_target_name(p, _fmodel(object())) == "mlhl"
  p.refine.refinement_target_name = "ls_wunit_k1"
  assert api.resolve_target_name(p, _fmodel(object())) == "ls_wunit_k1"
  p.refine.refinement_target_name = "mlhl"
  try: api.resolve_target_name(p, _fmodel(None))
  except Sorry: pass
  else: raise AssertionError

class _twin_fmodel(object):
  """Just the attribute api.refine looks at (mmtbx.f_model.manager)."""
  twin_law = "h,-h-k,-l"

def test_refine_refuses_twin():
  p = qr.get_default_params(); api.finalize_params(p, n_conformers=1)
  try: api.refine(None, p, fmodel=_twin_fmodel(), log=null_out())
  except Sorry as e: assert "twin" in str(e), str(e)
  else: raise AssertionError

def test_refine_run_needs_finalized():
  p = qr.get_default_params()
  try: refine.run(None, None, None, p, null_out(), None, None)
  except AssertionError as e:
    assert "finalize_params" in str(e)
  else: raise AssertionError

def test_n_conformers_required_and_checked():
  import iotbx.pdb, mmtbx.model
  from mmtbx.hydrogens.tst_xh_mapper import pdb_str
  p = qr.get_default_params()
  try: api.finalize_params(p, d_min=2.0, has_data=True)
  except TypeError: pass
  else: raise AssertionError
  api.finalize_params(p, n_conformers=1)
  assert p.n_conformers == 1
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str),
    log = null_out())
  from qrefine import conformers
  assert conformers.n_conformers(model.get_hierarchy()) > 1
  try: refine.run(model, None, None, p, null_out(), None, None)
  except AssertionError as e:
    assert "different conformer count" in str(e), str(e)
  else: raise AssertionError

# ---------------------------------------------------------------------------
# api.refine (fixture in.pdb, data.mtz written by run())
# ---------------------------------------------------------------------------

common_args = fast_args
cctbx_args = ["restraints=cctbx", "expansion=False", "clustering=False"]
# MOPAC: no weight search (each MOPAC evaluation costs seconds)
mopac_args = ["restraints=qm", "quantum.engine_name=mopac", "expansion=False",
  "clustering=False", "skip_weight_search=True"]

_root = [None] # the test folder (run() sets it)

def inputs():
  """Absolute paths of the fixture model and data (in the test folder)."""
  return [os.path.join(_root[0], f) for f in ["in.pdb", "data.mtz"]]

def parse_params(args):
  """qr.refine parameters from key=value strings, as the command line."""
  import iotbx.phil
  master = iotbx.phil.parse(qr.master_phil_str)
  interpreter = master.command_line_argument_interpreter(home_scope=None)
  sources = [interpreter.process(arg=a) for a in args]
  return master.fetch(sources=sources).extract()

def objects_like_cli(args, log=None, with_data=True):
  """
  The objects qr.refine builds from args (files and key=value strings):
  DataManager, qr.prepare_inputs, finalize_params, the checks and the
  refinement target. Returns prepare_inputs' group_args plus params.
  """
  from iotbx.data_manager import DataManager
  if(log is None): log = null_out()
  dm = DataManager(qr.Program.datatypes)
  phil_args = []
  for a in args:
    if(not os.path.isfile(a)):  phil_args.append(a)
    elif(a.endswith(".mtz")):
      if(with_data): dm.process_miller_array_file(a)
    else: dm.process_model_file(a)
  params = parse_params(phil_args)
  o = qr.prepare_inputs(data_manager=dm, params=params, log=log)
  has_map = o.map_data is not None
  api.finalize_params(params,
    n_conformers     = conformers.n_conformers(o.model.get_hierarchy()),
    d_min            = o.d_min,
    scattering_table = params.scattering_table,
    has_map          = has_map,
    has_data         = o.fmodel is not None or has_map,
    log              = log)
  api.check_engine(params)
  api.check_model(o.model, params, has_map=has_map)
  if(o.fmodel is not None):
    o.fmodel.set_target_name(
      target_name = api.resolve_target_name(params, o.fmodel))
  o.params = params
  return o

def reported_r_work(cli):
  """The last r_work qr.refine reports, as printed (%6.4f)."""
  key = "Best Rwork, Rfree (after refinement):"
  lines = [l for l in cli.std_out_text.splitlines() if l.startswith(key)]
  assert len(lines) > 0
  return lines[-1].split(":")[1].split()[0]

class in_folder(object):
  """with in_folder(name): run in a fresh sub-folder (created), then back."""
  def __init__(self, name): self.name = name
  def __enter__(self):
    self.cwd = os.getcwd()
    os.makedirs(self.name)
    os.chdir(self.name)
    return os.getcwd()
  def __exit__(self, *a):
    os.chdir(self.cwd)

def listing(folder):
  """Recursive listing (relative paths) of folder."""
  result = []
  for root, dirs, files in os.walk(folder):
    for n in dirs + files:
      result.append(os.path.relpath(os.path.join(root, n), folder))
  return sorted(result)

def api_vs_cli(args, label):
  """
  qr.refine on args in folder cli_<label>, api.refine on the same objects
  in api_<label> (working_folder api_ase): same r_work as reported, same
  sites (to the PDB rounding), and qr.write_output of the API result gives
  the command line's file byte for byte.
  """
  args = inputs() + args
  cli = run_qrefine_in(args, "in", folder="cli_"+label)
  with in_folder("api_"+label):
    o = objects_like_cli(args)
    r = api.refine(o.model, o.params, fmodel=o.fmodel,
      working_folder="api_ase", log=null_out())
    assert r.mode == "refine"
    assert "%6.4f" % r.r_work == reported_r_work(cli), (r.r_work,
      reported_r_work(cli))
    d = r.sites_cart.rms_difference(cli.model.get_sites_cart())
    assert d < 1.e-3, d # PDB coordinates have three decimals
    fn = qr.write_output(result=r, params=o.params, prefix=o.prefix,
      log=null_out())
    assert fn == os.path.join("pdb", "in_refined.pdb"), fn
    assert open(fn).read() == open(cli.pdb).read()
  return r

def test_api_equals_cli_cctbx():
  r = api_vs_cli(common_args + cctbx_args, "cctbx")
  assert r.n_cycles > 0 and r.n_fev > 0
  assert r.cc_mask is None and r.r_free is not None
  assert r.data_weight is not None and r.restraints_weight_scale is not None
  assert r.b_rmsd is not None and r.a_rmsd is not None
  # the X-H mapper is active (n_gaussian): reported sites are X-ray ones
  assert r.sites_cart.rms_difference(r.sites_cart_nuclear) > 1.e-3
  assert r.model.get_sites_cart().rms_difference(r.sites_cart_nuclear) == 0

def test_api_equals_cli_mopac():
  api_vs_cli(common_args + mopac_args, "mopac")

def test_inputs_untouched():
  """
  The caller's model and fmodel (here sharing one xray_structure) are not
  changed by api.refine.
  """
  import mmtbx.f_model
  with in_folder("untouched"):
    o = objects_like_cli(inputs() + common_args + cctbx_args)
    fmodel = mmtbx.f_model.manager(
      f_obs          = o.fmodel.f_obs(),
      r_free_flags   = o.fmodel.r_free_flags(),
      xray_structure = o.model.get_xray_structure())
    fmodel.update_all_scales(remove_outliers=False,
      refine_hd_scattering=False)
    assert fmodel.xray_structure is o.model.get_xray_structure()
    sites = o.model.get_sites_cart().deep_copy()
    pdb_str = o.model.get_hierarchy().as_pdb_string()
    r_work = fmodel.r_work()
    r = api.refine(o.model, o.params, fmodel=fmodel,
      working_folder="ase", log=null_out())
    assert r.sites_cart_nuclear.rms_difference(sites) > 1.e-3
    assert o.model.get_sites_cart().rms_difference(sites) == 0
    assert o.model.get_hierarchy().as_pdb_string() == pdb_str
    assert fmodel.r_work() == r_work
    assert r.model is not o.model

def test_no_writes_outside():
  """
  clustering=True (fragment and expansion files): only the working folder
  appears in the run folder; intermediate_folder=None writes no cycle file.
  """
  with in_folder("outside") as folder:
    o = objects_like_cli(inputs() + common_args +
      ["restraints=cctbx", "expansion=False", "clustering=True"])
    before = sorted(os.listdir(folder))
    api.refine(o.model, o.params, fmodel=o.fmodel, working_folder="wf",
      log=null_out())
    after = sorted(os.listdir(folder))
    assert after == sorted(before + ["wf"]), (before, after)
    files = listing("wf")
    # altloc fixture clustered as is (per-fragment conformer variants)
    assert "expansion.pdb" in files, files
    assert "A/expansion.pdb" not in files, files
    assert [f for f in files if f.endswith("_cycle.pdb")] == []

def test_intermediate_folder():
  """
  intermediate_folder="cycles" (absent before) receives the cycle models;
  nothing else appears in the run folder (cctbx without clustering writes
  no scratch, so the working folder ase is not created).
  """
  with in_folder("intermediate") as folder:
    o = objects_like_cli(inputs() + common_args + cctbx_args)
    assert not os.path.exists("cycles")
    before = sorted(os.listdir(folder))
    api.refine(o.model, o.params, fmodel=o.fmodel, working_folder="ase",
      intermediate_folder="cycles", log=null_out())
    after = sorted(os.listdir(folder))
    assert set(after) - set(before) <= {"ase", "cycles"}, (before, after)
    assert "cycles" in after and set(before) <= set(after), (before, after)
    files = sorted(os.listdir("cycles"))
    assert "0_refine_cycle.pdb" in files, files
    assert all(f.endswith("_refine_cycle.pdb") for f in files), files

def test_called_twice():
  with in_folder("twice"):
    o = objects_like_cli(inputs() + common_args + cctbx_args)
    r1 = api.refine(o.model, o.params, fmodel=o.fmodel, working_folder="ase",
      log=null_out())
    r2 = api.refine(o.model, o.params, fmodel=o.fmodel, working_folder="ase",
      log=null_out())
    assert r1.r_work == r2.r_work and r1.r_free == r2.r_free
    assert r1.sites_cart.rms_difference(r2.sites_cart) == 0

class _failing_restraints(object):
  """Fragment restraints that fail (a QM job that crashed)."""
  def __init__(self, file_name): self.file_name = file_name
  def __call__(self, args): raise RuntimeError("engine failed")

class _fragment_manager(object):
  """Just enough of fragment.fragments for cluster_restraints.from_cluster."""
  def __init__(self, hierarchy, cs):
    from libtbx import group_args
    from cctbx.array_family import flex
    self.pdb_hierarchy_super = hierarchy
    self.expansion = group_args(cs_box = cs)
    self.fragment_selections = [flex.bool(hierarchy.atoms_size(), True)]
  def update_xyz(self, sites_cart): pass
  def get_fragment_extracts(self, with_super_sphere_grm=False): return None

def test_cluster_error_copy():
  """
  A failure of the fragment jobs: with keep_error_copy (a folder the caller
  chose) the working folder is copied into <folder>/error and the Sorry
  names it; without (api.refine's temporary folder) nothing is copied and
  nothing appears next to the folder.
  """
  import iotbx.pdb
  from libtbx import group_args
  from qrefine import cluster_restraints
  pdb_inp = iotbx.pdb.input(file_name=inputs()[0])
  h = pdb_inp.construct_hierarchy()
  parallel = group_args(method="threading", nproc=1, qsub_command=None)
  with in_folder("error_copy") as folder:
    for keep in [True, False]:
      wf = os.path.join(folder, "wf_%s" % keep)
      os.makedirs(wf)
      open(os.path.join(wf, "scratch.txt"), "w").write("x")
      c = cluster_restraints.from_cluster(
        restraints_manager = _failing_restraints(
          os.path.join(wf, "tmp_ase.pdb")),
        fragment_manager   = _fragment_manager(h, pdb_inp.crystal_symmetry()),
        parallel_params    = parallel,
        keep_error_copy    = keep)
      try: c.target_and_gradients(sites_cart=h.atoms().extract_xyz())
      except Sorry as e: msg = str(e)
      else: raise AssertionError
      error = os.path.join(wf, "error")
      if(keep):
        assert os.path.isfile(os.path.join(error, "scratch.txt"))
        assert "copied to %s" % error in msg, msg
      else:
        assert not os.path.exists(error) and "copied" not in msg, msg
    assert sorted(os.listdir(folder)) == ["wf_False", "wf_True"]

def _no_files_at_all(label, clustering):
  """
  working_folder=None, intermediate_folder=None: nothing is written into
  the run folder (recursively) nor into the temporary folder (a private
  one for this test; in the system one no qrefine_* entry appears).
  """
  private_tmp = os.path.abspath(label + "_tmp")
  os.makedirs(private_tmp)
  system_tmp = tempfile.gettempdir()
  def qrefine_entries():
    return sorted([n for n in os.listdir(system_tmp)
      if n.startswith("qrefine")])
  with in_folder(label) as folder:
    o = objects_like_cli(inputs() + common_args + ["restraints=cctbx",
      "expansion=False", "clustering=%s" % clustering])
    assert o.params.cluster.save_clusters # default: still nothing written
    before = listing(folder)
    system_before = qrefine_entries()
    saved = tempfile.tempdir
    tempfile.tempdir = private_tmp
    try:
      r = api.refine(o.model, o.params, fmodel=o.fmodel, working_folder=None,
        intermediate_folder=None, log=null_out())
    finally:
      tempfile.tempdir = saved
    assert r.r_work is not None
    assert listing(folder) == before, (before, listing(folder))
  assert os.listdir(private_tmp) == [], os.listdir(private_tmp)
  assert qrefine_entries() == system_before

def test_no_files_at_all():
  _no_files_at_all("nofiles", clustering=False)

def test_no_files_at_all_clustering():
  _no_files_at_all("nofiles_clustering", clustering=True)

def test_mopac_tempdir_removed():
  """
  An external-program engine with working_folder=None runs in a temporary
  folder that is named in the log and removed on return.
  """
  with in_folder("mopac_tmp") as folder:
    o = objects_like_cli(inputs() + common_args + ["restraints=qm",
      "quantum.engine_name=mopac", "clustering=False", "expansion=True",
      "skip_weight_search=True"])
    before = listing(folder)
    log = StringIO()
    r = api.refine(o.model, o.params, fmodel=o.fmodel, working_folder=None,
      log=log)
    assert r.r_work is not None
    lines = [l for l in log.getvalue().splitlines()
      if l.startswith("scratch folder: ")]
    assert len(lines) == 1, lines
    assert lines[0].endswith(" (removed on return)"), lines[0]
    path = lines[0][len("scratch folder: "):-len(" (removed on return)")]
    assert os.path.basename(path).startswith("qrefine_"), path
    assert not os.path.exists(path), path
    assert listing(folder) == before

def test_opt_mode():
  with in_folder("opt"):
    o = objects_like_cli(inputs() + cctbx_args +
      ["number_of_micro_cycles=2"], with_data=False)
    assert o.fmodel is None and o.params.refine.mode == "opt"
    r = api.refine(o.model, o.params, working_folder="ase", log=null_out())
    assert r.mode == "opt"
    assert r.r_work is None and r.r_free is None and r.cc_mask is None
    assert r.sites_cart.rms_difference(r.sites_cart_nuclear) == 0
    assert r.n_cycles > 0 and r.n_fev > 0
    assert r.b_rmsd is not None and r.a_rmsd is not None

def test_real_space_mode():
  """
  Map refinement (map from the fixture without its D atom: the map-model CC
  uses the electron table, which has no D): mode real_space, CC, nothing
  written with working_folder=None.
  """
  from cctbx.array_family import flex
  from mmtbx.hydrogens.tst_xh_mapper import pdb_str
  with open("in_no_d.pdb", "w") as fo:
    fo.write("\n".join([l for l in pdb_str.splitlines()
      if not l.endswith(" D")]) + "\n")
  model_file = os.path.abspath("in_no_d.pdb")
  with in_folder("real_space"):
    o = objects_like_cli([model_file] + cctbx_args +
      ["number_of_weight_search_cycles=1", "number_of_refine_cycles=1"])
    xrs = o.model.get_xray_structure()
    fc = xrs.structure_factors(d_min=2.0).f_calc()
    map_data = fc.fft_map(resolution_factor=0.25).real_map_unpadded()
    map_data = map_data - flex.mean(map_data)
    map_data = map_data / map_data.sample_standard_deviation()
    params = parse_params(cctbx_args + ["number_of_weight_search_cycles=1",
      "number_of_refine_cycles=1"])
    api.finalize_params(params, n_conformers=o.params.n_conformers,
      has_map=True)
    import sys
    saved, sys.stdout = sys.stdout, StringIO()
    try:
      r = api.refine(o.model, params, map_data=map_data, working_folder=None,
        log=None)
      printed = sys.stdout.getvalue()
    finally:
      sys.stdout = saved
    assert printed == "", printed # log=None: silent
    assert r.mode == "real_space"
    assert r.cc_mask is not None and r.r_work is None
    assert r.n_cycles == 1 and r.data_weight is not None
    assert r.sites_cart.rms_difference(r.sites_cart_nuclear) == 0
    assert r.sites_cart_nuclear.rms_difference(o.model.get_sites_cart()) > 0
    assert os.listdir(".") == [], os.listdir(".")

def test_raw_params_refused():
  with in_folder("raw"):
    o = objects_like_cli(inputs() + common_args + cctbx_args)
    try:
      api.refine(o.model, qr.get_default_params(), fmodel=o.fmodel,
        log=null_out())
    except AssertionError as e:
      assert "finalize_params" in str(e)
    else: raise AssertionError

def test_auto_target_mlhl():
  """
  fmodel with experimental phases (HL coefficients made up from the model's
  phases): the Auto target resolves to mlhl and api.refine runs with it.
  """
  import mmtbx.f_model
  with in_folder("mlhl"):
    o = objects_like_cli(inputs() + common_args + cctbx_args)
    hl = o.fmodel.f_calc().make_up_hl_coeffs(k_blur=1, b_blur=10)
    fmodel = mmtbx.f_model.manager(
      f_obs          = o.fmodel.f_obs(),
      r_free_flags   = o.fmodel.r_free_flags(),
      abcd           = hl,
      xray_structure = o.model.get_xray_structure().deep_copy_scatterers())
    fmodel.update_all_scales(remove_outliers=False,
      refine_hd_scattering=False)
    assert fmodel.hl_coeffs() is not None
    target_name = api.resolve_target_name(o.params, fmodel)
    assert target_name == "mlhl", target_name
    fmodel.set_target_name(target_name=target_name)
    r = api.refine(o.model, o.params, fmodel=fmodel, working_folder="ase",
      log=null_out())
    assert r.mode == "refine" and r.r_work is not None

def setup_fixture():
  """The fixture in.pdb and data.mtz in the current folder, which becomes
  the folder inputs() reads them from (tst_56, tst_72, tst_73)."""
  from mmtbx.hydrogens.tst_xh_mapper import pdb_str
  with open("in.pdb", "w") as fo:
    fo.write(pdb_str)
  run_fmodel(args=["in.pdb", "high_res=2.5", "r_free=0.05"], prefix="data")
  _root[0] = os.getcwd()

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Parameter rules in api.finalize_params; Auto target; removed params;
  api.refine: same result as qr.refine, inputs untouched. Where api.refine
  writes: tst_72; no files, temporary folders, opt and real-space modes,
  mlhl target: tst_73.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  test_finalize_auto_cust()
  test_finalize_refuses()
  test_finalize_modes()
  test_rst_file_and_max_atoms_gone()
  test_resolve_target_name()
  test_refine_refuses_twin()
  test_refine_run_needs_finalized()
  test_n_conformers_required_and_checked()
  setup_fixture()
  test_api_equals_cli_cctbx()
  test_api_equals_cli_mopac()
  test_inputs_untouched()

if(__name__ == "__main__"):
  run()
  print("OK")
