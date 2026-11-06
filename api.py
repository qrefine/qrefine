"""
Callable interface of qr.refine (work in progress).

Public functions:
  refine            refine (or optimize) a copy of a model; returns the result
  finalize_params   apply every parameter rule once; marks params finalized
  resolve_target_name  resolve the Auto refinement target against an fmodel
  check_engine      Sorry if the QM engine is not available (ENGINES)
  check_model       Sorry if the model is not ready for the restraints

qr.py (the command-line program) is the first caller. After finalize_params
no function changes the parameters: refine.run asserts the flag.
"""
from __future__ import absolute_import, division, print_function
import shutil
import tempfile
from libtbx import Auto, group_args
from libtbx.utils import Sorry, null_out


def _pr(log, *args):
  if(log is not None):
    print(*args, file=log)

def _auto_cust(params, d_min, has_map, has_data, log):
  """
  Presets for engine_name=aimnet2 (as used in the tests for the paper).
  Data (d_min given) wins over a map, as it always did.
  """
  if(params.quantum.engine_name != "aimnet2"): return
  if(has_data and d_min is None):
    raise Sorry("d_min is required when has_data is True (aimnet2 presets "
      "depend on resolution).")
  msg = """
The following settings have been auto-set to match refine.mode=refine
and quantum.engine_name=aimnet2:
"""
  if params.debug: params.restraints = "cctbx"
  else:            params.restraints = "qm"
  params.cluster.clustering = False
  params.refine.number_of_weight_search_cycles = 10
  params.refine.number_of_refine_cycles = 5
  params.cluster.select_within_radius = 7
  params.refine.max_iterations_weight = 100
  if(d_min is not None):
    params.expansion = True
    params.refine.minimizer = "lbfgs"
    params.refine.gradient_only = True
    if(d_min < 1.2):
      params.refine.max_bond_rmsd = 0.025
      params.refine.max_angle_rmsd = 2.5
    elif(d_min >= 1.2 and d_min < 3):
      params.refine.max_bond_rmsd = 0.02
      params.refine.max_angle_rmsd = 2.0
    else:
      params.refine.max_bond_rmsd = 0.01
      params.refine.max_angle_rmsd = 1.7
  elif(has_map):
    params.refine.max_bond_rmsd = 0.009
    params.expansion = False
    params.refine.minimizer = "lbfgsb"
    params.refine.gradient_only = False
  else: # Optimization
    assert params.refine.mode == "opt"
    params.expansion = True
  if(log is None): return
  p = params
  print(msg, file=log)
  print("  restraints                              ", p.restraints, file=log)
  print("  cluster.clustering                      ", p.cluster.clustering, file=log)
  print("  refine.minimizer                        ", p.refine.minimizer, file=log)
  print("  refine.number_of_weight_search_cycles   ", p.refine.number_of_weight_search_cycles, file=log)
  print("  refine.number_of_refine_cycles          ", p.refine.number_of_refine_cycles, file=log)
  print("  expansion                               ", p.expansion, file=log)
  print("  cluster.select_within_radius            ", p.cluster.select_within_radius, file=log)
  print("  refine.max_iterations_weight            ", p.refine.max_iterations_weight, file=log)
  print("  refine.max_bond_rmsd                    ", p.refine.max_bond_rmsd, file=log)
  print("  refine.max_angle_rmsd                   ", p.refine.max_angle_rmsd, file=log)
  print("  refine.gradient_only                    ", p.refine.gradient_only, file=log)
  print("\n***", file=log)
  print("To disable auto-custing of parameters use auto_cust=False", file=log)
  print("***\n", file=log)

def set_qm_defaults(params, log):
  """Engine method and basis defaults (was refine.set_qm_defaults)."""
  outl = ''
  if params.quantum.engine_name == 'mopac':
    if params.quantum.method == Auto:
      params.quantum.method = 'PM7'
      outl += '  Setting QM method to PM7 (mopac)\n'
    if params.quantum.basis == Auto:
      params.quantum.basis = ''
  if params.quantum.engine_name == 'xtb':
    if params.quantum.method == Auto:
      params.quantum.method = ' --gfn 2 --etemp 500 --acc 0.1 --gbsa h2o'
      _pr(log, ' Default method for xtb is %s' % params.quantum.method)
      params.quantum.basis = ''
  if params.quantum.engine_name == 'aimnet2':
    _pr(log, "")
    _pr(log, "Setting aimnet2 specific defaults:")
    if params.quantum.method == Auto:
      params.quantum.method = "aimnet2-qr"
      _pr(log, "  Default Aimnet2 model: %s" % params.quantum.method)
      # The check that aimnet2calc is importable is check_engine's job.
  else:
    if params.quantum.method == Auto:
      params.quantum.method = 'HF'
      outl += '  Setting QM method to HF\n'
    if params.quantum.basis == Auto:
      params.quantum.basis = '6-31g'
      outl += '  Setting QM basis to 6-31g\n'
  if outl:
    _pr(log, ' Setting QM defaults')
    _pr(log, outl)
  if params.quantum.engine_name == 'mopac':
    if params.quantum.basis:
      _pr(log, '  Because engine is %s basis set %s ignored' % (
        params.quantum.engine_name, params.quantum.basis))
      params.quantum.basis = ''
    if params.quantum.method == 'hf': # default
      _pr(log, ' Default method set as PM7 (mopac)')
      params.quantum.method = 'PM7'

def finalize_params(params, *, n_conformers, d_min=None,
                    scattering_table=None, has_map=False, has_data=False,
                    log=None):
  """
  Apply every parameter rule, once, and mark params finalized. Mutates and
  returns params. Rules, in this order:
    0. refine.mode: opt without data and map, refine with a map (applied
       first because rule 2 depends on it);
    1. expansion and cluster.clustering are exclusive (Sorry); lbfgsb refuses
       gradient_only (Sorry); expansion or clustering force minimizer=lbfgs
       and gradient_only=True;
    2. aimnet2 presets (skipped with auto_cust=False); d_min is the
       resolution of the data (None without reflection data);
    3. QM engine method and basis defaults;
    4. altlocs.handling=native with restraints=qm is a Sorry; more than one
       conformer (n_conformers, with altlocs.handling=conformers) forces
       minimizer=lbfgs and gradient_only=True;
    5. (the Auto refinement target is resolved by resolve_target_name, once
       the fmodel exists);
    6. params.finalized = True and params.n_conformers = n_conformers
       (refine.run checks the latter against the model).
  n_conformers, the number of conformers of the model, is required.
  scattering_table is accepted for the callers' convenience; the value of
  params.scattering_table is not changed. log may be None.
  """
  if(getattr(params, "finalized", False)):
    raise Sorry("qrefine parameters are already finalized: "
      "api.finalize_params is applied once per parameter extract (make a "
      "fresh extract to change parameters).")
  if(has_map): params.refine.mode = "refine"
  if(not has_data and not has_map): params.refine.mode = "opt"
  # 1
  if(params.expansion and params.cluster.clustering):
    raise Sorry("expansion and cluster.clustering are mutually exclusive.")
  if(params.refine.minimizer == "lbfgsb" and params.refine.gradient_only):
    raise Sorry("gradient_only must be False for lbfgsb.")
  if(params.expansion or params.cluster.clustering):
    params.refine.minimizer = "lbfgs"
    params.refine.gradient_only = True
    _pr(log, "\n  expansion/clustering require minimizer=lbfgs and "
      "gradient_only=True \n")
  # 2
  if(params.auto_cust):
    _auto_cust(params, d_min=d_min, has_map=has_map, has_data=has_data,
      log=log)
  # 3
  set_qm_defaults(params, log)
  # 4
  if(n_conformers > 1 and params.altlocs.handling == "native" and
     params.restraints == "qm"):
    raise Sorry("altlocs.handling=native needs restraints=cctbx: a QM engine cannot handle alternative conformations. Use altlocs.handling=conformers.")
  if(n_conformers > 1 and params.altlocs.handling == "conformers"):
    params.refine.minimizer = "lbfgs"
    params.refine.gradient_only = True
    _pr(log, "  alternative conformations require minimizer=lbfgs and "
      "gradient_only=True")
  if(params.cluster.clustering):
    params.refine.gradient_only = True
  # 6
  params.__inject__("finalized", True)
  params.__inject__("n_conformers", n_conformers)
  return params

def resolve_target_name(params, fmodel):
  """
  Refinement target for fmodel: Auto gives mlhl when fmodel carries
  experimental phases (hl_coeffs), else ml; an explicit choice is kept
  (mlhl without hl_coeffs is a Sorry).
  """
  name = params.refine.refinement_target_name
  has_hl = fmodel.hl_coeffs() is not None
  if(name in (None, Auto, "Auto")):
    return "mlhl" if has_hl else "ml"
  if(name == "mlhl" and not has_hl):
    raise Sorry("refinement_target_name=mlhl needs Hendrickson-Lattman "
      "coefficients in the data.")
  return name

# ---------------------------------------------------------------------------
# Engine availability (check_engine)
# ---------------------------------------------------------------------------

def _command_found(make_calculator):
  """
  True if the plugin's get_command() returns a command. The plugin object is
  made here, at check time, never at import (constructing some reads the
  environment). Any exception (xtb raises without the program) is False;
  its repr is kept in lookup.error (None otherwise), which check_engine
  appends to the looked_for text (a broken plugin import is then visible).
  """
  def lookup(params):
    lookup.error = None
    try: return bool(make_calculator().get_command())
    except Exception as e:
      lookup.error = repr(e)
      return False
  lookup.error = None
  return lookup

def _mopac():
  from qrefine.plugin.ase.mopac_qr import Mopac
  return Mopac()

def _xtb():
  from qrefine.plugin.ase.xtb_qr import GFNxTB
  return GFNxTB()

def _orca():
  from qrefine.plugin.ase.orca_qr import Orca
  return Orca()

def _terachem():
  from qrefine.plugin.ase.terachem_qr import TeraChem
  return TeraChem()

def _turbomole():
  from qrefine.plugin.ase.turbomole_qr import Turbomole
  return Turbomole()

def _gaussian_found(params):
  import shutil
  return bool(shutil.which("g16") or shutil.which("g09"))

def _importable(package):
  """Lookup: True if package imports (aimnet2calc warnings silenced)."""
  def lookup(params):
    import importlib, warnings
    with warnings.catch_warnings():
      warnings.filterwarnings("ignore", module="aimnet2calc")
      warnings.filterwarnings("ignore", message=".*PySisiphus.*")
      try: importlib.import_module(package)
      except Exception: return False
    return True
  return lookup

def _server_reachable(params):
  """HTTP GET of quantum.server_url (3 s): any HTTP answer counts."""
  import urllib.request, urllib.error
  try:
    urllib.request.urlopen(params.quantum.server_url, timeout=3).close()
  except urllib.error.HTTPError:
    return True # the server answered
  except Exception:
    return False
  return True

def _engine(needs_files, lookup, looked_for, hint, elements=None):
  return group_args(needs_files=needs_files, lookup=lookup,
    looked_for=looked_for, hint=hint, elements=elements)

_vendor = "Obtain it from its vendor and set %s"
_aimnet2_hint = "Install AIMNet2 with config/update_phenix.sh aimnet2"

# Engines of the quantum.engine_name choice. needs_files: the engine is an
# external program that reads and writes files (needs a working folder);
# lookup(params) -> bool; looked_for: what lookup looks for (for the Sorry);
# hint: install hint; elements: elements the engine supports (None: no
# restriction declared; aimnet2 reports its species at check time).
ENGINES = {
  "mopac": _engine(True, _command_found(_mopac),
    "Phenix MOPAC, mopac in PATH or MOPAC_COMMAND",
    "Install MOPAC from conda-forge (conda install -c conda-forge mopac) or "
    "set MOPAC_COMMAND"),
  "xtb": _engine(True, _command_found(_xtb),
    "xtb in PATH (or XTBHOME)",
    "Install xTB from conda-forge (conda install -c conda-forge xtb)"),
  "orca": _engine(True, _command_found(_orca),
    "ORCA_COMMAND or Orca_COMMAND", _vendor % "ORCA_COMMAND"),
  "terachem": _engine(True, _command_found(_terachem),
    "TeraChem_COMMAND", _vendor % "TeraChem_COMMAND"),
  "turbomole": _engine(True, _command_found(_turbomole),
    "TM_COMMAND", _vendor % "TM_COMMAND"),
  "gaussian": _engine(True, _gaussian_found,
    "g16 or g09 in PATH", "Obtain Gaussian from its vendor and put g16 or "
    "g09 in PATH"),
  "aimnet2": _engine(False, _importable("aimnet2calc"),
    "import of aimnet2calc", _aimnet2_hint),
  "aimnet2-old": _engine(False, _importable("aimnet2calc"),
    "import of aimnet2calc", _aimnet2_hint),
  # pyscf, torchani: unverified: the plugin sets a cwd-relative chkfile /
  # coordinates path at construction; until tested in-process, treat as
  # writing files.
  "torchani": _engine(True, _importable("torchani"),
    "import of torchani", "Install torchani (pip install torchani)",
    elements = {"H", "C", "N", "O"}),
  "pyscf": _engine(True, _importable("pyscf"),
    "import of pyscf", "Install PySCF (pip install pyscf)"),
  "fairchem": _engine(False, _importable("fairchem"),
    "import of fairchem", "Install fairchem (pip install fairchem-core)"),
  "server": _engine(False, _server_reachable,
    "an HTTP answer from quantum.server_url", "Start the server or set "
    "quantum.server_url to its address"),
}

def check_engine(params):
  """
  Sorry if the QM engine of quantum.engine_name is not available (program,
  Python package or server, see ENGINES). Nothing for restraints=cctbx.
  Reads only restraints, quantum.engine_name and quantum.server_url, so
  params need not be finalized.
  """
  if(params.restraints != "qm"): return
  name = params.quantum.engine_name
  e = ENGINES.get(name)
  if(e is None):
    raise Sorry("Unknown QM engine '%s'." % name)
  looked_for = e.looked_for
  if(name == "server"):
    looked_for = "%s at %s" % (looked_for, params.quantum.server_url)
  if(not e.lookup(params)):
    error = getattr(e.lookup, "error", None) # _command_found: the exception
    if(error is not None): looked_for = "%s (%s)" % (looked_for, error)
    raise Sorry("QM engine '%s' is not available: %s. %s" % (
      name, looked_for, e.hint))

# ---------------------------------------------------------------------------
# Model checks (check_model)
# ---------------------------------------------------------------------------

def _listing(items, n_max=10):
  """At most n_max items, comma separated, then ', ...'."""
  result = ", ".join(items[:n_max])
  if(len(items) > n_max): result += ", ..."
  return result

def _residue_label(rg, resname):
  """'<chain> <resname> <resseq><icode>'"""
  return "%s %s %s%s" % (rg.parent().id.strip(), resname.strip(),
    rg.resseq.strip(), rg.icode.strip())

def _residues_by_key(hierarchy):
  """
  Map of the residue ids missing_atoms_by_residue reports (atom group
  id_str, with the altloc replaced by any conformer label of the residue
  group) -> (residue group, resname).
  """
  result = {}
  for rg in hierarchy.residue_groups():
    labels = set([ag.altloc or " " for ag in rg.atom_groups()] + [" "])
    for ag in rg.atom_groups():
      s = ag.id_str() # altloc (or blank) + resname + residue group id
      for label in labels:
        result.setdefault(label + s[1:], (rg, ag.resname))
  return result

def _engine_elements(params):
  """Elements the engine supports (set of symbols), None if unrestricted."""
  name = params.quantum.engine_name
  if(name == "aimnet2"):
    import warnings
    with warnings.catch_warnings():
      warnings.filterwarnings("ignore", module="aimnet2calc")
      warnings.filterwarnings("ignore", message=".*PySisiphus.*")
      from qrefine.plugin.ase.aimnet2_qr import AIMNet2Calculator
      method = params.quantum.method
      model = "aimnet2-qr" if method == "rhf" else method
      calculator = AIMNet2Calculator(model)
    species = getattr(calculator, "implemented_species", None)
    if(species is None):
      # aimnet2calc looks for implemented_species on its base calculator;
      # the aimnet2-qr model file stores the buffer as impemented_species.
      m = getattr(getattr(calculator, "base_calc", None), "model", None)
      for attr in ("implemented_species", "impemented_species"):
        species = getattr(m, attr, None)
        if(species is not None):
          species = species.cpu().numpy()
          break
    if(species is None): return None
    from ase.data import chemical_symbols
    return set(chemical_symbols[int(z)] for z in species)
  e = ENGINES.get(name)
  if(e is None): return None
  return e.elements

def check_model(model, params, has_map=False):
  """
  Sorry if model cannot be refined with params. Runs on the caller's model
  (processed with restraints), which is not changed. In this order:
    - altlocs.handling=conformers with more than one conformer: the model
      must be completed (conformers.check_completed);
    - with restraints=qm:
      a. no residue with missing heavy atoms unless it carries a capping H
         (mmtbx.reduce2 qmready=true capping=true);
      b. no missing H, every water with two H (mmtbx.reduce2 qmready=true);
      c. only elements the engine supports (aimnet2, torchani);
      d. without a map and without expansion or clustering: no symmetry
         neighbours within cluster.select_within_radius
         (refine.check_qm_sees_crystal_neighbours).
  Offenders are listed (at most ten). params must be finalized.
  """
  assert getattr(params, "finalized", False), "call api.finalize_params first"
  from qrefine import conformers
  import iotbx.pdb
  hierarchy = model.get_hierarchy()
  if(params.altlocs.handling == "conformers" and
     conformers.n_conformers(hierarchy) > 1):
    conformers.check_completed(hierarchy)
  if(params.restraints != "qm"): return
  # Atom order of the copy is that of the model (selections apply).
  h = hierarchy.deep_copy()
  h.atoms().reset_i_seq()
  by_key = _residues_by_key(h)
  def labels(keys):
    result = []
    for k in keys:
      rg_rn = by_key.get(k)
      if(rg_rn is None): label = " ".join(k.split())
      else:              label = _residue_label(*rg_rn)
      if(label not in result): result.append(label)
    return result
  # a. truncated residues without caps
  missing_heavy = model.missing_atoms_by_residue(include_hydrogens=False)
  capping = model.selection("capping_h")
  uncapped = []
  for k in missing_heavy.keys():
    rg_rn = by_key.get(k)
    if(rg_rn is not None):
      i_seqs = rg_rn[0].atoms().extract_i_seq()
      if(capping.select(i_seqs).count(True) > 0): continue
    uncapped.append(k)
  if(len(uncapped) > 0):
    l = labels(uncapped)
    raise Sorry("Truncated residues without caps (%d): %s. With "
      "restraints=qm run mmtbx.reduce2 qmready=true capping=true on the "
      "model first." % (len(l), _listing(l)))
  # b. missing H; water with fewer than two H
  missing_all = model.missing_atoms_by_residue(include_hydrogens=True)
  missing_h = [k for k, names in missing_all.items()
    if set(names) - set(missing_heavy.get(k, []))]
  get_class = iotbx.pdb.common_residue_names_get_class
  def n_h(atom_groups):
    return sum([[a.element.strip().upper() in ("H", "D")
      for a in ag.atoms()].count(True) for ag in atom_groups])
  for rg in h.residue_groups():
    ags = rg.atom_groups()
    if(not [ag for ag in ags if get_class(name=ag.resname)=="common_water"]):
      continue
    # Per conformer: the blank-altloc group plus that conformer's group (O
    # blank and H/D in altloc groups is usual in neutron models).
    blank = [ag for ag in ags if not ag.altloc.strip()]
    altlocs = [ag for ag in ags if ag.altloc.strip()]
    conformers_h = [n_h(blank + [ag]) for ag in altlocs] or [n_h(blank)]
    if(min(conformers_h) < 2):
      missing_h.append(" " + ags[0].id_str()[1:]) # residue group label
  if(len(missing_h) > 0):
    l = labels(missing_h)
    raise Sorry("Residues or waters with missing hydrogen atoms (%d): %s. "
      "With restraints=qm run mmtbx.reduce2 qmready=true on the model "
      "first." % (len(l), _listing(l)))
  # c. elements
  supported = _engine_elements(params)
  if(supported is not None):
    present = set()
    for a in hierarchy.atoms():
      e = a.element.strip()
      # blank element column (old PDB files): element implied by the name
      if(not e): e = (a.determine_chemical_element_simple() or "").strip()
      if(not e): continue
      e = e.capitalize()
      if(e == "D"): e = "H"
      present.add(e)
    unsupported = sorted(present - set(supported))
    if(len(unsupported) > 0):
      raise Sorry("QM engine '%s' does not support the elements: %s. Use "
        "another engine or restraints=cctbx." % (
        params.quantum.engine_name, _listing(unsupported)))
  # d. crystal neighbours
  from qrefine import refine
  refine.check_qm_sees_crystal_neighbours(model=model, params=params,
    has_map=has_map)

# ---------------------------------------------------------------------------
# Refinement (refine)
# ---------------------------------------------------------------------------

def pdb_interpretation_params():
  """
  The model processing of qr.refine: neutron X-H distances, CDL off, atoms
  in input order.
  """
  import mmtbx.model
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  return params

def dictionaries_d_min(model, fmodel=None, map_data=None):
  """
  Resolution of the scattering dictionaries, as in phenix.refine: of the
  data (fmodel), else of the map (d99), else 1 A.
  """
  if(fmodel is not None):
    return fmodel.f_obs().d_min()
  if(map_data is not None):
    from cctbx import maptbx
    return maptbx.d99(map = map_data,
      crystal_symmetry = model.crystal_symmetry()).result.d99
  return 1.0

def _model_copy(model, fmodel, map_data):
  """
  qrefine's own model: built from a deep copy of the hierarchy of model
  with its crystal symmetry, restraint objects and secondary structure
  annotation, processed as qr.refine processes models, scattering
  dictionaries of model's table at dictionaries_d_min, the sites of model's
  structure.
  """
  import mmtbx.model
  result = mmtbx.model.manager(
    model_input       = None,
    pdb_hierarchy     = model.get_hierarchy().deep_copy(),
    crystal_symmetry  = model.crystal_symmetry(),
    restraint_objects = model.get_restraint_objects(),
    log               = null_out())
  ss_annotation = model.get_ss_annotation()
  if(ss_annotation is not None):
    result.set_ss_annotation(ann = ss_annotation.deep_copy())
  result.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = pdb_interpretation_params())
  result.setup_scattering_dictionaries(
    scattering_table = model.get_scattering_table(),
    d_min            = dictionaries_d_min(model=model, fmodel=fmodel,
                                          map_data=map_data))
  # Processing leaves sites that differ from model's in the last bits
  # (Cartesian-fractional round trips); take model's sites as they are, in
  # the hierarchy and in the structure.
  result.set_sites_cart(model.get_sites_cart().deep_copy())
  result.get_xray_structure().set_sites_frac(
    model.get_xray_structure().sites_frac().deep_copy())
  return result

def _needs_scratch_folder(params):
  """
  True if refinement with params cannot run without files: restraints=qm
  with an external-program engine (ENGINES needs_files) or with a
  quantum.qm_addon (stand-alone programs).
  """
  if(params.restraints != "qm"): return False
  e = ENGINES.get(params.quantum.engine_name)
  if(e is not None and e.needs_files): return True
  return params.quantum.qm_addon is not None

def refine(model, params, fmodel=None, map_data=None, working_folder="ase",
           intermediate_folder=None, bulk_solvent_params=None, log=None):
  """
  Refine (fmodel or map_data given) or optimize (neither) the sites of
  model with qrefine restraints. Neither model nor fmodel is changed: the
  refinement runs on qrefine's own model (_model_copy) and a deep copy of
  fmodel holding a copy of that model's structure. fmodel's scatterers must
  be in model's atom order (Sorry otherwise); its target is the caller's
  (api.resolve_target_name). params: after finalize_params (asserted).
  The fmodel's target is used as the caller set it:
  params.refine.refinement_target_name is not consulted here (qr.py resolves
  it with resolve_target_name before calling). Twinned fmodels (twin_law
  set) are refused (Sorry).

  working_folder: engine scratch, fragment, expansion and conformer files
  (created by the writer that needs it). None: no file is written; an
  external-program engine (ENGINES needs_files) then runs in a temporary
  folder that is removed on return. intermediate_folder: per-cycle models
  unless None. bulk_solvent_params: passed to every update_all_scales call
  on the fmodel copy (None: cctbx defaults). log: every print (None:
  nothing).

  Returns group_args(sites_cart, sites_cart_nuclear, r_work, r_free,
  data_weight, restraints_weight_scale, b_rmsd, a_rmsd, n_fev, n_cycles,
  cc_mask, mode, model), see refine.run; None with refine.dry_run.
  sites_cart are in the convention of the scattering table (X-ray X-H
  lengths with an X-ray table), sites_cart_nuclear the nuclear ones; model
  is qrefine's model (nuclear sites, atom order of the input).
  """
  from qrefine import refine as refine_module
  if(log is None): log = null_out()
  assert getattr(params, "finalized", False), "call api.finalize_params first"
  if(fmodel is not None and getattr(fmodel, "twin_law", None) is not None):
    raise Sorry("qrefine does not support twinned refinement "
      "(fmodel.twin_law is set).")
  if(fmodel is not None):
    refine_module.check_model_fmodel_order(model = model, fmodel = fmodel)
  own_model = _model_copy(model=model, fmodel=fmodel, map_data=map_data)
  own_fmodel = None
  if(fmodel is not None):
    own_fmodel = fmodel.deep_copy()
    own_fmodel.update_xray_structure(
      xray_structure = own_model.get_xray_structure().deep_copy_scatterers(),
      update_f_calc  = True,
      update_f_mask  = True)
  scratch = None
  # The copy of the working folder on a failure of the fragment jobs is for
  # a folder the caller chose, never for the temporary one.
  keep_error_copy = working_folder is not None
  if(working_folder is None and _needs_scratch_folder(params)):
    scratch = tempfile.mkdtemp(prefix="qrefine_")
    print("scratch folder: %s (removed on return)" % scratch, file=log)
    working_folder = scratch
  try:
    return refine_module.run(
      model               = own_model,
      fmodel              = own_fmodel,
      map_data            = map_data,
      params              = params,
      log                 = log,
      working_folder      = working_folder,
      intermediate_folder = intermediate_folder,
      bulk_solvent_params = bulk_solvent_params,
      keep_error_copy     = keep_error_copy)
  finally:
    if(scratch is not None):
      # The one recursive delete in qrefine: the folder this call created.
      shutil.rmtree(scratch, ignore_errors=True)
