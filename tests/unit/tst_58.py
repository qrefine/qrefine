from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os, glob, shutil
import iotbx.pdb
import mmtbx.model
from libtbx import easy_run
from libtbx.utils import Sorry, null_out
from qrefine import qr, api

def finalized(**kw):
  """
  Default qr.refine parameters with the given settings (engine sets
  quantum.engine_name and restraints=qm), finalized for X-ray data at 2 A.
  """
  p = qr.get_default_params()
  if("engine" in kw):
    p.quantum.engine_name = kw.pop("engine")
    p.restraints = "qm"
  for k, v in kw.items():
    if(k == "clustering"): p.cluster.clustering = v
    else: setattr(p, k, v)
  return api.finalize_params(p, n_conformers=1, d_min=2.0, has_data=True,
    log=null_out())

def assert_sorry(f, *args, **kw):
  """Call f, expect a Sorry whose text contains every string in kw["has"]."""
  has = kw.pop("has")
  if(not isinstance(has, (list, tuple))): has = [has]
  try: f(*args, **kw)
  except Sorry as e:
    for s in has: assert str(e).find(s) > -1, (s, str(e))
    return str(e)
  raise AssertionError("Expected Sorry containing %s" % str(has))

def get_model(file_name=None, lines=None):
  """mmtbx.model.manager processed with restraints."""
  if(file_name is not None):
    pdb_inp = iotbx.pdb.input(file_name=file_name)
  else:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=lines)
  m = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  m.process(make_restraints=True)
  return m

def variant(extra, file_name):
  """
  capped.pdb with the HETATM lines extra appended, written to file_name with
  write_pdb_or_mmcif_file. Returns the file name.
  """
  lines = [l for l in open("capped.pdb").read().splitlines()
    if not l.startswith(("END", "TER"))] + extra
  ph = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  cs = iotbx.pdb.input(file_name="capped.pdb").crystal_symmetry()
  return ph.write_pdb_or_mmcif_file(target_filename=file_name,
    crystal_symmetry=cs)

def pdbtools(args, file_name):
  """phenix.pdbtools capped.pdb <args>, output written to file_name."""
  assert easy_run.call("phenix.pdbtools capped.pdb %s output.filename=%s "
    "overwrite=true > %s.log" % (args, file_name, file_name)) == 0
  return file_name

def test_check_engine():
  p = finalized(restraints="cctbx"); api.check_engine(p)
  p = finalized(engine="mopac"); api.check_engine(p)
  if(shutil.which("xtb") is None and "XTBHOME" not in os.environ):
    p = finalized(engine="xtb")
    assert_sorry(api.check_engine, p, has=["not available", "xtb"])
  try: import torchani
  except ImportError:
    p = finalized(engine="torchani")
    assert_sorry(api.check_engine, p,
      has="QM engine 'torchani' is not available: import of torchani.")
  p = finalized(engine="server"); p.quantum.server_url = "http://127.0.0.1:1"
  assert_sorry(api.check_engine, p, has="server")
  # every engine of the PHIL choice has an entry
  for name in ("mopac fairchem aimnet2 aimnet2-old torchani terachem "
               "turbomole pyscf orca gaussian xtb server").split():
    e = api.ENGINES[name]
    assert e.needs_files == (name in ("mopac", "xtb", "orca", "terachem",
      "turbomole", "gaussian", "pyscf", "torchani")), name
  # params need not be finalized
  p = qr.get_default_params(); p.restraints = "cctbx"
  assert api.check_engine(p) is None
  p = qr.get_default_params(); p.restraints = "qm"
  p.quantum.engine_name = "mopac"
  assert api.check_engine(p) is None

def test_check_model_capped_ok():
  api.check_model(get_model("capped.pdb"), finalized(engine="mopac"))

def test_check_model_uncapped():
  m = get_model(pdbtools("remove=capping_h", "uncapped.pdb"))
  assert_sorry(api.check_model, m, finalized(engine="mopac"),
    has=["A ARG 3", "without caps", "reduce2 qmready=true"])

def test_check_model_missing_h():
  m = get_model(pdbtools('remove="element H"', "no_h.pdb"))
  assert_sorry(api.check_model, m, finalized(engine="mopac"),
    has="reduce2 qmready=true")
  # H removed from SER 4 only: the cap of ARG 3 stays, missing H reported
  m = get_model(pdbtools('remove="element H and resseq 4"', "no_h_4.pdb"))
  e = assert_sorry(api.check_model, m, finalized(engine="mopac"),
    has=["A SER 4", "missing hydrogen", "reduce2 qmready=true"])
  assert e.find("without caps") < 0 and e.find("ARG") < 0, e

def test_water_no_h():
  f = variant([
    "HETATM   25  O   HOH S   1      35.000  35.000  35.000  1.00  8.00"
    "           O"], "water_no_h.pdb")
  assert_sorry(api.check_model, get_model(f), finalized(engine="mopac"),
    has=["S HOH 1", "reduce2 qmready=true"])

def test_water_one_h():
  f = variant([
    "HETATM   25  O   HOH S   1      35.000  35.000  35.000  1.00  8.00"
    "           O",
    "HETATM   26  H1  HOH S   1      35.957  35.000  35.000  1.00  8.00"
    "           H"], "water_one_h.pdb")
  assert_sorry(api.check_model, get_model(f), finalized(engine="mopac"),
    has=["S HOH 1", "reduce2 qmready=true"])

def water_atom(serial, name, altloc, x, element):
  """HETATM line of water S 1 (occupancy 0.5 with an altloc, else 1)."""
  occ = "0.50" if altloc.strip() else "1.00"
  return ("HETATM%5d %-4s%sHOH S   1    %8.3f  35.000  35.000  %s  8.00"
    "          %2s" % (serial, " "+name if len(name) < 4 else name, altloc,
    x, occ, element))

def test_water_altloc_h():
  """
  The water rule counts H/D per conformer (blank group plus the conformer's
  group): O blank with D1/D2 in A and in B passes; O blank with two H in A
  and one in B fails (three H in the residue group); O and one H in the
  blank group fails. (A B group without H cannot be built: a B group needs
  an atom and O in B next to a blank O is refused by processing.)
  """
  o = water_atom(25, "O", " ", 35.0, "O")
  f = variant([o,
    water_atom(26, "D1", "A", 35.957, "D"),
    water_atom(27, "D2", "A", 34.043, "D"),
    water_atom(28, "D1", "B", 35.900, "D"),
    water_atom(29, "D2", "B", 34.100, "D")], "water_altloc_d.pdb")
  api.check_model(get_model(f), finalized(engine="mopac"))
  f = variant([o,
    water_atom(26, "H1", "A", 35.957, "H"),
    water_atom(27, "H2", "A", 34.043, "H"),
    water_atom(28, "H1", "B", 35.900, "H")], "water_altloc_b_one_h.pdb")
  assert_sorry(api.check_model, get_model(f), finalized(engine="mopac"),
    has=["S HOH 1", "missing hydrogen"])
  f = variant([o, water_atom(26, "H1", " ", 35.957, "H")],
    "water_blank_one_h.pdb")
  assert_sorry(api.check_model, get_model(f), finalized(engine="mopac"),
    has=["S HOH 1", "missing hydrogen"])

def test_check_model_elements():
  try:
    import aimnet2calc
  except Exception:
    print("test_check_model_elements: aimnet2calc not importable, skipped")
    return
  f = variant([
    "HETATM   25 ZN    ZN Z   1      10.000  10.000  10.000  1.00  8.00"
    "          ZN"], "zn.pdb")
  p = finalized(engine="aimnet2")
  assert p.restraints == "qm"
  api.check_model(get_model("capped.pdb"), p)
  assert_sorry(api.check_model, get_model(f), p,
    has=["Zn", "restraints=cctbx"])
  # blank element column (old PDB files; processing fills it in, so blank
  # it on the processed model): the element comes from the atom name, passes
  m = get_model("capped.pdb")
  og = [a for a in m.get_hierarchy().atoms() if a.name == " OG "]
  assert len(og) == 1
  og[0].element = "  "
  api.check_model(m, p)

def test_altloc_completion():
  from qrefine.tests.unit import tst_50, altloc_scenarios
  altloc_scenarios.write(tst_50.build(), "altloc_in.pdb")
  m = get_model("altloc_in.pdb")
  from qrefine import conformers
  assert conformers.n_conformers(m.get_hierarchy()) > 1
  p = finalized(restraints="cctbx")
  assert p.altlocs.handling == "conformers"
  assert_sorry(api.check_model, m, p, has="Altloc completion")

def test_check_model_cctbx_lenient():
  m = get_model(pdbtools("remove=capping_h", "uncapped.pdb"))
  api.check_model(m, finalized(restraints="cctbx"))

def test_crystal_neighbours():
  from mmtbx.hydrogens.tst_xh_mapper import pdb_str
  lines = pdb_str.splitlines()
  assert lines[0].startswith("CRYST1")
  lines[0] = "CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 1"
  m = get_model(lines=lines)
  p = finalized(engine="mopac", expansion=False, clustering=False)
  assert_sorry(api.check_model, m, p, has="Symmetry-related")
  api.check_model(m, p, has_map=True)

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  api.check_engine (engine availability) and api.check_model (altloc
  completion; with restraints=qm: caps, missing H and water H, elements the
  engine supports, crystal neighbours) raise Sorry naming the offenders and
  the preparation command.
  Fixture: the tst_173 ARG/SER fragment (ARG 3 truncated at CB) through
  mmtbx.reduce2 qmready=true capping=true (capped.pdb, one cap HCG).
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  from qrefine.tests.unit.tst_53 import pdb_str
  open("in.pdb", "w").write(pdb_str)
  assert easy_run.call(
    "mmtbx.reduce2 qmready=true capping=true in.pdb > reduce2.log") == 0
  capped = glob.glob("inH_charge*.pdb")
  assert len(capped) == 1, capped
  shutil.copyfile(capped[0], "capped.pdb")
  test_check_engine()
  test_check_model_capped_ok()
  test_check_model_uncapped()
  test_check_model_missing_h()
  test_water_no_h()
  test_water_one_h()
  test_water_altloc_h()
  test_check_model_elements()
  test_altloc_completion()
  test_check_model_cctbx_lenient()
  test_crystal_neighbours()

if(__name__ == "__main__"):
  run()
  print("OK")
