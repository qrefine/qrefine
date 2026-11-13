from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os, glob
import iotbx.pdb
import mmtbx.model
from libtbx import easy_run
from libtbx.utils import null_out
from qrefine.tests.unit import run_qrefine_in, run_fmodel
from qrefine.tests.unit import fast_args, fast_args_no_weight_search

# mmtbx/hydrogens/tst_capping.py, pdb_arg: ARG 3 truncated at CB
pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A   3      27.911  23.777  26.214  1.00  5.99           N
ATOM      2  CA  ARG A   3      27.558  25.104  25.723  1.00  6.90           C
ATOM      3  C   ARG A   3      26.047  25.238  25.567  1.00  6.28           C
ATOM      4  O   ARG A   3      25.450  26.196  26.055  1.00  7.24           O
ATOM      5  CB  ARG A   3      28.246  25.375  24.384  1.00  6.79           C
ATOM      6  N   SER A   4      25.443  24.270  24.884  1.00  5.99           N
ATOM      7  CA  SER A   4      24.000  24.258  24.649  1.00  6.90           C
ATOM      8  C   SER A   4      23.547  25.461  23.816  1.00  6.28           C
ATOM      9  O   SER A   4      23.472  25.372  22.590  1.00  7.24           O
ATOM     10  OXT SER A   4      23.253  26.536  24.382  1.00  6.28           O
ATOM     11  CB  SER A   4      23.234  24.186  25.975  1.00  6.79           C
ATOM     12  OG  SER A   4      21.834  24.229  25.757  1.00  8.89           O
END
"""

def capping_atoms(file_name):
  """
  Sorted (resseq, name) of the atoms selected by capping_h.
  """
  m = mmtbx.model.manager(
    model_input = iotbx.pdb.input(file_name = file_name), log = null_out())
  atoms = m.get_hierarchy().atoms()
  return sorted([(atoms[i].parent().parent().resseq.strip(),
    atoms[i].name.strip()) for i in m.selection("capping_h").iselection()])

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine on a model capped by mmtbx.reduce2 qmready=true capping=true
  (ARG 3 truncated at CB: one cap, HCG; HB2 and HB3 are dictionary H).
  Asserted: the capped model has exactly the cap HCG on ARG 3; qr.refine
  loads and refines it with cctbx restraints and with MOPAC
  (clustering=False, so the total charge is computed on the whole capped
  model), without a Sorry; the refined models keep exactly the same caps,
  HCG at occupancy 0; the MOPAC input carries the total charge from the
  reduce2 file name; phenix.pdbtools remove=capping_h strips HCG and keeps
  HB2 and HB3.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  open("in.pdb", "w").write(pdb_str)
  assert easy_run.call(
    "mmtbx.reduce2 qmready=true capping=true in.pdb > reduce2.log") == 0
  # qmready=true adds the total charge to the default output file name
  capped = glob.glob("inH_charge*.pdb")
  assert len(capped) == 1, capped
  capped = capped[0]
  charge = int(capped.split("charge")[1].split(".pdb")[0])
  print("capped model:", capped, "total charge:", charge)
  caps_before = capping_atoms(capped)
  assert caps_before == [("3", "HCG")], caps_before
  r = run_fmodel(args=[capped, "high_res=2.5", "r_free=0.05"], prefix="data")
  for folder, engine_args in [
      ("cctbx", ["restraints=cctbx"] + fast_args),
      ("mopac", ["restraints=qm", "quantum.engine_name=mopac",
                 "clustering=False"] + fast_args_no_weight_search)]:
    out = run_qrefine_in(args=[capped, r.mtz] + engine_args, prefix=prefix,
      folder=folder)
    assert out.pdb is not None, folder
    caps = capping_atoms(out.pdb)
    assert caps == caps_before, (folder, caps)
    occ = [a.occ for a in out.model.get_hierarchy().atoms()
      if a.name.strip() == "HCG"]
    assert occ == [0.0], (folder, occ)
    if(folder == "mopac"):
      # plugin/ase/mopac_qr.py writes 'CHARGE= ' + str(charge)
      mop = glob.glob(os.path.join(out.folder, "ase", "*.mop"))
      assert len(mop) == 1, mop
      header = open(mop[0]).readline()
      print("MOPAC keywords:", header.strip())
      assert ("CHARGE= %d " % charge) in header, header
  cmd = " ".join(["phenix.pdbtools", out.pdb, "remove=capping_h",
    "output.file_name=stripped.pdb", "> pdbtools.log"])
  assert easy_run.call(cmd) == 0
  stripped = iotbx.pdb.input(file_name="stripped.pdb").construct_hierarchy()
  names = [(a.parent().parent().resseq.strip(), a.name.strip())
    for a in stripped.atoms()]
  n_in = len(iotbx.pdb.input(file_name=out.pdb).construct_hierarchy().atoms())
  assert ("3", "HCG") not in names, names
  assert ("3", "HB2") in names and ("3", "HB3") in names, names
  assert len(names) == n_in - 1, (len(names), n_in)

if(__name__ == "__main__"):
  run()
  print("OK")
