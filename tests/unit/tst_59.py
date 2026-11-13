from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os, shutil
from mmtbx.hydrogens.tst_xh_mapper import pdb_str
from qrefine.tests.unit import run_qrefine_in, run_fmodel, fast_args

common_args = ["restraints=cctbx", "expansion=False"] + fast_args

def inputs():
  """Absolute paths of the fixture model and data (in the test folder)."""
  return [os.path.abspath("in.pdb"), os.path.abspath("data.mtz")]

def test_no_cwd_writes():
  """
  A command-line run writes nothing into the folder it runs in except the
  pdb/ folder (intermediate and final models); <prefix>.log is the run's
  stdout captured by run_qrefine_in. A cctbx run without clustering writes
  no engine scratch, so the working folder ase/ is not created at all.
  """
  os.makedirs("run", exist_ok=True)
  for f in inputs():
    shutil.copy(f, "run")
  args = [os.path.abspath(os.path.join("run", os.path.basename(f)))
    for f in inputs()]
  r = run_qrefine_in(args + common_args + ["clustering=False"], "in",
    folder="run")
  listing = sorted(os.listdir(r.folder))
  expected = sorted(["in.pdb", "data.mtz", "in.log", "pdb"])
  assert listing == expected, listing
  assert r.pdb is not None
  assert os.path.isfile(os.path.join(r.folder, "pdb", "0_refine_cycle.pdb"))

def test_dump_gradients_sorry():
  """
  dump_gradients=g.pkl dumps the three gradient arrays into the working
  folder and stops with a Sorry (not a zero exit).
  """
  r = run_qrefine_in(inputs() + common_args + ["clustering=False",
    "dump_gradients=g.pkl"], "in", folder="dump", sorry_expected=True)
  assert "gradients dumped to" in r.std_out_text, r.std_out_text
  for s in ["dg", "rg", "g"]:
    fn = os.path.join(r.folder, "ase", "g.pkl_%s" % s)
    assert os.path.isfile(fn), fn
    assert not os.path.exists(os.path.join(r.folder, "g.pkl_%s" % s))

def test_expansion_file_in_working_folder():
  """
  clustering=True (cctbx fragments, no QM): the expansion file is written
  into the working folder, not into the current directory. The fixture has
  conformers A and B: the altloc model is clustered (per-fragment conformer
  variants), one expansion file ase/expansion.pdb.
  """
  r = run_qrefine_in(inputs() + common_args + ["clustering=True"], "in",
    folder="cluster")
  assert "per-fragment variants" in r.std_out_text
  assert "Conformer A:" not in r.std_out_text
  fn = os.path.join(r.folder, "ase", "expansion.pdb")
  assert os.path.isfile(fn), fn
  assert not os.path.exists(os.path.join(r.folder, "ase", "A"))
  assert not os.path.exists(os.path.join(r.folder, "expansion.pdb"))
  assert sorted(os.listdir(r.folder)) == ["ase", "in.log", "pdb"]

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Where qr.refine writes: only pdb/ and ase/ in the run folder; dump_gradients
  into ase/ and a Sorry; expansion.pdb into ase/.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  with open("in.pdb", "w") as fo:
    fo.write(pdb_str)
  run_fmodel(args=["in.pdb", "high_res=2.5", "r_free=0.05"], prefix="data")
  test_no_cwd_writes()
  test_dump_gradients_sorry()
  test_expansion_file_in_working_folder()

if(__name__ == "__main__"):
  run()
  print("OK")
