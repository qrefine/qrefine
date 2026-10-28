from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import shutil
import iotbx.pdb
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from qrefine import qr, refine, conformers
from qrefine.tests.unit.tst_31 import duplicate
from qrefine.tests.unit.tst_11 import get_model, process, qr_unit_tests
from qrefine.tests.unit import run_qrefine_in

def get_params(expansion, clustering):
  params = qr.get_default_params()
  params.restraints = "qm"
  params.quantum.engine_name = "mopac"
  refine.set_qm_defaults(params, log=null_out())
  params.expansion = expansion
  params.cluster.clustering = clustering
  params.cluster.maxnum_residues_in_cluster = 2
  params.parallel.nproc = 1
  return params

def exercise_opt_reclustering():
  """
  qr.refine mode=opt restraints=qm (MOPAC) with clustering on gly2_1.pdb
  (A/B conformers): re-clustering goes through the conformer fragment
  manager proxy, scratch folders ase/A, ase/B are used, the output model keeps
  its altlocs and atoms.
  """
  pdb = os.path.join(qr_unit_tests, "data_files", "gly2_1.pdb")
  r_qr = run_qrefine_in(folder="qrefine_opt", prefix="opt", args=[pdb,
    "mode=opt", "restraints=qm", "quantum.engine_name=mopac",
    "clustering=True", "maxnum_residues_in_cluster=1",
    "number_of_micro_cycles=2", "max_iterations_refine=3",
    "cluster.re_calculate_rmsd_tolerance=1e-5", "parallel.nproc=1"])
  assert r_qr.std_out_text.find("need to redo clustering") > -1
  for label in ["A", "B"]:
    assert os.path.isdir(os.path.join(r_qr.folder, "ase", label, "0"))
  assert r_qr.model is not None
  h_in = iotbx.pdb.input(file_name=pdb).construct_hierarchy()
  h = r_qr.model.get_hierarchy()
  assert h.atoms_size() == h_in.atoms_size()
  assert sorted(h.altloc_indices().keys()) == ["", "A", "B"]
  assert [a.parent().altloc for a in h.atoms()] == \
    [a.parent().altloc for a in h_in.atoms()]

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Duplicate-conformer identity with MOPAC: a single-conformer model
  duplicated into identical conformers gives the single-conformer gradients
  through plain, expansion and clustering (three labels: plain only).
  Conformer scratch folders ase/<label> are used.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  # Folder checks below must not pass on files of an earlier run.
  for folder in ["ase", "qrefine_opt"]:
    if(os.path.isdir(folder)): shutil.rmtree(folder)
  #
  single = get_model("m2_complete_box_large.pdb")
  for labels, routes in [
      (["A","B"],     [(False, False), (True, False), (False, True)]),
      (["A","B","C"], [(False, False)])]:
    dup, original = duplicate(model=single, labels=labels,
      residue_indices=[2,3,4])
    process(dup)
    assert conformers.conformer_labels(dup.get_hierarchy()) == labels
    for expansion, clustering in routes:
      params = get_params(expansion=expansion, clustering=clustering)
      rm_single = refine.create_restraints_manager(params=params, model=single)
      rm_dup = refine.create_restraints_manager(params=params, model=dup)
      assert not isinstance(rm_single, conformers.from_conformers)
      assert isinstance(rm_dup, conformers.from_conformers)
      assert len(rm_dup.conformers) == len(labels)
      g_single = rm_single.target_and_gradients(
        sites_cart=single.get_sites_cart())[1]
      g_dup = rm_dup.target_and_gradients(sites_cart=dup.get_sites_cart())[1]
      assert g_dup.size() == dup.get_number_of_atoms()
      assert approx_equal(g_dup, g_single.select(original), 1.e-4)
  #
  assert os.path.isfile(os.path.join("ase", "A", "tmp_ase.pdb"))
  assert os.path.isdir(os.path.join("ase", "B", "0")) # clustering fragments
  assert os.path.isfile(os.path.join("ase", "C", "tmp_ase.pdb"))
  #
  exercise_opt_reclustering()

if(__name__ == "__main__"):
  run()
