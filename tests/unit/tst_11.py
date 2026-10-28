from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import shutil
import iotbx.pdb
import libtbx.load_env
import mmtbx.model
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx.array_family import flex
from qrefine import qr, refine, conformers
from qrefine.tests.unit.tst_31 import duplicate

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests = os.path.join(qrefine, "tests","unit")

def process(model):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  params.pdb_interpretation.sort_atoms = False
  params.pdb_interpretation.const_shrink_donor_acceptor=0.6
  model.process(make_restraints=True, grm_normalization=False,
    pdb_interpretation_params = params)
  return model

def get_model(file_name):
  file_name = os.path.join(qr_unit_tests,"data_files",file_name)
  pdb_inp = iotbx.pdb.input(file_name)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  return process(model)

def get_restraints_manager(file_name, expansion, clustering, handling,
                           model=None):
  """
  Restraints manager and the (processed) model it is built for; model, if
  given, is used instead of reading file_name.
  """
  if(model is None): model = get_model(file_name=file_name)
  params = qr.get_default_params()
  params.restraints = "cctbx"
  params.expansion = expansion
  params.cluster.clustering = clustering
  params.altlocs.handling = handling
  params.cluster.maxnum_residues_in_cluster = 999
  params.parallel.nproc = 1
  return refine.create_restraints_manager(params=params, model=model), model

def exercise_routes_agree():
  """
  clustering vs expansion with altlocs, both handlings (old tst_25).
  """
  for handling in ["conformers", "native"]:
    for fn in ["gly2_1.pdb", "gly2_2.pdb", "altlocs2.pdb", "altlocs.pdb",
               "h_altconf_complete.pdb"]:
      rm1, model = get_restraints_manager(fn, expansion=False,
        clustering=True, handling=handling)
      rm2, _ = get_restraints_manager(fn, expansion=True,
        clustering=False, handling=handling)
      assert conformers.n_conformers(model.get_hierarchy()) == 2
      assert isinstance(rm1, conformers.from_conformers) == \
        (handling == "conformers")
      sites_cart = model.get_sites_cart()
      g1 = rm1.target_and_gradients(sites_cart=sites_cart)[1]
      g2 = rm2.target_and_gradients(sites_cart=sites_cart)[1]
      assert g1.size() == sites_cart.size()
      diff = flex.max(flex.abs(g1.as_double()-g2.as_double()))
      assert diff < 1.e-6, [handling, fn, diff]

def exercise_duplicate_routes():
  """
  A single-conformer model duplicated into identical conformers gives the
  single-conformer gradients through plain, expansion and clustering.
  """
  pdb_inp = iotbx.pdb.input(os.path.join(qr_unit_tests, "data_files",
    "h_altconf_complete.pdb"))
  ph = pdb_inp.construct_hierarchy()
  ph.remove_alt_confs(always_keep_one_conformer=True)
  single = mmtbx.model.manager(model_input=None, pdb_hierarchy=ph,
    crystal_symmetry=pdb_inp.crystal_symmetry(), log=null_out())
  process(single)
  for labels in [["A","B"], ["A","B","C"]]:
    dup, original = duplicate(model=single, labels=labels,
      residue_indices=[2,3,4])
    process(dup)
    assert conformers.conformer_labels(dup.get_hierarchy()) == labels
    for expansion, clustering in [(False, False), (True, False), (False, True)]:
      rm_single, _ = get_restraints_manager(None, expansion=expansion,
        clustering=clustering, handling="conformers", model=single)
      rm_dup, _ = get_restraints_manager(None, expansion=expansion,
        clustering=clustering, handling="conformers", model=dup)
      assert not isinstance(rm_single, conformers.from_conformers)
      assert isinstance(rm_dup, conformers.from_conformers)
      assert len(rm_dup.conformers) == len(labels)
      g_single = rm_single.target_and_gradients(
        sites_cart=single.get_sites_cart())[1]
      g_dup = rm_dup.target_and_gradients(sites_cart=dup.get_sites_cart())[1]
      assert approx_equal(g_dup, g_single.select(original))

def exercise_folders():
  """
  Conformer scratch folders ase/A and ase/B are created for clustering; a
  model without altlocs keeps ase/ (no conformer subfolders).
  """
  cwd = os.getcwd()
  for name in ["folders_altlocs", "folders_no_altlocs"]:
    shutil.rmtree(name, ignore_errors=True) # no files of an earlier run
    os.makedirs(name)
    os.chdir(name)
    if(name == "folders_altlocs"):
      rm, model = get_restraints_manager("h_altconf_complete.pdb",
        expansion=False, clustering=True, handling="conformers")
      assert isinstance(rm, conformers.from_conformers)
      assert rm.fragment_manager is not None
      rm.target_and_gradients(sites_cart=model.get_sites_cart())
      for label in ["A", "B"]:
        assert os.path.isdir(os.path.join("ase", label))
        assert os.path.isfile(os.path.join("ase", label, "tmp_ase.pdb"))
      assert not os.path.exists(os.path.join("ase", "tmp_ase.pdb"))
    else:
      rm, model = get_restraints_manager("m00_good.pdb",
        expansion=False, clustering=True, handling="conformers")
      assert not isinstance(rm, conformers.from_conformers)
      rm.target_and_gradients(sites_cart=model.get_sites_cart())
      assert os.path.isfile(os.path.join("ase", "tmp_ase.pdb"))
      assert sorted(os.listdir("ase")) == ["tmp_ase.pdb"]
    os.chdir(cwd)

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  create_restraints_manager with altlocs: conformers vs native handling,
  all routes, conformer scratch folders.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  exercise_folders()
  exercise_duplicate_routes()
  exercise_routes_agree()

if(__name__ == "__main__"):
  run()
