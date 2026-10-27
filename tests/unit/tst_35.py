from __future__ import division
from __future__ import absolute_import

import os
import glob
import iotbx.pdb
import libtbx.load_env
from libtbx import easy_pickle
from qrefine.command_line import granalyse
from qrefine.tests.unit import run_qrefine

qrefine = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = os.path.join(qrefine, "tests", "unit", "data_files")

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Gradients from clustered QM (mopac) match gradients from the QM calculation
  on the whole molecule. A 7-residue peptide in a large box (no crystal
  neighbours) is split into 4 clusters of at most 2 residues; each fragment is
  the cluster plus one buffer residue on either side, and the cut peptide bonds
  are capped with hydrogens.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  #
  pdb_in = os.path.join(qr_unit_tests_data, "m2_complete_box_large.pdb")
  n_atoms = iotbx.pdb.input(file_name = pdb_in).atoms().size()
  #
  common = [
    pdb_in,
    "mode=opt",
    "stpmax=0.2",
    "gradient_only=true",
    "minimizer=lbfgs",
    "restraints=qm",
    "quantum.engine_name=mopac",
    "debug=True",
    
    #"distance_buffer=true",    # This makes d less than 2
    #"distance_buffer_radius=7",# This makes d less than 2
    
    "number_of_micro_cycles=1",
    "max_iterations_refine=1",
    "parallel.nproc=4"]
  r_qr = run_qrefine(prefix = prefix, args=common+[
    "clustering=False", "expansion=False"])
  r_qr = run_qrefine(prefix = prefix, args=common+[
    "clustering=True", "maxnum_residues_in_cluster=2", "two_buffers=False"])
  #
  # The model was actually split into fragments and the cuts were capped.
  #
  fragment_dirs = sorted(glob.glob(os.path.join("ase", "[0-9]*")))
  assert len(fragment_dirs) == 4, fragment_dirs
  for d in fragment_dirs:
    i = os.path.basename(d)
    raw    = iotbx.pdb.input(file_name = os.path.join(d, "%s.pdb"%i))
    capped = iotbx.pdb.input(file_name = os.path.join(d, "%s_capping.pdb"%i))
    n_raw, n_capped = raw.atoms().size(), capped.atoms().size()
    print("    fragment %s: %d atoms, %d with caps"%(i, n_raw, n_capped))
    assert n_raw < n_atoms          # a real fragment, not the whole molecule
    assert n_capped > n_raw         # cut bonds were capped
    elements = [a.element.strip().upper() for a in capped.atoms()]
    assert elements[n_raw:].count("H") == n_capped-n_raw  # with hydrogens
  #
  # Gradients at the starting coordinates match. Later gradients are not
  # comparable: the two optimizations diverge after the first step.
  #
  g1 = easy_pickle.load("debug_optimization_clusteringFalse.pkl").gradients[0]
  g2 = easy_pickle.load("debug_optimization_clusteringTrue.pkl").gradients[0]
  assert g1.size() == n_atoms and g2.size() == n_atoms
  d = max(granalyse.get_grad_wdelta(ref = g1, g = g2))
  print("    get_grad_wdelta (cluster vs no cluster): %.2f"%d)
  assert d < 10, d

if(__name__ == "__main__"):
  run()
