from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import iotbx.pdb
import iotbx.mrcfile
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from qrefine import conformers
from qrefine.tests.unit import run_qrefine_in
from qrefine.tests.unit.tst_06 import write_inputs, get_model

def write_map(pdb_file, file_name):
  """
  Sigma-scaled 2 A Fcalc map of the model, CCP4 format.
  """
  xrs = iotbx.pdb.input(file_name=pdb_file).xray_structure_simple()
  fft_map = xrs.structure_factors(d_min=2.0).f_calc().fft_map(
    resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  iotbx.mrcfile.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = xrs.unit_cell(),
    space_group = xrs.space_group(),
    map_data    = fft_map.real_map_unpadded(),
    labels      = flex.std_string([""]))
  return os.path.abspath(file_name)

def cc_mask_values(text):
  result = []
  for line in text.splitlines():
    if(line.strip().startswith("CC_mask  :")):
      result.append(float(line.split()[-1]))
  return result

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  Real-space refinement (model + map) of a model with alternative
  conformations, cctbx restraints, altlocs.handling=conformers: runs, keeps
  altlocs and occupancies, map-model CC improves.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb_good, pdb_poor = write_inputs(prefix)
  map_file = write_map(pdb_good, "%s_good.ccp4" % prefix)
  # data_weight=10: with the default start weight (1) the restraints dominate
  # on this sigma-scaled map and the CC drifts down while the geometry is
  # idealized (0.959 -> 0.954); with 10 it goes 0.959 -> 0.962.
  r_qr = run_qrefine_in(prefix=prefix, args=[pdb_poor,
    map_file, "restraints=cctbx", "clustering=False",
    "number_of_weight_search_cycles=1", "data_weight=10"])
  # refine.run writes real_space_refined.pdb and returns before pdb/ output.
  assert r_qr.std_out_text.find(
    "alternative conformations require minimizer=lbfgs") > -1
  refined = os.path.join(r_qr.folder, "real_space_refined.pdb")
  assert os.path.isfile(refined)
  h = get_model(refined).get_hierarchy()
  h_poor = get_model(pdb_poor).get_hierarchy()
  assert conformers.conformer_labels(h_poor) == ["A", "B"]
  assert sorted(h.altloc_indices().keys()) == ["", "A", "B"]
  assert h.atoms_size() == h_poor.atoms_size()
  assert approx_equal(h.atoms().extract_occ(), h_poor.atoms().extract_occ())
  cc = cc_mask_values(r_qr.std_out_text)
  print("CC_mask:", cc)
  assert len(cc) >= 2
  assert cc[-1] > cc[0], cc

if(__name__ == "__main__"):
  run()
