from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
import time
import iotbx.pdb
from libtbx.utils import null_out
from qrefine.tests.unit import run_qrefine_in

# ALA-ALA, both residues with A and B conformers, B = A shifted by 4 A along
# z (as in cctbx iotbx/pdb/tst_hierarchy_reconcile_altlocs.py dipeptide()).
# Residue 2 has the wrong labels: its A conformer is the shifted one, so the
# C(1)-N(2) peptide bond is broken in both conformers.
pdb_str_swapped = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N  AALA A   1      -1.521   3.898   0.713  0.50 20.00           N
ATOM      1  CA AALA A   1      -1.879   2.678   1.370  0.50 20.00           C
ATOM      1  CB AALA A   1      -3.346   2.715   1.747  0.50 20.00           C
ATOM      1  C  AALA A   1      -1.551   1.486   0.458  0.50 20.00           C
ATOM      1  O  AALA A   1      -2.007   1.439  -0.662  0.50 20.00           O
ATOM      1  N  BALA A   1      -1.521   3.898   4.713  0.50 20.00           N
ATOM      1  CA BALA A   1      -1.879   2.678   5.370  0.50 20.00           C
ATOM      1  CB BALA A   1      -3.346   2.715   5.747  0.50 20.00           C
ATOM      1  C  BALA A   1      -1.551   1.486   4.458  0.50 20.00           C
ATOM      1  O  BALA A   1      -2.007   1.439   3.338  0.50 20.00           O
ATOM      1  N  BALA A   2      -0.744   0.404   0.952  0.50 20.00           N
ATOM      1  CA BALA A   2      -0.493  -0.749   0.119  0.50 20.00           C
ATOM      1  CB BALA A   2      -1.476  -1.885   0.331  0.50 20.00           C
ATOM      1  C  BALA A   2       0.972  -1.193   0.100  0.50 20.00           C
ATOM      1  O  BALA A   2       1.603  -1.231   1.129  0.50 20.00           O
ATOM      1  N  AALA A   2      -0.744   0.404   4.952  0.50 20.00           N
ATOM      1  CA AALA A   2      -0.493  -0.749   4.119  0.50 20.00           C
ATOM      1  CB AALA A   2      -1.476  -1.885   4.331  0.50 20.00           C
ATOM      1  C  AALA A   2       0.972  -1.193   4.100  0.50 20.00           C
ATOM      1  O  AALA A   2       1.603  -1.231   5.129  0.50 20.00           O
END
"""

def labels_of(pdb_file, resseq):
  """
  (altloc, z) of every atom of residue resseq, in atom order.
  """
  h = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy()
  result = []
  for rg in h.residue_groups():
    if(rg.resseq_as_int() != resseq): continue
    for ag in rg.atom_groups():
      for a in ag.atoms():
        result.append((ag.altloc, a.xyz[2]))
  return result

def check_labels(pdb_file, resseq, unshifted, shifted):
  """
  Atoms of the unshifted conformer (z < 3) carry label unshifted, those of
  the shifted one label shifted.
  """
  labels = labels_of(pdb_file, resseq)
  assert len(labels) == 10, labels
  for altloc, z in labels:
    if(z < 3): assert altloc == unshifted, labels
    else:      assert altloc == shifted, labels

def atom_names(pdb_file):
  """
  (resseq, name) of every atom, in file order.
  """
  h = iotbx.pdb.input(file_name=pdb_file).construct_hierarchy(
    sort_atoms=False)
  return [(a.parent().parent().resseq, a.name) for a in h.atoms()]

def run(prefix = "qrefine_"+os.path.basename(__file__).replace(".py","")):
  """
  qr.refine relabels altlocs of interacting residues before restraints are
  built (altlocs.reconcile=True, default); altlocs.reconcile=False keeps
  the input labels.
  """
  os.makedirs(prefix, exist_ok=True)
  os.chdir(prefix)
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_swapped)
  h = pdb_inp.construct_hierarchy()
  rg2 = list(h.residue_groups())[1]
  assert rg2.resseq_as_int() == 2
  r = h.reconcile_altlocs(crystal_symmetry=pdb_inp.crystal_symmetry(),
    check_only=True, log=null_out())
  assert r.swaps == [(rg2.id_str(), {"A":"B", "B":"A"})], r.swaps
  assert r.unresolved == []
  pdb_swapped = os.path.abspath("%s_swapped.pdb"%prefix)
  with open(pdb_swapped, "w") as fo:
    fo.write(pdb_str_swapped)
  check_labels(pdb_swapped, 2, unshifted="B", shifted="A")
  common = [pdb_swapped, "restraints=cctbx", "mode=opt", "clustering=False",
    "number_of_micro_cycles=1", "max_iterations_refine=1"]
  # reconcile on: residue 2 relabeled to the pairing of residue 1
  r_qr = run_qrefine_in(folder="on", prefix=prefix+"_on", args=common)
  assert r_qr.std_out_text.find("Reconciling altlocs") > -1
  assert r_qr.std_out_text.find("1 residues relabeled") > -1
  check_labels(r_qr.pdb, 1, unshifted="A", shifted="B")
  check_labels(r_qr.pdb, 2, unshifted="A", shifted="B")
  names_on = atom_names(r_qr.pdb)
  # reconcile off: labels as in the input
  r_qr = run_qrefine_in(folder="off", prefix=prefix+"_off",
    args=common+["altlocs.reconcile=False"])
  assert r_qr.std_out_text.find("Reconciling altlocs") == -1
  check_labels(r_qr.pdb, 1, unshifted="A", shifted="B")
  check_labels(r_qr.pdb, 2, unshifted="B", shifted="A")
  # The relabel does not change the atom order of the output (both keep the
  # input-file order).
  names_off = atom_names(r_qr.pdb)
  assert names_on == names_off, (names_on, names_off)
  assert names_off == atom_names(pdb_swapped)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
