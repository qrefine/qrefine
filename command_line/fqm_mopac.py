from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME qr.fqm_mopac
import sys
import time
import os.path
import libtbx
import iotbx.pdb
from qrefine.fragment import fragments
from qrefine.fragment import fragment_extracts
from qrefine.fragment import get_qm_file_name_and_pdb_hierarchy
from qrefine.fragment import charge
from qrefine.fragment import write_mm_charge_file
from qrefine.plugin.ase.mopac_qr import Mopac
from qrefine.restraints import ase_atoms_from_pdb_hierarchy

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_yoink_path =os.path.join(qrefine_path, "plugin","yoink","yoink")

log = sys.stdout

legend = """\
Cluster a system into many small pieces
"""

def run(pdb_file, log):
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  print >> log, '\n\tfragmenting "%s"' % pdb_file
  t0=time.time()
  fq = fragments(
    pdb_hierarchy=ph,
    crystal_symmetry=cs,
    charge_embedding=True,
    debug=True,
    qm_engine_name="mopac")
  print >> log, '\n\tfragmenting took %0.1f\n' % (time.time()-t0)
  print >> log, "Residue indices for each cluster:\n", fq.clusters
  fq_ext = fragment_extracts(fq)
  qm_engine = Mopac()
  for i in xrange(len(fq.clusters)):
      # add capping for the cluster and buffer
      print >> log, "capping frag:", i
      qm_pdb_file, ph = get_qm_file_name_and_pdb_hierarchy(
                          fragment_extracts=fq_ext,
                          index=i)
      qm_charge = charge(fragment_extracts=fq_ext,
                                      index=i)
      atoms = ase_atoms_from_pdb_hierarchy(ph)
      qm_engine.label = qm_pdb_file[:-4]
      qm_engine.run_qr(atoms,charge=qm_charge, pointcharges=None,
          coordinates=qm_pdb_file[:-4]+".xyz", define_str=None)
      energy = self.qm_engine.energy_free*unit_convert
      print >> log, " frag:", i, "  energy:", energy 


if (__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args[0], log)
  print >> log, "Time: %6.4f" % (time.time() - t0)
