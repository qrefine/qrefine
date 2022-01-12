from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME qr.fqm_mopac
import sys
import time
import os.path
import libtbx
import iotbx.pdb
import cPickle as pickle
from qrefine.fragment import fragments
from qrefine.fragment import fragment_extracts
from qrefine.fragment import get_qm_file_name_and_pdb_hierarchy
from qrefine.fragment import charge
from qrefine.fragment import write_mm_charge_file
from qrefine.plugin.ase.mopac_qr import Mopac
from qrefine.plugin.ase.orca_qr import Orca
from qrefine.restraints import ase_atoms_from_pdb_hierarchy

qrefine_path = libtbx.env.find_in_repositories("qrefine")
qr_yoink_path =os.path.join(qrefine_path, "plugin","yoink","yoink")

log = sys.stdout

legend = """\
Cluster a system into many small pieces
"""

def print_time(s):
  outl = ''
  secs=None
  if s>60:
    secs = s%60
    s = s//60
  if secs:
    return '%02dm:%0.1fs' % (int(s), secs)
  return '%0.1fs' % s

def get_qm_energy(qm_engine, fragments_extracted, index):
    qm_pdb_file, ph = get_qm_file_name_and_pdb_hierarchy(
      fragment_extracts=fragments_extracted,
      index=index)
    qm_charge = charge(fragment_extracts=fragments_extracted,
                       index=index)
    atoms = ase_atoms_from_pdb_hierarchy(ph)
    qm_engine.label = qm_pdb_file[:-4]
    print('LABEL',qm_engine.label)
    print('qm_charge',qm_charge)
    qm_engine.run_qr(atoms,
                     charge=qm_charge,
                     pointcharges=None,
                     coordinates=qm_pdb_file[:-4]+".xyz",
                     define_str=None)
    energy = qm_engine.energy_free #*unit_convert
    return energy

class qm_energy_manager (object) :
  def __init__ (self, qm_engine, fragments_extracted) :
    self.qm_engine = qm_engine
    self.fragments_extracted = fragments_extracted

  def __call__ (self, index) :
    return get_qm_energy(self.qm_engine,
                         self.fragments_extracted,
                         index,
    )

def run_all (qm_engine,
             fragments_extracted,
             indices,
             method='multiprocessing',
             processes=1,
             qsub_command=None,
             callback=None) :
  qm_engine_object = qm_energy_manager(qm_engine, fragments_extracted)
  from libtbx.easy_mp import parallel_map
  return parallel_map(
    func=qm_engine_object,
    iterable=indices,
    method=method,
    processes=processes,
    callback=callback,
    qsub_command=qsub_command)

def run(pdb_file, log):
  pdb_inp = iotbx.pdb.input(pdb_file)
  ph = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  print('\n\tfragmenting "%s"' % pdb_file, file=log)
  t0=time.time()
  #
  # use a pickle file of fragments so testing of parallel is snappy
  #
  fq_fn = '%s_fq.pickle' % pdb_file.replace('.pdb','')
  if not os.path.exists(fq_fn):
    fq = fragments(
      pdb_hierarchy=ph,
      crystal_symmetry=cs,
      charge_embedding=True,
      debug=False,
      qm_engine_name="orca")
    print('\n\tfragmenting took %0.1f\n' % (time.time()-t0), file=log)
    print("Residue indices for each cluster:\n", fq.clusters, file=log)
    fq_ext = fragment_extracts(fq)
    f=open(fq_fn, 'w')
    pickle.dump(fq_ext, f)
    f.close()
  else:
    f=open(fq_fn, 'r')
    fq_ext = pickle.load(f)
    f.close()
  #
  # get QM engine
  #
  qm_engine = Orca()
  if 0:
    #
    # use parallel_map
    #
    #indices = range(len(fq.clusters))
    indices = range(len(fq_ext.fragment_selections))
    rc = run_all(qm_engine,
                 fq_ext,
                 indices,
                 processes=6,
    )
    print(rc)
  elif 0:
    #
    # use multi_core_run
    #
    from libtbx import easy_mp
    nproc=6
    argss = []
    results = []
    for i in range(len(fq_ext.fragment_selections)):
      argss.append([qm_engine, fq_ext, i])
    for args, res, err_str in easy_mp.multi_core_run(get_qm_energy,
                                                      argss,
                                                      nproc,
                                                      ):
      print('%sTotal time: %6.2f (s)' % (' '*7, time.time()-t0))
      print(args[-1], res)
      if err_str:
        print('Error output from %s' % args)
        print(err_str)
        print('_'*80)
      results.append([args, res, err_str])
    print('-'*80)
    for args, res, err_str in results:
      print(args[-1], res)
  else:
    #
    # serial
    #
    for i in range(len(fq_ext.fragment_selections)):
      # add capping for the cluster and buffer
      print("capping frag:", i, file=log)
      energy = get_qm_energy(qm_engine, fq_ext, i)
      print(" frag:", i, "  energy:", energy, ' time:', print_time(time.time()-t0), file=log)
      time.sleep(10)

if (__name__ == "__main__"):
  t0 = time.time()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args[0], log)
  print("Time: %s" % print_time(time.time() - t0), file=log)
