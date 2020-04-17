from __future__ import print_function
import os, sys, shutil
from libtbx import easy_run
from multiprocessing import Pool
from StringIO import StringIO
import libtbx.load_env

qrefine_path = libtbx.env.find_in_repositories("qrefine")
pdb_dir='./tmp'

results = {}

from qrefine.tests.unit.skip import skip

def callback(args):
  #print "call back",args
  if args[0]!=0:
    results[args[1]]=args[0]
  return args

def generate_test_pdb_filenames(d):
  for filename in os.listdir(d):
    print(d, filename)
    if not filename.endswith(".pdb"): continue
    #if len(filename.split('.'))!=2: continue
    #if len(filename.split('.')[0])!=4: continue
    tf = "%s.pdb" % filename[:-4]
    print(tf)
    shutil.copyfile(os.path.join(d, filename), tf)
    if tf in skip:
      print('\n\tSKIPPING %s\n' % tf)
      continue
    yield tf

def _process_pdb_filename(pdb_file):
  print('_process_pdb_filename',pdb_file)
  complete_file=pdb_file[:-4]+"_complete.pdb"
  rc=0
  if ( pdb_file.endswith("pdb") and not os.path.exists(complete_file) ):
    cmd = "qr.finalise %s | tee %s" % (
      pdb_file,
      pdb_file.replace('.pdb', ".log"),
      )
    print('\n\t~> %s\n' % cmd)
    ero = easy_run.fully_buffered(command=cmd)
    err = StringIO()
    ero.show_stderr(out=err)
    rc = err.getvalue()
    if rc=='': rc=0
  return rc, pdb_file

def run(folder,
        nproc=8,
        only_code=None,
      ):
  filenames = os.listdir(os.getcwd())
  try: nproc=int(nproc)
  except: nproc=1
  pool = None
  if nproc>1:
    pool = Pool(processes=nproc)
  for pdb_file in generate_test_pdb_filenames(folder):
    print(pdb_file)
    if only_code is not None and only_code!=pdb_file.split('.')[0]: continue
    if nproc==1:
      rc, pdb_file = _process_pdb_filename(pdb_file)
      print('rc',rc)
      assert rc==0, 'return code != 0'
    else:
      rc = pool.apply_async(
        _process_pdb_filename,
        [pdb_file],
        callback=callback,
        )
    #break
  if pool:
    pool.close()
    pool.join()

  pdb_files = os.listdir(os.getcwd())
  ok_pdbs = []
  for pdb_file in pdb_files:
    if "complete.pdb" in pdb_file:
      ok_pdbs.append(pdb_file[:4])
  print(ok_pdbs)
  print(results)
  for pdb, rc in sorted(results.items()):
    print('-'*80)
    print(pdb)
    print(rc)
  return ok_pdbs

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
