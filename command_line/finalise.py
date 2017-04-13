# LIBTBX_SET_DISPATCHER_NAME qr.finalise

import os, sys, shutil
from libtbx import easy_run
from multiprocessing import Pool

pdb_dir='./tmp'

def callback(args):
  print "call back"
  return args
  
def _process_pdb_filename(pdb_file):
  os.chdir(pdb_dir)
  complete_file=pdb_file[:-4]+"_complete.pdb"
  if ( pdb_file.endswith("pdb")  and not os.path.exists(complete_file) ):
    cmd = "phenix.python ../../../qr-core/finalise.py %s" % pdb_file + "> "+ pdb_file[:-4]+".log"
    print '\n\t~> %s\n' % cmd
    easy_run.call(cmd)
  os.chdir("../")
  return None

def run(folder,
    nproc=8,
    only_code=None, 
    ):
  pdb_dir = './tmp/' 
  filenames = os.listdir(pdb_dir)
  assert not filenames, 'script must be run in a empty directory'
  cmd = "cp "+ folder + "/./*.pdb "+ pdb_dir
  os.system(cmd)
  try: nproc=int(nproc)
  except: nproc=1  
  pool = None
  if nproc>1:
    pool = Pool(processes=nproc)
  os.chdir(pdb_dir)
  pdb_dir="./"
  pdb_files = os.listdir(pdb_dir)
  for pdb_file in pdb_files:
    cmd='mv  '+ pdb_file +'  ' + pdb_file[:4]+'.pdb'
    os.system(cmd)
  pdb_files = os.listdir(pdb_dir)
  print "pdbs need to be completed:"
  print pdb_files
  for pdb_file in pdb_files:
    if not pdb_file.endswith(".pdb"): continue
    if len(pdb_file.split('.'))!=2: continue
    if len(pdb_file.split('.')[0])!=4: continue
    if only_code is not None and only_code!=pdb_file.split('.')[0]: continue
    if nproc==1:
      _process_pdb_filename(pdb_file)
    else:
      rc = pool.apply_async(
        _process_pdb_filename,
        [pdb_file],
        callback=callback,
        )
  if pool:
    pool.close()
    pool.join()
  pdb_files = os.listdir(pdb_dir)
  ok_pdbs = []
  for pdb_file in pdb_files:
    if "complete.pdb" in pdb_file:
      ok_pdbs.append(pdb_file[:4]) 
  print ok_pdbs
  os.chdir("../")
  return ok_pdbs
  
if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
