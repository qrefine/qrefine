# LIBTBX_SET_DISPATCHER_NAME qr.development.status
import os, sys
import time
from collections import OrderedDict

results = {}

file_status = OrderedDict([
  # ORCA
  ('ORCA SCF GRADIENT CALCULATION', 'Started Gradients'),
  ('***  Starting incremental Fock matrix formation  ***', 'Started SCF'),
  ('TOTAL RUN TIME:', 'ORCA DONE'),
  # MOPAC
  ('== MOPAC DONE ==', 'MOPAC DONE'),
  ])

def process_file(filename):
  print 'process_file',filename
  f=file(filename, 'rb')
  lines = f.read()
  f.close()
  done=False
  for line in lines.splitlines():
    #print line
    for lookup, answer in file_status.items():
      if line.find(lookup)>-1:
        done=answer
  print line
  return done

def process_dir(d, engine_name='orca'):
  for filename in os.listdir(d):
    if not filename.endswith('.out'): continue
    rc = process_file(os.path.join(d, filename))
    results.setdefault(os.path.join(d,filename), None)
    results[os.path.join(d,filename)] = rc

def check_output_pdbs(d):
  for filename in os.listdir(d):
    if filename.find('cycle.pdb')>-1:
      if filename.find('weight')>-1:
        results.setdefault('weight', [])
        results['weight'].append(os.path.join(d, filename))
      elif filename.find('refine')>-1:
        results.setdefault('refine', [])
        results['refine'].append(os.path.join(d, filename))
    elif filename.find('refined.pdb')>-1:
      results.setdefault('refined', [])
      results['refined'].append(os.path.join(d, filename))


def show_item(files):
  def _cmp_mtime(f1, f2):
    if os.stat(f1).st_mtime<os.stat(f2).st_mtime:
      return -1
    return 1

  files.sort(_cmp_mtime)
  for f in files:
    print '  %s : %s' % (f, time.asctime(time.localtime(os.stat(f).st_mtime)))

def run(cwd=None):
  if cwd is None:
    cwd=os.getcwd()
  else:
    os.chdir(cwd)
  for root, dirs, files in os.walk(cwd):
    if os.path.basename(root)=='ase_error':
      print 'ASE error directory exists'
    elif root.find('ase_error')>-1:
      pass
    elif root.find('ase')>-1:
      try:
        i = int(os.path.basename(root))
      except ValueError, e:
        i = None
      if i is not None:
        process_dir(root)

    if os.path.basename(root)=='pdb':
      check_output_pdbs(root)

  print '_'*80
  print '\nResults'
  print '_'*80
  print
  weight_dates = []
  refine_dates = []
  for key, item in sorted(results.items()):
    if key in ['weight', 'refine', 'refined']:
      print key
      show_item(item)
      for f in item:
        if key=='weight': weight_dates.append(os.stat(f).st_mtime)
        if key=='refine': refine_dates.append(os.stat(f).st_mtime)
    else:
      print '  %s : "%s"' % (key, item)
  if weight_dates and refine_dates:
    if max(weight_dates)>min(refine_dates):
      print '\n\n. *** refine output models are older than weight models ***\n\n'

if __name__=='__main__':
  run(*tuple(sys.argv[1:]))
